/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This class converts conservative to primitive variables
#ifndef PRIMITIVERECOVERY_HPP_
#define PRIMITIVERECOVERY_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "DefaultEoS.hpp"
#include "EoS.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "UsingNamespace.H"
#include "VarsTools.hpp"
#include "simd.hpp"

// template <class eos_t = DefaultEoS>
class PrimitiveRecovery : public EoS
{
  public:
    template <class data_t> struct Vars
    {
        Tensor<2, data_t> h;
        data_t chi, D, tau, rho, eps;
        Tensor<1, data_t> Sj, vi;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function);
    };

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto vars = current_cell.template load_vars<Vars>();

        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);

        data_t tolerance = 1e-8;
        data_t Wa2, Wa, xn, diff;

        data_t P_of_rho = 0.;

        data_t S2 = 0.;
        FOR(i, j) S2 += vars.chi * h_UU[i][j] * vars.Sj[i] * vars.Sj[j];

        data_t q = vars.tau / vars.D;
        data_t r = S2 / pow(vars.D, 2.);
        data_t xa = 1.5 * (1. + q);
        Wa = sqrt(pow(xa, 2.) / (pow(xa, 2.) - r));

        vars.rho = pow(vars.chi, 1.5) * vars.D / Wa;
        vars.eps = -1. + xa / Wa * (1. - Wa * Wa) + Wa * (1. + q);
        EoS::compute_eos(P_of_rho, vars);
        xn = Wa * (1. + vars.eps + P_of_rho / vars.rho);
        diff = abs(xn - xa);

        int i = 0;

        data_t empty;
        while (!simd_all_false(simd_compare_lt(tolerance, diff), empty))
        {
            i++;
            Wa = sqrt(pow(xa, 2.) / (pow(xa, 2.) - r));

            vars.rho = pow(vars.chi, 1.5) * vars.D / Wa;
            vars.eps = -1. + xa / Wa * (1. - Wa * Wa) + Wa * (1. + q);
            EoS::compute_eos(P_of_rho, vars);
            xn = Wa * (1. + vars.eps + P_of_rho / vars.rho);
            diff = abs(xn - xa);
            xa = xn;
            if (i >= 100)
                break;
        }
        data_t hh = 1. + vars.eps + P_of_rho / vars.rho;
        Tensor<1, data_t> vi_D;
        FOR(i) vi_D[i] = vars.Sj[i] / vars.rho / hh / Wa / Wa;

        FOR(i)
        {
            vars.vi[i] = 0.;
            FOR(j) vars.vi[i] += vars.chi * h_UU[i][j] * vi_D[j];
        }

        current_cell.store_vars(vars);
    }
};

template <class data_t>
template <typename mapping_function_t>
void PrimitiveRecovery::Vars<data_t>::enum_mapping(
    mapping_function_t mapping_function)
{
    VarsTools::define_symmetric_enum_mapping(mapping_function,
                                             GRInterval<c_h11, c_h33>(), h);
    VarsTools::define_enum_mapping(mapping_function, c_chi, chi);
    VarsTools::define_enum_mapping(mapping_function, c_D, D);
    VarsTools::define_enum_mapping(mapping_function, GRInterval<c_Sj1, c_Sj3>(),
                                   Sj);
    VarsTools::define_enum_mapping(mapping_function, c_tau, tau);
    VarsTools::define_enum_mapping(mapping_function, GRInterval<c_vi1, c_vi3>(),
                                   vi);
    VarsTools::define_enum_mapping(mapping_function, c_rho, rho);
    VarsTools::define_enum_mapping(mapping_function, c_eps, eps);
}

#endif /* PRIMITIVERECOVERY_HPP_ */

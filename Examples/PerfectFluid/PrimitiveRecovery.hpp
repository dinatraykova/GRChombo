/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This class converts conservative to primitive variables
#ifndef PRIMITIVERECOVERY_HPP_
#define PRIMITIVERECOVERY_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

class PrimitiveRecovery
{
  public:
    template <class data_t> struct Vars
    {
        Tensor<2, data_t> h;
        data_t chi, D, tau, rho, eps, nn, Jt;
        Tensor<1, data_t> Sj, vi;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function);
    };

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto vars = current_cell.template load_vars<Vars>();

        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);

        data_t E = vars.tau + vars.D;
        data_t S2 = 0.;
        FOR(i, j) S2 += vars.chi * h_UU[i][j] * vars.Sj[i] * vars.Sj[j];

        data_t E2 = E * E;
        data_t sqrt1 = sqrt(4. * E2 - 3. * S2);
        data_t sqrt2 = sqrt(2. * E2 - 3. * S2 + E * sqrt1);

        // eps
        vars.eps =
            sqrt2 / (2. * vars.D) - 1.; // Is it dangerous to divide by D?
        // rho
        vars.rho = (sqrt1 - E) / (1. + vars.eps);
        // vi_D
        Tensor<1, data_t> vi_D;
        FOR(i) vi_D[i] = 3. * vars.Sj[i] * (sqrt1 - E) / (sqrt2 * sqrt2);
        // vi
        FOR(i)
        {
            vars.vi[i] = 0.;
            FOR(j) vars.vi[i] += vars.chi * h_UU[i][j] * vi_D[j];
        }

        data_t v2 = 0.;
        FOR(i) v2 += vars.vi[i] * vi_D[i];

        vars.nn = vars.Jt / sqrt(1. - v2);

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
    VarsTools::define_enum_mapping(mapping_function, c_Jt, Jt);
    VarsTools::define_enum_mapping(mapping_function, GRInterval<c_Sj1, c_Sj3>(),
                                   Sj);
    VarsTools::define_enum_mapping(mapping_function, c_tau, tau);
    VarsTools::define_enum_mapping(mapping_function, GRInterval<c_vi1, c_vi3>(),
                                   vi);
    VarsTools::define_enum_mapping(mapping_function, c_rho, rho);
    VarsTools::define_enum_mapping(mapping_function, c_eps, eps);
    VarsTools::define_enum_mapping(mapping_function, c_nn, nn);
}

#endif /* PRIMITIVERECOVERY_HPP_ */

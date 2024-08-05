/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This compute class enforces the dominant energy conditions
#ifndef POSITIVEDENSITY_HPP_
#define POSITIVEDENSITY_HPP_

#include "Cell.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"

class PositiveDensity
{
    // Only variables needed are chi, D, tau, Sj, h
    template <class data_t> struct Vars
    {
        data_t chi;
        data_t D;
        data_t tau;
        Tensor<1, data_t> Sj;
        Tensor<2, data_t> h;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
            define_enum_mapping(mapping_function, c_D, D);
            define_enum_mapping(mapping_function, c_tau, tau);
            define_enum_mapping(mapping_function, GRInterval<c_Sj1, c_Sj3>(),
                                Sj);
            define_symmetric_enum_mapping(mapping_function,
                                          GRInterval<c_h11, c_h33>(), h);
        }
    };

  private:
    const double m_min_D;
    const double m_min_v;

  public:
    //! Constructor for class
    PositiveDensity(const double a_min_D = 1e-12, const double a_min_v = 0.)
        : m_min_D(a_min_D), m_min_v(a_min_v)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto vars = current_cell.template load_vars<Vars>();
        const data_t chi_regularised = simd_max(1e-6, vars.chi);

        // Take into account that conservative variables are conformally
        // rescaled

        // D >= min_D
        vars.D = simd_max(vars.D, m_min_D / pow(chi_regularised, 1.5));

        data_t S2_over_chi = 0.;
        FOR(i, j) S2_over_chi += vars.h[i][j] * vars.Sj[i] * vars.Sj[j];

        //  S^2 <= (D + tau)^2
        data_t factor =
            simd_min(1., vars.chi * (vars.D + vars.tau) / sqrt(S2_over_chi));
        FOR(i) vars.Sj[i] *= factor;

        current_cell.store_vars(vars.D, c_D);
        current_cell.store_vars(vars.Sj, GRInterval<c_Sj1, c_Sj3>());
    }
};

#endif /* POSITIVEDENSITY_HPP_ */

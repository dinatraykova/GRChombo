/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This compute class enforces the positive chi and alpha condition
#ifndef POSITIVEDENSITY_HPP_
#define POSITIVEDENSITY_HPP_

#include "Cell.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"

class PositiveDensity
{
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
        auto D = current_cell.load_vars(c_D);
        auto rho = current_cell.load_vars(c_rho);
        auto tau = current_cell.load_vars(c_tau);
        Tensor<1, data_t> vi, Sj;
        vi[0] = current_cell.load_vars(c_vi1);
        vi[1] = current_cell.load_vars(c_vi2);
        vi[2] = current_cell.load_vars(c_vi3);

        auto make_zero = simd_compare_lt(D, m_min_D);
        D = simd_conditional(make_zero, D, m_min_D);
        // tau = simd_conditional(make_zero, tau, 1e-4);
        FOR(i) vi[i] = simd_conditional(make_zero, vi[i], m_min_v);

        // current_cell.store_vars(D, c_D);
        //  current_cell.store_vars(rho, c_rho);
        // current_cell.store_vars(vi, GRInterval<c_vi1, c_vi3>());
        //  current_cell.store_vars(Sj, GRInterval<c_Sj1, c_Sj3>());
        //  current_cell.store_vars(tau, c_tau);
    }
};

#endif /* POSITIVEDENSITY_HPP_ */

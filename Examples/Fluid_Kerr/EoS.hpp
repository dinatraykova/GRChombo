/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EOS_HPP_
#define EOS_HPP_

#include "simd.hpp"

class EoS
{
  public:
    struct params_t
    {
        double eos_w;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    EoS(params_t a_params) : m_params(a_params) {}

    //! Set the pressure of the perfect fluid here
    template <class data_t, template <typename> class vars_t>
    void compute_eos(data_t &P_over_rho, const vars_t<data_t> &vars) const
    {
        // The pressure value in function of rho
        P_over_rho = m_params.eos_w * (1. + vars.eps);
    }
};

#endif /* EOS_HPP_ */

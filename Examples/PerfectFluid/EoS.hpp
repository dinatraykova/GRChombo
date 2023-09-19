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
    void compute_eos(data_t &P_of_rho, data_t &dPdrho,
                     const vars_t<data_t> &vars) const
    {
        double vareps = 0.;
        // The pressure value in function of rho
        P_of_rho = m_params.eos_w * vars.rho * (1. + vareps);

        // The pressure gradient wrt rho
        dPdrho = m_params.eos_w * (1. + vareps);
    }
};

#endif /* EOS_HPP_ */

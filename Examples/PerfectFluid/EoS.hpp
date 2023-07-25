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

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_eos(data_t &P_of_rho0, data_t &dPdrho0,
                           const vars_t<data_t> &vars) const
    {
      // The pressure value at rho0
      // w rho0
      P_of_rho0 = m_params.eos_w * vars.rho0;

      // The pressure gradient at rho0
      dPdrho0 = m_params.eos_w;
    }
};

#endif /* EOS_HPP_ */

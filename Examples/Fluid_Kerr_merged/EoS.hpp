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
    //! The constructor
    EoS() {}

    //! Set the pressure of Gamma-2 polytrope
    template <class data_t, template <typename> class vars_t>
    void compute_eos(data_t &P_of_rho, const vars_t<data_t> &vars) const
    {
        // The pressure value as a function of rho
        const double K = 1. / 3.; // 100.;
        const double n = 1.;
        const double Gamma = 1.; // + 1. / n;
        P_of_rho = K * pow(vars.rho, Gamma);
    }
};

#endif /* EOS_HPP_ */

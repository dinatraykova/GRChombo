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

    //! Set the pressure of the perfect fluid here
    template <class data_t, template <typename> class vars_t>
    void compute_eos(data_t &P_over_rho, const vars_t<data_t> &vars) const
    {
        // The pressure value as a function of rho
        P_over_rho = (1. + vars.eps) / 3.;
    }
};

#endif /* EOS_HPP_ */

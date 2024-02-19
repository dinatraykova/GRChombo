/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DEFAULTEOS_HPP_
#define DEFAULTEOS_HPP_

#include "Tensor.hpp"
#include "simd.hpp"

class DefaultEoS
{
  public:
    //! The constructor
    DefaultEoS() {}

    //! Set the pressure of the perfect fluid here to zero
    template <class data_t, template <typename> class vars_t>
    void compute_eos(data_t &P_over_rho, const vars_t<data_t> &vars) const
    {
        // The pressure value in function of rho
        P_over_rho = 1. / 3.;
    }
};

#endif /* DEFAULTEOS_HPP_ */

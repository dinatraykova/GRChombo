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

    //! Set the potential function for the scalar field here to zero
    template <class data_t, template <typename> class vars_t>
    void compute_eos(data_t &P_of_rho0, data_t &dPdrho0,
                           const vars_t<data_t> &vars) const
    {
        // The pressure value at rho0
        P_of_rho0 = 0.0;

	// The pressure gradient at rho0
	dPdrho0 = 0.0;
    }
};

#endif /* DEFAULTEOS_HPP_ */

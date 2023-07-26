/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALFLUIDDATA_HPP_
#define INITIALFLUIDDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "PerfectFluid.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "PrimitiveRecovery.hpp"

//! Class which sets the initial scalar field matter config
class InitialFluidData
{
  public:
    //! A structure for the input params for fluid properties and initial
    //! conditions
    struct params_t
    {
        double amplitude; //
        std::array<double, CH_SPACEDIM>
            center;   //!< Centre of perturbation in initial SF bubble
        double delta;
        double width; //!< Width of bump in initial SF bubble
    };

    //! The constructor
    InitialFluidData(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        data_t rr = coords.get_radius();
        data_t rr2 = rr * rr;
	
        // calculate the field value
        data_t rho = m_params.amplitude *
	                     (exp(-pow(rr / m_params.width, 2.0))) + m_params.delta;
	data_t v2 = 0.;
	data_t eps = 0.;
	data_t P = rho * (1. + eps) / 3.;
        data_t WW = 1./(1. - v2);
        data_t hh = 1. + eps + P / rho;

	data_t D = rho * sqrt(WW);
	data_t tau = rho * hh * WW - P - D;
        // store the vars
        current_cell.store_vars(D, c_D);
	current_cell.store_vars(tau, c_tau);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALFLUIDDATA_HPP_ */

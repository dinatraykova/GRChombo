/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALFLUIDDATA_HPP_
#define INITIALFLUIDDATA_HPP_

#include "ADMConformalVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FluidCCZ4RHS.hpp"
#include "PerfectFluid.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "PrimitiveRecovery.hpp"

//! Class which sets the initial fluid matter config
// aka Kevin data
class InitialFluidData
{
  // Now the non grid ADM vars
  template <class data_t>
  using MetricVars = ADMConformalVars::VarsWithGauge<data_t>;

  public:
    //! A structure for the input params for fluid properties and initial
    //! conditions
    struct params_t
    {
        std::array<double, CH_SPACEDIM>
	    center;   //!< Centre of perturbation in initial SF bubble
        double rho0;
        double uflow;
        double amplitude;
        double awidth;
        double sigma;
        std::array<double, 2>
	    ycenter;
    };

    //! The constructor
  InitialFluidData(params_t a_params, double a_dx)
        :  m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

	MetricVars<data_t> metric_vars;
	
	data_t x = coords.x;
	double y = coords.y;
	double z = coords.z;
	
	data_t vx, vy, vz, eps;
	data_t nn;

	vx = m_params.uflow*(tanh((y-m_params.ycenter[0])/m_params.awidth)
		    - tanh((y-m_params.ycenter[1])/m_params.awidth) - 1.);
	vy = m_params.amplitude * sin(2.*M_PI*x)
	  *(exp(-pow((y-m_params.ycenter[0])/m_params.sigma,2))
	    + exp(-pow((y-m_params.ycenter[1])/m_params.sigma,2)));
	vz = 0.;

	//DT: Thought eps was meant to be 0?
	nn = 1. + 0.5*(tanh((y-m_params.ycenter[0])/m_params.awidth)
		       - tanh((y-m_params.ycenter[1])/m_params.awidth));
	eps = 0.;
	//gamma_fac = 1./sqrt(1. - vx*vx - vy*vy - vz*vz);
	//ux = gamma_fac * vx;
	//uy = gamma_fac * vy;
	//uz = gamma_fac * vz;

        data_t rho = m_params.rho0;
	data_t v2 = vx*vx + vy*vy + vz*vz;
	data_t P = rho * (1. + eps) / 3.;
	data_t WW = 1./(1. - v2);
        data_t hh = 1. + eps + P / rho;
	
	data_t D = rho * sqrt(WW);
	data_t tau = rho * hh * WW - P - D;
	data_t Sj1 = rho * hh * WW * vx;
	data_t Sj2 = rho * hh * WW * vy;
	data_t Sj3 = rho * hh * WW * vz;
	
        // store the vars
	current_cell.store_vars(rho, c_rho);
	current_cell.store_vars(vx, c_vi1);
	current_cell.store_vars(vy, c_vi2);
	current_cell.store_vars(vz, c_vi3);
        current_cell.store_vars(D, c_D);
	current_cell.store_vars(Sj1, c_Sj1);
	current_cell.store_vars(Sj2, c_Sj2);
	current_cell.store_vars(Sj3, c_Sj3);
	current_cell.store_vars(tau, c_tau);
    }
protected:
  double m_dx;
  const params_t m_params;
};

#endif /* INITIALFLUIDDATA_HPP_ */

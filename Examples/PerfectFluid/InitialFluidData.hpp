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

//! Class which sets the initial fluid matter config
// aka Kevin data
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
        double rho0;
        double uflow;
        double amp;
        double awidth;
        double sigma;
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
	data_t x = coords.x;
	data_t y = coords.y;
	data_t z = coords.z;
	
	double xx, yy;
	double vx, vy, nn;
	double gamma_fac;
	double ux, uy;

	//DT: What is L? also ycenter?
        const double L[2] = {1.,2.};
	//double uflow = 1./(4.*sqrt(3));
	double ycenter[2] = {-0.5,0.5};
	//KevinHelmholtz(1.,1./(4.*sqrt(3)),0.01,0.05,0.2,ycenter);

	vx = uflow*(tanh((y-ycenter[0])/awidth)
		    - tanh((y-ycenter[1])/awidth) - 1.);
	vy = amp * sin(2.*M_PI*x/L[0])*(exp(-pow((y-ycenter[0])/sigma,2))
				       + exp(-pow((y-ycenter[1])/sigma,2)));
	vz = 0.;
	//DT: Thought eps was meant to be 0?
	esp = 1. + 0.5*(tanh((y-ycenter[0])/awidth) - tanh((y-ycenter[1])/awidth));
	gamma_fac = 1./sqrt(1-vx*vx-vy*vy - vz*vz);
	ux = gamma_fac * vx;
	uy = gamma_fac * vy;
	uz = gamma_fac * vz;

        data_t rho = rho0;
	data_t v2 = ux*ux + uy*uy + uz*uz;
	data_t P = rho * (1. + eps) / 3.;
        data_t WW = 1./(1. - v2);
        data_t hh = 1. + eps + P / rho;

	data_t D = rho * sqrt(WW);
	data_t tau = rho * hh * WW - P - D;
	data_t Sj1 = rho * hh * WW * ux;
	data_t Sj2 = rho * hh * WW * uy;
	data_t Sj3 = rho * hh * WW * uz;

        // store the vars
	current_cell.store_vars(rho, c_rho);
	current_cell.store_vars(u1, c_vi1);
	current_cell.store_vars(u2, c_vi2);
	current_cell.store_vars(u3, c_vi3);
	current_cell.store_vars(eps, c_eps);
        current_cell.store_vars(D, c_D);
	current_cell.store_vars(Sj3, c_Sj1);
	current_cell.store_vars(Sj2, c_Sj2);
	current_cell.store_vars(Sj1, c_Sj3);
	current_cell.store_vars(tau, c_tau);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALFLUIDDATA_HPP_ */

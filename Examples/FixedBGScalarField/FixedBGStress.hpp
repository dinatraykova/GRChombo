/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGSTRESS_HPP_
#define FIXEDBGSTRESS_HPP_

#include "ADMFixedBGVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "CoordinateTransformations.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the momentum flux S_i with type matter_t and writes it to the grid
template <class matter_t, class background_t> class FixedBGStress
{
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const matter_t m_matter;                        //!< The matter object
    const double m_dx;                              //!< The grid spacing
    const background_t m_background;                //!< The metric background
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center

  public:
    FixedBGStress(matter_t a_matter, background_t a_background, double a_dx,
                   std::array<double, CH_SPACEDIM> a_center)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and derivs
        const auto vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);

        // get the metric vars from the background
        MetricVars<data_t> metric_vars;
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        m_background.compute_metric_background(metric_vars, coords);

        using namespace TensorAlgebra;
	using namespace CoordinateTransformations;
	//	const auto gamma = metric_vars.gamma;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, metric_vars, d1, gamma_UU, chris_phys.ULL);
	const auto lapse = metric_vars.lapse;
	const auto shift = metric_vars.shift;
	const data_t det_gamma = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);

        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t R = coords.get_radius();
	const auto gamma_spher = cartesian_to_spherical_LL(metric_vars.gamma, x, y, z);

	Tensor<1, data_t> si;
	si[0] = x/R;
	si[1] = y/R;
	si[2] = z/R;
	
	data_t si_norm = 0.0;
	FOR2(j, k)
	  {
	    si_norm += si[j]*si[k]*metric_vars.gamma[j][k];
	  }

	FOR1(i)
	{
	  si[i] = si[i]/sqrt(si_norm);
	}

	data_t Stress = 0.0;
	FOR1(i)
	{
	  Stress += si[i]*emtensor.Sij[i][0];
	}
	Tensor<1, data_t> si_spher;
        si_spher[0] = 1.0/sqrt(gamma_spher[0][0]);
        si_spher[1] = 0.0;
        si_spher[2] = 0.0;
	
	Tensor<2, data_t> Proj_spher;

	FOR2(i, j)
	  {
            Proj_spher[i][j] = delta(i, j);
	    FOR1(k)
            {
                Proj_spher[i][j] += -gamma_spher[i][k] * si_spher[k] * si_spher[j];
            }
	  }


        Tensor<2, data_t> Sigma;
	FOR2(i, j)
	  {
            Sigma[i][j] = 0.0;
	    FOR2(m, n)
	      {
                Sigma[i][j] +=
		  Proj_spher[i][m] * Proj_spher[j][n] * gamma_spher[m][n];
	      }
	  }

	Tensor<1, data_t> Ni;
        Ni[0] = coords.x / R;
        Ni[1] = coords.y / R;
        Ni[2] = coords.z / R;

        // The integrand for the x-momentum flux out of a radial
        // shell at the current position
        data_t Mdot = 0;
        FOR1(i)
        {
	  Mdot += metric_vars.lapse * emtensor.Sij[0][i] * Ni[i];
        }
        Mdot *= sqrt(det_gamma);

        // assign values of conserved density in output box,
        // including factors of lapse and det_gamma
        //current_cell.store_vars(Mdot, c_Stress);

	//const data_t detSigma = Sigma[1][1] * Sigma[2][2] - Sigma[1][2] * Sigma[2][1];
	const data_t dArea = sqrt(Sigma[1][1] * Sigma[2][2] - Sigma[1][2] * Sigma[2][1]);
	//pout()<< "Lapse in stress vars" << metric_vars.lapse <<endl;
	//pout()<< "Stress" << Stress <<endl;
	//pout()<< "dArea" << dArea;
        // assign values of Momentum flux in output box
	current_cell.store_vars(emtensor.rho, c_rho);
        current_cell.store_vars(Stress, c_Stress);
	current_cell.store_vars(dArea, c_dArea);
    }
};

#endif /* FIXEDBGSTRESS_HPP_ */

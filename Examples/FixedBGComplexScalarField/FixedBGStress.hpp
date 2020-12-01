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
	const data_t det_gamma = compute_determinant_sym(metric_vars.gamma);

        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t R = coords.get_radius();
	const auto gamma_spher = cartesian_to_spherical_LL(metric_vars.gamma, x, y, z);
	
        // the unit vector in the radial direction
	Tensor<1, data_t> si_L;
	si_L[0] = x/R;
	si_L[1] = y/R;
	si_L[2] = z/R;

	data_t si_norm = 0.0;
	FOR2(i, j)
	{
	  si_norm += si_L[i] * si_L[j] * (gamma_UU[i][j]);
	  //	  pout()<<"si_norm = "<<si_norm<<endl;
	}
	//	pout()<<"si_norm = "<<si_norm<<endl; 

	FOR1(i)	{ si_L[i] = si_L[i]/sqrt(si_norm);}

	data_t Stress = 0.0;
	FOR1(i)
	{
          Stress += -shift[i] * si_L[i] * emtensor.Si[0];
	  FOR1(j)
	  {
	    Stress += lapse * si_L[i] * gamma_UU[i][j] * emtensor.Sij[0][j];
	    //pout()<<"Stress = "<<Stress<<endl;
	  }
	}
	
	const auto dArea = area_element_sphere(gamma_spher);
	
	data_t Xmom = -emtensor.Si[0] * sqrt(det_gamma);

	data_t S1 = -emtensor.rho * metric_vars.d1_lapse[0];
	data_t Source = S1;

	data_t S2 = 0;
	data_t S3 = 0;
	FOR1(i)
	{
	  S2 += emtensor.Si[i] * metric_vars.d1_shift[i][0];
	  Source += emtensor.Si[i] * metric_vars.d1_shift[i][0];
	  FOR2(j,k)
	    {
	      S3 += lapse * gamma_UU[i][k]*emtensor.Sij[k][j] *
                chris_phys.ULL[j][i][0];
	      Source += lapse * gamma_UU[i][k]*emtensor.Sij[k][j] *
                chris_phys.ULL[j][i][0];
	    }
	}

	S1 = S1 * sqrt(det_gamma);
	S2 = S2 * sqrt(det_gamma);
	S3 = S3 * sqrt(det_gamma);
	Source = Source * sqrt(det_gamma);
	data_t rho = emtensor.rho;
	auto cuto = simd_compare_gt(R, 200.);
	auto cuti = simd_compare_lt(R, 3.);

	rho = simd_conditional(cuto, 0.0, rho);
	rho = simd_conditional(cuti, 0.0, rho);
	Source = simd_conditional(cuto, 0.0, Source);
	Source = simd_conditional(cuti, 0.0, Source);
	S1 = simd_conditional(cuto, 0.0, S1);
        S1 = simd_conditional(cuti, 0.0, S1);
	S2 = simd_conditional(cuto, 0.0, S2);
        S2 = simd_conditional(cuti, 0.0, S2);
	S3 = simd_conditional(cuto, 0.0, S3);
        S3 = simd_conditional(cuti, 0.0, S3);
        Xmom = simd_conditional(cuto, 0.0, Xmom);
	Xmom = simd_conditional(cuti, 0.0, Xmom);

	current_cell.store_vars(emtensor.rho, c_rhofull);
	current_cell.store_vars(rho, c_rho);
	current_cell.store_vars(Source, c_Source);
	current_cell.store_vars(S1, c_S1);
	current_cell.store_vars(S2, c_S2);
	current_cell.store_vars(S3, c_S3);
	current_cell.store_vars(Xmom, c_Xmom);
        current_cell.store_vars(Stress, c_Stress);
	current_cell.store_vars(dArea, c_dArea);
    }
};

#endif /* FIXEDBGSTRESS_HPP_ */
/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WENODERIVATIVES_HPP_
#define WENODERIVATIVES_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include <array>

class WENODerivatives
{
  public:
    enum
    {
      LEFT_MINUS,
      LEFT_PLUS,
      RIGHT_MINUS,
      RIGHT_PLUS,
    };

  public:
    WENODerivatives(double a_eW=1.0): m_eW(a_eW){}

    template <class data_t>
    ALWAYS_INLINE data_t get_Pface(const double *in_ptr, const int idx,
				   const int stride, int dir_switch) const
    {
        double beta[3] = {0.,0.,0.};
        const double dd[3] = {3./10.,3./5.,1./10.};
	//const double eW = 1.;
        double alpha[3] = {0.,0.,0.};
        double sum_alpha;
        double weights[3] = {0.,0.,0.};
        double v[3] = {0.,0.,0.};
        double u[3] = {0.,0.,0.};
        double pim2, pim1, pi0, pip1, pip2;

	auto in = SIMDIFY<data_t>(in_ptr);

	
	if (dir_switch == 0) // LEFT_MINUS
	  {
	    pim2 = in[idx - 3.*stride];
	    pim1 = in[idx - 2.*stride];
	    pi0  = in[idx - stride];
	    pip1 = in[idx];
	    pip2 = in[idx + stride];

	    // ENO polynomials
	    v[0] = (2.*pim2 - 7.*pim1 + 11.*pi0 )/6.;
	    v[1] = (  -pim1 + 5.*pi0  +  2.*pip1)/6.;
	    v[2] = (2.*pi0  + 5.*pip1 -     pip2)/6.;

	    // Smoothness indicators
	    beta[0] = (13./12.)*pow(pim2-2.*pim1+pi0,2)
	      + (1./4.)*pow(pim2  -4.*pim1+3.*pi0,2);
	    beta[1] = (13./12.)*pow(pim1-2.*pi0 +pip1,2)
	      + (1./4.)*pow(pim1  -pip1,2);
	    beta[2] = (13./12.)*pow(pi0 -2.*pip1+pip2,2)
	      + (1./4.)*pow(3.*pi0-4.*pip1+pip2,2);

	    // Weights
	    sum_alpha = 0.;
	    FOR(j)
	    {
	      alpha[j] = dd[j]/((a_eW+beta[j])*(a_eW+beta[j]));
	      sum_alpha += alpha[j];
	    }
	    FOR(j){ weights[j] = alpha[j]/sum_alpha; }

	    // primitive variables on the left side of the
	    // cell interface i-1/2: p^{L-}_{i-1/2}
	    double p_LeftMinus = 0.;
	    FOR(j) p_LeftMinus += weights[j]*v[j];
	    return p_LeftMinus;
	  }

	if (dir_switch == 3) // LEFT_PLUS
	  {
	    // Choose negative fluxes to compute primitive variables
            // at the right cell boundary i+1/2: p^{L+}_{i-1/2}
            pim2 = [idx - 2.*stride];
            pim1 = [idx - stride];
            pi0  = [idx];
            pip1 = [idx + stride];
            pip2 = [idx + 2.*stride];

            // ENO polynomials
            u[0] = (   -pim2 + 5.*pim1 + 2.*pi0 )/6.;
            u[1] = ( 2.*pim1 + 5.*pi0  -    pip1)/6.;
            u[2] = (11.*pi0  - 7.*pip1 + 2.*pip2)/6.;

            // Smoothness return p_RightMinus;indicators
            beta[0] = (13./12.)*pow(pim2-2.*pim1+pi0,2)
              + (1./4.)*pow(pim2  -4.*pim1+3.*pi0,2);
            beta[1] = (13./12.)*pow(pim1-2.*pi0 +pip1,2)
              + (1./4.)*pow(pim1  -pip1,2);
            beta[2] = (13./12.)*pow(pi0 -2.*pip1+pip2,2)
              + (1./4.)*pow(3.*pi0-4.*pip1+pip2,2);

            // Weights
            sum_alpha = 0.;
            FOR(j)
            {
              alpha[j] = dd[j]/((a_eW+beta[j])*(a_eW+beta[j]));
              sum_alpha += alpha[j];
            }
            FOR(j)
            {
              weights[j] = alpha[j]/sum_alpha;
            }

            // primitive variables on the right side of the
            // cell interface i+1/2: p^{L+}_{i+1/2}
            double p_LeftPlus = 0.;
            FOR(j) {p_LeftPlus += weights[j]*v[j];}
            return p_LeftPlus;
	  }
	
	if (dir_switch == 2) // RIGHT_MINUS
	  {
	    pim2 = in[idx - 2.*stride];
	    pim1 = in[idx - stride];
	    pi0  = in[idx];
	    pip1 = in[idx + stride];
	    pip2 = in[idx + 2.*stride];

	    // ENO polynomials
	    v[0] = (2.*pim2 - 7.*pim1 + 11.*pi0 )/6.;
	    v[1] = (  -pim1 + 5.*pi0  +  2.*pip1)/6.;
	    v[2] = (2.*pi0  + 5.*pip1 -     pip2)/6.;

	    // Smoothness indicators
	    beta[0] = (13./12.)*pow(pim2-2.*pim1+pi0,2)
	      + (1./4.)*pow(pim2  -4.*pim1+3.*pi0,2);
	    beta[1] = (13./12.)*pow(pim1-2.*pi0 +pip1,2)
	      + (1./4.)*pow(pim1  -pip1,2);
	    beta[2] = (13./12.)*pow(pi0 -2.*pip1+pip2,2)
	      + (1./4.)*pow(3.*pi0-4.*pip1+pip2,2);

	    // Weights
	    sum_alpha = 0.;
	    FOR(j)
	    {
	      alpha[j] = dd[j]/((a_eW+beta[j])*(a_eW+beta[j]));
	      sum_alpha += alpha[j];
	    }
	    FOR(j){ weights[j] = alpha[j]/sum_alpha; }

	    // primitive variables on the left side of the
	    // cell interface i+1/2: p^{R-}_{i+1/2}
	    double p_RightMinus = 0.;
	    FOR(j) {p_RightMinus += weights[j]*v[j];}
	    return p_RightMinus;
	  }
	
	if (dir_switch == 3) // RIGHT_PLUS
	  {
	    // Choose negative fluxes to compute primitive variables
	    // at the right cell boundary i+1/2: p^{R+}_{i+1/2}
	    pim2 = [idx - stride];
	    pim1 = [idx];
	    pi0  = [idx + stride];
	    pip1 = [idx + 2.*stride];
	    pip2 = [idx + 3.*stride];

	    // ENO polynomials
	    u[0] = (   -pim2 + 5.*pim1 + 2.*pi0 )/6.;
	    u[1] = ( 2.*pim1 + 5.*pi0  -    pip1)/6.;
	    u[2] = (11.*pi0  - 7.*pip1 + 2.*pip2)/6.;
	    
	    // Smoothness return p_RightMinus;indicators
	    beta[0] = (13./12.)*pow(pim2-2.*pim1+pi0,2)
	      + (1./4.)*pow(pim2  -4.*pim1+3.*pi0,2);
	    beta[1] = (13./12.)*pow(pim1-2.*pi0 +pip1,2)
	      + (1./4.)*pow(pim1  -pip1,2);
	    beta[2] = (13./12.)*pow(pi0 -2.*pip1+pip2,2)
	      + (1./4.)*pow(3.*pi0-4.*pip1+pip2,2);
	    
	    // Weights
	    sum_alpha = 0.;
	    FOR(j)
	    {
	      alpha[j] = dd[j]/((a_eW+beta[j])*(a_eW+beta[j]));
	      sum_alpha += alpha[j];
	    }
	    FOR(j)
	    {
	      weights[j] = alpha[j]/sum_alpha;
	    }
	    
	    // primitive variables on the right side of the
            // cell interface i+1/2: p^{R+}_{i+1/2}
            double p_RightPlus = 0.;
            FOR(j) {p_RightPlus += weights[j]*v[j];}
            return p_RightPlus;
	  }
    }

    template <class data_t>
    void get_Pface(Tensor<1, <data_t> &diff_value, const Cell<data_t> &current_cell,
		   int direction, int ivar, int dir_switch) const
    {
        const int stride =
            current_cell.get_box_pointers().m_in_stride[direction];
        const int in_index = current_cell.get_in_index();
        diff_value[direction] = get_Pface<data_t>(
	    current_cell.get_box_pointers().m_in_ptr[ivar], in_index, stride, dir_switch);
    }

};

#endif /* WENODERIVATIVES_HPP_ */

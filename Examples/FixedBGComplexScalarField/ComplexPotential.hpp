/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXPOTENTIAL_HPP_
#define COMPLEXPOTENTIAL_HPP_

#include "simd.hpp"

class ComplexPotential
{
protected:
  const double m_mu;
  const double m_lambda;
  const double m_kappa;

public:
    //! The constructor
  ComplexPotential(const double a_mu, const double a_lambda, const double a_kappa) : m_mu(a_mu), m_lambda(a_lambda), m_kappa(a_kappa) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi_re,
                           data_t &dVdphi_im, const vars_t<data_t> &vars) const
    {
      //        const double m = m_params.scalar_mass;
        // The potential value at phi
        // 1/2 m^2 phi^2
        V_of_phi = 0.5 * m_mu * m_mu * vars.phi_Re * vars.phi_Re +
	           0.5 * m_mu * m_mu * vars.phi_Im * vars.phi_Im + 
	           0.25 * m_lambda * pow(vars.phi_Re * vars.phi_Re + vars.phi_Im * vars.phi_Im, 2.) +
	           1./6. * m_kappa * pow(vars.phi_Re * vars.phi_Re + vars.phi_Im * vars.phi_Im, 3.);

        // The potential gradient at phi
        // m^2 phi
	//      dVdphi_re = m_mu * m_mu * vars.phi_Re + m_lambda * vars.phi_Re * (vars.phi_Re * vars.phi_Re + vars.phi_Im * vars.phi_Im);
	dVdphi_re = vars.phi_Re * (m_mu * m_mu + m_lambda * (vars.phi_Re * vars.phi_Re + vars.phi_Im * vars.phi_Im) 
				   + m_kappa * pow(vars.phi_Re * vars.phi_Re + vars.phi_Im * vars.phi_Im,2.));
	//	dVdphi_im = m_mu * m_mu * vars.phi_Im + m_lambda * vars.phi_Im * (vars.phi_Re * vars.phi_Re + vars.phi_Im * vars.phi_Im);
	dVdphi_im = vars.phi_Im * (m_mu * m_mu + m_lambda * (vars.phi_Re * vars.phi_Re + vars.phi_Im * vars.phi_Im)
                                   + m_kappa * pow(vars.phi_Re * vars.phi_Re + vars.phi_Im * vars.phi_Im,2.));
	//pout()<< "Scalar mass in potential" << m_mu <<endl;
    }
};

#endif /* COMPLEXPOTENTIAL_HPP_ */
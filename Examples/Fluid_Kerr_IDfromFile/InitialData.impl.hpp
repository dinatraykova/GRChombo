/* GRChombﬁo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(INITIALDATA_HPP_)
#error "This file should only be included through InitialData.hpp"
#endif

#ifndef INITIALDATA_IMPL_HPP_
#define INITIALDATA_IMPL_HPP_

#include "DimensionDefinitions.hpp"

template <class data_t>
void InitialData::compute(Cell<data_t> current_cell) const
{
    // set up vars for the metric and extrinsic curvature, shift and lapse in
    // spherical coords
    Tensor<2, data_t> spherical_g;
    Tensor<2, data_t> spherical_K;
    Tensor<1, data_t> spherical_shift;

    // The cartesian variables and coords
    auto metric_vars = current_cell.template load_vars<MetricVars>();
    auto matter_vars = current_cell.template load_vars<Vars>();
    Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

    // Compute the components in spherical coords as per 1401.1548
    compute_spherical(spherical_g, spherical_K, spherical_shift, coords);

    // work out where we are on the grid
    data_t x = coords.x;
    double y = coords.y;
    double z = coords.z;
    data_t rr = coords.get_radius();
    data_t rr2 = rr * rr;

    int ind_L = static_cast<int>(floor(rr / m_params.spacing));
    int ind_H = static_cast<int>(ceil(rr / m_params.spacing));
    double rho_L = *(m_params.rho_1D + ind_L);
    double rho_H = *(m_params.rho_1D + ind_H);
    double rho_interp =
        rho_L + (rr / m_params.spacing - ind_L) * (rho_H - rho_L);

    double eps_L = *(m_params.eps_1D + ind_L);
    double eps_H = *(m_params.eps_1D + ind_H);
    double eps_interp =
        eps_L + (rr / m_params.spacing - ind_L) * (eps_H - eps_L);

    double lapse_L = *(m_params.lapse_1D + ind_L);
    double lapse_H = *(m_params.lapse_1D + ind_H);
    double lapse_interp =
        lapse_L + (rr / m_params.spacing - lapse_L) * (lapse_H - lapse_L);

    double phi_L = *(m_params.phi_1D + ind_L);
    double phi_H = *(m_params.phi_1D + ind_H);
    double phi_interp =
        phi_L + (rr / m_params.spacing - ind_L) * (phi_H - phi_L);
    
    using namespace CoordinateTransformations;
    // Convert spherical components to cartesian components using coordinate
    // transforms
    metric_vars.h = spherical_to_cartesian_LL(spherical_g, x, y, z);
    metric_vars.A = spherical_to_cartesian_LL(spherical_K, x, y, z);
    metric_vars.shift = spherical_to_cartesian_U(spherical_shift, x, y, z);

    using namespace TensorAlgebra;
    // Convert to BSSN vars
    data_t deth = compute_determinant(metric_vars.h);
    auto h_UU = compute_inverse_sym(metric_vars.h);
    metric_vars.chi = exp(-4.*phi_interp); //pow(deth, -1. / 3.);

    // transform extrinsic curvature into A and TrK - note h is still non
    // conformal version which is what we need here
    metric_vars.K = compute_trace(metric_vars.A, h_UU);
    make_trace_free(metric_vars.A, metric_vars.h, h_UU);

    // Make conformal
    FOR(i, j)
    {
        metric_vars.h[i][j] *= metric_vars.chi;
        metric_vars.A[i][j] *= metric_vars.chi;
    }

    // use a pre collapsed lapse, could also use analytic one
    metric_vars.lapse = lapse_interp;
    //metric_vars.lapse = pow(metric_vars.chi, 0.5);

    // Populate the variables on the grid
    // NB We stil need to set Gamma^i which is NON ZERO
    // but we do this via a separate class/compute function
    // as we need the gradients of the metric which are not yet available
    current_cell.store_vars(metric_vars);

    // load vars
    // auto metric_vars = current_cell.template load_vars<MetricVars>();
    // auto matter_vars = current_cell.template load_vars<Vars>();
    VarsTools::assign(matter_vars, 0.);

    data_t chi_regularised = simd_max(metric_vars.chi, 1e-6);

    // calculate the field value
    matter_vars.rho = rho_interp;
      //m_params.rho0 * (exp(-pow(rr / 2. / m_params.awidth, 2.0))) +
                       // m_params.delta;
    matter_vars.eps = eps_interp;
    data_t v2 = 0.;
    FOR(i, j)
    v2 += metric_vars.h[i][j] * matter_vars.vi[i] * matter_vars.vi[j] /
          chi_regularised;

    data_t P_of_rho = 0.;
    EoS::compute_eos(P_of_rho, matter_vars);

    data_t WW = 1. / (1. - v2);
    data_t hh = 1. + matter_vars.eps + P_of_rho / matter_vars.rho;

    data_t rho_conformal = matter_vars.rho / pow(chi_regularised, 1.5);
    data_t P_conformal = P_of_rho / pow(chi_regularised, 1.5);

    matter_vars.D = rho_conformal * sqrt(WW);
    matter_vars.tau = rho_conformal * hh * WW - P_conformal - matter_vars.D;

    FOR(i)
    {
        matter_vars.Sj[i] = 0.;
        FOR(j)
        matter_vars.Sj[i] += rho_conformal * hh * WW * metric_vars.h[i][j] *
                             matter_vars.vi[j] / chi_regularised;
    }

    // store the vars
    current_cell.store_vars(matter_vars);
}

template <class data_t>
void InitialData::compute_spherical(Tensor<2, data_t> &spherical_g,
                               Tensor<2, data_t> &spherical_K,
                               Tensor<1, data_t> &spherical_shift,
                               const Coordinates<data_t> coords) const
{
    // Kerr black hole params - mass M and spin a
    double M = m_params.mass;
    double a = m_params.spin;

    // work out where we are on the grid
    data_t x = coords.x;
    double y = coords.y;
    double z = coords.z;

    // the radius, subject to a floor
    data_t r = coords.get_radius();
    data_t r2 = r * r;

    // the radius in xy plane, subject to a floor
    data_t rho2 = simd_max(x * x + y * y, 1e-12);
    data_t rho = sqrt(rho2);

    // calculate useful position quantities
    data_t cos_theta = z / r;
    data_t sin_theta = rho / r;
    data_t cos_theta2 = cos_theta * cos_theta;
    data_t sin_theta2 = sin_theta * sin_theta;

    // calculate useful metric quantities
    //double r_plus = M + sqrt(M * M - a * a);
    //double r_minus = M - sqrt(M * M - a * a);

    // The Boyer-Lindquist coordinate
    //data_t r_BL = r * pow(1.0 + 0.25 * r_plus / r, 2.0);

    // Other useful quantities per 1001.4077
    //data_t Sigma = r_BL * r_BL + a * a * cos_theta2;
    // data_t Delta = r_BL * r_BL - 2.0 * M * r_BL + a * a;
    // In the paper this is just 'A', but not to be confused with A_ij
    //data_t AA = pow(r_BL * r_BL + a * a, 2.0) - Delta * a * a * sin_theta2;
    // The rr component of the conformal spatial matric
    //data_t gamma_rr = 1;
    //        Sigma * pow(r + 0.25 * r_plus, 2.0) / (r * r2 * (r_BL - r_minus));

    // Metric in semi isotropic Kerr-Schild coordinates, r, theta (t or th), phi
    // (p)
    FOR(i, j) { spherical_g[i][j] = 0.0; }
    spherical_g[0][0] = 1.;                   // gamma_rr
    spherical_g[1][1] = r2;                   // gamma_tt
    spherical_g[2][2] = r2 * sin_theta2;      // gamma_pp

    // Extrinsic curvature
    FOR(i, j) { spherical_K[i][j] = 0.0; }

    // set non zero elements of Krtp - K_rp, K_tp
    /*spherical_K[0][2] = 
        a * M * sin_theta2 / (Sigma * sqrt(AA * Sigma)) *
        (3.0 * pow(r_BL, 4.0) + 2 * a * a * r_BL * r_BL - pow(a, 4.0) -
         a * a * (r_BL * r_BL - a * a) * sin_theta2) *
        (1.0 + 0.25 * r_plus / r) / sqrt(r * r_BL - r * r_minus);
    spherical_K[2][0] = spherical_K[0][2];
    spherical_K[2][1] = -2.0 * pow(a, 3.0) * M * r_BL * cos_theta * sin_theta *
                        sin_theta2 / (Sigma * sqrt(AA * Sigma)) *
                        (r - 0.25 * r_plus) * sqrt(r_BL / r - r_minus / r);
			spherical_K[1][2] = spherical_K[2][1];*/
    // set the analytic lapse
    //kerr_lapse = sqrt(Delta * Sigma / AA);

    // set the shift (only the phi component is non zero)
    spherical_shift[0] = 0.0;
    spherical_shift[1] = 0.0;
    spherical_shift[2] = 0.0; //-2.0 * M * a * r_BL / AA;
}

#endif /* INITIALDATA_IMPL_HPP_ */

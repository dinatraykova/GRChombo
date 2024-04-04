/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(INITIALFLUIDDATA_HPP_)
#error "This file should only be included through PrimitiveRecovery.hpp"
#endif

#ifndef INITIALFLUIDDATA_IMPL_HPP_
#define INITIALFLUIDDATA_IMPL_HPP_

inline InitialFluidData::InitialFluidData(params_t a_params, double a_dx)
    : m_dx(a_dx), m_params(a_params)
{
}

template <class data_t>
void InitialFluidData::compute(Cell<data_t> current_cell) const
{
    // load vars
    // FluidCCZ4RHS<PerfectFluid<>>::Vars<data_t> vars;
    // VarsTools::assign(vars, 0.);
    const auto metric_vars = current_cell.template load_vars<MetricVars>();

    // where am i?
    Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
    data_t rr = coords.get_radius();
    data_t rr2 = rr * rr;

    data_t x = coords.x;
    double y = coords.y;
    double z = coords.z;

    Tensor<1, data_t> vi, Sj;
    data_t chi_regularised = simd_max(metric_vars.chi, 1e-6);

    FOR(i) { vi[i] = 0.; }
    data_t D = 0.;
    data_t tau = 0.;
    data_t rho = 0.;
    data_t eps = 0.;

    // calculate the field value
    rho =
        m_params.rho0 * (exp(-pow(rr / m_params.awidth, 2.0))) + m_params.delta;
    data_t v2 = 0.;
    FOR(i, j) v2 += metric_vars.h[i][j] * vi[i] * vi[j] / chi_regularised;

    data_t P_over_rho = (1. + eps) / 3.;
    // data_t P_over_rho = 0.;
    // EoS::compute_eos(P_over_rho, vars);

    data_t WW = 1. / (1. - v2);
    data_t hh = 1. + eps + P_over_rho;

    D = rho * sqrt(WW);
    tau = rho * (hh * WW - P_over_rho) - D;
    FOR(i)
    {
        Sj[i] = 0.;
        FOR(j)
        Sj[i] += rho * hh * WW * metric_vars.h[i][j] * vi[j] / chi_regularised;
    }

    // store the vars
    // current_cell.store_vars(vars);
    current_cell.store_vars(rho, c_rho);
    current_cell.store_vars(vi, GRInterval<c_vi1, c_vi3>());
    current_cell.store_vars(D, c_D);
    current_cell.store_vars(Sj, GRInterval<c_Sj1, c_Sj3>());
    current_cell.store_vars(tau, c_tau);
}

#endif /* INITIALFLUIDDATA_IMPL_HPP_ */

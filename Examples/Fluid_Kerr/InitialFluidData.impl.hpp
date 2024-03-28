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
    // Set up EoS
    // data_t P_over_rho = 0.;
    // my_eos.compute_eos(P_over_rho, vars);
    // load vars
    FluidCCZ4RHS<PerfectFluid<>>::Vars<data_t> vars;
    VarsTools::assign(vars, 0.);
    // where am i?
    Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
    data_t rr = coords.get_radius();
    data_t rr2 = rr * rr;

    data_t x = coords.x;
    double y = coords.y;
    double z = coords.z;

    // Tensor<1, data_t> vi, Sj;
    data_t chi_regularised = simd_max(vars.chi, 1e-6);

    // vi[0] = 0.;
    // vi[1] = 0.;
    // vi[2] = 0.;

    // calculate the field value
    vars.rho =
        m_params.rho0 * (exp(-pow(rr / m_params.awidth, 2.0))) + m_params.delta;
    data_t v2 = 0.;
    FOR(i, j) v2 += vars.h[i][j] * vars.vi[i] * vars.vi[j] / chi_regularised;
    vars.eps = 0.;
    data_t P_over_rho = vars.rho * (1. + vars.eps) / 3.;
    data_t WW = 1. / (1. - v2);
    data_t hh = 1. + vars.eps + P_over_rho;

    vars.D = vars.rho * sqrt(WW);
    vars.tau = vars.rho * (hh * WW - P_over_rho) - vars.D;
    FOR(i)
    {
        vars.Sj[i] = 0.;
        FOR(j)
        vars.Sj[i] +=
            vars.rho * hh * WW * vars.h[i][j] * vars.vi[j] / chi_regularised;
    }

    // store the vars
    current_cell.store_vars(vars);
}

#endif /* INITIALFLUIDDATA_IMPL_HPP_ */

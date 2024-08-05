/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(INITIALFLUIDDATA_HPP_)
#error "This file should only be included through PrimitiveRecovery.hpp"
#endif

#ifndef INITIALFLUIDDATA_IMPL_HPP_
#define INITIALFLUIDDATA_IMPL_HPP_

template <class data_t>
void InitialFluidData::compute(Cell<data_t> current_cell) const
{
    // load vars
    auto metric_vars = current_cell.template load_vars<MetricVars>();
    auto matter_vars = current_cell.template load_vars<Vars>();
    VarsTools::assign(matter_vars, 0.);

    // where am i?
    Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
    data_t rr = coords.get_radius();
    data_t rr2 = rr * rr;

    data_t x = coords.x;
    double y = coords.y;
    double z = coords.z;

    data_t chi_regularised = simd_max(metric_vars.chi, 1e-6);

    // calculate the field value
    matter_vars.rho =
        m_params.rho0; // * (exp(-pow(rr / 2. / m_params.awidth, 2.0))) +
                       // m_params.delta;
    data_t v2 = 0.;
    FOR(i, j)
    v2 += metric_vars.h[i][j] * matter_vars.vi[i] * matter_vars.vi[j] /
          chi_regularised;

    data_t P_of_rho = 0.;
    EoS::compute_eos(P_of_rho, matter_vars);

    data_t WW = 1. / (1. - v2);
    data_t hh = 1. + matter_vars.eps + P_of_rho / matter_vars.rho;

    data_t rho_conformal = matter_vars.rho / pow(chi_regularised, 1.5);

    matter_vars.D = rho_conformal * sqrt(WW);
    matter_vars.tau = rho_conformal * hh * WW - P_of_rho - matter_vars.D;

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

#endif /* INITIALFLUIDDATA_IMPL_HPP_ */
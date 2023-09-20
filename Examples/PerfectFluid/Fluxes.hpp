/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FLUXES_HPP_
#define FLUXES_HPP_

#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

namespace Fluxes
{
template <class data_t, template <typename> class vars_t>
vars_t<data_t> compute_flux(const vars_t<data_t> &vars, const int idir)
{
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    vars_t<data_t> out;
    out.D = (vars.lapse * vars.vi[idir] - vars.shift[idir]) * vars.D;
    Tensor<1, data_t> Sj_U;
    FOR(i)
    {
        Sj_U[i] = 0.;
        FOR(j) Sj_U[i] += vars.chi * h_UU[i][j] * vars.Sj[j];
    }

    data_t chi_regularised = simd_max(1e-6, vars.chi);
    Tensor<1, data_t> vi_D;
    FOR(i)
    {
        vi_D[i] = 0.;
        FOR(j) { vi_D[i] += vars.h[i][j] * vars.vi[j] / chi_regularised; }
    }

    data_t v2 = 0.;
    FOR(i) v2 += vars.vi[i] * vi_D[i];

    data_t P =
        vars.rho * (1. + vars.eps) / 3.; // for now we assume a conformal fluid
    data_t WW = 1. / (1. - v2);
    data_t hh = 1. + vars.eps + P / vars.rho;

    FOR(j)
    {
        out.Sj[j] = vars.lapse * vars.rho * hh * WW * vars.vi[idir] * vi_D[j] -
                    vars.shift[idir] * vars.Sj[j];
        FOR(k) out.Sj[j] += vars.lapse * P * h_UU[idir][k] * vars.h[j][k];
    }

    out.tau = vars.lapse * (Sj_U[idir] - vars.D * vars.vi[idir]) -
              vars.shift[idir] * vars.tau;

    out.Jt = vars.nn * vars.vi[idir] * sqrt(WW);

    return out;
}
template <class data_t, template <typename> class vars_t>
vars_t<data_t> compute_num_flux(const vars_t<data_t> &vars, const int idir,
                                const double lambda, const int sign)
{
    vars_t<data_t> out = compute_flux(vars, idir);
    out.D += sign * lambda * vars.D;
    FOR(j) out.Sj[j] += sign * lambda * vars.Sj[j];
    out.tau += sign * lambda * vars.tau;
    out.Jt += sign * lambda * vars.Jt;
    return out;
}

} // namespace Fluxes
#endif /* FLUXES_HPP_ */

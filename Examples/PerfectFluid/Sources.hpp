/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SOURCES_HPP_
#define SOURCES_HPP_

#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

namespace Sources
{
// Primitive to conservative variables
template <class data_t, template <typename> class vars_t>
vars_t<data_t> compute_source(const vars_t<data_t> &vars,
                              const vars_t<Tensor<1, data_t>> &d1)
{
    data_t chi_regularised = simd_max(1e-6, vars.chi);
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    vars_t<data_t> out;

    data_t v2 = 0.;
    FOR(i, j) v2 += vars.h[i][j] * vars.vi[j] * vars.vi[i] / chi_regularised;

    data_t P =
        vars.rho * (1. + vars.eps) / 3.; // for now we assume a conformal fluid
    data_t WW = 1. / (1. - v2);
    data_t hh = 1. + vars.eps + P / vars.rho;

    out.D = 0.;
    FOR(j)
    {
        out.Sj[j] = -(vars.tau + vars.D) * d1.lapse[j];
        FOR(i)
        {
            out.Sj[j] += vars.Sj[i] * d1.shift[i][j];
            FOR(k)
            {
                out.Sj[j] += vars.lapse / 2. *
                             (d1.h[i][k][j] -
                              vars.h[i][k] * d1.chi[j] / chi_regularised) *
                             (vars.rho * hh * WW * vars.vi[i] * vars.vi[k] /
                                  chi_regularised +
                              P * h_UU[i][k]);
            }
        }
    }
    out.tau = 0.;
    FOR(i, j)
    {
        out.tau += vars.lapse * (vars.A[i][j] + vars.h[i][j] / 3. * vars.K) *
                       (vars.rho * hh * WW * vars.vi[i] * vars.vi[j] /
                            chi_regularised +
                        P * h_UU[i][j]) -
                   vars.chi * h_UU[i][j] * vars.Sj[i] * d1.lapse[j];
    }

    out.Jt = 0.;

    return out;
}

} // namespace Sources
#endif /* SOURCES_HPP_ */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CONSERVATIVERECOVERY_HPP_
#define CONSERVATIVERECOVERY_HPP_

#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

namespace ConservativeRecovery
{
// Primitive to conservative variables
template <class data_t, template <typename> class vars_t>
void PtoC(vars_t<data_t> &vars)
{
    data_t chi_regularised = simd_max(1e-6, vars.chi);
    Tensor<1, data_t> vi_D;
    FOR(i)
    {
        vi_D[i] = 0.;
        FOR(j) { vi_D[i] += vars.h[i][j] * vars.vi[j] / chi_regularised; }
    }

    data_t v2 = 0.;
    FOR(i) v2 += vars.vi[i] * vi_D[i];

    data_t Poverrho =
        (1. + vars.eps) / 3.; // for now we assume a conformal fluid
    data_t WW = 1. / (1. - v2);
    data_t hh = 1. + vars.eps + Poverrho;

    vars.D = vars.rho * sqrt(WW);
    vars.tau = vars.rho * hh * WW - vars.rho * Poverrho - vars.D;

    // S_j (note lower index) = - n^a T_ai
    FOR(i) { vars.Sj[i] = vars.rho * hh * WW * vi_D[i]; }
}

/*// Conservative to primitive variables
template <class data_t, template <typename> class vars_t>
void CtoP(vars_t<data_t> &vars)
{
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    data_t E = vars.tau + vars.D;
    data_t S2 = 0.;
    FOR(i, j) S2 += vars.chi * h_UU[i][j] * vars.Sj[i] * vars.Sj[j];

    data_t E2 = E * E;
    data_t sqrt1 = sqrt(4. * E2 - 3. * S2);
    data_t sqrt2 = sqrt(2. * E2 - 3. * S2 + E * sqrt1);

    // eps
    vars.eps =
        sqrt2 / (2. * vars.D) - 1.; // Is it dangerous to divide by D?
    // rho
    vars.rho = (sqrt1 - E) / (1. + vars.eps);
    // vi_D
    Tensor<1, data_t> vi_D;
    FOR(i) vi_D[i] = 3. * vars.Sj[i] * (sqrt1 - E) / (sqrt2 * sqrt2);
    // vi
    FOR(i)
    {
        vars.vi[i] = 0.;
        FOR(j) vars.vi[i] += vars.chi * h_UU[i][j] * vi_D[j];
    }
}*/
} // namespace ConservativeRecovery
#endif /* CONSERVATIVERECOVERY_HPP_ */

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
    vars.Jt = vars.nn * sqrt(1. - v2);
}

} // namespace ConservativeRecovery
#endif /* CONSERVATIVERECOVERY_HPP_ */

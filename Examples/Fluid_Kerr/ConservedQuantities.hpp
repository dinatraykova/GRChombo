/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CONSERVEDQUANTITIES_HPP_
#define CONSERVEDQUANTITIES_HPP_

#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

namespace ConservedQuantities
{
// Primitive to conservative variables
template <class data_t, template <typename> class vars_t>
void PtoC(const data_t P_over_rho, vars_t<data_t> &vars)
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

    data_t WW = 1. / (1. - v2);
    data_t hh = 1. + vars.eps + P_over_rho;

    data_t rho_conformal = vars.rho / pow(vars.chi, 1.5);

    vars.D = rho_conformal * sqrt(WW);
    vars.tau = rho_conformal * (hh * WW - P_over_rho) - vars.D;

    // S_j (note lower index) = - n^a T_ai
    FOR(i) { vars.Sj[i] = rho_conformal * hh * WW * vi_D[i]; }
}

} // namespace ConservedQuantities
#endif /* CONSERVEDQUANTITIES_HPP_ */

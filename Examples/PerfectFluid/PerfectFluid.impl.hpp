/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(PERFECTFLUID_HPP_)
#error "This file should only be included through PerfectFluid.hpp"
#endif

#ifndef PERFECTFLUID_IMPL_HPP_
#define PERFECTFLUID_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class eos_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> PerfectFluid<eos_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;
    
    // Set up EoS
    data_t P_of_rho = 0.0;
    data_t dPdrho = 0.0;
    my_eos.compute_eos(P_of_rho, dPdrho, vars);
  
    // Useful quantities
    data_t chi_regularised = simd_max(1e-6, vars.chi);
    data_t v2 = 0.0;
    FOR(i, j)
    {
        v2 +=  vars.h[i][j] * vars.vi[i] * vars.vi[j] / chi_regularised;
    }
    data_t WW = 1. / (1. - v2);
    data_t hh = 1. + vars.eps + P_of_rho / vars.rho;

    Tensor<1, data_t> vi_D;
    FOR(i)
    {
        vi_D[i] = 0;
        FOR(j)
        {
          vi_D[i] += vars.h[i][j] * vars.vi[j] / chi_regularised;
        }
    }
    
    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR(i, j)
    {
      out.Sij[i][j] = vars.rho * hh * WW * vi_D[i] * vi_D[j]
	  + vars.h[i][j] * P_of_rho / chi_regularised;
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // S_i (note lower index) = - n^a T_ai
    FOR(i) { out.Si[i] = vars.rho * hh * WW * vi_D[i]; }

    // rho = n^a n^b T_ab
    out.rho = vars.rho * hh * WW - P_of_rho;

    return out;
}

template <class eos_t>
template <class data_t, template <typename> class vars_t,
	  //          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void PerfectFluid<eos_t>::add_matter_rhs(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &lm,
    const vars_t<Tensor<1, data_t>> &lp,
    const vars_t<Tensor<1, data_t>> &rm,
    const vars_t<Tensor<1, data_t>> &rp) const
//const vars_t<Tensor<1, data_t>> &d1,
//   const vars_t<data_t> &advec) const
{
    using namespace TensorAlgebra;

    //const auto h_UU = compute_inverse_sym(vars.h);
    //    const auto chris = compute_christoffel(d1.h, h_UU);

    data_t P_of_rho = 0.0;
    data_t dPdrho = 0.0;
    my_eos.compute_eos(P_of_rho, dPdrho, vars);

    vars_t<data_t> source = Sources::compute_source(vars, lp);
    // evolution equations for the fluid conservative variables
    rhs.D = source.D;
    FOR(i) rhs.Sj[i] = source.Sj[i];
    rhs.tau = source.tau;

    FOR(idir) {
        vars_t<data_t> vars_right_p;
        vars_right_p.rho = rp.rho[idir];
        vars_right_p.eps = rp.eps[idir];
        FOR(j) vars_right_p.vi[j] = rp.vi[j][idir];
        PrimitiveRecovery::PtoC(vars_right_p);
        vars_t<data_t> flux_right_p = Fluxes::compute_flux(vars_right_p, idir, m_lambda);

        vars_t<data_t> vars_right_m;
        vars_right_m.rho = rm.rho[idir];
        vars_right_m.eps = rm.eps[idir];
        FOR(j) vars_right_m.vi[j] = rm.vi[j][idir];
        PrimitiveRecovery::PtoC(vars_right_m);
        vars_t<data_t> flux_right_m = Fluxes::compute_flux(vars_right_m, idir, m_lambda);

        rhs.D -= 1. / (2. * m_dx) * (flux_right_p.D + flux_right_m.D);
        FOR(j) rhs.Sj[j] -= 1. / (2. * m_dx) * (flux_right_p.Sj[j] + flux_right_m.Sj[j]);
        rhs.tau -= 1. / (2. * m_dx) * (flux_right_p.tau + flux_right_m.tau);
    }

    FOR(idir) {
        vars_t<data_t> vars_left_p;
        vars_left_p.rho = lp.rho[idir];
        vars_left_p.eps = lp.eps[idir];
        FOR(j) vars_left_p.vi[j] = lp.vi[j][idir];
        PrimitiveRecovery::PtoC(vars_left_p);
        vars_t<data_t> flux_left_p = Fluxes::compute_flux(vars_left_p, idir, m_lambda);

        vars_t<data_t> vars_left_m;
        vars_left_m.rho = lm.rho[idir];
        vars_left_m.eps = lm.eps[idir];
        FOR(j) vars_left_m.vi[j] = lm.vi[j][idir];
        PrimitiveRecovery::PtoC(vars_left_m);
        vars_t<data_t> flux_left_m = Fluxes::compute_flux(vars_left_m, idir, m_lambda);

        rhs.D += 1. / (2. * m_dx) * (flux_left_p.D + flux_left_m.D);
        FOR(j) rhs.Sj[j] += 1. / (2. * m_dx) * (flux_left_p.Sj[j] + flux_left_m.Sj[j]); 
        rhs.tau += 1. / (2. * m_dx) * (flux_left_p.tau + flux_left_m.tau);
    }

    rhs.rho = 0.;
    rhs.eps = 0.;
    FOR(i) rhs.vi[i] = 0.;

}

#endif /* SCALARFIELD_IMPL_HPP_ */

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
    data_t P_of_rho0 = 0.0;
    data_t dPdrho0 = 0.0;
    my_eos.compute_eos(P_of_rho0, dPdrho0, vars);
  
    // Useful quantities
    data_t v2 = 0.0;
    FOR(i, j)
    {
        v2 +=  vars.h[i][j] * vars.vi[i] * vars.vi[j] / vars.chi;
    }
    data_t WW = 1./(1. - v2);
    data_t hh = 1. + vars.eps + P_of_rho0/vars.rho0;

    Tensor<1, data_t> vi_D;
    FOR(i)
    {
        vi_D[i] = 0;
        FOR(j)
        {
          vi_D[i] += vars.h[i][j] * vars.vi[j] / vars.chi;
        }
    }
    
    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR(i, j)
    {
        out.Sij[i][j] = vars.rho0 * hh * WW * vi_D[i] * vi_D[j]
	  + vars.h[i][j] * P_of_rho0 / vars.chi;
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // S_i (note lower index) = - n^a T_ai
    FOR(i) { out.Si[i] = vars.rho0 * hh * WW * vi_D[i]; }

    // rho = n^a n^b T_ab
    out.rho = vars.rho0 * hh * WW - P_of_rho0;
}

template <class eos_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void PerfectFluid<eos_t>::add_matter_rhs(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    data_t P_of_rho0 = 0.0;
    data_t dPdrho0 = 0.0;
    my_eos.compute_eos(P_of_rho0, dPdrho0, vars);

    // evolution equations for the fluid conservative variables
    rhs.D = 0.0;
    rhs.Ec = 0.0;
    rhs.rho0 = 0.0;
    rhs.eps = 0.0;

    FOR(i)
    {
      rhs.Sj[i] = 0.0;
      rhs.vi[i] = 0.0;
    }
}

#endif /* SCALARFIELD_IMPL_HPP_ */

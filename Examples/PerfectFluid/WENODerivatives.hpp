/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WENODERIVATIVES_HPP_
#define WENODERIVATIVES_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include <array>

class WENODerivatives
{
  public:
    enum
    {
        LEFT_MINUS,
        LEFT_PLUS,
        RIGHT_MINUS,
        RIGHT_PLUS,
    };

    // public:
    //  WENODerivatives(double a_eW) : m_eW(a_eW){}
    WENODerivatives() {}
    template <class data_t>
    ALWAYS_INLINE data_t get_Pface(const double *in_ptr, const int idx,
                                   const int stride, int dir_switch) const
    {
        data_t beta[3] = {0., 0., 0.};
        const double dd[3] = {3. / 10., 3. / 5., 1. / 10.};
        const double eW = 1.;
        data_t alpha[3] = {0., 0., 0.};
        data_t sum_alpha;
        data_t weights[3] = {0., 0., 0.};
        data_t v[3] = {0., 0., 0.};
        // double u[3] = {0.,0.,0.};
        data_t pim2, pim1, pi0, pip1, pip2;

        auto in = SIMDIFY<data_t>(in_ptr);

        if (dir_switch == 0) // LEFT_MINUS
        // Negative fluxes to compute primitive variables
        // at the left cell boundary i-1/2: p^{L-}_{i-1/2}
        {
            pim2 = in[idx - 3 * stride];
            pim1 = in[idx - 2 * stride];
            pi0 = in[idx - stride];
            pip1 = in[idx];
            pip2 = in[idx + stride];

            // ENO polynomials
            v[0] = (2. * pim2 - 7. * pim1 + 11. * pi0) / 6.;
            v[1] = (-pim1 + 5. * pi0 + 2. * pip1) / 6.;
            v[2] = (2. * pi0 + 5. * pip1 - pip2) / 6.;
        }

        if (dir_switch == 1) // LEFT_PLUS
        {
            // Positive fluxes to compute primitive variables
            // at the left cell boundary i-1/2: p^{L+}_{i-1/2}
            pim2 = in[idx - 2 * stride];
            pim1 = in[idx - stride];
            pi0 = in[idx];
            pip1 = in[idx + stride];
            pip2 = in[idx + 2 * stride];

            // ENO polynomials
            v[0] = (-pim2 + 5. * pim1 + 2. * pi0) / 6.;
            v[1] = (2. * pim1 + 5. * pi0 - pip1) / 6.;
            v[2] = (11. * pi0 - 7. * pip1 + 2. * pip2) / 6.;
        }

        if (dir_switch == 2) // RIGHT_MINUS
        // Negative fluxes to compute primitive variables
        // at the right cell boundary i+1/2: p^{R-}_{i+1/2}
        {
            pim2 = in[idx - 2 * stride];
            pim1 = in[idx - stride];
            pi0 = in[idx];
            pip1 = in[idx + stride];
            pip2 = in[idx + 2 * stride];

            // ENO polynomials
            v[0] = (2. * pim2 - 7. * pim1 + 11. * pi0) / 6.;
            v[1] = (-pim1 + 5. * pi0 + 2. * pip1) / 6.;
            v[2] = (2. * pi0 + 5. * pip1 - pip2) / 6.;
        }

        if (dir_switch == 3) // RIGHT_PLUS
        {
            // Positive fluxes to compute primitive variables
            // at the right cell boundary i+1/2: p^{R+}_{i+1/2}
            pim2 = in[idx - stride];
            pim1 = in[idx];
            pi0 = in[idx + stride];
            pip1 = in[idx + 2 * stride];
            pip2 = in[idx + 3 * stride];

            // ENO polynomials
            v[0] = (-pim2 + 5. * pim1 + 2. * pi0) / 6.;
            v[1] = (2. * pim1 + 5. * pi0 - pip1) / 6.;
            v[2] = (11. * pi0 - 7. * pip1 + 2. * pip2) / 6.;
        }

        // Smoothness indicators
        beta[0] = (13. / 12.) * pow(pim2 - 2. * pim1 + pi0, 2) +
                  (1. / 4.) * pow(pim2 - 4. * pim1 + 3. * pi0, 2);
        beta[1] = (13. / 12.) * pow(pim1 - 2. * pi0 + pip1, 2) +
                  (1. / 4.) * pow(pim1 - pip1, 2);
        beta[2] = (13. / 12.) * pow(pi0 - 2. * pip1 + pip2, 2) +
                  (1. / 4.) * pow(3. * pi0 - 4. * pip1 + pip2, 2);

        // Weights
        sum_alpha = 0.;
        FOR(j)
        {
            alpha[j] = dd[j] / ((eW + beta[j]) * (eW + beta[j]));
            sum_alpha += alpha[j];
        }
        FOR(j) { weights[j] = alpha[j] / sum_alpha; }

        // primitive variables on the cell interface
        data_t pFace = 0.;
        FOR(j) { pFace += weights[j] * v[j]; }
        return pFace;
    }

    template <template <typename> class vars_t, class data_t>
    auto get_Pface(const Cell<data_t> &current_cell, int dir_switch) const
    {
        const auto in_index = current_cell.get_in_index();
        const auto strides = current_cell.get_box_pointers().m_in_stride;
        vars_t<Tensor<1, data_t>> d1;
        d1.enum_mapping(
            [&](const int &ivar, Tensor<1, data_t> &var)
            {
                FOR(idir)
                {
                    var[idir] = get_Pface<data_t>(
                        current_cell.get_box_pointers().m_in_ptr[ivar],
                        in_index, strides[idir], dir_switch);
                }
            });
        return d1;
    }

    template <class data_t>
    void get_Pface(Tensor<1, data_t> &diff_value,
                   const Cell<data_t> &current_cell, int direction, int ivar,
                   int dir_switch) const
    {
        const int stride =
            current_cell.get_box_pointers().m_in_stride[direction];
        const int in_index = current_cell.get_in_index();
        diff_value[direction] =
            get_Pface<data_t>(current_cell.get_box_pointers().m_in_ptr[ivar],
                              in_index, stride, dir_switch);
    }
};

#endif /* WENODERIVATIVES_HPP_ */

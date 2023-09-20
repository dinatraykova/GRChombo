/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALFLUIDDATA_HPP_
#define INITIALFLUIDDATA_HPP_

#include "ADMConformalVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FluidCCZ4RHS.hpp"
#include "PerfectFluid.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which sets the initial fluid matter config
// aka Kevin data
class InitialFluidData
{
    // Now the non grid ADM vars
    template <class data_t>
    using MetricVars = ADMConformalVars::VarsWithGauge<data_t>;

  public:
    //! A structure for the input params for fluid properties and initial
    //! conditions
    struct params_t
    {
        std::array<double, CH_SPACEDIM>
            center; //!< Centre of perturbation in initial SF bubble
        double L;
        double rho0;
        double uflow;
        double amplitude;
        double awidth;
        double asigma;
        std::array<double, 2> ycenter;
    };

    //! The constructor
    InitialFluidData(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        const auto metric_vars = current_cell.template load_vars<MetricVars>();

        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;

        Tensor<1, data_t> vi, Sj; // Flow velocity, covariant momentum
        data_t chi_regularised = simd_max(metric_vars.chi, 1e-6);

        vi[0] = m_params.uflow *
                (tanh((y - m_params.ycenter[0]) / m_params.awidth) -
                 tanh((y - m_params.ycenter[1]) / m_params.awidth) - 1.);
        vi[1] = m_params.amplitude * sin(2. * M_PI * x / m_params.L) *
                (exp(-pow((y - m_params.ycenter[0]) / m_params.asigma, 2)) +
                 exp(-pow((y - m_params.ycenter[1]) / m_params.asigma, 2)));
        vi[2] = 0.;

        // number density
        data_t nn =
            1. + 0.5 * (tanh((y - m_params.ycenter[0]) / m_params.awidth) -
                        tanh((y - m_params.ycenter[1]) / m_params.awidth));

        data_t eps = 0.;

        data_t rho = m_params.rho0;
        data_t v2 = 0.;
        FOR(i, j) v2 += metric_vars.h[i][j] * vi[i] * vi[j] / chi_regularised;
        data_t P = rho * (1. + eps) / 3.;
        data_t W2 = 1. / (1. - v2);
        data_t hh = 1. + eps + P / rho;

        data_t D = rho * sqrt(W2);
        data_t tau = rho * hh * W2 - P - D;
        FOR(i)
        {
            Sj[i] = 0.;
            FOR(j)
            Sj[i] +=
                rho * hh * W2 * metric_vars.h[i][j] * vi[j] / chi_regularised;
        }

        data_t Jt = nn * sqrt(1. + v2);

        // store the vars

        current_cell.store_vars(rho, c_rho);
        current_cell.store_vars(vi, GRInterval<c_vi1, c_vi3>());
        current_cell.store_vars(D, c_D);
        current_cell.store_vars(Sj, GRInterval<c_Sj1, c_Sj3>());
        current_cell.store_vars(tau, c_tau);
        current_cell.store_vars(nn, c_nn);
        current_cell.store_vars(Jt, c_Jt);
    }

  protected:
    double m_dx;
    const params_t m_params;
};

#endif /* INITIALFLUIDDATA_HPP_ */

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
  public:
    //! A structure for the input params for fluid properties and initial
    //! conditions
    struct params_t
    {
        std::array<double, CH_SPACEDIM>
            center; //!< Centre of perturbation in initial SF bubble
        double rho0;
        double awidth;
        double delta;
    };

    //! The constructor
    InitialFluidData(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // load vars
        FluidCCZ4RHS<PerfectFluid<>>::Vars<data_t> vars;
        VarsTools::assign(vars, 0.);
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        data_t rr = coords.get_radius();
        data_t rr2 = rr * rr;

        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;

        // Tensor<1, data_t> vi, Sj;
        data_t chi_regularised = simd_max(vars.chi, 1e-6);

        // vi[0] = 0.;
        // vi[1] = 0.;
        // vi[2] = 0.;

        // calculate the field value
        vars.rho = m_params.rho0 * (exp(-pow(rr / m_params.awidth, 2.0))) +
                   m_params.delta;
        data_t v2 = 0.;
        FOR(i, j)
        v2 += vars.h[i][j] * vars.vi[i] * vars.vi[j] / chi_regularised;
        vars.eps = 0.;
        data_t P = vars.rho * (1. + vars.eps) / 3.;
        data_t WW = 1. / (1. - v2);
        data_t hh = 1. + vars.eps + P / vars.rho;

        vars.D = vars.rho * sqrt(WW);
        vars.tau = vars.rho * hh * WW - P - vars.D;
        FOR(i)
        {
            vars.Sj[i] = 0.;
            FOR(j)
            vars.Sj[i] += vars.rho * hh * WW * vars.h[i][j] * vars.vi[j] /
                          chi_regularised;
        }

        // store the vars
        current_cell.store_vars(vars);
    }

  protected:
    double m_dx;
    const params_t m_params;
};

#endif /* INITIALFLUIDDATA_HPP_ */

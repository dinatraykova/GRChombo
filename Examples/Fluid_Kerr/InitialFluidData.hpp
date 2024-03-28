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
    InitialFluidData(params_t a_params, double a_dx);
    //              : m_dx(a_dx), m_params(a_params)
    //{
    //}

    //! Function to compute the value of all the initial vars on the grid
    //! The constructor
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    double m_dx;
    const params_t m_params;
};

#include "InitialFluidData.impl.hpp"

#endif /* INITIALFLUIDDATA_HPP_ */

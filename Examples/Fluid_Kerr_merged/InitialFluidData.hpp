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
#include <cmath>

//! Class which sets the initial fluid matter config
class InitialFluidData : public EoS
{
    template <class data_t>
    using MetricVars = ADMConformalVars::VarsWithGauge<data_t>;

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
        double spacing;
        double *rho_1D;
    };

    //! The constructor
    InitialFluidData(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }
    //! Function to compute the value of all the initial vars on the grid
    //! The constructor
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        data_t D;
        Tensor<1, data_t> Sj;
        data_t tau;
        //        data_t Jt;
        Tensor<1, data_t> vi;
        data_t rho;
        data_t eps;
        //        data_t nn;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
            // VarsTools
            define_enum_mapping(mapping_function, c_D, D);
            define_enum_mapping(mapping_function, GRInterval<c_Sj1, c_Sj3>(),
                                Sj);
            define_enum_mapping(mapping_function, c_tau, tau);
            //            define_enum_mapping(mapping_function, c_Jt, Jt);
            define_enum_mapping(mapping_function, GRInterval<c_vi1, c_vi3>(),
                                vi);
            define_enum_mapping(mapping_function, c_rho, rho);
            define_enum_mapping(mapping_function, c_eps, eps);
            //           define_enum_mapping(mapping_function, c_nn, nn);
        }
    };

  protected:
    double m_dx;
    const params_t m_params;
};

#include "InitialFluidData.impl.hpp"

#endif /* INITIALFLUIDDATA_HPP_ */

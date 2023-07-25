/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PERFECTFLUID_HPP_
#define PERFECTFLUID_HPP_

#include "CCZ4Geometry.hpp"
#include "DefaultEoS.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "WENODerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

//!  Calculates the matter type specific elements such as the EMTensor and
//   matter evolution
/*!
     This class is an example of a matter_t object which calculates the
     matter type specific elements for the RHS update and the evaluation
     of the constraints. This includes the Energy Momentum Tensor, and
     the matter evolution terms. In this case, a scalar field,
     the matter elements are phi and (minus) its conjugate momentum, Pi.
     It is templated over a potential function potential_t which the
     user must specify in a class, although a default is provided which
     sets dVdphi and V_of_phi to zero.
     It assumes minimal coupling of the field to gravity.
     \sa MatterCCZ4(), ConstraintsMatter()
*/
template <class eos_t = DefaultEoS> class PerfectFluid
{
  protected:
    //! The local copy of the potential
    eos_t my_eos;

  public:
    //!  Constructor of class ScalarField, inputs are the matter parameters.
    PerfectFluid(const eos_t a_eos) : my_eos(a_eos) {}

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        data_t D;
        Tensor<1, data_t> Sj;
        data_t Ec;
        Tensor<1, data_t> vi;
        data_t rho0;
        data_t eps;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
	    using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_D, D);
            define_enum_mapping(mapping_function,
                                GRInterval<c_Sj1, c_Sj3>(), Sj);
	    define_enum_mapping(mapping_function, c_Ec, Ec);
            define_enum_mapping(mapping_function,
                                GRInterval<c_vi1, c_vi3>(), vi);
	    define_enum_mapping(mapping_function, c_rho0, rho0);
	    define_enum_mapping(mapping_function, c_eps, eps);
        }
    };

    //! Structure containing the rhs variables for the matter fields requiring
    //!  2nd derivs
    template <class data_t> struct Diff2Vars
    {
        data_t rho0;

        /// Defines the mapping between members of Vars and Chombo grid
        ///  variables (enum in User_Variables)
      template <typename mapping_function_t>
      void enum_mapping(mapping_function_t mapping_function)
      {
          VarsTools::define_enum_mapping(mapping_function, c_rho0, rho0);
  //        VarsTools::define_enum_mapping(
  //              mapping_function, GRInterval<c_Avec1, c_Avec3>(), Avec);
      }
  };

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices)
        const Tensor<3, data_t> &chris_ULL)
        const; //!< the conformal christoffel symbol

    //! The function which adds in the RHS for the matter field vars
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &rhs,             //!< value of the RHS for all vars
        const vars_t<data_t> &vars,          //!< value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< value of the 2nd derivs
        const vars_t<data_t> &advec)
        const; //!< the value of the advection terms

};

#include "PerfectFluid.impl.hpp"

#endif /* PERFECTFLUID_HPP_ */

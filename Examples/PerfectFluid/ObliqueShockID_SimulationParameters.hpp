/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "EoS.hpp"
#include "InitialFluidData.hpp"
#include "Minkowski.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // Initial scalar field data
        initial_params.center =
            center; // already read in SimulationParametersBase
        initial_params.L = L;
        pp.load("G_Newton", G_Newton,
                0.0); // for now the example neglects backreaction
        pp.load("fluid_emax", initial_params.emax, 3.);
        pp.load("fluid_emin", initial_params.emin, 2.97);
        pp.load("fluid_vx_in", initial_params.vx_in, 0.97);
        pp.load("fluid_vy_in", initial_params.vy_in, 0.97);
        pp.load("fluid_nn_in", initial_params.nn_in, 0.2);
        pp.load("lambda", lambda, 1.); // eigenvalue for numerical flux
        pp.load("eos_w", eos_params.eos_w, 1. / 3.);

        // Initial Kerr data
        pp.load("kerr_mass", kerr_params.mass, 1.0);
        pp.load("kerr_spin", kerr_params.spin, 0.0);
        pp.load("kerr_center", kerr_params.center, center);
    }

    void check_params()
    {
        //        warn_parameter("scalar_mass", potential_params.scalar_mass,
        //                potential_params.scalar_mass <
        //                   0.2 / coarsest_dx / dt_multiplier,
        //               "oscillations of scalar field do not appear to be "
        //               "resolved on coarsest level");
        // warn_parameter("scalar_width", initial_params.width,
        //               initial_params.width < 0.5 * L,
        //               "is greater than half the domain size");
        warn_parameter("kerr_mass", kerr_params.mass, kerr_params.mass >= 0.0,
                       "should be >= 0.0");
        check_parameter("kerr_spin", kerr_params.spin,
                        std::abs(kerr_params.spin) <= kerr_params.mass,
                        "must satisfy |a| <= M = " +
                            std::to_string(kerr_params.mass));
        FOR(idir)
        {
            std::string name = "kerr_center[" + std::to_string(idir) + "]";
            warn_parameter(
                name, kerr_params.center[idir],
                (kerr_params.center[idir] >= 0) &&
                    (kerr_params.center[idir] <= (ivN[idir] + 1) * coarsest_dx),
                "should be within the computational domain");
        }
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    double lambda;
    InitialFluidData::params_t initial_params;
    EoS::params_t eos_params;
    Minkowski::params_t kerr_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
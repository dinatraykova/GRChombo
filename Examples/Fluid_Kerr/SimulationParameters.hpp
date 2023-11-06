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
#include "KerrBH.hpp"

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
        // Initial fluid data
        initial_params.center =
            center; // already read in SimulationParametersBase
        pp.load("G_Newton", G_Newton,
                0.0); // for now the example neglects backreaction
        pp.load("fluid_rho0", initial_params.rho0, 1.0);
        pp.load("fluid_width", initial_params.awidth, 0.1);
        pp.load("fluid_delta", initial_params.delta, 0.2);
        pp.load("lambda", lambda, 1.); // eigenvalue for numerical flux
        pp.load("eos_w", eos_params.eos_w, 1. / 3.);

        // Initial Kerr data
        pp.load("kerr_mass", kerr_params.mass);
        pp.load("kerr_spin", kerr_params.spin);
        pp.load("kerr_center", kerr_params.center, center);
    }

    void check_params()
    {
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
    KerrBH::params_t kerr_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

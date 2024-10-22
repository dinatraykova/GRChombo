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
#include "InitialData.hpp"
// #include "KerrBH.hpp"

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
	pp.load("min_D", min_D, 1.e-6);

        // Reading data
        pp.load("spacing", initial_params.spacing);
        pp.load("lines", lines);
        double rho_1D[lines];
        double eps_1D[lines];
        double lapse_1D[lines];
        double phi_1D[lines];

        double tmp_data;
        ifstream read_file_rho("data_files/rho_vs_r.txt");
        ifstream read_file_eps("data_files/eps_vs_r.txt");
        ifstream read_file_lapse("data_files/lapse_vs_r.txt");
        ifstream read_file_phi("data_files/phi_vs_r.txt");

        for (int i = 0; i < lines; ++i)
        {
            read_file_rho >> tmp_data;
            rho_1D[i] = tmp_data;

            read_file_eps >> tmp_data;
            eps_1D[i] = tmp_data;

            read_file_lapse >> tmp_data;
            lapse_1D[i] = tmp_data;

            read_file_phi >> tmp_data;
            phi_1D[i] = tmp_data;
        }
        initial_params.rho_1D = rho_1D;
        initial_params.eps_1D = eps_1D;
        initial_params.lapse_1D = lapse_1D;
        initial_params.phi_1D = phi_1D;

        // Initial Kerr data
        pp.load("kerr_mass", initial_params.mass);
        pp.load("kerr_spin", initial_params.spin);
        pp.load("kerr_center", initial_params.center, center);
    }

    void check_params()
    {
        warn_parameter("kerr_mass", initial_params.mass,
                       initial_params.mass >= 0.0, "should be >= 0.0");
        check_parameter("kerr_spin", initial_params.spin,
                        std::abs(initial_params.spin) <= initial_params.mass,
                        "must satisfy |a| <= M = " +
                            std::to_string(initial_params.mass));
        FOR(idir)
        {
            std::string name = "kerr_center[" + std::to_string(idir) + "]";
            warn_parameter(name, initial_params.center[idir],
                           (initial_params.center[idir] >= 0) &&
                               (initial_params.center[idir] <=
                                (ivN[idir] + 1) * coarsest_dx),
                           "should be within the computational domain");
        }
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    double lambda;
    double min_D;
    int lines;
    InitialData::params_t initial_params;
    // KerrBH::params_t kerr_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

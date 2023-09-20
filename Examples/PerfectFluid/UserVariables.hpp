/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4UserVariables.hpp"
#include "DiagnosticVariables.hpp"

// assign an enum to each variable
enum
{
    // todo: Here add conservative variables
    // Note that it is important that the first enum value is set to 1 more than
    // the last CCZ4 var enum
    c_D = NUM_CCZ4_VARS, // matter field added
    c_Sj1,
    c_Sj2,
    c_Sj3,
    c_tau,
    c_Jt,
    c_vi1,
    c_vi2,
    c_vi3,
    c_rho,
    c_eps,
    c_nn,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_names = {"D",   "Sj1", "Sj2", "Sj3", "tau", "Jt",
                           "vi1", "vi2", "vi3", "rho", "eps", "nn"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */

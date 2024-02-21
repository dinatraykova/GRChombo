/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "PerfectFluidLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "PrimitiveRecovery.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"
#include "WENODerivatives.hpp"

// For RHS update
#include "FluidCCZ4RHS.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"

// For tag cells
// #include "FixedGridsTaggingCriterion.hpp"
#include "ChiTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "EoS.hpp"
#include "GammaCalculator.hpp"
#include "InitialFluidData.hpp"
#include "KerrBH.hpp"
#include "PerfectFluid.hpp"
#include "PositiveDensity.hpp"
#include "SetValue.hpp"

// Things to do at each advance step, after the RK4 is calculated
void PerfectFluidLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha(),
                                     PositiveDensity()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void PerfectFluidLevel::initialData()
{
    CH_TIME("PerfectFluidLevel::initialData");
    if (m_verbosity)
        pout() << "PerfectFluidLevel::initialData " << m_level << endl;

    // First set everything to zero then initial conditions for scalar field -
    // here a Kerr BH and a scalar field profile
    //    EoS eos(m_p.eos_params);
    // PerfectFluidEoS perfect_fluid(eos);
    BoxLoops::loop(
        make_compute_pack(SetValue(0.), KerrBH(m_p.kerr_params, m_dx),
                          InitialFluidData(m_p.initial_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

#ifdef CH_USE_HDF5
// Things to do before outputting a checkpoint file
void PerfectFluidLevel::prePlotLevel()
{
    fillAllGhosts();
    EoS eos(m_p.eos_params);
    PerfectFluidEoS perfect_fluid(m_dx, m_p.lambda, eos);
    BoxLoops::loop(MatterConstraints<PerfectFluidEoS>(perfect_fluid, m_dx,
                                                      m_p.G_Newton, c_Ham,
                                                      Interval(c_Mom1, c_Mom3)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}
#endif

// Things to do in RHS update, at each RK4 step
void PerfectFluidLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                        const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveDensity(),
                                     PositiveChiAndAlpha(),
                                     PrimitiveRecovery()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS/*, disable_simd()*/);

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    EoS eos(m_p.eos_params);
    PerfectFluidEoS perfect_fluid(m_dx, m_p.lambda, eos);
    if (m_p.max_spatial_derivative_order == 4)
    {
        FluidCCZ4RHS<PerfectFluidEoS, MovingPunctureGauge,
                     FourthOrderDerivatives>
            // FourthOrderDerivatives>
            my_ccz4_matter(perfect_fluid, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        FluidCCZ4RHS<PerfectFluidEoS, MovingPunctureGauge,
                     SixthOrderDerivatives>
            // FourthOrderDerivatives>
            my_ccz4_matter(perfect_fluid, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void PerfectFluidLevel::specificUpdateODE(GRLevelData &a_soln,
                                          const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveDensity(),
                                     PositiveChiAndAlpha(),
                                     PrimitiveRecovery()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS, disable_simd());
}

void PerfectFluidLevel::preTagCells()
{
    // We only use chi in the tagging criterion so only fill the ghosts for chi
    fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
}

// void PerfectFluidLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
//                                                 const FArrayBox
//                                                 &current_state)
//{
//     BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.L,
//     m_p.center),
//                    current_state, tagging_criterion);
// }

void PerfectFluidLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    BoxLoops::loop(ChiTaggingCriterion(m_dx), current_state, tagging_criterion);
}

#include "PomLevel.hpp"
#include "initialise.hpp"

using namespace amrex;

int      PomLevel::verbose         = 0;
Real     PomLevel::cfl             = 0.9;
int      PomLevel::do_reflux       = 1;

int      PomLevel::NUM_STATE       = 5;  // Number of variables in the state
int      PomLevel::NUM_GROW        = 2;  // number of ghost cells

// Default constructor - builds invalid objetc
PomLevel::PomLevel()
{
    flux_register = 0;
}

// The basic constructor.
PomLevel(amrex::Amr& papa,
        int lev,
        const amrex::Geometry& level_geom,
        const amrex::BoxArray& bl,
        const amrex::DistributionMapping& dm,
        amrex::Real time)
    :
    AmrLevel(papa, lev, level_geom, bl, dm, time)
{
    flux_reg = 0;
    if (level > 0 && do_reflux)
        flux_register = new FluxRegister(grids, dmap, crse_ratio, level, NUM_STATE);
}

// The destructor
PomLevel::~PomLevel()
{
    delete flux_register;
}

// Initialise grid data at problem start-up.
void PomLevel::initData()
{
    // Loop over grids, call function to init with data.
    MultiFab& S_new = get_new_data(Phi_Type);

    if (verbose) {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.validbox();
        amrex::FArrayBox& fab = S_new[mfi];

        setGeometrySodX(box, fab, geom, 1.4);
    }

    if (verbose) {
	amrex::Print() << "Done initializing the level " << level 
                   << " data " << std::endl;
    }
}


// Initialise data on this level from another PomLevel (during regrid).
//
void PomLevel::init(AmrLevel &old)
{
    PomLevel* oldlev = (PomLevel*) &old;

    // Create new grid data by fillpatching from old.
    Real dt_new = parent->dtLevel(level);
    Real cur_time  = oldlev->state[Phi_Type].curTime();
    Real prev_time = oldlev->state[Phi_Type].prevTime();
    Real dt_old = cur_time - prev_time;

    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(Phi_Type);

    FillPatch(old, S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

// Initialise data on this level after regridding if old level did not previously exist
void AmrLevelAdv::init()
{
    Real dt = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

// Advance grids at this level in time.
Real
AmrLevelAdv::advance (Real time,
                      Real dt,
                      int iteration,
                      int ncycle)
{
   for (int k = 0; k < NUM_STATE_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

        MultiFab& S_new = get_new_data(Phi_Type);

    const Real prev_time = state[Phi_Type].prevTime();
    const Real cur_time = state[Phi_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    // Get pointers to Flux registers, or set pointer to zero if not there.
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;
    
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level) {
        fine = &getFluxRegister(level+1);
        fine->setVal(0.0);
    }

    if (do_reflux && level > 0) {
	    current = &getFluxRegister(level);
    }

    MultiFab fluxes[BL_SPACEDIM];

    if (do_reflux) {
        for (int j = 0; j < BL_SPACEDIM; j++) {
            BoxArray ba = S_new.boxArray();
            ba.surroundingNodes(j);
            fluxes[j].define(ba, dmap, NUM_STATE, 0);
        }
    }

    // State with ghost cells
    MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW);
    FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

}
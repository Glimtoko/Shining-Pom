#include <filesystem>
#include <iostream>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>

#include "PomLevel.hpp"
#include "initialise.hpp"
#include "shiningpom.hpp"

#include "hydro/hydro.hpp"

using namespace amrex;

int      PomLevel::NUM_STATE       = 5;  // Number of variables in the state
int      PomLevel::NUM_GROW        = 2;  // number of ghost cells

// Default constructor - builds invalid objetc
PomLevel::PomLevel()
{
    // std::cout << "Construct\n" << std::flush;
    flux_register = 0;
}

// The basic constructor.
PomLevel::PomLevel(amrex::Amr& papa,
        int lev,
        const amrex::Geometry& level_geom,
        const amrex::BoxArray& bl,
        const amrex::DistributionMapping& dm,
        amrex::Real time)
    :
    AmrLevel(papa, lev, level_geom, bl, dm, time)
{
    // Read_Inputs();
    verbose = true;

    flux_register = 0;
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
    // Print() << "Oh\n";
    // Loop over grids, call function to init with data.
    MultiFab& S_new = get_new_data(Phi_Type);

    if (verbose) {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.validbox();
        amrex::FArrayBox& fab = S_new[mfi];
        amrex::Array4<amrex::Real> const& a = fab.array();

        setGeometrySodX(box, a, geom, 1.4);
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
    // Print() << "In (old)\n";
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
void PomLevel::init()
{
    // std::cout << "In\n" << std::flush;


    Real dt = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time, dt_old, dt);
    MultiFab& S_new = get_new_data(Phi_Type);
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

// Advance grids at this level in time.
Real
PomLevel::advance (amrex::Real time,
                   amrex::Real dt,
                   int iteration,
                   int ncycle)
{   
    for (int k = 0; k < NUM_STATE_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }


    MultiFab& S_new = get_new_data(Phi_Type);
    const Real* dx = geom.CellSize();

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
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif

    {
	FArrayBox flux[BL_SPACEDIM];

	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    FArrayBox& statein = Sborder[mfi];
	    FArrayBox& stateout = S_new[mfi];

        // Allocate fabs for fluxes and MUSCL-Hancock construction.
        for (int i = 0; i < BL_SPACEDIM ; i++) {
            const Box& bxtmp = amrex::surroundingNodes(bx,i);
            flux[i].resize(bxtmp, S_new.nComp());
        }

        // MUSCL-Hancock/Godunov update
        amrex::Array4<amrex::Real> const& sOld = statein.array();
        amrex::Array4<amrex::Real> const& sNew = stateout.array();
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Hydro::MUSCLHancock2D(
                sOld, i, j, k, 0,
                eos_gamma, dt, dx[0], dx[1],
                sNew
            );
        });

        if (do_reflux) {
            for (int i = 0; i < BL_SPACEDIM ; i++)
                fluxes[i][mfi].copy<RunOn::Host>(flux[i],mfi.nodaltilebox(i));
        }
	}
    }

    // if (do_reflux) {
    //     if (current) {
    //         for (int i = 0; i < BL_SPACEDIM ; i++)
    //         current->FineAdd(fluxes[i],i,0,0,NUM_STATE,1.);
    //     }
    //     if (fine) {
    //         for (int i = 0; i < BL_SPACEDIM ; i++)
    //         fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
    //     }
    // }

    return dt;
}


void
PomLevel::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    // Get options, set phys_bc - This needs to be moved...
    // Read_Inputs();

    desc_lst.addDescriptor(
        Phi_Type,
        IndexType::TheCellType(),
        StateDescriptor::Point,
        NUM_GROW,
        NUM_STATE,
        &cell_cons_interp
    );

    // This needs to be updated
    // ------------------------------------------------------------
    int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        lo_bc[i] = hi_bc[i] = BCType::reflect_even;   // reflective boundaries
    }

    // BCRec bc(lo_bc, hi_bc);
    // ------------------------------------------------------------

    amrex::Vector<std::string> variables {"Rho", "MomU", "MomV", "E", "dt"};

    for (int j=0; j<NUM_STATE; j++) {
        BCRec bc(lo_bc, hi_bc);
        desc_lst.setComponent(
            Phi_Type, 
            j, 
            variables[j], 
            bc,
            StateDescriptor::BndryFunc(pomBoundary)
        );
    }
}

// Cleanup data descriptors at end of run.
void
PomLevel::variableCleanUp () 
{
    desc_lst.clear();
}

amrex::Real 
PomLevel::estTimeStep (amrex::Real dt_old)
{
    BL_PROFILE("AmrCorePom::EstTimeStep()");

    amrex::Real dt_est = std::numeric_limits<Real>::max();

    const amrex::Real* dx = geom.CellSize();
    // const amrex::Real* prob_lo = geom.ProbLo();
    // const amrex::Real cur_time = t_new;
    amrex::MultiFab& S_new = get_new_data(Phi_Type);

	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
        const amrex::Box& box = mfi.validbox();

        amrex::FArrayBox& fOld = S_new[mfi];
        amrex::Array4<amrex::Real> const& sOld = fOld.array();

        amrex::ParallelFor(box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Hydro::getCellTimestep(
                sOld,
                i, j, k,
                eos_gamma, dx[0], dx[1], cfl, 0.1
                );
        });
	}
    dt_est = S_new.min(QUANT_DT);

    ParallelDescriptor::ReduceRealMin(dt_est);

    if (verbose) {
        amrex::Print() << "PomLevel::estTimeStep at level " << level 
                       << ":  dt_est = " << dt_est << std::endl;
    }

    return dt_est;
}

amrex::Real
PomLevel::initialTimeStep()
{
    return estTimeStep(0.0);
}

void 
PomLevel::computeInitialDt (int finest_level,
                            int sub_cycle,
                            amrex::Vector<int>& n_cycle,
                            const amrex::Vector<amrex::IntVect>& ref_ratio,
                            amrex::Vector<amrex::Real>& dt_level,
                            amrex::Real end_time)
{
    // Grids have been constructed, compute dt for all levels.
    if (level > 0)
        return;

    amrex::Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    // Limit dt by stop_time.
    const amrex::Real eps = 0.001*dt_0;
    amrex::Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}


void 
PomLevel::computeNewDt (int finest_level,
                        int sub_cycle,
                        amrex::Vector<int>& n_cycle,
                        const amrex::Vector<amrex::IntVect>& ref_ratio,
                        amrex::Vector<amrex::Real>& dt_min,
                        amrex::Vector<amrex::Real>& dt_level,
                        amrex::Real end_time,
                        int post_regrid_flag)
{
    if (level > 0)
        return;

    // Get timestep on each level using the estTimeStep function
    for (int i = 0; i <= finest_level; i++)
    {
        PomLevel& pom_level = getLevel(i);
        dt_min[i] = pom_level.estTimeStep(dt_level[i]);
    }

    // dt Limiter - either pre-regrid dt or max change
    if (post_regrid_flag == 1)
    {
        // Limit dt by pre-regrid dt
        for (int i = 0; i <= finest_level; i++)
        {
            dt_min[i] = std::min(dt_min[i], dt_level[i]);
        }
    } else {
        // Limit dt by max_dt_change * old dt
        for (int i = 0; i <= finest_level; i++)
        {
            dt_min[i] = std::min(dt_min[i],max_dt_change*dt_level[i]);
        }
    }

    // Find the minimum over all levels
    amrex::Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    // Limit dt by stop_time.
    const amrex::Real eps = 0.001*dt_0;
    amrex::Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

// Do work after timestep()
void
PomLevel::post_timestep (int iteration)
{
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    int finest_level = parent->finestLevel();

    // if (do_reflux && level < finest_level)
    //     reflux();

    if (level < finest_level)
        avgDown();
}

// Do work after regrid()
void
PomLevel::post_regrid (int lbase, int new_finest) 
{
    // Do nothing
}

// Do work after a restart()
void
PomLevel::post_restart() 
{
    // Do nothing
}

// Do work after init().
void
PomLevel::post_init (Real stop_time)
{
    if (level > 0)
        return;

    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();
}

void
PomLevel::reflux ()
{
    // Do nothing
}

void
PomLevel::avgDown ()
{
    if (level == parent->finestLevel()) return;
    avgDown(Phi_Type);
}

void
PomLevel::avgDown (int state_indx)
{
    if (level == parent->finestLevel()) return;

    PomLevel& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  S_crse   = get_new_data(state_indx);
    
    amrex::average_down(S_fine, S_crse,
                        fine_lev.geom, geom,
                        0, S_fine.nComp(), parent->refRatio(level));
}


void PomLevel::Read_Inputs()
{
    // Core problem definition
    {
        ParmParse pp("pom");
        pp.query("problem", problem);
        pp.query("end_time", stop_time);
        pp.query("max_step", max_step);
        pp.query("cfl", cfl);
        pp.query("do_reflux", do_reflux);
    }

    // AMR setup - note a lot of the AMR input is read automatically
    // by AMReX
    {
        ParmParse pp("amr");
        pp.query("regrid_int", regrid_int);

        pp.queryarr("e_refine", e_refine);
        pp.queryarr("r_refine", r_refine);
        pp.queryarr("u_refine", u_refine);
        
    }

    // Output
    {
        ParmParse pp("output");
        pp.query("plot_int", plot_int);
        pp.query("plot_file", plot_file);

        pp.query("chk_int", chk_int);
        pp.query("chk_file", chk_file);
        pp.query("restart", restart_chkfile);

        std::string subdir {""};
        pp.query("subdir", subdir);
        if (subdir != "") {
            if (ParallelDescriptor::IOProcessor()) {
                std::filesystem::create_directory(subdir);
            }

            plot_file = subdir + "/" + plot_file;
            chk_file = subdir + "/" + chk_file;
        }
    }



    // amrex::Print() << "Problem: " << problem << std::endl;
}

void 
PomLevel::errorEst (amrex::TagBoxArray& tags,
                    int clearval,
                    int tagval,
                    amrex::Real time,
                    int n_error_buf,
                    int ngrow)
{
    // const int clearval = TagBox::CLEAR;
    // const int tagval = TagBox::SET;

    // const Real* dx      = geom.CellSize();
    // const Real* prob_lo = geom.ProbLo();

    MultiFab& state = get_new_data(Phi_Type);

    Vector<int> itags;
    for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
	{
	    const Box& tilebox  = mfi.tilebox();
        TagBox& tagfab  = tags[mfi];
        tagfab.get_itags(itags, tilebox);

        const auto lo = lbound(tilebox);
	    const auto hi = ubound(tilebox);

        amrex::FArrayBox& fab = state[mfi];
        amrex::Array4<amrex::Real> const& a = fab.array();

        tagfab.get_itags(itags, tilebox);

        int idx = 0;
        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (a(i, j, k, QUANT_MOMU) > u_refine[level]) {
                        itags[idx] = 1;
                        idx++;
                    }
                }
            }
        }

        tagfab.tags_and_untags(itags, tilebox);
    }
}

extern "C"
{
void AMREX_ATTRIBUTE_WEAK amrex_probinit (const int* init,
                     const int* name,
                     const int* namelen,
                     const amrex_real* problo,
                     const amrex_real* probhi)
{
}
}

void pomBoundary(
    Box const& bx, FArrayBox& data,
    const int dcomp, const int numcomp,
    Geometry const& geom, const Real time,
    const Vector<BCRec>& bcr, const int bcomp,
    const int scomp)
{
    /*
        AMRLevel doesn't apply boundary conditions, so hard code them here
    */

    const Box& domain = geom.Domain();

    Array4<Real> const& a = data.array();
    Dim3 lo = lbound(a);
    Dim3 hi = ubound(a);
   
    if (lo.x == -2) {
        for (int j = lo.y; j <= hi.y; ++j) {
            a(-2, j, 0, dcomp) = a(1, j, 0, dcomp);
            a(-1, j, 0, dcomp) = a(0, j, 0, dcomp);
        }
    }

    if (hi.x > domain.bigEnd(0)) {
        int idx = domain.bigEnd(0);
        for (int j = lo.y; j <= hi.y; ++j) {
            a(idx+1, j, 0, dcomp) = a(idx, j, 0, dcomp);
            a(idx+2, j, 0, dcomp) = a(idx-1, j, 0, dcomp);
        }
    }

    if (lo.y == -2) {
        for (int i = lo.x; i <= hi.x; ++i) {
            a(i, -2, 0, dcomp) = a(i, 1, 0, dcomp);
            a(i, -1, 0, dcomp) = a(i, 0, 0, dcomp);
        }
    }

    if (hi.y > domain.bigEnd(1)) {
        int idx = domain.bigEnd(1);
        for (int i = lo.x; i <= hi.x; ++i) {
            a(i, idx+1, 0, dcomp) = a(i, idx, 0, dcomp);
            a(i, idx+2, 0, dcomp) = a(i, idx-1, 0, dcomp);
        }
    }


}
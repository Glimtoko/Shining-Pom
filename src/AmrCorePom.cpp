#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>

#include "AmrCorePom.hpp"
#include "initialise.hpp"
#include "shiningpom.hpp"
#include "hydro/hydro.hpp"


using namespace amrex; 

AmrCorePom::AmrCorePom()
{
    int nlevs_max = max_level + 1;

    // Update istep and nsubsteps
    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
	    nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    // Update timestep arrays
    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    // Update data and BC arrays
    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);
    bcs.resize(5);

    // TODO: This should come from problem spec
    auto xLeft = amrex::BCType::reflect_even;
    auto xRight = amrex::BCType::reflect_even;
    auto yDown = amrex::BCType::reflect_even;
    auto yUp = amrex::BCType::reflect_even;

    Read_Inputs();
    // ParmParse pp("pom");
    // pp.query("problem", problem);
    // amrex::Print() << "Problem: " << problem << std::endl;

    if (problem == 3) {
        xRight = amrex::BCType::reflect_even;
        yDown = amrex::BCType::reflect_even;
        yUp = amrex::BCType::reflect_even;
    }

    // X boundaries
    bcs[QUANT_RHO].setLo(0, amrex::BCType::reflect_even);
    bcs[QUANT_RHO].setHi(0, amrex::BCType::reflect_even);
    bcs[QUANT_MOMU].setLo(0, xLeft);
    bcs[QUANT_MOMU].setHi(0, xRight);
    bcs[QUANT_MOMV].setLo(0, xLeft);
    bcs[QUANT_MOMV].setHi(0, xRight);
    bcs[QUANT_E].setLo(0, amrex::BCType::reflect_even);
    bcs[QUANT_E].setHi(0, amrex::BCType::reflect_even);
    bcs[QUANT_DT].setLo(0, amrex::BCType::reflect_even);
    bcs[QUANT_DT].setHi(0, amrex::BCType::reflect_even);

    // Y boundaries
    bcs[QUANT_RHO].setLo(1, amrex::BCType::reflect_even);
    bcs[QUANT_RHO].setHi(1, amrex::BCType::reflect_even);
    bcs[QUANT_MOMU].setLo(1, yDown);
    bcs[QUANT_MOMU].setHi(1, yUp);
    bcs[QUANT_MOMV].setLo(1, yDown);
    bcs[QUANT_MOMV].setHi(1, yUp);
    bcs[QUANT_E].setLo(1, amrex::BCType::reflect_even);
    bcs[QUANT_E].setHi(1, amrex::BCType::reflect_even);
    bcs[QUANT_DT].setLo(1, amrex::BCType::reflect_even);
    bcs[QUANT_DT].setHi(1, amrex::BCType::reflect_even);

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] is never actually used in the reflux operation
    flux_reg.resize(nlevs_max+1);
}


// Destructor - do nothing
AmrCorePom::~AmrCorePom()
{
}

void AmrCorePom::Read_Inputs()
{
    {
        ParmParse pp("pom");
        pp.query("problem", problem);
        pp.query("end_time", stop_time);
        pp.query("max_step", max_step);
        pp.query("cfl", cfl);
        pp.query("do_reflux", do_reflux);
    }

    {
        ParmParse pp("amr");
        pp.query("regrid_int", regrid_int);

        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);

        pp.query("chk_int", chk_int);
        pp.query("chk_file", chk_file);
        pp.query("restart", restart_chkfile);

        pp.queryarr("e_refine", e_refine);
        pp.queryarr("r_refine", r_refine);
        pp.queryarr("u_refine", u_refine);
        
    }



    amrex::Print() << "Problem: " << problem << std::endl;
}


// Advance solution to final time
void AmrCorePom::Evolve()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;
    int last_chkpt_file_step = 0;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step) {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt();

        int lev = 0;
        int iteration = 1;
        timeStep(lev, cur_time, iteration);

        cur_time += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                        << " DT = " << dt[0]  << std::endl;

        // sync up time
        for (lev = 0; lev <= finest_level; ++lev) {
            t_new[lev] = cur_time;
        }

        if (plot_int > 0 && (step+1) % plot_int == 0) {
            last_plot_file_step = step+1;
            WritePlotFile();
        }

        if (chk_int > 0 && (step+1) % chk_int == 0) {
            last_chkpt_file_step = step+1;
            WriteCheckpointFile();
        }

        if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
	    WritePlotFile();
    }

    if (chk_int > 0 && istep[0] > last_chkpt_file_step) {
	    WriteCheckpointFile();
    }
}

// Initialise
void AmrCorePom::InitData()
{
    // For now, we always initialise from scratch
    if (restart_chkfile == "") {
        const Real time = 0.0;
        InitFromScratch(time);
        AverageDown();

        if (chk_int > 0) {
            WriteCheckpointFile();
        }
    } else {
        ReadCheckpointFile();
    }

    if (plot_int > 0) WritePlotFile();
}


// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void AmrCorePom::MakeNewLevelFromCoarse (
    int lev,
    Real time,
    const BoxArray& ba,
    const DistributionMapping& dm
)
{
    const int ncomp = phi_new[lev-1].nComp();
    const int nghost = phi_new[lev-1].nGrow();

    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	    flux_reg[lev].reset(
            new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp)
        );
    }

    FillCoarsePatch(lev, time, phi_new[lev], 0, ncomp);
}


// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void AmrCorePom::RemakeLevel (
    int lev,
    Real time,
    const BoxArray& ba,
    const DistributionMapping& dm
)
{
    const int ncomp = phi_new[lev].nComp();
    const int nghost = phi_new[lev].nGrow();

    MultiFab new_state(ba, dm, ncomp, nghost);
    MultiFab old_state(ba, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	    flux_reg[lev].reset(
            new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp)
        );
    }
}


// Delete level data
// overrides the pure virtual function in AmrCore
void AmrCorePom::ClearLevel (int lev)
{
    phi_new[lev].clear();
    phi_old[lev].clear();
    flux_reg[lev].reset(nullptr);
}


// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void AmrCorePom::MakeNewLevelFromScratch(
    int lev,
    Real time,
    const BoxArray& ba,
    const DistributionMapping& dm
)
{
    const int ncomp = 5;
    const int nghost = 2;

    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	    flux_reg[lev].reset(
            new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp)
        );
    }

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real cur_time = t_new[lev];

    MultiFab& state = phi_new[lev];

    // TODO: This needs to account for problem number
    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.validbox();
        amrex::FArrayBox& fab = state[mfi];
        amrex::Array4<amrex::Real> const& a = fab.array();

        switch(problem) {
            case 1:
                setGeometrySodX(box, a, geom[lev], 1.4);
                break;
            case 2:
                setGeometrySodY(box, a, geom[lev], 1.4);
                break;
            case 3:
                setGeometryTriple(box, a, geom[lev], 1.4);
                break;
        }
    }
}


// Tag cells for refinement
// overrides the pure virtual function in AmrCore
// TODO: Re-write this completely
void AmrCorePom::ErrorEst(
    int lev,
    TagBoxArray& tags,
    Real time,
    int ngrow
)
{
    if (lev >= e_refine.size()) return;

    const int clearval = TagBox::CLEAR;
    const int tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    MultiFab& state = phi_new[lev];

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
                    if (a(i, j, k, QUANT_MOMU) > u_refine[lev]) {
                        itags[idx] = 1;
                        idx++;
                    }
                }
            }
        }

        tagfab.tags_and_untags(itags, tilebox);
    }
}


// Set covered coarse cells to be the average of overlying fine cells
void
AmrCorePom::AverageDown()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
	amrex::average_down(
        phi_new[lev+1], phi_new[lev],
        geom[lev+1], geom[lev],
        0, phi_new[lev].nComp(), refRatio(lev)
    );
    }
}


// More flexible version of AverageDown() that lets you average down across multiple levels
void
AmrCorePom::AverageDownTo(int crse_lev)
{
    amrex::average_down(
        phi_new[crse_lev+1], phi_new[crse_lev],
        geom[crse_lev+1], geom[crse_lev],
        0, phi_new[crse_lev].nComp(), refRatio(crse_lev))
    ;
}


// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void AmrCorePom::FillPatch(
    int lev,
    Real time,
    MultiFab& mf,
    int icomp,
    int ncomp
)
{
    if (lev == 0) {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime);

        CpuBndryFuncFab bfunc;
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        
        amrex::FillPatchSingleLevel(
            mf,
            time,
            smf,
            stime,
            0,
            icomp,
            ncomp,
            geom[lev],
            physbc,
            0
        );
    } else {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev-1, time, cmf, ctime);
        GetData(lev, time, fmf, ftime);

        CpuBndryFuncFab bfunc;
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(
            mf,
            time,
            cmf,
            ctime,
            fmf,
            ftime,
            0, icomp, ncomp, geom[lev-1], geom[lev],
            cphysbc, 0, fphysbc, 0,
            refRatio(lev-1), mapper, bcs, 0
        );
    }
}


// Fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void AmrCorePom::FillCoarsePatch (
    int lev,
    Real time,
    MultiFab& mf,
    int icomp,
    int ncomp
)
{
    BL_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime);

    if (cmf.size() != 1) {
	    amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    CpuBndryFuncFab bfunc;
    PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
    PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

    Interpolater* mapper = &cell_cons_interp;

    amrex::InterpFromCoarseLevel(
        mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
        cphysbc, 0, fphysbc, 0, refRatio(lev-1),
        mapper, bcs, 0
    );
}


// Utility to copy in data from phi_old and/or phi_new into another multifab
void AmrCorePom::GetData (
    int lev,
    Real time,
    Vector<MultiFab*>& data,
    Vector<Real>& datatime
)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)     {
        data.push_back(&phi_new[lev]);
        datatime.push_back(t_new[lev]);
    } else if (time > t_old[lev] - teps && time < t_old[lev] + teps)  {
        data.push_back(&phi_old[lev]);
        datatime.push_back(t_old[lev]);
    } else {
        data.push_back(&phi_old[lev]);
        data.push_back(&phi_new[lev]);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

// advance a level by dt
// includes a recursive call for finer levels
void AmrCorePom::timeStep (int lev, Real time, int iteration)
{
    // We may need to regrid
    if (regrid_int > 0)  { 
        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static Vector<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if
        // it was taken care of during a coarser regrid
        if (lev < max_level && istep[lev] > last_regrid_step[lev]) {
            if (istep[lev] % regrid_int == 0) {
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
		        int old_finest = finest_level;
		        regrid(lev, time);

                // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) {
                    last_regrid_step[k] = istep[k];
                }

                // if there are newly created levels, set the time step
		        for (int k = old_finest+1; k <= finest_level; ++k) {
		            dt[k] = dt[k-1] / MaxRefRatio(k-1);
		        }
	        }
	    }
    }

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
        amrex::Print() << "ADVANCE with time = " << t_new[lev]
                       << " dt = " << dt[lev] << std::endl;
    }

    // Advance a single level for a single time step, updates flux registers
    Advance(lev, time, dt[lev], iteration, nsubsteps[lev]);

    // Increment step
    ++istep[lev];

    if (lev < finest_level) {
        // recursive call for next-finer level
	    for (int i = 1; i <= nsubsteps[lev+1]; ++i) {
            timeStep(lev+1, time+(i-1)*dt[lev+1], i);
        }

        if (do_reflux) {
            // update lev based on coarse-fine flux mismatch
            flux_reg[lev+1]->Reflux(phi_new[lev], 1.0, 0, 0, phi_new[lev].nComp(), geom[lev]);
        }

	    AverageDownTo(lev); // average lev+1 down to lev
    }
}


// advance a single level for a single time step, updates flux registers
void AmrCorePom::Advance(
    int lev,
    Real time,
    Real dt_lev,
    int iteration,
    int ncycle
)
{
    constexpr int num_grow = 2;

    std::swap(phi_old[lev], phi_new[lev]);
    t_old[lev] = t_new[lev];
    t_new[lev] += dt_lev;

    MultiFab& S_new = phi_new[lev];

    const Real old_time = t_old[lev];
    const Real new_time = t_new[lev];
    const Real ctr_time = 0.5*(old_time+new_time);

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    // Storage for fluxes
    MultiFab fluxes[BL_SPACEDIM];
    if (do_reflux) {
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            BoxArray ba = grids[lev];
            ba.surroundingNodes(i);
            fluxes[i].define(ba, dmap[lev], S_new.nComp(), 0);
        }
    }

    // State with ghost cells
    MultiFab Sborder(grids[lev], dmap[lev], S_new.nComp(), num_grow);
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());
    {
	    FArrayBox flux[BL_SPACEDIM];

        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            FArrayBox& statein = Sborder[mfi];
            FArrayBox& stateout      =   S_new[mfi];

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
                    gamma, dt_lev, dx[0], dx[1],
                    sNew
                );
            });

            // TODO: Fluxes for reflux
            // if (do_reflux) {
            //     for (int i = 0; i < BL_SPACEDIM ; i++) {
            //         fluxes[i][mfi].copy<RunOn::Host>(flux[i],mfi.nodaltilebox(i));
            //     }
            // }
        }
    }
}


void AmrCorePom::ComputeDt()
{
    Vector<Real> dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev) {
	    dt_tmp[lev] = EstTimeStep(lev, true);
    }
    ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

    constexpr Real change_max = 1.1;
    Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = std::min(dt_tmp[lev], change_max*dt[lev]);
        n_factor *= nsubsteps[lev];
        dt_0 = std::min(dt_0, n_factor*dt_tmp[lev]);
    }

    // Limit dts by the value of stop_time.
    const Real eps = 1.e-3*dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) {
	    dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
	    dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
}


// Compute dt from CFL considerations
Real AmrCorePom::EstTimeStep(int lev, bool local)
{
    BL_PROFILE("AmrCorePom::EstTimeStep()");

    Real dt_est = std::numeric_limits<Real>::max();

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    const Real cur_time = t_new[lev];

    MultiFab& S_new = phi_new[lev];

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
                gamma, dx[0], dx[1], cfl, 0.1
                );
        });
	}
    dt_est = S_new.min(QUANT_DT);

    if (!local) {
	    ParallelDescriptor::ReduceRealMin(dt_est);
    }

    return dt_est;
}


// Get plotfile name
std::string AmrCorePom::PlotFileName(int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}


// Put together an array of multifabs for writing
Vector<const MultiFab*> AmrCorePom::PlotFileMF() const
{
    Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
	    r.push_back(&phi_new[i]);
    }
    return r;
}


// set plotfile variable names
Vector<std::string> AmrCorePom::PlotFileVarNames() const
{
    return {"Rho", "MomU", "MomV", "E", "dt"};
}

// write plotfile to disk
void AmrCorePom::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();

    amrex::Print() 
        << "Writing plotfile " 
        << plotfilename 
        << " at (course) step " << istep[0]
        << ", time = " << t_new[0] << "s"
        << std::endl;

    amrex::WriteMultiLevelPlotfileHDF5(
        plotfilename,
        finest_level + 1,
        mf,
        varnames,
		Geom(),
        t_new[0],
        istep,
        refRatio()
    );
}

// This routine stolen from the AMR advection example...
void AmrCorePom::WriteCheckpointFile () const
{

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save
    //                     (e.g., finest_level, t_new, etc.) and also the BoxArrays 
    //                     at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each 
    //                     level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(chk_file,istep[0]);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       std::ofstream HeaderFile;
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
       HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
		                               std::ofstream::trunc |
                                               std::ofstream::binary);
       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       HeaderFile.precision(17);

       // write out title line
       HeaderFile << "Checkpoint file for AmrCoreAdv\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write out array of istep
       for (int i = 0; i < istep.size(); ++i) {
           HeaderFile << istep[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of dt
       for (int i = 0; i < dt.size(); ++i) {
           HeaderFile << dt[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of t_new
       for (int i = 0; i < t_new.size(); ++i) {
           HeaderFile << t_new[i] << " ";
       }
       HeaderFile << "\n";

       // write the BoxArray at each level
       for (int lev = 0; lev <= finest_level; ++lev) {
           boxArray(lev).writeOn(HeaderFile);
           HeaderFile << '\n';
       }
   }

   // write the MultiFab data to, e.g., chk00010/Level_0/
   for (int lev = 0; lev <= finest_level; ++lev) {
       VisMF::Write(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "phi"));
   }

}


namespace {
// utility to skip to next line in Header
void GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}
}


// This routine stolen from the AMR advection example...
void AmrCorePom::ReadCheckpointFile ()
{

    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab and FluxRegister data
        int ncomp = 1;
        int nghost = 0;
        phi_old[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        phi_new[lev].define(grids[lev], dmap[lev], ncomp, nghost);

        if (lev > 0 && do_reflux) {
            flux_reg[lev].reset(new FluxRegister(grids[lev], dmap[lev], refRatio(lev-1), lev, ncomp));
        }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "phi"));
    }

}
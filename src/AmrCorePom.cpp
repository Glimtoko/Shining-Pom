#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>

#include "AmrCorePom.hpp"
#include "shiningpom.hpp"


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

    // X boundaries
    bcs[QUANT_RHO].setLo(0, amrex::BCType::reflect_even);
    bcs[QUANT_RHO].setHi(0, amrex::BCType::reflect_even);
    bcs[QUANT_MOMU].setLo(0, xLeft);
    bcs[QUANT_MOMU].setHi(0, xRight);
    bcs[QUANT_MOMV].setLo(0, xLeft);
    bcs[QUANT_MOMV].setHi(0, xRight);
    bcs[QUANT_E].setLo(0, amrex::BCType::reflect_even);
    bcs[QUANT_E].setHi(0, amrex::BCType::reflect_even);

    // Y boundaries
    bcs[QUANT_RHO].setLo(1, amrex::BCType::reflect_even);
    bcs[QUANT_RHO].setHi(1, amrex::BCType::reflect_even);
    bcs[QUANT_MOMU].setLo(1, yDown);
    bcs[QUANT_MOMU].setHi(1, yUp);
    bcs[QUANT_MOMV].setLo(1, yDown);
    bcs[QUANT_MOMV].setHi(1, yUp);
    bcs[QUANT_E].setLo(1, amrex::BCType::reflect_even);
    bcs[QUANT_E].setHi(1, amrex::BCType::reflect_even);

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


// Initialise
AmrCorePom::InitData()
{
    // For now, we always initialise from scratch
    InitFromScratch();

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
    const int ncomp = 1;
    const int nghost = 0;

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

        setGeometrySodX(box, fab, geom, 1.4);
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
    static Vector<Real> e_refine {2.5, 2.7, 3.0, 3.5};

    if (lev >= e_refine.size()) return;

    const int clearval = TagBox::CLEAR;
    const int tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    const MultiFab& state = phi_new[lev];

    Vector<int>  itags;
    for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
	{
	    const Box& tilebox  = mfi.tilebox();
        TagBox& tagfab  = tags[mfi];
        tagfab.get_itags(itags, tilebox);

        const int* lo = tilebox.loVect();
	    const int* hi = tilebox.hiVect();

        amrex::FArrayBox& fab = state[mfi];

        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                }
                if (fab(i, j, k, QUANT_E) > e_refine[lev]) {
                    tagfab(i, j, k) = 1;
                }
            }
        }
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
void AmrCoreAdv::FillPatch(
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

        BndryFuncArray bfunc(phifill);
        PhysBCFunct<BndryFuncArray> physbc(geom[lev], bcs, bfunc);
        
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

        BndryFuncArray bfunc(phifill);
        PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ],bcs,bfunc);

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

    BndryFuncArray bfunc(phifill);
    PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1],bcs,bfunc);
    PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ],bcs,bfunc);

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
	    FArrayBox& flux[BL_SPACEDIM];
        FArrayBox& mhLR[BL_SPACEDIM];
        FArrayBox& mhUD[BL_SPACEDIM];

        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            const FArrayBox& statein = Sborder[mfi];
            FArrayBox& stateout      =   S_new[mfi];

            // Allocate fabs for fluxes and MUSCL-Hancock construction.
            for (int i = 0; i < BL_SPACEDIM ; i++) {
                const Box& bxtmp = amrex::surroundingNodes(bx,i);
                flux[i].resize(bxtmp, S_new.nComp());
                mhLR[i].resize(bxtmp, S_new.nComp());
                mhUD[i].resize(bxtmp, S_new.nComp());
            }

            // MUSCL-Hancock Reconstruction
            amrex::Array4<amrex::Real> const& sOld = statein.array();
            amrex::Array4<amrex::Real> const& lMH = mhLR[0].array();
            amrex::Array4<amrex::Real> const& rMH = mhLR[1].array();
            amrex::Array4<amrex::Real> const& uMH = mhUD[0].array();
            amrex::Array4<amrex::Real> const& dMH = mhUD[1].array();

            amrex::ParallelFor(box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Hydro::MUSCLHancock2D_Reconstruct(
                    sOld,
                    i, j, k,
                    gamma, dt, dx[0], dx[1],
                    lMH, rMH, dMH, uMH
                );
            });

            // Calculate fluxes
            amrex::Array4<amrex::Real> const& fX = flux[0].array();
            amrex::Array4<amrex::Real> const& fY = flux[0].array();

            amrex::ParallelFor(box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Hydro::calculateFluxes(
                    lMH, rMH, dMH, uMH,
                    i, j, k,
                    gamma,
                    fX, fY
                );
            });

            // compute velocities on faces (prescribed function of space and time)
            // get_face_velocity(&lev, &ctr_time,
            //         AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
            //             BL_TO_FORTRAN(uface[1]),
            //             BL_TO_FORTRAN(uface[2])),
            //         dx, prob_lo);

                // compute new state (stateout) and fluxes.
            //     advect(&time, bx.loVect(), bx.hiVect(),
            // BL_TO_FORTRAN_3D(statein),
            // BL_TO_FORTRAN_3D(stateout),
            // AMREX_D_DECL(BL_TO_FORTRAN_3D(uface[0]),
            //     BL_TO_FORTRAN_3D(uface[1]),
            //     BL_TO_FORTRAN_3D(uface[2])),
            // AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]),
            //     BL_TO_FORTRAN_3D(flux[1]),
            //     BL_TO_FORTRAN_3D(flux[2])),
            // dx, &dt_lev);

            if (do_reflux) {
                for (int i = 0; i < BL_SPACEDIM ; i++) {
                    fluxes[i][mfi].copy<RunOn::Host>(flux[i],mfi.nodaltilebox(i));
                }
            }
        }
    }
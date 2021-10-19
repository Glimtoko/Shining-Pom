#ifndef _AmrCorePom_H_
#define _AMRCorePom_H_

#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_BCRec.H>

#include "rapidjson/document.h"

using namespace amrex;

class AmrCorePom: public amrex::AmrCore
{
public:
    // Constructor and destructor
    AmrCorePom(
        Geometry const& level_0_gome, 
        Vector<BCRec> bcs,
        AmrInfo const& amr_info
    );
    virtual ~AmrCorePom();

    // Advance solution to final time
    void Evolve();

    // User Input
    void Read_Inputs();

    // Initialise data
    void InitData();

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm
    ) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm
    ) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel(int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm
    ) override;

    void SetGeometry(    
        amrex::Box const& bx,
        amrex::Array4<amrex::Real> const& a,
        amrex::Geometry const& geom,
        double gamma
    );

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst(
        int lev,
        amrex::TagBoxArray& tags,
        amrex::Real time,
        int ngrow
    ) override;

private:
    // Set covered coarse cells to be the average of overlying fine cells
    void AverageDown();

    // More flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo(int crse_lev);

    // Compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating
    // from coarse)
    void FillPatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mf,
        int icomp,
        int ncomp
    );

    // Fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
     void FillCoarsePatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mf,
        int icomp,
        int ncomp
    );

    // Utility to copy in data from phi_old and/or phi_new into another multifab
    void GetData (
        int lev,
        amrex::Real time,
        amrex::Vector<amrex::MultiFab*>& data,
        amrex::Vector<amrex::Real>& datatime
    );

    // Advance a level by dt
    // Includes a recursive call for finer levels
    void timeStep (
        int lev,
        amrex::Real time,
        int iteration
    );

    // Advance a single level for a single time step, updates flux registers
    void Advance (
        int lev,
        amrex::Real time,
        amrex::Real dt_lev,
        int iteration,
        int ncycle
    );

    // Wrapper for EstTimeStep()
    void ComputeDt ();

    // Timestep calculation
    Real EstTimeStep (int lev, bool local);

    // Compute dt from CFL considerations
    Real EstTimeStep (
        int lev,
        bool local=false
    ) const;

    // Output
    // ======

    // Get plotfile name
    std::string PlotFileName(int lev) const;

    // Put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF() const;

    // Write plotfile to disk
    void WritePlotFile() const;

    // Write checkpoint file
    void WriteCheckpointFile () const;

    // Restart from checkpoint file
    void ReadCheckpointFile ();


    // Private data members
    // ====================

    amrex::Vector<std::string> variables {"Rho", "MomU", "MomV", "E", "dt"};
    int nghost {2};

    amrex::Vector<int> istep;      // Which step?
    amrex::Vector<int> nsubsteps;  // How many substeps on each level?

    // Keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // Array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap".
    amrex::Vector<amrex::MultiFab> phi_new;
    amrex::Vector<amrex::MultiFab> phi_old;

    // BCRec is essentially a 2*DIM integer array storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<BCRec> bcs;  // 1-component

    // Stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[nlevs_max] are never actually
    // used in the reflux operation
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg;


    // Runtime parameters
    // ==================

    // Problem ID
    int problem = 3;

    // Maximum number of steps and stop time
    int max_step = 10000; //std::numeric_limits<int>::max();
    amrex::Real stop_time = 40; //std::numeric_limits<amrex::Real>::max();

    // CFL number - dt = CFL*dx/umax
    amrex::Real cfl = 0.7;

    // Gamma
    amrex::Real gamma = 1.4;

    // Refinement criteria
    amrex::Vector<amrex::Real> e_refine {1.0, 2.2, 2.3, 3.5};
    amrex::Vector<amrex::Real> u_refine {0.1, 0.15, 0.2, 3.5};
    amrex::Vector<amrex::Real> r_refine {0.8, 1.0};

    // How often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // Hyperbolic refluxing as part of multilevel synchronization
    int do_reflux = 0;

    // Plotfile prefix and frequency
    std::string plot_file {"output/plt"};
    int plot_int = -1; // Default -> Disable plots
    amrex::Real plot_dt = -1.0; // Default -> Disable plots

    // Checkpoint prefex and frequency
    std::string chk_file {"output/chkpt"};
    int chk_int = 25;

    std::string restart_chkfile {""};

    amrex::Real omega = 0.0;
};

#endif

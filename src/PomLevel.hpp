#ifndef _PomLevel_H_
#define _PomLevel_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_FluxRegister.H>

enum StateType {
    Phi_Type = 0,
    NUM_STATE_TYPE
};

using namespace amrex;

class PomLevel: public amrex::AmrLevel
{
public:
    // Default constructor.  Builds invalid object.
    PomLevel();

    // The basic constructor.
    PomLevel(amrex::Amr& papa,
	         int lev,
             const amrex::Geometry& level_geom,
             const amrex::BoxArray& bl,
             const amrex::DistributionMapping& dm,
             amrex::Real time);

    // The destructor.
    virtual ~PomLevel() override;

    // Read input parameters
    void Read_Inputs();

    // Data at problem start-up
    virtual void initData () override;

    // Initialize data on this level from another AmrLevelAdv (during regrid).
    virtual void init (amrex::AmrLevel& old) override;

    // Initialize data on this level after regridding if old level did not previously exist
    virtual void init () override;

    // Define data descriptors.
    static void variableSetUp ();

    // Cleanup data descriptors at end of run.
    static void variableCleanUp ();

    // Advance grids at this level in time.
    virtual amrex::Real advance (amrex::Real time,
                                 amrex::Real dt,
                                 int iteration,
                                 int ncycle) override;

    // Estimate timestep on a level
    amrex::Real estTimeStep (amrex::Real dt_old);

    // Compute initial time step.
    amrex::Real initialTimeStep ();

    // Compute initial `dt'.
    virtual void computeInitialDt (int finest_level,
                                   int sub_cycle,
                                   amrex::Vector<int>& n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>& dt_level,
                                   amrex::Real end_time) override;

    // Compute new `dt'.
    virtual void computeNewDt (int finest_level,
                               int sub_cycle,
                               amrex::Vector<int>& n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>& dt_min,
                               amrex::Vector<amrex::Real>& dt_level,
                               amrex::Real end_time,
                               int post_regrid_flag) override;

    // Do work after timestep().
    virtual void post_timestep (int iteration) override;

    // Do work after regrid().
    virtual void post_regrid (int lbase, int new_finest) override;

    // Do work after init().
    virtual void post_init (amrex::Real stop_time) override;

    // Do work after restart()
    virtual void post_restart () override;

    // Error estimation for regridding.
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int clearval,
                           int tagval,
                           amrex::Real time,
			               int n_error_buf = 0,
                           int ngrow = 0) override;

    static int  NUM_STATE;
    static int  NUM_GROW;

protected:

    // Inlined functions
    PomLevel& getLevel(int level);
    amrex::FluxRegister& getFluxRegister();
    amrex::FluxRegister& getFluxRegister(int level);

    // Other utility functions
    void reflux ();
    void avgDown ();
    void avgDown (int state_indx);

    // The data.
    amrex::FluxRegister* flux_register;

    // Static data members.

// private:

    amrex::Vector<std::string> variables {"Rho", "MomU", "MomV", "E", "dt"};

    amrex::Real max_dt_change = 1.1;

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
    amrex::Real eos_gamma = 1.4;

    // Refinement criteria
    amrex::Vector<amrex::Real> e_refine {1.0, 2.2, 2.3, 3.5};
    amrex::Vector<amrex::Real> u_refine {0.1, 0.15, 0.2, 3.5};
    amrex::Vector<amrex::Real> r_refine {0.8, 1.0};

    // How often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // Hyperbolic refluxing as part of multilevel synchronization
    int do_reflux = 0;

    int verbose = 0;

    // Plotfile prefix and frequency
    std::string plot_file {"output/plt"};
    int plot_int = 10;

    // Checkpoint prefex and frequency
    std::string chk_file {"output/chkpt"};
    int chk_int = 25;

    std::string restart_chkfile {""};

};

inline PomLevel& PomLevel::getLevel(int level)
{
    return *(PomLevel *) &parent->getLevel(level);
}

inline amrex::FluxRegister& PomLevel::getFluxRegister()
{
    BL_ASSERT(flux_register);
    return *flux_register;
}

inline amrex::FluxRegister& PomLevel::getFluxRegister(int level)
{
    return getLevel(level).getFluxRegister();
}

extern "C" {
void boundnull(Real* data, AMREX_ARLIM_P(lo), AMREX_ARLIM_P(hi),
              const int* dom_lo, const int* dom_hi,
              const Real* dx, const Real* grd_lo,
              const Real* time, const int* bc);
}
#endif
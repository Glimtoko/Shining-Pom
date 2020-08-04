#ifndef _PomLevel_H_
#define _PomLevel_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>

enum StateType {
    Phi_Type = 0,
    NUM_STATE_TYPE
};

class PomLevel: public amrex::AmrLevel
{
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

    // Data at problem start-up
    virtual void initData () override;

    // Initialize data on this level from another AmrLevelAdv (during regrid).
    virtual void init (amrex::AmrLevel& old) override;

    // Initialize data on this level after regridding if old level did not previously exist
    virtual void init () override;

    // Advance grids at this level in time.
    virtual amrex::Real advance (amrex::Real time,
                                 amrex::Real dt,
                                 int iteration,
                                 int ncycle) override;

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

    // The data.
    amrex::FluxRegister* flux_register;

    // Static data members.
    static int verbose;
    static amrex::Real cfl;
    static int do_reflux;
}

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
    return getLevel(level).getFluxReg();
}

#endif
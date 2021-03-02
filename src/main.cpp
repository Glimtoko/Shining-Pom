
#include <iostream>

#include <AMReX.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>

#include "PomLevel.hpp"

#include <fenv.h>

using namespace amrex;

int main(int argc, char* argv[])
{
    // feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);
    amrex::Initialize(argc,argv);

    Print() << "Starting\n";

    // wallclock time
    auto dRunTime1 = amrex::second();

    int  max_step;
    Real strt_time = 0.0;
    Real stop_time;

    ParmParse pp("pom");
    pp.query("end_time", stop_time);
    pp.query("max_step", max_step);

	Amr amr;

    Print() << "Init\n";
	amr.init(strt_time, stop_time);
    Print() << "Init done\n";

    amr.writePlotFile();

	while ( amr.okToContinue() &&
  	       (amr.levelSteps(0) < max_step || max_step < 0) &&
	       (amr.cumTime() < stop_time || stop_time < 0.0) )

	{
	    //
	    // Do a coarse timestep.  Recursively calls timeStep()
	    //
        Print() << "Timestep\n";
	    amr.coarseTimeStep(stop_time);
	}

	// Write final checkpoint and plotfile
	if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
	    amr.checkPoint();
	}

	if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
	    amr.writePlotFile();
	}

    auto dRunTime2 = amrex::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run time = " << dRunTime2 << std::endl;

    amrex::Finalize();
}
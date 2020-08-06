
#include <iostream>

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include "AmrCorePom.hpp"

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    // wallclock time
    const Real strt_total = amrex::second();

    {
        // constructor - reads in parameters from inputs file
        //             - sizes multilevel arrays and data structures
        Print() << "Starting" << std::endl;
        AmrCorePom amr_core_adv;
	
            // initialize AMR data
        amr_core_adv.InitData();

            // advance solution to final time
        amr_core_adv.Evolve();
        
            // wallclock time
        Real end_total = amrex::second() - strt_total;
        
            // print wallclock time
        ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
        if (amr_core_adv.Verbose()) {
                amrex::Print() << "\nTotal Time: " << end_total << '\n';
        }
    }

    amrex::Finalize();
}
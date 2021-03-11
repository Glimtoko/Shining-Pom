
#include <iostream>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

#include "AmrCorePom.hpp"
#include "utilities.hpp"

#include <fenv.h>

using namespace amrex;

int main(int argc, char* argv[])
{
    // feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW);
    // feenableexcept(FE_ALL_EXCEPT);
    amrex::Initialize(argc,argv);

    // wallclock time
    const Real strt_total = amrex::second();

    {
        // Read and store AMR input
        AmrInfo amr_info = pom::GetAmrInfo();

        // Read mesh data
        std::string meshJSONFile {""};
        std::string meshJSONEmbed {""};
        char *meshJSON;
        int JSONlen;
        ParmParse pp("mesh");
        pp.query("file", meshJSONFile);
        pp.query("embed", meshJSONEmbed);
        
        if (meshJSONEmbed != "") {
            Print() << "Using embedded mesh description" << std::endl;

            // Force lower case
            std::transform(meshJSONEmbed.begin(), meshJSONEmbed.end(), meshJSONEmbed.begin(),
                                [](unsigned char c){ return std::tolower(c); });

            // RapidJSON requires double quoted strings, but AMReX input format
            // forced single quotes, so convert these
            std::replace(meshJSONEmbed.begin(), meshJSONEmbed.end(), '\'', '"');

            JSONlen = meshJSONEmbed.size() + 1;
            const char *meshJSONConst = meshJSONEmbed.c_str();

            meshJSON = new char[JSONlen];
            strcpy(meshJSON, meshJSONConst);
        } else {
            Print() << "Using mesh description from file " << meshJSONFile << std::endl;
            if (ParallelDescriptor::IOProcessor()) {

                std::ifstream t(meshJSONFile);
                std::stringstream buffer;
                buffer << t.rdbuf();
                std::string meshJSONStr = buffer.str();

                std::transform(meshJSONStr.begin(), meshJSONStr.end(), meshJSONStr.begin(),
                    [](unsigned char c){ return std::tolower(c); });

                JSONlen = meshJSONStr.size() + 1;
                const char *meshJSONConst = meshJSONStr.c_str();

                meshJSON = new char[JSONlen];

                if (ParallelDescriptor::IOProcessor()) {
                    strcpy(meshJSON, meshJSONConst);
                }
            }

            ParallelDescriptor::Bcast<int>(&JSONlen, 1);
            if (!ParallelDescriptor::IOProcessor()) {
                meshJSON = new char[JSONlen];
            }
            ParallelDescriptor::Bcast<char>(meshJSON, JSONlen);
        }


        amrex::Vector<std::string> variables {"Rho", "MomU", "MomV", "E", "dt"};
        pom::PomMesh mesh = pom::GetPomMesh(meshJSON, variables);


        Print() << "Starting" << std::endl;
        AmrCorePom amr_core_adv(mesh.geom, mesh.bcs, amr_info);
	
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
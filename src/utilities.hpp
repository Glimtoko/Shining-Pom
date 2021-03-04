#ifndef _PomUtilities_H_
#define _PomUtilities_H_

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>


namespace pom {
struct PomMesh {
    amrex::Geometry geom;
    amrex::Vector<amrex::BCRec> bcs;
};

PomMesh GetPomMesh(char *meshJSON, amrex::Vector<std::string> variables);
amrex::AmrInfo GetAmrInfo ();
}

#endif
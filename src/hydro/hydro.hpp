#ifndef FLUXH
#define FLUXH

#include "AMReX_FabArray.H"

// #include "mesh2d.hpp"

namespace Hydro {
    struct Flux {
        double rho;
        double momU;
        double momV;
        double E;
    };

    Flux getFluxHLLC(
        double uL, double vL, double rhoL, double pL,
        double uR, double vR, double rhoR, double pR,
        double gamma
    );

    void MUSCLHancock1D(
        double* rho, double* E, double* momN, double* momT,
        int ni, int iUpper, double gamma, double dt, double dx
    );

    void MUSCLHancock2D(
        amrex::Array4<amrex::Real> const& stateOld,
        int iIndex, int jIndex, int kIndex, int niGhosts,
        double gamma, double dt, double dx, double dy,
        amrex::Array4<amrex::Real> const& stateNew
    );

    double getTimestep(
        double *rho, double *momU, double *momV, double *E,
        double dx, double dy,
        int nCells,
        double gamma, double cfl, double dtmax
    );

    void getCellTimestep(
    amrex::Array4<amrex::Real> const& stateOld,
    int i, int j, int k,
    double gamma, double dx, double dy, double cfl, double dtmax,
    amrex::Array4<amrex::Real> const& stateNew
    );
}

#endif

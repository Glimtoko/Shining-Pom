#ifndef FLUXH
#define FLUXH

#define RHO 0
#define MOMU 1
#define MOMV 2
#define ENERGY 3
#define DT 4
#define GET(A, D, I, J) A(I, J, 0, D)

#include "AMReX_FabArray.H"

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

    void SetLimiter(int limiter);

    void MUSCLHancock2D(
        amrex::Array4<amrex::Real> const& stateOld,
        int iIndex, int jIndex, int kIndex, int niGhosts,
        double gamma, double dt, double dx, double dy,
        amrex::Array4<amrex::Real> const& stateNew
    );

    void update(
        amrex::Array4<amrex::Real> const& stateOld,
        amrex::Array4<amrex::Real> const& fluxX,
        amrex::Array4<amrex::Real> const& fluxY,
        int i, int j, int k,
        double dt, double dx, double dy
    );

    void getCellTimestep(
    amrex::Array4<amrex::Real> const& stateOld,
    int i, int j, int k,
    double gamma, double dx, double dy, double cfl, double dtmax
    );
}

#endif

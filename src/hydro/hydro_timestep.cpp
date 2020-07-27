#include "hydro.hpp"
#include <math.h>
#include <iostream>

#define RHO 0
#define MOMU 1
#define MOMV 2
#define ENERGY 3
#define DT 4

#define GETOLD(D, I, J) stateOld(I, J, 0, D)
#define GETNEW(D, I, J) stateNew(I, J, 0, D)

double Hydro::getTimestep(
    double *rho, double *momU, double *momV, double *E,
    double dx, double dy,
    int nCells,
    double gamma, double cfl, double dtmax
) {
    double Sx = 0.0, Sy = 0.0;
    for (int i=0; i<nCells; i++) {
        double u = momU[i]/rho[i];
        double v = momV[i]/rho[i];
        double p = (gamma - 1.0)*(E[i] - 0.5*rho[i]*u*u - 0.5*rho[i]*u*u);
        double a = sqrt((gamma*p)/rho[i]);

        Sx = std::max(Sx, a + abs(u));
        Sy = std::max(Sy, a + abs(v));
    }

    return std::min(dtmax, std::min(cfl*dx/Sx, cfl*dy/Sy));
}


void Hydro::getCellTimestep(
    amrex::Array4<amrex::Real> const& stateOld,
    int i, int j, int k,
    double gamma, double dx, double dy, double cfl, double dtmax,
    amrex::Array4<amrex::Real> const& stateNew
) {
    double rho = GETOLD(RHO, i, j);
    double momU = GETOLD(MOMU, i, j);
    double momV = GETOLD(MOMV, i, j);
    double E = GETOLD(ENERGY, i, j);

    double u = momU/rho;
    double v = momV/rho;

    double p = (gamma - 1.0)*(E - 0.5*rho*u*u - 0.5*rho*u*u);
    double a = sqrt((gamma*p)/rho);

    double Sx = a + abs(u);
    double Sy = a + abs(v);

    GETNEW(DT, i, j) = std::min(dtmax, std::min(cfl*dx/Sx, cfl*dy/Sy));

//     std::cout << i << " " << j << " " << GETNEW(DT, i, j) <<
}

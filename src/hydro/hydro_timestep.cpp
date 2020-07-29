#include "hydro.hpp"
#include <math.h>
#include <iostream>

#define RHO 0
#define MOMU 1
#define MOMV 2
#define ENERGY 3
#define DT 4

#define GETOLD(D, I, J) stateOld(I, J, 0, D)

void Hydro::getCellTimestep(
    amrex::Array4<amrex::Real> const& stateOld,
    int i, int j, int k,
    double gamma, double dx, double dy, double cfl, double dtmax
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

    GETOLD(DT, i, j) = std::min(dtmax, std::min(cfl*dx/Sx, cfl*dy/Sy));
}

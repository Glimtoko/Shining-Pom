#include "hydro.hpp"

#include <math.h>
#include <iostream>
#include "AMReX_FabArray.H"

double getXSlope(amrex::Array4<amrex::Real> stateOld, int U, int i, int j, double omega);
double getYSlope(amrex::Array4<amrex::Real> stateOld, int U, int i, int j, double omega);
double getLimiter(double di1, double di2, double omega);

inline double calcFluxRho(double rho, double u) {
    return rho*u;
}


inline double calcFluxMom(double rho, double u, double v, double p) {
    return rho*u*v + p;
}


inline double calcFluxE(double u, double E, double p) {
    return u*(E + p);
}

void Hydro::MUSCLHancock2D_Reconstruct(
    amrex::Array4<amrex::Real> const& stateOld,
    int i, int j, int k,
    double gamma, double dt, double dx, double dy,
    amrex::Array4<amrex::Real> const& left,
    amrex::Array4<amrex::Real> const& right,
    amrex::Array4<amrex::Real> const& down,
    amrex::Array4<amrex::Real> const& up
)
{
    double omega = 0.0;

    // Initial slope construction
    double di = 0.5*getXSlope(stateOld, RHO, i, j, omega);
    double rhoL = GET(stateOld, RHO, i, j) - di;
    double rhoR = GET(stateOld, RHO, i, j) + di;

    di = 0.5*getXSlope(stateOld, MOMU, i, j, omega);
    double momUL = GET(stateOld, MOMU, i, j) - di;
    double momUR = GET(stateOld, MOMU, i, j) + di;

    di = 0.5*getXSlope(stateOld, MOMV, i, j, omega);
    double momVL = GET(stateOld, MOMV, i, j) - di;
    double momVR = GET(stateOld, MOMV, i, j) + di;

    di = 0.5*getXSlope(stateOld, ENERGY, i, j, omega);
    double EL = GET(stateOld, ENERGY, i, j) - di;
    double ER = GET(stateOld, ENERGY, i, j) + di;

    di = 0.5*getYSlope(stateOld, RHO, i, j, omega);
    double rhoD = GET(stateOld, RHO, i, j) - di;
    double rhoU = GET(stateOld, RHO, i, j) + di;

    di = 0.5*getYSlope(stateOld, MOMU, i, j, omega);
    double momUD = GET(stateOld, MOMU, i, j) - di;
    double momUU = GET(stateOld, MOMU, i, j) + di;

    di = 0.5*getYSlope(stateOld, MOMV, i, j, omega);
    double momVD = GET(stateOld, MOMV, i, j) - di;
    double momVU = GET(stateOld, MOMV, i, j) + di;

    di = 0.5*getYSlope(stateOld, ENERGY, i, j, omega);
    double ED = GET(stateOld, ENERGY, i, j) - di;
    double EU = GET(stateOld, ENERGY, i, j) + di;

    // Set half timestep factors
    double fx = 0.5*dt/dx;
    double fy = 0.5*dt/dy;

    // Advance to half timestep
    double uL = momUL/rhoL;
    double uR = momUR/rhoR;

    double vL = momVL/rhoL;
    double vR = momVR/rhoR;

    double pL = (gamma - 1.0)*(
        EL - 0.5*rhoL*uL*uL
           - 0.5*rhoL*vL*vL
    );

    double pR = (gamma - 1.0)*(
        ER - 0.5*rhoR*uR*uR
           - 0.5*rhoR*vR*vR
    );


    double dFx_rho = fx*(calcFluxRho(rhoL, uL) -
                         calcFluxRho(rhoR, uR));

    double dFx_momN = fx*(calcFluxMom(rhoL, uL, uL, pL) -
                          calcFluxMom(rhoR, uR, uR, pR));

    double dFx_momT = fx*(calcFluxMom(rhoL, uL, vL, 0.0) -
                          calcFluxMom(rhoR, uR, vR, 0.0));

    double dFx_E = fx*(calcFluxE(uL, EL, pL) -
                       calcFluxE(uR, ER, pR));


    double uD = momUD/rhoD;
    double uU = momUU/rhoU;

    double vD = momVD/rhoD;
    double vU = momVU/rhoU;

    double pD = (gamma - 1.0)*(
        ED - 0.5*rhoD*uD*uD
           - 0.5*rhoD*vD*vD
    );

    double pU = (gamma - 1.0)*(
        EU - 0.5*rhoU*uU*uU
           - 0.5*rhoU*vU*vU
    );


    double dFy_rho = fy*(calcFluxRho(rhoD, vD) -
                         calcFluxRho(rhoU, vU));

    double dFy_momN = fy*(calcFluxMom(rhoD, vD, vD, pD) -
                          calcFluxMom(rhoU, vU, vU, pU));

    double dFy_momT = fy*(calcFluxMom(rhoD, uD, vD, 0.0) -
                          calcFluxMom(rhoU, uU, vU, 0.0));

    double dFy_E = fy*(calcFluxE(vD, ED, pD) -
                       calcFluxE(vU, EU, pU));


    rhoL += dFx_rho + dFy_rho;
    rhoR += dFx_rho + dFy_rho;
    momUL += dFx_momN + dFy_momT;
    momUR += dFx_momN + dFy_momT;
    momVL += dFx_momT + dFy_momN;
    momVR += dFx_momT + dFy_momN;
    EL += dFx_E + dFy_E;
    ER += dFx_E + dFy_E;

    rhoD += dFx_rho + dFy_rho;
    rhoU += dFx_rho + dFy_rho;
    momUD += dFx_momN + dFy_momT;
    momUU += dFx_momN + dFy_momT;
    momVD += dFx_momT + dFy_momN;
    momVU += dFx_momT + dFy_momN;
    ED += dFx_E + dFy_E;
    EU += dFx_E + dFy_E;

    // Store results
    GET(left, RHO, i, j) = rhoL;
    GET(right, RHO, i, j) = rhoR;
    GET(up, RHO, i, j) = rhoU;
    GET(down, RHO, i, j) = rhoD;

    GET(left, MOMU, i, j) = momUL;
    GET(right, MOMU, i, j) = momUR;
    GET(up, MOMU, i, j) = momUU;
    GET(down, MOMU, i, j) = momUD;

    GET(left, MOMV, i, j) = momVL;
    GET(right, MOMV, i, j) = momVR;
    GET(up, MOMV, i, j) = momVU;
    GET(down, MOMV, i, j) = momVD;

    GET(left, ENERGY, i, j) = EL;
    GET(right, ENERGY, i, j) = ER;
    GET(up, ENERGY, i, j) = EU;
    GET(down, ENERGY, i, j) = ED;
}

double getXSlope(amrex::Array4<amrex::Real> stateOld, int U, int i, int j, double omega) {
    double di1 = GET(stateOld, U, i, j) - GET(stateOld, U, i-1, j);
    double di2 = GET(stateOld, U, i+1, j) - GET(stateOld, U, i, j);

    double diU = 0.5*(1.0 + omega)*di1 + 0.5*(1.0 - omega)*di2;

    diU *= getLimiter(di1, di2, omega);

    return diU;
}

double getYSlope(amrex::Array4<amrex::Real> stateOld, int U, int i, int j, double omega) {
    double di1 = GET(stateOld, U, i, j) - GET(stateOld, U, i, j-1);
    double di2 = GET(stateOld, U, i, j+1) - GET(stateOld, U, i, j);

    double diU = 0.5*(1.0 + omega)*di1 + 0.5*(1.0 - omega)*di2;

    diU *= getLimiter(di1, di2, omega);

    return diU;
}

double getLimiter(double di1, double di2, double omega) {
    double xi;
    // Slope limiter - Van Leer
    if (di2 == 0) {
        xi = 0.0;
    } else {
        double r = di1/di2;
        if (r <= 0.0) {
            xi = 0.0;
        } else {
            double xiR = 2.0/(1.0 - omega + (1 + omega)*r);
            xi = std::min(2*r/(1+r), xiR);
        }
    }
    return xi;
}

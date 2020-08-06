#include "hydro.hpp"

#include <math.h>
#include <iostream>
#include "AMReX_FabArray.H"

double getLimiter(double di1, double di2, double omega);
double getSlopeX(amrex::Array4<amrex::Real> stateOld, int U, int i, int j, int niGhosts, double omega);
double getSlopeY(amrex::Array4<amrex::Real> stateOld, int U, int i, int j, int niGhosts, double omega);

double calcFluxRho(double rho, double u);
double calcFluxMom(double rho, double u, double v, double p);
double calcFluxE(double u, double E, double p);

#define RHO 0
#define MOMU 1
#define MOMV 2
#define ENERGY 3
#define DT 4

#define GETOLD(D, I, J) stateOld(I, J, 0, D)
#define GETNEW(D, I, J) stateNew(I, J, 0, D)

void Hydro::MUSCLHancock2D(
    amrex::Array4<amrex::Real> const& stateOld,
    int iIndex, int jIndex, int kIndex, int niGhosts,
    double gamma, double dt, double dx, double dy,
    amrex::Array4<amrex::Real> const& stateNew
) {
    double omega = 0.0;

//     std::cout << iIndex << " " << jIndex << std::endl;

    // 5 Point stencil
    int stencil[5][2] =
    {
        {-1,0}, {0,0}, {1,0},
        {0,-1}, {0,1}
    };


    // Data reconstruction in X
    double rhoL[5];
    double rhoR[5];
    double momUL[5];
    double momUR[5];
    double momVL[5];
    double momVR[5];
    double EL[5];
    double ER[5];
    double rhoD[5];
    double rhoU[5];
    double momUD[5];
    double momUU[5];
    double momVD[5];
    double momVU[5];
    double ED[5];
    double EU[5];

    for (int n=0; n<5; n++) {
        int i = iIndex + stencil[n][0];
        int j = jIndex + stencil[n][1];

        double di = 0.5*getSlopeX(stateOld, RHO, i, j, niGhosts, omega);
        rhoL[n] = GETOLD(RHO, i, j) - di;
        rhoR[n] = GETOLD(RHO, i, j) + di;

        di = 0.5*getSlopeX(stateOld, MOMU, i, j, niGhosts, omega);
        momUL[n] = GETOLD(MOMU, i, j) - di;
        momUR[n] = GETOLD(MOMU, i, j) + di;

        di = 0.5*getSlopeX(stateOld, MOMV, i, j, niGhosts, omega);
        momVL[n] = GETOLD(MOMV, i, j) - di;
        momVR[n] = GETOLD(MOMV, i, j) + di;

        di = 0.5*getSlopeX(stateOld, ENERGY, i, j, niGhosts, omega);
        EL[n] = GETOLD(ENERGY, i, j) - di;
        ER[n] = GETOLD(ENERGY, i, j) + di;

        di = 0.5*getSlopeY(stateOld, RHO, i, j, niGhosts, omega);
        rhoD[n] = GETOLD(RHO, i, j) - di;
        rhoU[n] = GETOLD(RHO, i, j) + di;

        di = 0.5*getSlopeY(stateOld, MOMU, i, j, niGhosts, omega);
        momUD[n] = GETOLD(MOMU, i, j) - di;
        momUU[n] = GETOLD(MOMU, i, j) + di;

        di = 0.5*getSlopeY(stateOld, MOMV, i, j, niGhosts, omega);
        momVD[n] = GETOLD(MOMV, i, j) - di;
        momVU[n] = GETOLD(MOMV, i, j) + di;

        di = 0.5*getSlopeY(stateOld, ENERGY, i, j, niGhosts, omega);
        ED[n] = GETOLD(ENERGY, i, j) - di;
        EU[n] = GETOLD(ENERGY, i, j) + di;
    }

    double fx = 0.5*dt/dx;
    double fy = 0.5*dt/dy;
    for (int i=0; i<5; i++) {
        double uL = momUL[i]/rhoL[i];
        double uR = momUR[i]/rhoR[i];

        double vL = momVL[i]/rhoL[i];
        double vR = momVR[i]/rhoR[i];

        double pL = (gamma - 1.0)*(
            EL[i] - 0.5*rhoL[i]*uL*uL
                  - 0.5*rhoL[i]*vL*vL
        );

        double pR = (gamma - 1.0)*(
            ER[i] - 0.5*rhoR[i]*uR*uR
                  - 0.5*rhoR[i]*vR*vR
        );


        double dFx_rho = fx*(calcFluxRho(rhoL[i], uL) -
                             calcFluxRho(rhoR[i], uR));

        double dFx_momN = fx*(calcFluxMom(rhoL[i], uL, uL, pL) -
                              calcFluxMom(rhoR[i], uR, uR, pR));

        double dFx_momT = fx*(calcFluxMom(rhoL[i], uL, vL, 0.0) -
                              calcFluxMom(rhoR[i], uR, vR, 0.0));

        double dFx_E = fx*(calcFluxE(uL, EL[i], pL) -
                           calcFluxE(uR, ER[i], pR));


        double uD = momUD[i]/rhoD[i];
        double uU = momUU[i]/rhoU[i];

        double vD = momVD[i]/rhoD[i];
        double vU = momVU[i]/rhoU[i];

        double pD = (gamma - 1.0)*(
            ED[i] - 0.5*rhoD[i]*uD*uD
                  - 0.5*rhoD[i]*vD*vD
        );

        double pU = (gamma - 1.0)*(
            EU[i] - 0.5*rhoU[i]*uU*uU
                  - 0.5*rhoU[i]*vU*vU
        );


        double dFy_rho = fy*(calcFluxRho(rhoD[i], vD) -
                             calcFluxRho(rhoU[i], vU));

        double dFy_momN = fy*(calcFluxMom(rhoD[i], vD, vD, pD) -
                              calcFluxMom(rhoU[i], vU, vU, pU));

        double dFy_momT = fy*(calcFluxMom(rhoD[i], uD, vD, 0.0) -
                              calcFluxMom(rhoU[i], uU, vU, 0.0));

        double dFy_E = fy*(calcFluxE(vD, ED[i], pD) -
                           calcFluxE(vU, EU[i], pU));


        rhoL[i] += dFx_rho + dFy_rho;
        rhoR[i] += dFx_rho + dFy_rho;
        momUL[i] += dFx_momN + dFy_momT;
        momUR[i] += dFx_momN + dFy_momT;
        momVL[i] += dFx_momT + dFy_momN;
        momVR[i] += dFx_momT + dFy_momN;
        EL[i] += dFx_E + dFy_E;
        ER[i] += dFx_E + dFy_E;

        rhoD[i] += dFx_rho + dFy_rho;
        rhoU[i] += dFx_rho + dFy_rho;
        momUD[i] += dFx_momN + dFy_momT;
        momUU[i] += dFx_momN + dFy_momT;
        momVD[i] += dFx_momT + dFy_momN;
        momVU[i] += dFx_momT + dFy_momN;
        ED[i] += dFx_E + dFy_E;
        EU[i] += dFx_E + dFy_E;

    }


    Flux fluxL, fluxR, fluxU, fluxD;
    double rhoX1 = rhoR[0];
    double uX1 = momUR[0]/rhoR[0];
    double vX1 = momVR[0]/rhoR[0];
    double pX1 = (gamma - 1.0)*(
        ER[0] - 0.5*rhoR[0]*uX1*uX1
              - 0.5*rhoR[0]*vX1*vX1
    );

    double rhoX2 = rhoL[1];
    double uX2 = momUL[1]/rhoL[1];
    double vX2 = momVL[1]/rhoL[1];
    double pX2 = (gamma - 1.0)*(
        EL[1] - 0.5*rhoL[1]*uX2*uX2
              - 0.5*rhoL[1]*vX2*vX2
    );

    double rhoX3 = rhoR[1];
    double uX3 = momUR[1]/rhoR[1];
    double vX3 = momVR[1]/rhoR[1];
    double pX3 = (gamma - 1.0)*(
        ER[1] - 0.5*rhoR[1]*uX3*uX3
              - 0.5*rhoR[1]*vX3*vX3
    );

    double rhoX4 = rhoL[2];
    double uX4 = momUL[2]/rhoL[2];
    double vX4 = momVL[2]/rhoL[2];
    double pX4 = (gamma - 1.0)*(
        EL[2] - 0.5*rhoL[2]*uX4*uX4
              - 0.5*rhoL[2]*vX4*vX4
    );

    fluxL = Hydro::getFluxHLLC(
        uX1, vX1, rhoX1, pX1,
        uX2, vX2, rhoX2, pX2,
        gamma);

    fluxR = Hydro::getFluxHLLC(
        uX3, vX3, rhoX3, pX3,
        uX4, vX4, rhoX4, pX4,
        gamma);

    // Y
    double rhoY1 = rhoU[3];
    double uY1 = momUU[3]/rhoU[3];
    double vY1 = momVU[3]/rhoU[3];
    double pY1 = (gamma - 1.0)*(
        EU[3] - 0.5*rhoU[3]*uY1*uY1
              - 0.5*rhoU[3]*vY1*vY1
    );

    double rhoY2 = rhoD[1];
    double uY2 = momUD[1]/rhoD[1];
    double vY2 = momVD[1]/rhoD[1];
    double pY2 = (gamma - 1.0)*(
        ED[1] - 0.5*rhoD[1]*uY2*uY2
              - 0.5*rhoD[1]*vY2*vY2
    );

    double rhoY3 = rhoU[1];
    double uY3 = momUU[1]/rhoU[1];
    double vY3 = momVU[1]/rhoU[1];
    double pY3 = (gamma - 1.0)*(
        EU[1] - 0.5*rhoU[1]*uY3*uY3
              - 0.5*rhoU[1]*vY3*vY3
    );

    double rhoY4 = rhoD[4];
    double uY4 = momUD[4]/rhoD[4];
    double vY4 = momVD[4]/rhoD[4];
    double pY4 = (gamma - 1.0)*(
        ED[4] - 0.5*rhoD[4]*uY4*uY4
              - 0.5*rhoD[4]*vY4*vY4
    );

    fluxD = Hydro::getFluxHLLC(
        vY1, uY1, rhoY1, pY1,
        vY2, uY2, rhoY2, pY2,
        gamma);

    fluxU = Hydro::getFluxHLLC(
        vY3, uY3, rhoY3, pY3,
        vY4, uY4, rhoY4, pY4,
        gamma);

    fx = dt/dx;
    fy = dt/dy;

    GETNEW(RHO, iIndex, jIndex) = GETOLD(RHO, iIndex, jIndex) + fx*(fluxL.rho - fluxR.rho) + fy*(fluxD.rho - fluxU.rho);
    GETNEW(MOMU, iIndex, jIndex) = GETOLD(MOMU, iIndex, jIndex) + fx*(fluxL.momU - fluxR.momU) + fy*(fluxD.momV - fluxU.momV);
    GETNEW(MOMV, iIndex, jIndex) = GETOLD(MOMV, iIndex, jIndex) + fx*(fluxL.momV - fluxR.momV) + fy*(fluxD.momU - fluxU.momU);
    GETNEW(ENERGY, iIndex, jIndex) = GETOLD(ENERGY, iIndex, jIndex) + fx*(fluxL.E - fluxR.E) + fy*(fluxD.E - fluxU.E);

}

double getSlopeX(amrex::Array4<amrex::Real> stateOld, int U, int i, int j, int niGhosts, double omega) {
    double di1 = GETOLD(U, i, j) - GETOLD(U, i-1, j);
    double di2 = GETOLD(U, i+1, j) - GETOLD(U, i, j);

    double diU = 0.5*(1.0 + omega)*di1 + 0.5*(1.0 - omega)*di2;

    diU *= getLimiter(di1, di2, omega);

    return diU;
}

double getSlopeY(amrex::Array4<amrex::Real> stateOld, int U, int i, int j, int niGhosts, double omega) {
    double di1 = GETOLD(U, i, j) - GETOLD(U, i, j-1);
    double di2 = GETOLD(U, i, j+1) - GETOLD(U, i, j);

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

inline double calcFluxRho(double rho, double u) {
    return rho*u;
}


inline double calcFluxMom(double rho, double u, double v, double p) {
    return rho*u*v + p;
}


inline double calcFluxE(double u, double E, double p) {
    return u*(E + p);
}

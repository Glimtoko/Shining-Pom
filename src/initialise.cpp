#include "AMReX.H"
#include "AMReX_Geometry.H"
#include "AMReX_FabArray.H"
#include <AMReX_ArrayLim.H>

#include "shiningpom.hpp"

void setGeometrySodX (
    amrex::Box const& bx,
    amrex::Array4<amrex::Real> const& a,
    amrex::Geometry const& geom,
    double gamma
)
{
    // Physical parameters for Sod
    const double uL = 0.0;
    const double uR = 0.0;
    const double rhoL = 1.0;
    const double rhoR = 0.125;
    const double pL = 1.0;
    const double pR = 0.1;
    const double xInt = 0.5;

    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    const auto x0 = geom.ProbLo(0) + dx[0]/2.0;

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                auto x = x0 + dx[0]*i;
                double e = x < xInt ? pL/((gamma - 1.0)*rhoL) : pR/((gamma - 1.0)*rhoR);

                a(i, j, k, QUANT_RHO) = x < xInt ? rhoL : rhoR;
                a(i, j, k, QUANT_MOMU) = x < xInt ? rhoL * uL : rhoR * uR;
                a(i, j, k, QUANT_MOMV) = 0.0;
                a(i, j, k, QUANT_E) = x < xInt ? rhoL*(0.5*uL*uL + e) : rhoR*(0.5*uR*uR + e);
                a(i, j, k, QUANT_DT) = 100.0;
            }
        }
    }
}


void setGeometrySodY (
    amrex::Box const& bx,
    amrex::Array4<amrex::Real> const& a,
    amrex::Geometry const& geom,
    double gamma
)
{
    // Physical parameters for Sod
    const double uL = 0.0;
    const double uR = 0.0;
    const double rhoL = 1.0;
    const double rhoR = 0.125;
    const double pL = 1.0;
    const double pR = 0.1;
    const double yInt = 0.5;

    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    const auto y0 = geom.ProbLo(1) + dx[1]/2.0;

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                auto y = y0 + dx[1]*j;
                double e = y < yInt ? pL/((gamma - 1.0)*rhoL) : pR/((gamma - 1.0)*rhoR);

                a(i, j, k, QUANT_RHO) = y < yInt ? rhoL : rhoR;
                a(i, j, k, QUANT_MOMV) = y < yInt ? rhoL * uL : rhoR * uR;
                a(i, j, k, QUANT_MOMU) = 0.0;
                a(i, j, k, QUANT_E) = y < yInt ? rhoL*(0.5*uL*uL + e) : rhoR*(0.5*uR*uR + e);
                a(i, j, k, QUANT_DT) = 100.0;
            }
        }
    }
}

void setGeometryTriple (
    amrex::Box const& bx,
    amrex::Array4<amrex::Real> const& a,
    amrex::Geometry const& geom,
    double gamma
)
{
    double xInt = 0, yInt = 5, x1 = 10.0, x2 = 20.0, y1 = 3.0, y2 = 7.0;
    double p;

    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    const auto x0 = geom.ProbLo(0) + dx[0]/2.0;
    const auto y0 = geom.ProbLo(1) + dx[1]/2.0;

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                auto x = x0 + dx[0]*i;
                auto y = y0 + dx[1]*j;
                double r = sqrt(pow(y - yInt, 2) + pow(x - xInt, 2));
                if (r <= 2.0) {
                    a(i, j, k, QUANT_RHO) = 1.0;
                    p = 1.0;
                } else {
                    a(i, j, k, QUANT_RHO) = 0.125;
                    p = 0.1;
                }

                r = sqrt(pow(y - yInt, 2) + pow(x - x1, 2));
                if (r <= 2.0) a(i, j, k, QUANT_RHO) = 0.5;

                r = sqrt(pow(y - y1, 2) + pow(x - x2, 2));
                if (r <= 1.5) a(i, j, k, QUANT_RHO) = 0.5;

                r = sqrt(pow(y - y2, 2) + pow(x - x2, 2));
                if (r <= 1.5) a(i, j, k, QUANT_RHO) = 0.5;

                double e = p/((gamma - 1.0)*a(i, j, k, QUANT_RHO));
                a(i, j, k, QUANT_E) = e*a(i, j, k, QUANT_RHO);
                a(i, j, k, QUANT_MOMU) = 0.0;
                a(i, j, k, QUANT_MOMV) = 0.0;
                a(i, j, k, QUANT_DT) = 100.0;
            }
        }
    }

}

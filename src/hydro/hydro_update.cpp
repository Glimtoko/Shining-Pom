#include "hydro.hpp"

#include "AMReX_FabArray.H"

void Hydro::update(
    amrex::Array4<amrex::Real> const& stateOld,
    amrex::Array4<amrex::Real> const& fluxX,
    amrex::Array4<amrex::Real> const& fluxY,
    int i, int j, int k,
    double dt, double dx, double dy
) {
    double fx = dt/dx;
    double fy = dt/dy;

    GET(stateOld, RHO, i, j) = GET(stateOld, RHO, i, j)
        + fx*(GET(fluxX, RHO, i-1, j) - GET(fluxX, RHO, i, j))
        + fy*(GET(fluxY, RHO, i, j-1) - GET(fluxY, RHO, i, j));

    GET(stateOld, MOMU, i, j) = GET(stateOld, MOMU, i, j)
        + fx*(GET(fluxX, MOMU, i-1, j) - GET(fluxX, MOMU, i, j))
        + fy*(GET(fluxY, MOMV, i, j-1) - GET(fluxY, MOMV, i, j));

    GET(stateOld, MOMV, i, j) = GET(stateOld, MOMV, i, j)
        + fx*(GET(fluxX, MOMV, i-1, j) - GET(fluxX, MOMV, i, j))
        + fy*(GET(fluxY, MOMU, i, j-1) - GET(fluxY, MOMU, i, j));

    GET(stateOld, ENERGY, i, j) = GET(stateOld, ENERGY, i, j)
        + fx*(GET(fluxX, ENERGY, i-1, j) - GET(fluxX, ENERGY, i, j))
        + fy*(GET(fluxY, ENERGY, i, j-1) - GET(fluxY, ENERGY, i, j));

}


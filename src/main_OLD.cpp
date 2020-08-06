#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Print.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Geometry.H"
#include "AMReX_FabArray.H"
#include "AMReX_MultiFab.H"
#include "AMReX_PlotFileUtil.H"
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>

#include "mpi.h"

#include "hydro/hydro.hpp"
#include "initialise.hpp"
#include "shiningpom.hpp"

#include <iostream>
#include <vector>
#include <fenv.h>

int main(int argc, char* argv[]) {
//     feenableexcept(FE_INVALID | FE_OVERFLOW);
    amrex::Initialize(argc, argv, MPI::COMM_WORLD);
    {

    // TODO: AMR

    // Test basic parallel information
    int myRank = amrex::ParallelDescriptor::MyProc();
    std::cout << "Rank: " << myRank << std::endl;

    // Parse user input
    amrex::ParmParse pp;

    // Gamma
    double gamma;
    if (!pp.query("gamma", gamma)) {
        gamma = 1.4;
    }

    // Max timesteps
    int maxSteps = -1;
    pp.query("nsteps", maxSteps);

    double tend = 0.25;
    pp.query("tend", tend);

    double dtOut = 2.0;
    pp.query("dtOut", dtOut);

    // Grid resolution at base layer
    amrex::Vector<int> nCells;
    pp.getarr("ncells", nCells);

    // X and Y physical ranges
    amrex::Vector<amrex::Real> xr {0.0, 1.0};
    if (!pp.queryarr("xrange", xr)) {
        amrex::Print() << "Cannot find xrange in input, "
                       << "so the default {0.0, 1.0} will be used\n";
    }

    amrex::Vector<amrex::Real> yr {0.0, 1.0};
    if (!pp.queryarr("yrange", yr)) {
        amrex::Print() << "Cannot find yrange in input, "
                       << "so the default {0.0, 1.0} will be used\n";
    }


    amrex::Print() << "maxSteps = " << maxSteps << std::endl;
    amrex::Print() << AMREX_SPACEDIM << std::endl;

    // Default boundary conditions
    auto xLeft = amrex::BCType::reflect_even;
    auto xRight = amrex::BCType::reflect_even;
    auto yDown = amrex::BCType::reflect_even;
    auto yUp = amrex::BCType::reflect_even;

    // Function pointer to set up
    void (*initFunction) (amrex::Box const&,
    amrex::Array4<amrex::Real> const&,
    amrex::Geometry const&,
    double);


    // Set problem details
    int problem;
    pp.get("problem", problem);

    if (problem == 1) {
        initFunction = &setGeometrySodX;
    } else if (problem == 2) {
        initFunction = &setGeometrySodY;
    } else if (problem == 4) {
        initFunction = &setGeometryTriple;
        xRight = amrex::BCType::reflect_odd;
        yDown = amrex::BCType::reflect_odd;
        yUp = amrex::BCType::reflect_odd;
    } else {
        // Default to X-Sod
        initFunction = &setGeometrySodX;
    }


    // Construct base layer box
    amrex::Box domainLogical(amrex::IntVect{AMREX_D_DECL(0, 0, 0)},
                             amrex::IntVect{AMREX_D_DECL(nCells[0]-1, nCells[1]-1, 0)});
    amrex::Print() << domainLogical << std::endl;

    // Problem geometry
    amrex::RealBox domainPhysical(AMREX_D_DECL(xr[0], yr[0], 0),
                                  AMREX_D_DECL(xr[1], yr[1], 0));
    amrex::Print() << domainPhysical << std::endl;

    // Coordinate system
    int coordType = 0;  // Cartesian geometry

    // Periodic boundaries?
    amrex::Array<int,AMREX_SPACEDIM> isPeriodic {AMREX_D_DECL(0, 0, 0)};

    // Box Array
    amrex::BoxArray boxArray;
    boxArray.define(domainLogical);
    boxArray.maxSize(250);
    amrex::Print() << boxArray << std::endl;


    // Full problem geometry
    amrex::Geometry geom(domainLogical, domainPhysical, coordType, isPeriodic);

    // Distribute
    amrex::DistributionMapping distMap {boxArray};
    amrex::Print() << distMap << std::endl;

    // Set data storage (5 quantities, 2 ghost layers)
    amrex::MultiFab stateOld(boxArray, distMap, 5, 2);

    // Set MH reconstruction storage (4 quants, due to no DT)
    amrex::MultiFab leftMH(boxArray, distMap, 4, 2);
    amrex::MultiFab rightMH(boxArray, distMap, 4, 2);
    amrex::MultiFab upMH(boxArray, distMap, 4, 2);
    amrex::MultiFab downMH(boxArray, distMap, 4, 2);

    // Set flux storage
    amrex::MultiFab fluxX(boxArray, distMap, 4, 2);
    amrex::MultiFab fluxY(boxArray, distMap, 4, 2);

    // Set boundary conditions
    amrex::Vector<amrex::BCRec> boundaries(stateOld.nComp());

    // X boundaries
    boundaries[QUANT_RHO].setLo(0, amrex::BCType::reflect_even);
    boundaries[QUANT_RHO].setHi(0, amrex::BCType::reflect_even);
    boundaries[QUANT_MOMU].setLo(0, xLeft);
    boundaries[QUANT_MOMU].setHi(0, xRight);
    boundaries[QUANT_MOMV].setLo(0, xLeft);
    boundaries[QUANT_MOMV].setHi(0, xRight);
    boundaries[QUANT_E].setLo(0, amrex::BCType::reflect_even);
    boundaries[QUANT_E].setHi(0, amrex::BCType::reflect_even);

    // Y boundaries
    boundaries[QUANT_RHO].setLo(1, amrex::BCType::reflect_even);
    boundaries[QUANT_RHO].setHi(1, amrex::BCType::reflect_even);
    boundaries[QUANT_MOMU].setLo(1, yDown);
    boundaries[QUANT_MOMU].setHi(1, yUp);
    boundaries[QUANT_MOMV].setLo(1, yDown);
    boundaries[QUANT_MOMV].setHi(1, yUp);
    boundaries[QUANT_E].setLo(1, amrex::BCType::reflect_even);
    boundaries[QUANT_E].setHi(1, amrex::BCType::reflect_even);

    // Problem initialisation
    for (amrex::MFIter mfi(stateOld); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.validbox();
        amrex::FArrayBox& fab = stateOld[mfi];
        amrex::Array4<amrex::Real> const& a = fab.array();

        initFunction(box, a, geom, gamma);
    }

    // Update boundaries
    stateOld.FillBoundary(geom.periodicity());
    amrex::FillDomainBoundary(stateOld, geom, boundaries);

    // Initialise reconstruction and flux "states"
    amrex::MultiFab::Copy(leftMH, stateOld, 0, 0, 4, 2);
    amrex::MultiFab::Copy(rightMH, stateOld, 0, 0, 4, 2);
    amrex::MultiFab::Copy(upMH, stateOld, 0, 0, 4, 2);
    amrex::MultiFab::Copy(downMH, stateOld, 0, 0, 4, 2);

    amrex::MultiFab::Copy(fluxX, stateOld, 0, 0, 4, 2);
    amrex::MultiFab::Copy(fluxY, stateOld, 0, 0, 4, 2);

    const auto dx = geom.CellSizeArray();

    amrex::Vector<std::string> varNames = {
        "Rho", "MomU", "MomV", "E", "DT"
    };

    double t = 0;
    double outNext = t + dtOut;
    for (int step=1; step<maxSteps; step++) {
        for (amrex::MFIter mfi(stateOld); mfi.isValid(); ++mfi) // Loop over grids
        {
            const amrex::Box& box = mfi.validbox();

            amrex::FArrayBox& fOld = stateOld[mfi];
            amrex::Array4<amrex::Real> const& sOld = fOld.array();

            amrex::ParallelFor(box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Hydro::getCellTimestep(
                    sOld,
                    i, j, k,
                    1.4, dx[0], dx[1], 0.6, 0.1
                    );
            });
        }
        double dt = stateOld.min(QUANT_DT);
        dt = std::min(dt, outNext - t);
        t += dt;
        amrex::Print() << "Step " << step
                        << ", t = " << t
                        << ", dt = " << dt
                        << std::endl;


        // MUSCL-Hancock data reconstruction
        for (amrex::MFIter mfi(stateOld); mfi.isValid(); ++mfi) // Loop over grids
        {
            const amrex::Box& box = mfi.growntilebox(1);

            amrex::FArrayBox& fOld = stateOld[mfi];
            amrex::FArrayBox& fLeft = leftMH[mfi];
            amrex::FArrayBox& fRight = rightMH[mfi];
            amrex::FArrayBox& fUp = upMH[mfi];
            amrex::FArrayBox& fDown = downMH[mfi];

            amrex::Array4<amrex::Real> const& sOld = fOld.array();
            amrex::Array4<amrex::Real> const& lMH = fLeft.array();
            amrex::Array4<amrex::Real> const& rMH = fRight.array();
            amrex::Array4<amrex::Real> const& uMH = fUp.array();
            amrex::Array4<amrex::Real> const& dMH = fDown.array();

            amrex::ParallelFor(box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Hydro::MUSCLHancock2D_Reconstruct(
                    sOld,
                    i, j, k,
                    gamma, dt, dx[0], dx[1],
                    lMH,
                    rMH,
                    dMH,
                    uMH
                );
            });
        }


        // Calculate and store fluxes
        for (amrex::MFIter mfi(stateOld); mfi.isValid(); ++mfi) // Loop over grids
        {
            const amrex::Box& box = mfi.growntilebox(1);

            amrex::FArrayBox& fLeft = leftMH[mfi];
            amrex::FArrayBox& fRight = rightMH[mfi];
            amrex::FArrayBox& fUp = upMH[mfi];
            amrex::FArrayBox& fDown = downMH[mfi];

            amrex::FArrayBox& ffX = fluxX[mfi];
            amrex::FArrayBox& ffY = fluxY[mfi];

            amrex::Array4<amrex::Real> const& lMH = fLeft.array();
            amrex::Array4<amrex::Real> const& rMH = fRight.array();
            amrex::Array4<amrex::Real> const& uMH = fUp.array();
            amrex::Array4<amrex::Real> const& dMH = fDown.array();

            amrex::Array4<amrex::Real> const& fX = ffX.array();
            amrex::Array4<amrex::Real> const& fY = ffY.array();

            amrex::ParallelFor(box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Hydro::calculateFluxes(
                    lMH, rMH, dMH, uMH,
                    i, j, k,
                    gamma,
                    fX, fY
                );
            });
        }


        // Update to new state
        for (amrex::MFIter mfi(stateOld); mfi.isValid(); ++mfi) // Loop over grids
        {
            const amrex::Box& box = mfi.validbox();

            amrex::FArrayBox& fOld = stateOld[mfi];
            amrex::Array4<amrex::Real> const& sOld = fOld.array();

            amrex::FArrayBox& ffX = fluxX[mfi];
            amrex::FArrayBox& ffY = fluxY[mfi];

            amrex::Array4<amrex::Real> const& fX = ffX.array();
            amrex::Array4<amrex::Real> const& fY = ffY.array();

            amrex::ParallelFor(box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Hydro::update(
                    sOld, fX, fY,
                    i, j, k,
                    dt, dx[0], dx[1]
                );
            });

        }

        // Update boundaries
        stateOld.FillBoundary(geom.periodicity());
        amrex::FillDomainBoundary(stateOld, geom, boundaries);

        if (t >= outNext) {
            outNext += dtOut;
            const std::string& pfname = amrex::Concatenate("output/sod",step);

            amrex::WriteSingleLevelPlotfile(pfname, stateOld, varNames, geom, t, 1);
        }

        if (t >= tend) break;
    }

    amrex::Print() << "DONE" << std::endl;
    }
    amrex::Finalize();
    return 0;
}

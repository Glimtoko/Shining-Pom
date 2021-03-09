#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_AmrCore.H>

#include "rapidjson/document.h"

#include "utilities.hpp"

using namespace amrex;

amrex::AmrInfo pom::GetAmrInfo ()
{
    AmrInfo amr;
    ParmParse pp("amr");

    pp.query("v",amr.verbose);

    pp.get("max_level", amr.max_level);

    int nlev = amr.max_level + 1;

    // Make the default ref_ratio = 2 for all levels.
    amr.ref_ratio.resize(amr.max_level);
    for (int i = 0; i < amr.max_level; ++i)
    {
      amr.ref_ratio[i] = 2 * IntVect::TheUnitVector();
    }

    pp.query("n_proper",amr.n_proper);
    pp.query("grid_eff",amr.grid_eff);
    int cnt = pp.countval("n_error_buf");
    if (cnt > 0) {
        Vector<int> neb;
        pp.getarr("n_error_buf",neb);
        int n = std::min(cnt, amr.max_level+1);
        for (int i = 0; i < n; ++i) {
            amr.n_error_buf[i] = IntVect(neb[i]);
        }
        for (int i = n; i <= amr.max_level; ++i) {
            amr.n_error_buf[i] = IntVect(neb[cnt-1]);
        }
    }

    cnt = pp.countval("n_error_buf_x");
    if (cnt > 0) {
        int idim = 0;
        Vector<int> neb;
        pp.getarr("n_error_buf_x",neb);
        int n = std::min(cnt, amr.max_level+1);
        for (int i = 0; i < n; ++i) {
            amr.n_error_buf[i][idim] = neb[i];
        }
        for (int i = n; i <= amr.max_level; ++i) {
            amr.n_error_buf[i][idim] = neb[n-1];
        }
    }

#if (AMREX_SPACEDIM > 1)
    cnt = pp.countval("n_error_buf_y");
    if (cnt > 0) {
        int idim = 1;
        Vector<int> neb;
        pp.getarr("n_error_buf_y",neb);
        int n = std::min(cnt, amr.max_level+1);
        for (int i = 0; i < n; ++i) {
            amr.n_error_buf[i][idim] = neb[i];
        }
        for (int i = n; i <= amr.max_level; ++i) {
            amr.n_error_buf[i][idim] = neb[n-1];
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    cnt = pp.countval("n_error_buf_z");
    if (cnt > 0) {
        int idim = 2;
        Vector<int> neb;
        pp.getarr("n_error_buf_z",neb);
        int n = std::min(cnt, amr.max_level+1);
        for (int i = 0; i < n; ++i) {
            amr.n_error_buf[i][idim] = neb[i];
        }
        for (int i = n; i <= max_level; ++i) {
            amr.n_error_buf[i][idim] = neb[n-1];
        }
    }
#endif

    // Read in the refinement ratio IntVects as integer AMREX_SPACEDIM-tuples.
    if (amr.max_level > 0)
    {
        const int nratios_vect = amr.max_level*AMREX_SPACEDIM;

        Vector<int> ratios_vect(nratios_vect);

        int got_vect = pp.queryarr("ref_ratio_vect",ratios_vect,0,nratios_vect);

        Vector<int> ratios;

        const int got_int = pp.queryarr("ref_ratio",ratios);

        if (got_int == 1 && got_vect == 1)
        {
            amrex::Abort("Only input *either* ref_ratio or ref_ratio_vect");
        }
        else if (got_vect == 1)
        {
            int k = 0;
            for (int i = 0; i < amr.max_level; i++)
            {
                for (int n = 0; n < AMREX_SPACEDIM; n++,k++)
                    amr.ref_ratio[i][n] = ratios_vect[k];
            }
        }
        else if (got_int == 1)
        {
            const int ncnt = ratios.size();
            for (int i = 0; i < ncnt && i < amr.max_level; ++i)
            {
                for (int n = 0; n < AMREX_SPACEDIM; n++) {
                    amr.ref_ratio[i][n] = ratios[i];
                }
            }
            for (int i = ncnt; i < amr.max_level; ++i)
            {
                for (int n = 0; n < AMREX_SPACEDIM; n++) {
                    amr.ref_ratio[i][n] = ratios.back();
                }
            }
        }
        else
        {
            if (amr.verbose) {
                amrex::Print() << "Using default ref_ratio = 2 at all levels\n";
            }
        }
    }

    // Read in max_grid_size.  Use defaults if not explicitly defined.
    cnt = pp.countval("max_grid_size");
    if (cnt > 0) {
        Vector<int> mgs;
        pp.getarr("max_grid_size",mgs);
        int last_mgs = mgs.back();
        mgs.resize(amr.max_level+1,last_mgs);

        amr.max_grid_size.resize(amr.max_level+1);
        for (int i = 0; i <= amr.max_level; ++i) {
            amr.max_grid_size[i] = IntVect{AMREX_D_DECL(mgs[i],mgs[i],mgs[i])};
        }
    }

    cnt = pp.countval("max_grid_size_x");
    if (cnt > 0) {
        int idim = 0;
        Vector<int> mgs;
        pp.getarr("max_grid_size_x",mgs);
        int n = std::min(cnt, amr.max_level+1);
        for (int i = 0; i < n; ++i) {
            amr.max_grid_size[i][idim] = mgs[i];
        }
        for (int i = n; i <= amr.max_level; ++i) {
            amr.max_grid_size[i][idim] = mgs[n-1];
        }
    }

#if (AMREX_SPACEDIM > 1)
    cnt = pp.countval("max_grid_size_y");
    if (cnt > 0) {
        int idim = 1;
        Vector<int> mgs;
        pp.getarr("max_grid_size_y",mgs);
        int n = std::min(cnt, amr.max_level+1);
        for (int i = 0; i < n; ++i) {
            amr.max_grid_size[i][idim] = mgs[i];
        }
        for (int i = n; i <= amr.max_level; ++i) {
            amr.max_grid_size[i][idim] = mgs[n-1];
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    cnt = pp.countval("max_grid_size_z");
    if (cnt > 0) {
        int idim = 2;
        Vector<int> mgs;
        pp.getarr("max_grid_size_z",mgs);
        int n = std::min(cnt, amr.max_level+1);
        for (int i = 0; i < n; ++i) {
            amr.max_grid_size[i][idim] = mgs[i];
        }
        for (int i = n; i <= amr.max_level; ++i) {
            amr.max_grid_size[i][idim] = mgs[n-1];
        }
    }
#endif

    // Read in the blocking_factors.  Use defaults if not explicitly defined.
    cnt = pp.countval("blocking_factor");
    if (cnt > 0) {
        Vector<int> bf;
        pp.getarr("blocking_factor",bf);
        int last_bf = bf.back();
        bf.resize(amr.max_level+1,last_bf);

        amr.blocking_factor.resize(amr.max_level+1);
        for (int i = 0; i <= amr.max_level; ++i) {
            amr.blocking_factor[i] = IntVect{AMREX_D_DECL(bf[i],bf[i],bf[i])};
        }
    }

    cnt = pp.countval("blocking_factor_x");
    if (cnt > 0) {
        int idim = 0;
        Vector<int> bf;
        pp.getarr("blocking_factor_x",bf);
        int n = std::min(cnt, amr.max_level+1);
        for (int i = 0; i < n; ++i) {
            amr.blocking_factor[i][idim] = bf[i];
        }
        for (int i = n; i <= amr.max_level; ++i) {
            amr.blocking_factor[i][idim] = bf[n-1];
        }
    }

#if (AMREX_SPACEDIM > 1)
    cnt = pp.countval("blocking_factor_y");
    if (cnt > 0) {
        int idim = 1;
        Vector<int> bf;
        pp.getarr("blocking_factor_y",bf);
        int n = std::min(cnt, amr.max_level+1);
        for (int i = 0; i < n; ++i) {
            amr.blocking_factor[i][idim] = bf[i];
        }
        for (int i = n; i <= amr.max_level; ++i) {
            amr.blocking_factor[i][idim] = bf[n-1];
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    cnt = pp.countval("blocking_factor_z");
    if (cnt > 0) {
        int idim = 2;
        Vector<int> bf;
        pp.getarr("blocking_factor_z",bf);
        int n = std::min(cnt, amr.max_level+1);
        for (int i = 0; i < n; ++i) {
            amr.blocking_factor[i][idim] = bf[i];
        }
        for (int i = n; i <= amr.max_level; ++i) {
            amr.blocking_factor[i][idim] = bf[n-1];
        }
    }
#endif

    return amr;
}

amrex::BCType::mathematicalBndryTypes GetBoundary(
    rapidjson::Document &mesh, 
    std::string axis, 
    std::string variable,
    int idx
)
{
    if (mesh.HasMember("boundaries")) {
        if(mesh["boundaries"].HasMember(axis.c_str())) {
            rapidjson::Value &b =  mesh["boundaries"][axis.c_str()];
            if (b.HasMember(variable.c_str())) {
                rapidjson::Value &bc = b[variable.c_str()];
                std::string boundary{bc[idx].GetString()};
                if (boundary == "reflect_even")
                    return amrex::BCType::reflect_even;

                else if (boundary == "reflect_odd")
                    return amrex::BCType::reflect_odd;

                else if (boundary == "int_dir")
                    return amrex::BCType::int_dir;

                else if (boundary == "ext_dir")
                    return amrex::BCType::ext_dir;

                else if (boundary == "foextrap")
                    return amrex::BCType::foextrap;

                else return amrex::BCType::bogus;

            }
            return amrex::BCType::bogus;
        }
        return amrex::BCType::bogus;
    }
    return amrex::BCType::bogus;
}


pom::PomMesh pom::GetPomMesh(char *meshJSON, amrex::Vector<std::string> variables)
{
    pom::PomMesh pmesh;

    rapidjson::Document mesh;
    mesh.Parse(meshJSON);

    // Set boundary conditions
    pmesh.bcs.resize(variables.size());
    for (int i=0; i < variables.size(); i++)
    {
        std::string var = variables[i];
        std::transform(var.begin(), var.end(), var.begin(),
            [](unsigned char c){ return std::tolower(c); });

        pmesh.bcs[i].setLo(0, GetBoundary(mesh, "x", var, 0));
        pmesh.bcs[i].setHi(0, GetBoundary(mesh, "x", var, 1));
        pmesh.bcs[i].setLo(1, GetBoundary(mesh, "y", var, 3));
        pmesh.bcs[i].setHi(1, GetBoundary(mesh, "y", var, 2));
    }

    float x_low = 0.0;
    float x_high = 10.0;
    float y_low = 0.0;
    float y_high = 10.0;

    int x_cell = 10;
    int y_cell = 10;

    if (mesh.HasMember("mesh")) {
        rapidjson::Value &mesh_container = mesh["mesh"];

        if (mesh_container.HasMember("domain")) {
            rapidjson::Value &domain = mesh_container["domain"];
            if (domain.HasMember("bl")) {
                rapidjson::Value &bottom_left = domain["bl"];
                x_low = bottom_left[0].GetFloat();
                y_low = bottom_left[1].GetFloat();
            }

            if (domain.HasMember("tr")) {
                rapidjson::Value &top_right = domain["tr"];
                x_high = top_right[0].GetFloat();
                y_high = top_right[1].GetFloat();
            }
        }

        if (mesh_container.HasMember("cells")) {
            rapidjson::Value &cells = mesh_container["cells"];
            x_cell = cells[0].GetInt();
            y_cell = cells[1].GetInt();
        }
    }

    Print() << x_cell << " " << y_cell << std::endl;

        // This defines a Box with n_cell cells in each direction.
    Box domain(IntVect{AMREX_D_DECL(       0,        0,        0)},
            IntVect{AMREX_D_DECL(x_cell-1, y_cell-1, 0)});

    // This defines the physical box, [0,1] in each direction.
    RealBox real_box({AMREX_D_DECL(x_low, y_low, 0.0)},
                    {AMREX_D_DECL( x_high, y_high, 0.0)});

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to be non-periodic
    Array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(0,0,0)};

    // This defines a Geometry object
    Geometry geom(domain, real_box, coord, is_periodic);
    pmesh.geom = geom;

    return pmesh;
}

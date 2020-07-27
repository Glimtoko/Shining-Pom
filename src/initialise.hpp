void setGeometrySodX (
    amrex::Box const&,
    amrex::Array4<amrex::Real> const&,
    amrex::Geometry const&,
    double
);


void setGeometrySodY (
    amrex::Box const&,
    amrex::Array4<amrex::Real> const&,
    amrex::Geometry const&,
    double
);

void setGeometryTriple (
    amrex::Box const& bx,
    amrex::Array4<amrex::Real> const& a,
    amrex::Geometry const& geom,
    double gamma
);

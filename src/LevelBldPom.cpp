
#include <AMReX_LevelBld.H>
#include "PomLevel.hpp"

using namespace amrex;

class LevelBldPom
    :
    public LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual AmrLevel *operator() () override;
    virtual AmrLevel *operator() (
        Amr&            papa,
        int             lev,
        const Geometry& level_geom,
        const BoxArray& ba,
        const DistributionMapping& dm,
        Real            time
    ) override;
};

LevelBldPom Pom_bld;

LevelBld*
getLevelBld ()
{
    return &Pom_bld;
}

void
LevelBldPom::variableSetUp ()
{
    PomLevel::variableSetUp();
}

void
LevelBldPom::variableCleanUp ()
{
    PomLevel::variableCleanUp();
}

AmrLevel*
LevelBldPom::operator() ()
{
    return new PomLevel;
}

AmrLevel*
LevelBldPom::operator() (
    Amr&            papa,
    int             lev,
    const Geometry& level_geom,
    const BoxArray& ba,
    const DistributionMapping& dm,
    Real            time
)
{
    return new PomLevel(papa, lev, level_geom, ba, dm, time);
}

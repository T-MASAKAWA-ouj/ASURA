#include "config.h"

static bool SphericalAbsorptionBoundary(const int index){

    double R = NORM(Pbody[index]->Pos);

    if(R > BOUNDARY_CONDITION_SPHERICAL_SHELL_EDGE/Pall.UnitLength){
        return true;
    } else {
        return false;
    }
}

void ImposeBoundaryCondition(const int mode){

#ifdef USE_BOUNDARY_CONDITION //{
    switch(mode){
        case 0: // Periodic boundary
            break;
        case 1: // Box boundary with reflecting walls.
            break;
        case 2: // Absorption boundary with a box shape.
            break;
        case 3: // Absorption boundary with a spherical shell.
            ParticleRemoverArbitraryCriteria(SphericalAbsorptionBoundary);
            break;
    }
#endif // USE_BOUNDARY_CONDITION //}
    return ;
}


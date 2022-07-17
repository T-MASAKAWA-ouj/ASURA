#include "config.h"
#include "NeighborSearch.h"

/*
 * This file includes routines which check the members of structrures.
 * These are usuful for debugging.
 */

/*
 * This is the core part of the subroutine CheckHydroStructures(const int mode).
 */
static void CheckHydroStructuresEngine(const int Index){

    for(int k=0;k<3;k++){
        if((fpclassify(PhydroBody(Index)->Pos[k]) == FP_INFINITE)||
           (fpclassify(PhydroBody(Index)->Pos[k]) == FP_NAN)){
            fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s[%d] is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                    PhydroBody(Index)->GlobalID,Index,"Pos",k,PhydroBody(Index)->Pos[k]);
            fflush(NULL);
            assert(fpclassify(PhydroBody(Index)->Pos[k]) != FP_INFINITE);
            assert(fpclassify(PhydroBody(Index)->Pos[k]) != FP_NAN);
        }

        if((fpclassify(PhydroBody(Index)->Vel[k]) == FP_INFINITE)||
           (fpclassify(PhydroBody(Index)->Vel[k]) == FP_NAN)){
            fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s[%d] is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                    PhydroBody(Index)->GlobalID,Index,"Vel",k,PhydroBody(Index)->Vel[k]);
            fflush(NULL);
            assert(fpclassify(PhydroBody(Index)->Vel[k]) != FP_INFINITE);
            assert(fpclassify(PhydroBody(Index)->Vel[k]) != FP_NAN);
        }

        if((fpclassify(Phydro[Index]->HydroAcc[k]) == FP_INFINITE)||
           (fpclassify(Phydro[Index]->HydroAcc[k]) == FP_NAN)){
            fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s[%d] is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                    PhydroBody(Index)->GlobalID,Index,"HydroAcc",k,Phydro[Index]->HydroAcc[k]);
            fflush(NULL);
            assert(fpclassify(Phydro[Index]->HydroAcc[k]) != FP_INFINITE);
            assert(fpclassify(Phydro[Index]->HydroAcc[k]) != FP_NAN);
        }

        if((fpclassify(Phydro[Index]->Rot[k]) == FP_INFINITE)||
           (fpclassify(Phydro[Index]->Rot[k]) == FP_NAN)){
            fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s[%d] is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                    PhydroBody(Index)->GlobalID,Index,"Rot",k,Phydro[Index]->Rot[k]);
            fflush(NULL);
            assert(fpclassify(Phydro[Index]->Rot[k]) != FP_INFINITE);
            assert(fpclassify(Phydro[Index]->Rot[k]) != FP_NAN);
        }
    }

    if((fpclassify(Phydro[Index]->Rho) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Rho) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Rho",Phydro[Index]->Rho);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Rho) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Rho) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->Kernel) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Kernel) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Kernel",Phydro[Index]->Kernel);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Kernel) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Kernel) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->EnergyDensity) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->EnergyDensity) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"EnergyDensity",Phydro[Index]->EnergyDensity);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->EnergyDensity) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->EnergyDensity) != FP_NAN);
    }

#ifdef USE_GRAD_H //{
    if((fpclassify(Phydro[Index]->Gradh) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Gradh) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Gradh",Phydro[Index]->Gradh);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Gradh) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Gradh) != FP_NAN);
    }
#ifdef USE_GRAD_N //{
    if((fpclassify(Phydro[Index]->NumberDensity) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->NumberDensity) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"NumberDensity",Phydro[Index]->NumberDensity);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->NumberDensity) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->NumberDensity) != FP_NAN);
    }
    if((fpclassify(Phydro[Index]->GradN) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->GradN) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"GradN",Phydro[Index]->GradN);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->GradN) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->GradN) != FP_NAN);
    }
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
    
    if((fpclassify(Phydro[Index]->Div) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Div) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Div",Phydro[Index]->Div);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Div) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Div) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->F) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->F) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"F",Phydro[Index]->F);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->F) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->F) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->U) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->U) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"U",Phydro[Index]->U);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->U) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->U) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->Du) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Du) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Du",Phydro[Index]->Du);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Du) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Du) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->DuCooling) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->DuCooling) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"DuCooling",Phydro[Index]->DuCooling);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->DuCooling) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->DuCooling) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->DQheat) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->DQheat) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"DQheat",Phydro[Index]->DQheat);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->DQheat) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->DQheat) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->Z) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Z) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Z",Phydro[Index]->Z);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Z) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Z) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->Vsig) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Vsig) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Vsig",Phydro[Index]->Vsig);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Vsig) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Vsig) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->Alpha) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Alpha) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Alpha",Phydro[Index]->Alpha);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Alpha) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Alpha) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->DAlpha) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->DAlpha) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"DAlpha",Phydro[Index]->DAlpha);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->DAlpha) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->DAlpha) != FP_NAN);
    }

#ifdef USE_CELIB //{
#endif // USE_CELIB //}
    return ;
}

/*
 * This function checks members in the hydro structure. If there are inf or nan,
 * this function prints the result and calls ``assert''.
 */
void CheckHydroStructures(const int mode){

    if(mode == NONE){ // all hydro structures.
        for(int i=0;i<Pall.Nhydro;i++){
            CheckHydroStructuresEngine(i);
        }
    } else { // a special hydro structure.
        CheckHydroStructuresEngine(mode);
    }

    return ;
}



#include "config.h"
#include "KernelFunctions.h"

void ClearHydroData(void){

    for(int i=0;i<Pall.Nhydro;i++){
        if(PhydroActive(i)){
#ifdef EVALUATE_KERNEL_BY_ITERATION
            Phydro[i]->Rho = 
#endif // EVALUATE_KERNEL_BY_ITERATION
#ifdef USE_DISPH //{
            Phydro[i]->EnergyDensity = 
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
            Phydro[i]->PseudoDensity = 
#endif // USE_SPSPH //}
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
            Phydro[i]->SmoothedNumber = 
            Phydro[i]->DZw = 
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
#ifdef USE_GRAD_H //{
            Phydro[i]->Gradh = 
#ifdef USE_GRAD_N //{
            Phydro[i]->NumberDensity = 
            Phydro[i]->GradN = 
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
            Phydro[i]->DZdiff = 
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}
            Phydro[i]->Div = 
            Phydro[i]->Rot[0] = 
            Phydro[i]->Rot[1] = 
            Phydro[i]->Rot[2] = 
            Phydro[i]->Du = 
            Phydro[i]->HydroAcc[0] = 
            Phydro[i]->HydroAcc[1] = 
            Phydro[i]->HydroAcc[2] = 0.e0;
        }
    }
    return;
}

void InitializeHydroInteractionFlags(void){

    NumberofBufferHydroInteractionFlagsAllocated = 0;

    return ;
}

void CheckSizeofBufferHydroInteractionFlags(const int NumberofBufferHydroInteractionFlags){

    static bool first = true;
    if(first){
        InitializeHydroInteractionFlags();
        first = false;
    }

    int NProcs = MPIGetNumProcs();
    int BufferSize = NumberofBufferHydroInteractionFlags*(NProcs-1);
    if(NumberofBufferHydroInteractionFlagsAllocated < BufferSize){
        BufferSize = (int)(ForAngelsShare*BufferSize);
        if(NumberofBufferHydroInteractionFlagsAllocated == 0){
            BufferSize = MAX(FirstAllocationSize,BufferSize);
            BufferHydroInteractionFlags = malloc(sizeof(bool)*BufferSize);
        } else {
            BufferHydroInteractionFlags = realloc(BufferHydroInteractionFlags,sizeof(bool)*BufferSize);
        }
        NumberofBufferHydroInteractionFlagsAllocated = BufferSize;
    }
    return;
}

void ReleaseBufferHydroInteractionFlags(void){
    free(BufferHydroInteractionFlags);
    return;
}

void WriteKernelShape(void){

    if(MPIGetMyID() == MPI_ROOT_RANK){
        FILE *fp;
        char fname[MaxCharactersInLine];
#if KERNEL_TYPE == 0
        sprintf(fname,"./Kernel.sp");
#elif KERNEL_TYPE == 1 // Wendland C2
        sprintf(fname,"./Kernel.WC2");
#elif KERNEL_TYPE == 2 // Wendland C4
        sprintf(fname,"./Kernel.WC4");
#elif KERNEL_TYPE == 3 // Wendland C6
        sprintf(fname,"./Kernel.WC6");
#endif
        FileOpen(fp,fname,"w");

        double dr = 2.0/100.0;
        double h = 1.0;
        for(int i=0;i<101;i++){
            double r = dr*i;
            double w = SPHKernel(r,1.0/h);
            double dw = dSPHKernel(r,1.0/h);
            fprintf(fp,"%g %g %g\n",r,w,dw);
        }
        fclose(fp);
    }

    return ;
}


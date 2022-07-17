#include "config.h"

void Kick1Drift(void){ // First Kick and Drift.

    double TimingResultThisRoutine = GetElapsedTime();

    for(int i=0;i<Pall.Ntotal;i++){
        if((Pbody[i]->Active)&&(Pbody[i]->Type != TypeHydro)){
            double dt_half_gravity = 0.5*Pbody[i]->dt;
            Pbody[i]->Velh[0] = Pbody[i]->Vel[0]+dt_half_gravity*Pbody[i]->Acc[0];
            Pbody[i]->Velh[1] = Pbody[i]->Vel[1]+dt_half_gravity*Pbody[i]->Acc[1];
            Pbody[i]->Velh[2] = Pbody[i]->Vel[2]+dt_half_gravity*Pbody[i]->Acc[2];

            Pbody[i]->Pos[0] += Pbody[i]->dt*Pbody[i]->Velh[0];
            Pbody[i]->Pos[1] += Pbody[i]->dt*Pbody[i]->Velh[1];
            Pbody[i]->Pos[2] += Pbody[i]->dt*Pbody[i]->Velh[2];

            Pbody[i]->AccOld[0] = Pbody[i]->Acc[0];
            Pbody[i]->AccOld[1] = Pbody[i]->Acc[1];
            Pbody[i]->AccOld[2] = Pbody[i]->Acc[2];
        }
    }

    const static double OneThird = 1.e0/3.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            double dt_half_hydro = 0.5*Phydro[i]->dt_hydro;
            PhydroBody(i)->Velh[0] = PhydroBody(i)->Vel[0] + dt_half_hydro*Phydro[i]->HydroAcc[0];
            PhydroBody(i)->Velh[1] = PhydroBody(i)->Vel[1] + dt_half_hydro*Phydro[i]->HydroAcc[1];
            PhydroBody(i)->Velh[2] = PhydroBody(i)->Vel[2] + dt_half_hydro*Phydro[i]->HydroAcc[2];
            if(Phydro[i]->GravityKickFlag){
                double dt_half_gravity = 0.5*PhydroBody(i)->dt;
                PhydroBody(i)->Velh[0] += dt_half_gravity*PhydroBody(i)->Acc[0];
                PhydroBody(i)->Velh[1] += dt_half_gravity*PhydroBody(i)->Acc[1];
                PhydroBody(i)->Velh[2] += dt_half_gravity*PhydroBody(i)->Acc[2];
                PhydroBody(i)->AccOld[0] = PhydroBody(i)->Acc[0];
                PhydroBody(i)->AccOld[1] = PhydroBody(i)->Acc[1];
                PhydroBody(i)->AccOld[2] = PhydroBody(i)->Acc[2];
            }
            PhydroBody(i)->Pos[0] += Phydro[i]->dt_hydro*PhydroBody(i)->Velh[0];
            PhydroBody(i)->Pos[1] += Phydro[i]->dt_hydro*PhydroBody(i)->Velh[1];
            PhydroBody(i)->Pos[2] += Phydro[i]->dt_hydro*PhydroBody(i)->Velh[2];

            Phydro[i]->Kernel *= exp(OneThird*Phydro[i]->dt_hydro*Phydro[i]->Div);
#ifdef USE_MINIMUM_KERNEL_SIZE
            if(2.0*Phydro[i]->Kernel<=PhydroBody(i)->Eps*Pall.AdaptiveSofteningFactor)
                Phydro[i]->Kernel = 0.5*PhydroBody(i)->Eps*Pall.AdaptiveSofteningFactor;
#endif

#ifdef TASK_M2_COLLAPSE //{
            Phydro[i]->Du = 0.0;
#endif // TASK_M2_COLLAPSE //{

#if !(defined(ISOTHERMAL_EOS_RUN)||defined(BAROTROPIC_EOS_RUN))
            //Phydro[i]->U += dt_half_hydro*Phydro[i]->Du;

#ifdef COOLING_RUN //{
#ifdef USE_COOLING_ENERGYLOSS_LIMIT //{
            if(-Phydro[i]->U>4.0*dt_half_hydro*Phydro[i]->Du){
                Phydro[i]->U *= 0.75;
            } else 
#endif // USE_COOLING_ENERGYLOSS_LIMIT // }
#ifdef USE_MAXIMUM_TEMPERATURE //{
            if(Pall.ConvertUtoT*(Phydro[i]->U+dt_half_hydro*Phydro[i]->Du) > MAXIMUM_TEMPERATURE){
                Phydro[i]->U = Pall.ConvertTtoU*MAXIMUM_TEMPERATURE;
            } else 
#endif // USE_MAXIMUM_TEMPERATURE // }
            {
                Phydro[i]->U += dt_half_hydro*Phydro[i]->Du;
            }
#else // COOLING_RUN //}//{
            Phydro[i]->U += dt_half_hydro*Phydro[i]->Du;
#endif // COOLING_RUN //}
#endif


#ifdef USE_SPSPH //{
            Phydro[i]->Zw += dt_half_hydro*Phydro[i]->DZw;
#endif // USE_SPSPH //}

        }
    }    

    TimingResults.IntegralThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
}

void BuildPredictors(void){ // This routine must call just after calling the routine "KickFDrift()".

    double TimingResultThisRoutine = GetElapsedTime();

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type != TypeHydro){
            if(Pbody[i]->Active){
                Pbody[i]->PosP[0] = Pbody[i]->Pos[0];
                Pbody[i]->PosP[1] = Pbody[i]->Pos[1];
                Pbody[i]->PosP[2] = Pbody[i]->Pos[2];
            } else {
                double TDelta = Pall.EraLocal - Pbody[i]->EraLocal;
                double Velh[3]; // make velhalf
                double dt_half_gravity = 0.5*Pbody[i]->dt;
                Velh[0] = Pbody[i]->Vel[0] + dt_half_gravity*Pbody[i]->Acc[0];
                Velh[1] = Pbody[i]->Vel[1] + dt_half_gravity*Pbody[i]->Acc[1];
                Velh[2] = Pbody[i]->Vel[2] + dt_half_gravity*Pbody[i]->Acc[2];
                Pbody[i]->PosP[0] = Pbody[i]->Pos[0]+TDelta*Velh[0];
                Pbody[i]->PosP[1] = Pbody[i]->Pos[1]+TDelta*Velh[1];
                Pbody[i]->PosP[2] = Pbody[i]->Pos[2]+TDelta*Velh[2];
            }
        }
    }

    const static double OneThird = 1.e0/3.e0;
    //const static double TwoThird = 2.e0/3.e0;
    for(int i=0;i<Pall.Nhydro;i++){

        if(Phydro[i]->Active){
            PhydroBody(i)->PosP[0] = PhydroBody(i)->Pos[0];
            PhydroBody(i)->PosP[1] = PhydroBody(i)->Pos[1];
            PhydroBody(i)->PosP[2] = PhydroBody(i)->Pos[2];

            Phydro[i]->PosP[0] = PhydroBody(i)->PosP[0];
            Phydro[i]->PosP[1] = PhydroBody(i)->PosP[1];
            Phydro[i]->PosP[2] = PhydroBody(i)->PosP[2];

            // There would be strange....
            double dt_half_hydro = 0.5*Phydro[i]->dt_hydro;
            Phydro[i]->VelP[0] = PhydroBody(i)->Velh[0]+dt_half_hydro*Phydro[i]->HydroAcc[0];
            Phydro[i]->VelP[1] = PhydroBody(i)->Velh[1]+dt_half_hydro*Phydro[i]->HydroAcc[1];
            Phydro[i]->VelP[2] = PhydroBody(i)->Velh[2]+dt_half_hydro*Phydro[i]->HydroAcc[2];
            //if(Phydro[i]->GravityKickFlag){
            if(PhydroBody(i)->Active){
                double dt_half_gravity = 0.5*PhydroBody(i)->dt;
                Phydro[i]->VelP[0] += dt_half_gravity*PhydroBody(i)->Acc[0];
                Phydro[i]->VelP[1] += dt_half_gravity*PhydroBody(i)->Acc[1];
                Phydro[i]->VelP[2] += dt_half_gravity*PhydroBody(i)->Acc[2];
            }

            Phydro[i]->KernelPred = Phydro[i]->Kernel;
#if !(defined(ISOTHERMAL_EOS_RUN)||defined(BAROTROPIC_EOS_RUN)) //{
            //Phydro[i]->UPred = Phydro[i]->U+dt_half_hydro*Phydro[i]->Du;
#ifdef COOLING_RUN //{

#ifdef TASK_M2_COLLAPSE //{
            Phydro[i]->Du = 0.0;
#endif // TASK_M2_COLLAPSE //{

#ifdef USE_COOLING_ENERGYLOSS_LIMIT //{
            if(-Phydro[i]->U>4.0*dt_half_hydro*Phydro[i]->Du){
                Phydro[i]->UPred = 0.66667*Phydro[i]->U;
            } else 
#endif // USE_COOLING_ENERGYLOSS_LIMIT // }
#ifdef USE_MAXIMUM_TEMPERATURE //{
            if(Pall.ConvertUtoT*(Phydro[i]->U+dt_half_hydro*Phydro[i]->Du) > MAXIMUM_TEMPERATURE){
                Phydro[i]->UPred = Pall.ConvertTtoU*MAXIMUM_TEMPERATURE;
            } else
#endif // USE_MAXIMUM_TEMPERATURE // }
            {
                Phydro[i]->UPred = Phydro[i]->U+dt_half_hydro*Phydro[i]->Du;
            }
#else // COOLING_RUN //}//{
            Phydro[i]->UPred = Phydro[i]->U+dt_half_hydro*Phydro[i]->Du;
#endif // COOLING_RUN //}

#if 0
#if defined(COOLING_RUN)&&defined(USE_COOLING_ENERGYLOSS_LIMIT) //{
            //Phydro[i]->UPred = Phydro[i]->U+dt_half_hydro*Phydro[i]->Du;
            if(-Phydro[i]->U>4.0*dt_half_hydro*Phydro[i]->Du){
                Phydro[i]->UPred = 0.66667*Phydro[i]->U;
            } else {
                Phydro[i]->UPred = Phydro[i]->U+dt_half_hydro*Phydro[i]->Du;
            }
#else // defined(COOLING_RUN)&&defined(USE_COOLING_ENERGYLOSS_LIMIT) //}//{
            Phydro[i]->UPred = Phydro[i]->U+dt_half_hydro*Phydro[i]->Du;
#endif // defined(COOLING_RUN)&&defined(USE_COOLING_ENERGYLOSS_LIMIT) //}
#endif

            Phydro[i]->DuPrev = Phydro[i]->Du;
#endif //}
            if(Phydro[i]->UPred < 0.e0){
                fprintf(stderr,"Predictor error %g U[%d] %g %s;%d\n",Phydro[i]->UPred,i,Phydro[i]->U,__FUNCTION__,__LINE__);
                fprintf(stderr,"UPred, U, dt_half_hydro, Du %g %g %g %g\n",
                        Phydro[i]->UPred,Phydro[i]->U,dt_half_hydro,Phydro[i]->Du);
                fprintf(stderr,"[%02d] ID, GlobalID %d %ld\n",MPIGetMyID(),i,PhydroBody(i)->GlobalID);
                
                StructureReportPhydro(i);

                fflush(NULL);
                abort();
            }
#ifdef USE_SPSPH //{
            Phydro[i]->ZwPred = Phydro[i]->Zw + dt_half_hydro*Phydro[i]->DZw;
#endif // USE_SPSPH //}
        } else { 

            double TDelta = Pall.EraLocal - Phydro[i]->EraLocal_hydro;

            double Velh[3]; // make velhalf
            double dt_half_hydro = 0.5*Phydro[i]->dt_hydro;
            Velh[0] = PhydroBody(i)->Vel[0]+dt_half_hydro*Phydro[i]->HydroAcc[0];
            Velh[1] = PhydroBody(i)->Vel[1]+dt_half_hydro*Phydro[i]->HydroAcc[1];
            Velh[2] = PhydroBody(i)->Vel[2]+dt_half_hydro*Phydro[i]->HydroAcc[2];
            if(Phydro[i]->GravityKickFlag){
                double dt_half_gravity = 0.5*PhydroBody(i)->dt;
                Velh[0] += dt_half_gravity*PhydroBody(i)->Acc[0];
                Velh[1] += dt_half_gravity*PhydroBody(i)->Acc[1];
                Velh[2] += dt_half_gravity*PhydroBody(i)->Acc[2];
            }

            PhydroBody(i)->PosP[0] = PhydroBody(i)->Pos[0]+TDelta*Velh[0];
            PhydroBody(i)->PosP[1] = PhydroBody(i)->Pos[1]+TDelta*Velh[1];
            PhydroBody(i)->PosP[2] = PhydroBody(i)->Pos[2]+TDelta*Velh[2];
            Phydro[i]->PosP[0] = PhydroBody(i)->PosP[0];
            Phydro[i]->PosP[1] = PhydroBody(i)->PosP[1];
            Phydro[i]->PosP[2] = PhydroBody(i)->PosP[2];

            Phydro[i]->VelP[0] = PhydroBody(i)->Vel[0]+TDelta*Phydro[i]->HydroAcc[0];
            Phydro[i]->VelP[1] = PhydroBody(i)->Vel[1]+TDelta*Phydro[i]->HydroAcc[1];
            Phydro[i]->VelP[2] = PhydroBody(i)->Vel[2]+TDelta*Phydro[i]->HydroAcc[2];
            if(Phydro[i]->GravityKickFlag){
                Phydro[i]->VelP[0] += TDelta*PhydroBody(i)->Acc[0];
                Phydro[i]->VelP[1] += TDelta*PhydroBody(i)->Acc[1];
                Phydro[i]->VelP[2] += TDelta*PhydroBody(i)->Acc[2];
            }

// #ifdef USE_DENSITYKERNEL_PREDICTOR //{
            Phydro[i]->RhoPred = Phydro[i]->Rho*exp(-TDelta*Phydro[i]->Div);
            Phydro[i]->KernelPred = Phydro[i]->Kernel*exp(OneThird*TDelta*Phydro[i]->Div);
// #else // USE_DENSITYKERNEL_PREDICTOR //}//{
            // Phydro[i]->RhoPred = Phydro[i]->Rho;
            // Phydro[i]->KernelPred = Phydro[i]->Kernel;
// #endif // USE_DENSITYKERNEL_PREDICTOR //}
#ifdef USE_MINIMUM_KERNEL_SIZE
            if(2.0*Phydro[i]->KernelPred<=PhydroBody(i)->Eps*Pall.AdaptiveSofteningFactor)
                Phydro[i]->Kernel = 0.5*PhydroBody(i)->Eps*Pall.AdaptiveSofteningFactor;
#endif // USE_MINIMUM_KERNEL_SIZE


#ifdef USE_DISPH
            Phydro[i]->EnergyDensityPred = Phydro[i]->EnergyDensity*exp(-TDelta*Phydro[i]->Div);
#endif //USE_DISPH
#ifdef USE_SPSPH //{
            // Need prediction of pseudo-density
            Phydro[i]->PseudoDensityPred = Phydro[i]->PseudoDensity*exp(-TDelta*Phydro[i]->Div);
            Phydro[i]->ZwPred = Phydro[i]->Zw + TDelta*Phydro[i]->DZw;
#endif // USE_SPSPH //}

#ifdef TASK_M2_COLLAPSE //{
            Phydro[i]->Du = 0.0;
#endif // TASK_M2_COLLAPSE //{
#if (defined(ISOTHERMAL_EOS_RUN)||defined(BAROTROPIC_EOS_RUN))
            Phydro[i]->UPred = Phydro[i]->U;
#else 
#ifdef COOLING_RUN //{
#ifdef USE_COOLING_ENERGYLOSS_LIMIT //{
            if(-Phydro[i]->U>2.0*TDelta*Phydro[i]->Du){
                Phydro[i]->UPred = 0.5*Phydro[i]->U;
            } else 
#endif // USE_COOLING_ENERGYLOSS_LIMIT //}
#ifdef USE_MAXIMUM_TEMPERATURE //{
            if(Pall.ConvertUtoT*(Phydro[i]->U+TDelta*Phydro[i]->Du) > MAXIMUM_TEMPERATURE){
                Phydro[i]->UPred = Pall.ConvertTtoU*MAXIMUM_TEMPERATURE;
            } else 
#endif // USE_MAXIMUM_TEMPERATURE // }
            {
                Phydro[i]->UPred = Phydro[i]->U+TDelta*Phydro[i]->Du;
            }
#else // COOLING_RUN //}//{
            Phydro[i]->UPred = Phydro[i]->U+TDelta*Phydro[i]->Du;
#endif // COOLING_RUN //}
#endif
            if(Phydro[i]->UPred < 0.e0){ // here!
                fprintf(stderr,"Predictor[%ld] error T %g %g %s;%d\n",PhydroBody(i)->GlobalID,Pall.ConvertUtoT*Phydro[i]->UPred,
                        Pall.ConvertNumberDensityToCGS*Phydro[i]->Rho,__FUNCTION__,__LINE__);
                fprintf(stderr,"Predictor[%ld] error %g %s;%d\n",PhydroBody(i)->GlobalID,Phydro[i]->UPred,__FUNCTION__,__LINE__);
                fprintf(stderr,"Upred %g Uold %g dT %g Du %g %s;%d\n",
                        Phydro[i]->UPred,Phydro[i]->U,TDelta,Phydro[i]->Du,__FUNCTION__,__LINE__);
                
                StructureReportPhydro(i);

                fflush(NULL);
                abort();
            }
        }
    }

#ifdef PERIODIC_RUN //{
    PeriodicWrapping();
#endif // PERIODIC_RUN //}


    TimingResults.IntegralThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
}

void Kick2(void){ // Second Kick.

    double TimingResultThisRoutine = GetElapsedTime();

    // The second kick by hydro acc for new stars/sinks are added in the
    // formation routine.
    for(int i=0;i<Pall.Ntotal;i++){
        if((Pbody[i]->Active)&&(Pbody[i]->Type != TypeHydro)){
            double dt_half = 0.5*Pbody[i]->dt;
            Pbody[i]->Vel[0] = Pbody[i]->Velh[0]+dt_half*Pbody[i]->Acc[0];
            Pbody[i]->Vel[1] = Pbody[i]->Velh[1]+dt_half*Pbody[i]->Acc[1];
            Pbody[i]->Vel[2] = Pbody[i]->Velh[2]+dt_half*Pbody[i]->Acc[2];
        }
    }

    double EnergyLossThisStep = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            double dt_half_hydro = 0.5*Phydro[i]->dt_hydro;
            PhydroBody(i)->Vel[0] = PhydroBody(i)->Velh[0]+dt_half_hydro*Phydro[i]->HydroAcc[0];
            PhydroBody(i)->Vel[1] = PhydroBody(i)->Velh[1]+dt_half_hydro*Phydro[i]->HydroAcc[1];
            PhydroBody(i)->Vel[2] = PhydroBody(i)->Velh[2]+dt_half_hydro*Phydro[i]->HydroAcc[2];
#ifdef GRAVITY_RUN
            if(PhydroBody(i)->Active){
                double dt_half_gravity = 0.5*PhydroBody(i)->dt;
                PhydroBody(i)->Vel[0] += dt_half_gravity*PhydroBody(i)->Acc[0];
                PhydroBody(i)->Vel[1] += dt_half_gravity*PhydroBody(i)->Acc[1];
                PhydroBody(i)->Vel[2] += dt_half_gravity*PhydroBody(i)->Acc[2];
            }
#endif // GRAVITY_RUN


#ifdef TASK_M2_COLLAPSE //{
            Phydro[i]->Du = 0.0;
#endif // TASK_M2_COLLAPSE //{

#if !(defined(ISOTHERMAL_EOS_RUN)||defined(BAROTROPIC_EOS_RUN))
#ifdef COOLING_RUN //{
#ifdef USE_COOLING_ENERGYLOSS_LIMIT //{
            if(-Phydro[i]->U>4.0*dt_half_hydro*Phydro[i]->Du){
                Phydro[i]->U *= 0.75;
            } else 
#endif // USE_COOLING_ENERGYLOSS_LIMIT // }
#ifdef USE_MAXIMUM_TEMPERATURE //{
            if(Pall.ConvertUtoT*(Phydro[i]->U+dt_half_hydro*Phydro[i]->Du) > MAXIMUM_TEMPERATURE){
                Phydro[i]->U = Pall.ConvertTtoU*MAXIMUM_TEMPERATURE;
            } else 
#endif // USE_MAXIMUM_TEMPERATURE // }
            {
                Phydro[i]->U += dt_half_hydro*Phydro[i]->Du;
            }
#else // COOLING_RUN //}//{
            Phydro[i]->U += dt_half_hydro*Phydro[i]->Du;
#endif // COOLING_RUN //}
            if(Phydro[i]->U < 0.e0){
                fprintf(stderr,"%ld Uold %g U %g, Du %g, dt %g %g %g %s:%d\n",PhydroBody(i)->GlobalID,
                        Phydro[i]->U-2.0*dt_half_hydro*Phydro[i]->Du,
                        Phydro[i]->U,Phydro[i]->Du,2.0*dt_half_hydro,
                        Phydro[i]->dt_hydro,PhydroBody(i)->dt,
                        __FILE__,__LINE__);
                StructureReportPhydro(i);

                fflush(NULL);
                abort();
            }
#endif
#ifdef USE_SPSPH //{
            Phydro[i]->Zw += dt_half_hydro*Phydro[i]->DZw;
#endif // USE_SPSPH //}
#ifdef USE_VARIABLE_ALPHA
            Phydro[i]->Alpha += Phydro[i]->DAlpha*Phydro[i]->dt_hydro;
            Phydro[i]->Alpha = fmax(fmin(Phydro[i]->Alpha,Pall.ViscousAlphaMax),Pall.ViscousAlphaMin);
#endif // USE_VARIABLE_ALPHA
        }
    }
    Pall.CoolingEnergyLoss += EnergyLossThisStep;

    TimingResults.IntegralThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return ;
}

void UpdateGravityKickFlag(void){ // Gravity first kick for Hydro

    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            if(PhydroBody(i)->Active){
                Phydro[i]->GravityKickFlag = true;
            } else {
                Phydro[i]->GravityKickFlag = false;
            }
        }
    }
    return;
}


void SetInitialHydroVelh(void){

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->VelP[0] = PhydroBody(i)->Vel[0];
        Phydro[i]->VelP[1] = PhydroBody(i)->Vel[1];
        Phydro[i]->VelP[2] = PhydroBody(i)->Vel[2];
    }
    return ;
}

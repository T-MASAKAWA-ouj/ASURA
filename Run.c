#include "config.h"
#include "PreDecomposition.h"
#include "Decomposition.h"
#include "ForceMisc.h"
#include "ForceGRAPE.h"
#include "ForceParallelTreeGRAPE.h"
#include "ForceFromExternalPotentials.h"
#include "PlantGravityTree.h"
#include "PlantHydroTree.h"
#include "Integral.h"
#include "TimeStep.h"
#include "HydroDensity.h"
#include "HydroAcc.h"
#include "HydroMisc.h"
#include "HydroKernel.h"
#include "HydroExtraOperations.h"
#include "TimeStep.h"
#include "Cooling.h"
#include "Heating.h"
#include "StarFormation.h"
#include "Delayed.h"
#include "StellarFeedback.h"
#include "HIIregion.h"
#include "SetUpTestRun.h"
#include "SinkParticle.h"
#include "SizeDetermination.h"
#include "ThermalConductivity.h"
#include "Logs.h"
#include "RunLogs.h"

#include <mcheck.h>


static int DecompFrequencyCounter = 0;

#if 1
void Run(void){

#if 0
    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./dt.%03d",MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){

        double ParticleEpsSize = Pbody[i]->Eps*Pall.AdaptiveSofteningFactor;
        double Acc[3] = {Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]};
        double normAcc = NORM(Acc);
        double DiffAcc[3] = {Acc[0]-Pbody[i]->AccOld[0],
                             Acc[1]-Pbody[i]->AccOld[1],
                             Acc[2]-Pbody[i]->AccOld[2]};
        double normDiffAcc = NORM(DiffAcc);
        double dt_diffa = TFactorDiffAcc*normAcc/normDiffAcc*Pbody[i]->dt;
        double dt_a = TFactorAcc*sqrt(ParticleEpsSize/normAcc);

        double dt_h = 0.01/Pall.HubbleZ;


        if(Pbody[i]->Type == TypeHydro){
            fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,Pbody[i]->dt,
                    dt_diffa,dt_a,dt_h,Pbody[i]->Mass,Pbody[i]->Eps,ParticleEpsSize,
                    PbodyHydro(i)->dt_hydro,PbodyHydro(i)->Kernel,PbodyHydro(i)->Vsig,
                    PbodyHydro(i)->Kernel/PbodyHydro(i)->Vsig);
        } else {
            fprintf(fp,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,Pbody[i]->dt,
                    dt_diffa,dt_a,dt_h,Pbody[i]->Mass,Pbody[i]->Eps,ParticleEpsSize);
        }
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat ./dt.??? | sort -n > ./dt.dat");
        fflush(NULL);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    // if(MPIGetMyID() == MPI_ROOT_RANK){
        // system("rm -rf ./dt.???");
        // fflush(NULL);
    // }

#endif

    while(Pall.TCurrent < Pall.TEnd){

        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"%ld : ",Pall.TStepTotal);
            fflush(NULL);
        }

#if 0
        gprintlmpi(Pall.AdaptiveSofteningFactor);
        gprintlmpi(Pall.Redshift);

    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(1);
#endif


        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();
        
        UpdateCosmologicalParameters();

        if(Pall.TCurrent >= Pall.Era){
            //PreDomainDecomposition();
            // InspectParticleFluctuation();
            PreDomainDecompositionNew(0);
            DomainDecompositionOnlyDataExchange();

            // InspectParticleFluctuation();

            SortAllStructures();

            FirstTimeStep();
            BuildHierarchicalTimeStep();
        } else {
            if(DecompFrequencyCounter > DECOMPOSITION_FREQUENCY){
                if (10*Pall.NActives_t > Pall.Ntotal_t){
                    PreDomainDecompositionNew(1);
                    DomainDecompositionOnlyDataExchange();
                    DecompFrequencyCounter = 0;
                }
            } else {
                DecompFrequencyCounter ++;
            }
            // InspectParticleFluctuation();
            BuildNewTimeStep();
        }

        RaiseActiveFlags();

        Kick1Drift(); 
        BuildPredictors();

#ifdef GRAVITY_RUN //{
        if(Pall.NActives_t>0){ // Gravity
            PlantGravityTree();

            ClearGravitationalForce();
            ForceParallelTreeGRAPE(); // Split this function? insert PlantHydroTree and ClearHydroData?
            ForceFromExternalPotentials();
        } else {
            if(MPIGetMyID() == MPI_ROOT_RANK)
                fprintf(stderr,"Skip Gravity\n");
        }
#endif // GRAVITY_RUN //}

#ifdef HYDRO_RUN //{
        if(Pall.NActivesHydro_t>0){ // Hydro
            PlantHydroTree();

            ClearHydroData();

#ifdef EVALUATE_SIZES_ALL_TOGETHER //{
            CalcSize();
#else // EVALUATE_SIZES_ALL_TOGETHER //}//{
            CalcKernel();
#endif // EVALUATE_SIZES_ALL_TOGETHER //}
            CalcDensityDivRot(); 
            CalcDuDtAcc();
        }
#endif // HYDRO_RUN //}

        Kick2();

#ifdef HYDRO_RUN //{
        if(Pall.NActivesHydro_t>0){ // Hydro
#ifdef USE_CELIB //{
            StellarFeedback();
#else
            DelayedSNe();
#endif  // USE_CELIB //}
            CalcCooling();
            HIIregions();

            CalcDuDtAccEnergyDensityForCorrection();

            StarFormation();
            SinkParticles();
        }
#endif // HYDRO_RUN //}

        UpdateGravityKickFlag();

        // Calc New Time Step or Out Put Logs
        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();
        //LogsThisTimeStep();

        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            if(Pall.Nhydro_t>0)
                LogStarFormationRate();
            OutPutAllParticlesInASCIIFormat();

            if(Pall.Nhydro_t>0){
                CountNeighborNumber();
            }

            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;

            // DataDump();
            FileOutPutConstantInterval();

            LogOutPutEnergyMomentumAngularMomentum();
            LogOutPutElapsedTime();
            LogTotalMass();
            ReportAllocatedMemorySizes();
            EndLogs();
            fflush(NULL);
            MPI_Barrier(MPI_COMM_WORLD);

            // DataFullDump();

            //ImposeBoundaryCondition(BOUNDARY_CONDITION_TYPE);

        }else{ 
            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;
        }

        DataFullDump();
    }

    return ;
}
#else 
//if(MPIGetMyID() == MPI_ROOT_RANK){ fprintf(stderr,"%s:%d\n",__func__,__LINE__); fflush(NULL); } MPI_Barrier(MPI_COMM_WORLD);
enum {
    TimeCosmo,
    TimePreDecomp,
    TimeDecomp,
    TimeSort,
    TimeFirstT,
    TimeBuildT,
    TimeNewT,
    TimeRaise,
    TimePredict,
    TimeKick1,
    TimeDrift,
    TimeKick2,
    TimeGravityAll,
    TimeGravKick,
    TimeGravity,
    TimeGravityTree,
    TimeGravityExt,
    TimeHydroAll,
    TimeHydroKernel,
    TimeHydroDensity,
    TimeHydroAcc,
    TimeHydroAccCor,
    TimeHydroTree,
    TimeOthersAll,
    TimeCooling,
    TimeSF,
    TimeFB,
    TimeHII,
    TimeSink,
    TimeIO,
    TimeAll,
    NTime,
};


#define Start(_x) { \
    MPI_Barrier(MPI_COMM_WORLD); \
    Times[_x] = GetElapsedTime();\
    }

#define End(_x) {\
    MPI_Barrier(MPI_COMM_WORLD); \
    Times[_x] = GetElapsedTime()-Times[_x];\
    }

#define Add(_x) {\
    TimesAll[_x] += Times[_x];\
    }


#define ShowResults {\
    MPI_Barrier(MPI_COMM_WORLD); \
    if(MPIGetMyID() == MPI_ROOT_RANK){ \
        fprintf(stderr,"All %g\n",Times[TimeAll]); \
        fprintf(stderr,"Cosm %g, PreDec %g, Dec %g, Sort %g\n",Times[TimeCosmo],Times[TimePreDecomp],Times[TimeDecomp],Times[TimeSort]); \
        fprintf(stderr,"First %g, Build %g, New %g, Rais %g, Pred %g\n",Times[TimeFirstT],Times[TimeBuildT],Times[TimeNewT],Times[TimeRaise]); \
        fprintf(stderr,"Kick1 %g, Kick2 %g, GravKick %g, Pred %g\n",Times[TimeFirstT],Times[TimeBuildT],Times[TimeNewT],Times[TimeRaise]); \
        fprintf(stderr,"GravTree %g, Grav %g, GravExt %g\n",Times[TimeGravityTree],Times[TimeGravity],Times[TimeGravityExt]); \
        fprintf(stderr,"HTree %g, HK %g, HD %g, HA %g, HAC %g\n",Times[TimeHydroTree],Times[TimeHydroKernel],Times[TimeHydroDensity],Times[TimeHydroAcc],Times[TimeHydroAccCor]); \
        fprintf(stderr,"Cooling %g, SF %g, FB %g, HII %g, Sink %g IO %g\n",Times[TimeCooling],Times[TimeSF],Times[TimeFB],Times[TimeHII],Times[TimeSink],Times[TimeIO]); \
        fflush(stderr); \
    }\
    }

void Run(void){

    // Start
    // End(TimeGravityTree)

    double Times[NTime];
    double TimesAll[NTime];
    for(int i=0;i<NTime;i++)
        TimesAll[i] = 0.e0;


    int counter = 0;
    while(Pall.TCurrent < Pall.TEnd){


        for(int i=0;i<NTime;i++)
            Times[i] = 0.e0;

        Start(TimeAll)

        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"%ld : ",Pall.TStepTotal);
            fflush(NULL);
        }

        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();
        
        Start(TimeCosmo)
        UpdateCosmologicalParameters();
        End(TimeCosmo)

        if(Pall.TCurrent >= Pall.Era){
            Start(TimePreDecomp)
            //PreDomainDecomposition();
            PreDomainDecompositionNew(0);
            End(TimePreDecomp)

            Start(TimeDecomp)
            DomainDecompositionOnlyDataExchange();
            End(TimeDecomp)

            Start(TimeSort)
            SortAllStructures();
            End(TimeSort)

            Start(TimeFirstT)
            FirstTimeStep();
            End(TimeFirstT)

            Start(TimeBuildT)
            BuildHierarchicalTimeStep();
            End(TimeBuildT)

            HydroRoot.LifeGauge = 0;
        } else {

            if (10*Pall.NActives_t > Pall.Ntotal_t){
                Start(TimePreDecomp)
                //PreDomainDecomposition();
                PreDomainDecompositionNew(1);
                End(TimePreDecomp)

                Start(TimeDecomp)
                DomainDecompositionOnlyDataExchange();
                End(TimeDecomp)
                HydroRoot.LifeGauge = 0;
            }

            Start(TimeNewT)
            BuildNewTimeStep();
            End(TimeNewT)
        }

        Start(TimeRaise)
        RaiseActiveFlags();
        End(TimeRaise)

        Start(TimeKick1)
        Kick1Drift(); 
        End(TimeKick1)

        Start(TimePredict)
        BuildPredictors();
        End(TimePredict)

#ifdef GRAVITY_RUN //{
        Start(TimeGravityAll)
        if(Pall.NActives_t>0){ // Gravity

            Start(TimeGravityTree)
            PlantGravityTree();
            End(TimeGravityTree)

            ClearGravitationalForce();

            Start(TimeGravity)
            ForceParallelTreeGRAPE();
            End(TimeGravity)

            Start(TimeGravityExt)
            ForceFromExternalPotentials();
            End(TimeGravityExt)

        } else {
            if(MPIGetMyID() == MPI_ROOT_RANK)
                fprintf(stderr,"Skip Gravity\n");
        }
        End(TimeGravityAll)
        Add(TimeGravityAll)
#endif // GRAVITY_RUN //}


#ifdef HYDRO_RUN //{
        Start(TimeHydroAll)
        if(Pall.NActivesHydro_t>0){ // Hydro
            Start(TimeHydroTree)
            PlantHydroTree();
            End(TimeHydroTree)

            ClearHydroData();
            Start(TimeHydroKernel)
#ifdef EVALUATE_SIZES_ALL_TOGETHER //{
            CalcSize();
#else // EVALUATE_SIZES_ALL_TOGETHER //}//{
            CalcKernel();
#endif // EVALUATE_SIZES_ALL_TOGETHER //}
            End(TimeHydroKernel)

            Start(TimeHydroDensity)
            CalcDensityDivRot(); 
            End(TimeHydroDensity)

            Start(TimeHydroAcc)
            CalcDuDtAcc();
            End(TimeHydroAcc)
        }
        End(TimeHydroAll)
        Add(TimeHydroAll)
#endif // HYDRO_RUN //}

        Start(TimeKick2)
        Kick2();
        End(TimeKick2)

        Start(TimeOthersAll)
#ifdef HYDRO_RUN //{
        if(Pall.NActivesHydro_t>0){ // Hydro
            Start(TimeFB)
#ifdef USE_CELIB //{
            StellarFeedback();
#else
            DelayedSNe();
#endif  // USE_CELIB //}
            End(TimeFB)

            Start(TimeCooling)
            CalcCooling();
            End(TimeCooling)

            Start(TimeHII)
            HIIregions();
            End(TimeHII)

            Start(TimeHydroAccCor)
            CalcDuDtAccEnergyDensityForCorrection();
            End(TimeHydroAccCor)

            Start(TimeSF)
            StarFormation();
            End(TimeSF)

            Start(TimeSink)
            SinkParticles();
            End(TimeSink)
        }
#endif // HYDRO_RUN //}
        End(TimeOthersAll)
        Add(TimeOthersAll)

        Start(TimeGravKick)
        UpdateGravityKickFlag();
        End(TimeGravKick)

        /// Post output
        // WriteCurrentActiveParticles(Pall.TStepTotal,"Post");

        // Calc New Time Step or Out Put Logs
        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();
        LogsThisTimeStep();

        Start(TimeIO)
        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            if(Pall.Nhydro_t>0)
                LogStarFormationRate();
            // OutPutAllParticlesInASCIIFormat();

            if(Pall.Nhydro_t>0){
                CountNeighborNumber();
            }

            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;

            // DataDump();
            // FileOutPutConstantInterval();

            // LogOutPutEnergyMomentumAngularMomentum();
            // LogOutPutElapsedTime();
            // LogTotalMass();
            // ReportAllocatedMemorySizes();
            // EndLogs();
            // fflush(NULL);
            // MPI_Barrier(MPI_COMM_WORLD);

            // DataFullDump();

            ImposeBoundaryCondition(BOUNDARY_CONDITION_TYPE);

        }else{ 
            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;
        }

        // DataFullDump();

        End(TimeIO)
        End(TimeAll)
        Add(TimeAll)

        ShowResults;

        if(counter > 30){

            if(MPIGetMyID() == MPI_ROOT_RANK){
                FILE *fp;
                char fname[MaxCharactersInLine];
                Snprintf(fname,"Bench.%04d",MPIGetNumProcs());
                FileOpen(fp,fname,"w");
                fprintf(fp,"%04d %g %g %g %g %g\n",MPIGetNumProcs(),TimesAll[TimeAll],
                TimesAll[TimeGravityAll]+TimesAll[TimeHydroAll]+TimesAll[TimeOthersAll],
                TimesAll[TimeGravityAll],TimesAll[TimeHydroAll],TimesAll[TimeOthersAll]);
                fflush(NULL);
            }


            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
        counter ++;
    }

    return ;
}


#endif

#if 0
enum {
    TimeCosmo,
    TimePreDecomp,
    TimeDecomp,
    TimeSort,
    TimeFirstT,
    TimeBuildT,
    TimeNewT,
    TimeRaise,
    TimePredict,
    TimeKick1,
    TimeDrift,
    TimeKick2,
    TimeGravKick,
    TimeGravity,
    TimeGravityTree,
    TimeGravityExt,
    TimeHydroKernel,
    TimeHydroDensity,
    TimeHydroAcc,
    TimeHydroAccCor,
    TimeHydroTree,
    TimeCooling,
    TimeSF,
    TimeFB,
    TimeHII,
    TimeSink,
    TimeIO,
    TimeAll,
    NTime,
};

#define Start(_x) { \
    MPI_Barrier(MPI_COMM_WORLD); \
    Times[_x] = GetElapsedTime();\
    }

#define End(_x) {\
    MPI_Barrier(MPI_COMM_WORLD); \
    Times[_x] = GetElapsedTime()-Times[_x];\
    }

#define ShowResults {\
    MPI_Barrier(MPI_COMM_WORLD); \
    if(MPIGetMyID() == MPI_ROOT_RANK){ \
        fprintf(stderr,"All %g\n",Times[TimeAll]); \
        fprintf(stderr,"Cosm %g, PreDec %g, Dec %g, Sort %g\n",Times[TimeCosmo],Times[TimePreDecomp],Times[TimeDecomp],Times[TimeSort]); \
        fprintf(stderr,"First %g, Build %g, New %g, Rais %g, Pred %g\n",Times[TimeFirstT],Times[TimeBuildT],Times[TimeNewT],Times[TimeRaise]); \
        fprintf(stderr,"Kick1 %g, Kick2 %g, GravKick %g, Pred %g\n",Times[TimeFirstT],Times[TimeBuildT],Times[TimeNewT],Times[TimeRaise]); \
        fprintf(stderr,"GravTree %g, Grav %g, GravExt %g\n",Times[TimeGravityTree],Times[TimeGravity],Times[TimeGravityExt]); \
        fprintf(stderr,"HTree %g, HK %g, HD %g, HA %g, HAC %g\n",Times[TimeHydroTree],Times[TimeHydroKernel],Times[TimeHydroDensity],Times[TimeHydroAcc],Times[TimeHydroAccCor]); \
        fprintf(stderr,"Cooling %g, SF %g, FB %g, HII %g, Sink %g IO %g\n",Times[TimeCooling],Times[TimeSF],Times[TimeFB],Times[TimeHII],Times[TimeSink],Times[TimeIO]); \
        fflush(stderr); \
    }\
    }

void Run(void){

    // Start
    // End(TimeGravityTree)

    double Times[NTime];
    double ___t;

    int counter = 0;
    while(Pall.TCurrent < Pall.TEnd){

        for(int i=0;i<NTime;i++)
            Times[i] = 0.e0;

        Start(TimeAll)

        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"%ld : ",Pall.TStepTotal);
            fflush(NULL);
        }

        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();
        
        Start(TimeCosmo)
        UpdateCosmologicalParameters();
        End(TimeCosmo)

        if(Pall.TCurrent >= Pall.Era){
            Start(TimePreDecomp)
            PreDomainDecomposition();
            End(TimePreDecomp)

            Start(TimeDecomp)
            DomainDecompositionOnlyDataExchange();
            End(TimeDecomp)

            Start(TimeSort)
            SortAllStructures();
            End(TimeSort)

            Start(TimeFirstT)
            FirstTimeStep();
            End(TimeFirstT)

            Start(TimeBuildT)
            BuildHierarchicalTimeStep();
            End(TimeBuildT)

            HydroRoot.LifeGauge = 0;
        } else {

            if (10*Pall.NActives_t > Pall.Ntotal_t){
                Start(TimePreDecomp)
                PreDomainDecomposition();
                End(TimePreDecomp)

                Start(TimeDecomp)
                DomainDecompositionOnlyDataExchange();
                End(TimeDecomp)
                HydroRoot.LifeGauge = 0;
            }

            Start(TimeNewT)
            BuildNewTimeStep();
            End(TimeNewT)
        }

        Start(TimeRaise)
        RaiseActiveFlags();
        End(TimeRaise)

        Start(TimeKick1)
        Kick1Drift(); 
        End(TimeKick1)

        Start(TimePredict)
        BuildPredictors();
        End(TimePredict)

        /// Pre output
        // WriteCurrentActiveParticles(Pall.TStepTotal,"Pre");


#ifdef GRAVITY_RUN //{
        if(Pall.NActives_t>0){ // Gravity


            Start(TimeGravityTree)
            PlantGravityTree();
            End(TimeGravityTree)

            ClearGravitationalForce();

            Start(TimeGravity)
            ForceParallelTreeGRAPE();
            End(TimeGravity)

            Start(TimeGravityExt)
            ForceFromExternalPotentials();
            End(TimeGravityExt)

        } else {
            if(MPIGetMyID() == MPI_ROOT_RANK)
                fprintf(stderr,"Skip Gravity\n");
        }
#endif // GRAVITY_RUN //}


#ifdef HYDRO_RUN //{
        if(Pall.NActivesHydro_t>0){ // Hydro
            Start(TimeHydroTree)
            PlantHydroTree();
            End(TimeHydroTree)

            ClearHydroData();
            Start(TimeHydroKernel)
            CalcKernel();
            End(TimeHydroKernel)

            Start(TimeHydroDensity)
            CalcDensityDivRot(); 
            End(TimeHydroDensity)

            Start(TimeHydroAcc)
            CalcDuDtAcc();
            End(TimeHydroAcc)
        }
#endif // HYDRO_RUN //}

        Start(TimeKick2)
        Kick2();
        End(TimeKick2)

#ifdef HYDRO_RUN //{
        if(Pall.NActivesHydro_t>0){ // Hydro
            Start(TimeFB)
#ifdef USE_CELIB //{
            StellarFeedback();
#else
            DelayedSNe();
#endif  // USE_CELIB //}
            End(TimeFB)

            Start(TimeCooling)
            CalcCooling();
            End(TimeCooling)

            Start(TimeHII)
            HIIregions();
            End(TimeHII)

            Start(TimeHydroAccCor)
            CalcDuDtAccEnergyDensityForCorrection();
            End(TimeHydroAccCor)

            Start(TimeSF)
            StarFormation();
            End(TimeSF)

            Start(TimeSink)
            SinkParticles();
            End(TimeSink)
        }
#endif // HYDRO_RUN //}

        Start(TimeGravKick)
        UpdateGravityKickFlag();
        End(TimeGravKick)

        /// Post output
        // WriteCurrentActiveParticles(Pall.TStepTotal,"Post");

        // Calc New Time Step or Out Put Logs
        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();
        //LogsThisTimeStep();

        Start(TimeIO)
        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            if(Pall.Nhydro_t>0)
                LogStarFormationRate();
            OutPutAllParticlesInASCIIFormat();

            if(Pall.Nhydro_t>0){
                CountNeighborNumber();
            }

            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;

            // DataDump();
            FileOutPutConstantInterval();

            LogOutPutEnergyMomentumAngularMomentum();
            LogOutPutElapsedTime();
            LogTotalMass();
            ReportAllocatedMemorySizes();
            EndLogs();
            fflush(NULL);
            MPI_Barrier(MPI_COMM_WORLD);

            // DataFullDump();

            ImposeBoundaryCondition(BOUNDARY_CONDITION_TYPE);

        }else{ 
            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;
        }

        DataFullDump();

        End(TimeIO)
        End(TimeAll)

        ShowResults;

        if(counter > 10){
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
        counter ++;
    }

    return ;
}
#endif 

#ifdef USE_SPSPH //{
static void ResetInitialZ(void){
    
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->Zw = Phydro[i]->ZwPred = PhydroBody(i)->Mass/Phydro[i]->Rho;
        Phydro[i]->PseudoDensity = Phydro[i]->PseudoDensityPred = 1.e0;
    }

    CalcDensityDivRot();

    return ;
}
#endif // USE_SPSPH //}

#if VISCOSITY_TYPE == 1 //{
void InitCopyHydroPredictors(void){

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->RhoPred = Phydro[i]->Rho;
        Phydro[i]->KernelPred = Phydro[i]->Kernel;
#ifdef USE_DISPH //{
        Phydro[i]->EnergyDensityPred = Phydro[i]->EnergyDensity;
#endif //USE_DISPH //}
#ifdef USE_SPSPH //{
        Phydro[i]->PseudoDensityPred;
        Phydro[i]->ZwPred = Phydro[i]->Zw;
#endif // USE_SPSPH //}
    }
    return ;
}
#endif // VISCOSITY_TYPE == 1 //}

//int NNN;
void InitializeRun(void){

    //Kick1Drift();
    BuildPredictors(); // Pos -> PosP/ Vel -> VelP
    SetInitialHydroVelh();
    

    InitializeDecomposition();
    DomainDecomposition();

    // First force / hydro calculation.
    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

#ifdef GRAVITY_RUN
    PlantGravityTree();
    ClearGravitationalForce();
    ForceParallelTreeGRAPE();
    ForceFromExternalPotentials();
#endif // GRAVITY_RUN

    if(Pall.Nhydro_t > 0){
        InitializeRootForHydro();
#ifdef COOLING_RUN //{
        InitializeCoolingTable();
#endif // COOLING_RUN //}
        InitializeFarUltraVioletField();
#ifdef USE_THERMAL_CONDUCTIVITY //{ 
        InitThermalConductivity();
#endif // USE_THERMAL_CONDUCTIVITY //}

        InitializeStarFormation();
#ifdef USE_CELIB
        InitializeStellarFeedback();
#endif
#ifdef DELAYED_FEEDBACK
        InitializeDelayedSNII();
#endif

        ClearHydroData();
        PlantHydroTree();
#ifdef EVALUATE_SIZES_ALL_TOGETHER //{
#if VISCOSITY_TYPE == 1 //{  
        CalcKernel();
        CalcDensityDivRot();
        InitCopyHydroPredictors();
#endif // VISCOSITY_TYPE == 1 //}
        CalcSize();
#else // EVALUATE_SIZES_ALL_TOGETHER //}//{
        CalcKernel();
#endif // EVALUATE_SIZES_ALL_TOGETHER //}
        //CalcKernel();
        CalcDensityDivRot();
#ifdef USE_SPSPH //{
        ResetInitialZ();
#endif // USE_SPSPH //}
        CalcDuDtAcc();
    }

#ifdef USE_SYMMETRIZED_SOFTENING
    //CalcSymmetrizedPotential();
#endif
    LogOutPutEnergyMomentumAngularMomentum();

    FirstTimeStep();
    if(Pall.Nhydro_t > 0){
        PlantHydroTree();
#ifdef HYDRO_TIMESTEP_LIMITER 
        for(int i=0;i<Pall.Nhydro;i++){
            Phydro[i]->k_hydro_localmin = MaximumTimeHierarchy;
            Phydro[i]->NextUpdateEra = Pall.TEnd;
        }
        BuildHierarchicalTimeStep();
        CalcDensityDivRot();
        for(int i=0;i<Pall.Nhydro;i++){
            Phydro[i]->k_hydro = MIN(Phydro[i]->k_hydro,Phydro[i]->k_hydro_localmin+MAX_K_LOCAL);
            Phydro[i]->dt_hydro = Pall.dtmin*exp2(Phydro[i]->k_hydro);
#ifndef USE_FAST_SCHEME 
            PhydroBody(i)->k = Phydro[i]->k_hydro;
            PhydroBody(i)->dt = Phydro[i]->dt_hydro;
#else
            if(Phydro[i]->k_hydro > PhydroBody(i)->k){
                Phydro[i]->k_hydro = PhydroBody(i)->k;
                Phydro[i]->dt_hydro = PhydroBody(i)->dt;
            }
#endif
        }
        for(int i=0;i<Pall.Nhydro;i++){
            Phydro[i]->k_hydro_localmin = MaximumTimeHierarchy;
            Phydro[i]->NextUpdateEra = Pall.TEnd;
        }
#endif
    }


    UpdateGravityKickFlag();

    FileOutPutConstantInterval();
    // Pall.OutPutFileNumber ++;

#ifdef WRITE_KERNEL_SHAPE //{
    WriteKernelShape();
#endif //WRITE_KERNEL_SHAPE //}

    return ;
}

void RestartRun(void){

    // Decomposition 
    InitializeDecomposition();
    // PreDomainDecomposition();

    // First force / hydro calculation.
#ifdef GRAVITY_RUN
    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

#endif // GRAVITY_RUN

    if(Pall.Nhydro_t > 0){
        InitializeRootForHydro();
        InitializeCoolingTable();
        InitializeFarUltraVioletField();
        InitializeStarFormation();
#ifdef USE_CELIB
        InitializeStellarFeedback();
#else
        InitializeDelayedSNII();
#endif
    }

    LogOutPutEnergyMomentumAngularMomentum();
    UpdateGravityKickFlag();

    ImposeBoundaryCondition(BOUNDARY_CONDITION_TYPE);

    return ;
}


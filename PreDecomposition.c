#include "config.h"
#include "Decomposition.h"

#define EnsembleMeanDecomposition   (ON)
#define InteractionListBase         (OFF)
#define ParticleNumberBase          (OFF)
#define Hybrid                      (ON)
#define BaseInteractionLength   (1.e+0)

struct StructDecomposition{
    int     Ntotal;     
    int     Allocated;     
    int     *InteractionList; 
    double  (*Pos)[3]; 
} InfoDecomposition;

struct StructPreDecompostion{
    double  Pos[3];
    int     InteractionList;
} *PreDecompositionSend,*PreDecompositionRecv;

struct StructPreDecompositionQsort{
    int InteractionList;
    double Pos[3];
};

static void BiSectionQsort(void);
static void BiSectionQsortFreeDirection(void);
static void BiSectionQsortFreeDirectionNoUpdate(void);
static void BiSectionQsortInteractionList(void);
static void BiSectionQsortHybrid(void);
static int ParticleSampling(void);
static int ParticleSampling_(const int NSample);
static int ParticleSamplingEnsemble(const int Element);
static int ParticleSamplingPositionOnly(const int Element, double Pos[restrict][3]);
static int ParticleSamplingWithInteractionList(const int Element, double Pos[restrict][3], int InteractionList[]);

static void WriteDomainBisectionInfo(char *name);
static void DomainOutPut(char *name);
static void DomainSamplesOutPut(char *name);

void InitializeDecomposition(void){

    //int NProcs = MPIGetNumProcs();

#if (EnsembleMeanDecomposition)
    int NSample = MAX(DecompNSample,MPIGetNumProcs()*PREDECOMPOSITION_SAMPLING_NUMBER);
    InfoDecomposition.Allocated = NSample;
    //InfoDecomposition.Allocated = DecompNSample;
    InfoDecomposition.Pos = malloc(sizeof(double)*3*InfoDecomposition.Allocated);
    InfoDecomposition.InteractionList = malloc(sizeof(int)*InfoDecomposition.Allocated);
#else
    if(MPIGetMyID() == MPI_ROOT_RANK){
        InfoDecomposition.Allocated = NProcs*DecompNSample;
        InfoDecomposition.Pos = malloc(sizeof(double)*3*InfoDecomposition.Allocated);
        InfoDecomposition.InteractionList = malloc(sizeof(int)*InfoDecomposition.Allocated);
    }else{
        InfoDecomposition.Allocated = (int)(ForAngelsShare*(DecompNSample/NProcs));
        InfoDecomposition.Pos = malloc(sizeof(double)*3*InfoDecomposition.Allocated);
        InfoDecomposition.InteractionList = malloc(sizeof(int)*InfoDecomposition.Allocated);
    }
#endif

    int SizeBisection = (1<<MPIGetNumProcsPower())-1;
    //dprintlmpi(SizeBisection);
    InfoBiSection = malloc(sizeof(struct StructInfoBiSection)*(SizeBisection+1));

    return;
}

#if (EnsembleMeanDecomposition)
void PreDomainDecomposition(void){

    double TimingResultThisRoutine = GetElapsedTime();

    if(MPIGetNumProcs() == 1) 
        return;
#ifdef PREDECOMPOSITION_INTERVAL_STEPS
    if((Pall.TStepTotal)%PREDECOMPOSITION_INTERVAL_STEPS!=0) 
        return;
#endif

    int NProcs = MPIGetNumProcs();
    //char fname[MaxCharactersInLine];
    MPI_Status mpi_status;

    //int NSample = DecompNSample;
    int NSample = MAX(DecompNSample,MPIGetNumProcs()*PREDECOMPOSITION_SAMPLING_NUMBER);
    int Element = (int)(((double)NSample/(double)Pall.Ntotal_t)*Pall.Ntotal);
    while(Element*MPIGetNumProcs() > Pall.Ntotal_t)
        Element --;
    //MPI_Barrier(MPI_COMM_WORLD);
    int NElements[NProcs],NElementsAll; 
    NElements[NProcs-1] = Element;

    NElementsAll = Element;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(&Element,1,MPI_INT,CommunicationTable[i].SendRank,TAG_PREDECOMPOSITION_PRECOMM,
            NElements+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_PREDECOMPOSITION_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NElementsAll += NElements[i];
    }
    assert(NElementsAll <= NSample);
    InfoDecomposition.Ntotal = NElementsAll;

    int NElementsMax = NElements[0];
    for(int i=0;i<NProcs-1;i++)
        NElementsMax = MAX(NElementsMax,NElements[i]);

    if(InfoDecomposition.Allocated < NElementsAll){
        InfoDecomposition.Allocated = (int)(ForAngelsShare*NElementsAll);
        free(InfoDecomposition.Pos);
        free(InfoDecomposition.InteractionList);
        InfoDecomposition.Pos = malloc(sizeof(double)*3*InfoDecomposition.Allocated);
        InfoDecomposition.InteractionList = malloc(sizeof(int)*InfoDecomposition.Allocated);
    }

    PreDecompositionSend = malloc(sizeof(struct StructPreDecompostion)*Element);
    PreDecompositionRecv = malloc(sizeof(struct StructPreDecompostion)*NElementsMax+1);

    int NSampleThisTime = ParticleSamplingEnsemble(Element);
    assert(NSampleThisTime == Element);
    for(int k=0;k<Element;k++){
        InfoDecomposition.Pos[k][0] = PreDecompositionSend[k].Pos[0];
        InfoDecomposition.Pos[k][1] = PreDecompositionSend[k].Pos[1];
        InfoDecomposition.Pos[k][2] = PreDecompositionSend[k].Pos[2];
        InfoDecomposition.InteractionList[k] = PreDecompositionSend[k].InteractionList+1;
    }

    int Noffset = Element;
    for(int i=0;i<NProcs-1;i++){
        NSampleThisTime = ParticleSamplingEnsemble(Element);
        assert(NSampleThisTime==Element);
        MPI_Sendrecv(PreDecompositionSend,NSampleThisTime*sizeof(struct StructPreDecompostion),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_PREDECOMPOSITION_PRECOMM,
                     PreDecompositionRecv,NElements[i]*sizeof(struct StructPreDecompostion),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_PREDECOMPOSITION_PRECOMM,
                    MPI_COMM_WORLD,&mpi_status);
        for(int k=0;k<NElements[i];k++){
            InfoDecomposition.Pos[k+Noffset][0] = PreDecompositionRecv[k].Pos[0];
            InfoDecomposition.Pos[k+Noffset][1] = PreDecompositionRecv[k].Pos[1];
            InfoDecomposition.Pos[k+Noffset][2] = PreDecompositionRecv[k].Pos[2];
            InfoDecomposition.InteractionList[k+Noffset] = PreDecompositionRecv[k].InteractionList+1;
        }
        Noffset += NElements[i];
    }
    free(PreDecompositionSend);
    free(PreDecompositionRecv);

#if (InteractionListBase)
    BiSectionQsortInteractionList();
#endif
#if (Hybrid)
    BiSectionQsortHybrid();
#endif
#if (ParticleNumberBase)
    BiSectionQsort();
#endif

    // here ensemble mean
    int NDomain = (1<<MPIGetNumProcsPower())-1;
    double Pos[NDomain],EnsemblePos[NDomain];
    for(int i=0;i<NDomain;i++)
        Pos[i] = InfoBiSection[i].Pos;
    MPI_Allreduce(Pos,EnsemblePos,NDomain,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    double iNProcs = 1.0/NProcs;
    for(int i=0;i<NDomain;i++)
        InfoBiSection[i].Pos = EnsemblePos[i]*iNProcs;

    /*
    char fname[MaxCharactersInLine];
    sprintf(fname,"Domain.%02d",MPIGetMyID());
    DomainOutPut(fname);
    sprintf(fname,"DomainSample.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    DomainSamplesOutPut(fname);
    */

    // char fname[MaxCharactersInLine];
    // sprintf(fname,"Bisection.%02d",MPIGetMyID());
    // WriteDomainBisectionInfo(fname);

    TimingResults.DecompositionThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
}
#else
void PreDomainDecomposition(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    char fname[MaxCharactersInLine];
    MPI_Status mpi_status;

    int NSample = DecompNSample;
    int Element = (((double)NSample/(double)Pall.Ntotal_t)*Pall.Ntotal);

    if(InfoDecomposition.Allocated < Element){
        while(InfoDecomposition.Allocated < Element)
            InfoDecomposition.Allocated = (int)(ForAngelsShare*InfoDecomposition.Allocated);
        Delete(InfoDecomposition.Pos);
        Delete(InfoDecomposition.InteractionList);
        InfoDecomposition.Pos = malloc(sizeof(double)*3*InfoDecomposition.Allocated);
        InfoDecomposition.InteractionList = malloc(sizeof(int)*InfoDecomposition.Allocated);
    }

    int NSampleThisTime = ParticleSampling_(Element);
    //int NSampleThisTime = ParticleSampling();

    int NSampleRecv[NProcs];
    MPI_Gather(&NSampleThisTime,1,MPI_INT,NSampleRecv,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    //for(int i=0;i<NProcs;i++)
        //fprintf(stderr,"[%02d] %d\n",MPIGetMyID(),NSampleRecv[i]);

    int NSampleRecvAll = 0;
    for(int i=0;i<NProcs;i++)
        NSampleRecvAll += NSampleRecv[i];

    if(MyID == MPI_ROOT_RANK){
        assert(NSampleRecvAll <= NSample);
        while(InfoDecomposition.Allocated < NSampleRecvAll)
            InfoDecomposition.Allocated = (int)(ForAngelsShare*InfoDecomposition.Allocated);
        InfoDecomposition.Pos = realloc(InfoDecomposition.Pos,sizeof(double)*3*InfoDecomposition.Allocated);
        InfoDecomposition.InteractionList = 
            realloc(InfoDecomposition.InteractionList,sizeof(int)*InfoDecomposition.Allocated);
        assert(NSampleRecvAll <= InfoDecomposition.Allocated);
    }


    if(MyID == MPI_ROOT_RANK){
        MPI_Request mpi_request_pos[NProcs-1],mpi_request_list[NProcs-1];
        int NOffset = NSampleRecv[0];
        for(int i=1;i<NProcs;i++){
            MPI_Irecv(InfoDecomposition.Pos+NOffset,3*NSampleRecv[i],MPI_DOUBLE,i,
                TAG_PREDECOMPOSITION_POS,MPI_COMM_WORLD,mpi_request_pos+(i-1));
            MPI_Irecv(InfoDecomposition.InteractionList+NOffset,NSampleRecv[i],MPI_INT,i,
                TAG_PREDECOMPOSITION_LIST,MPI_COMM_WORLD,mpi_request_list+(i-1));
            NOffset += NSampleRecv[i];
        }
        for(int i=1;i<NProcs;i++){
            MPI_Wait(mpi_request_pos+(i-1),&mpi_status);
            MPI_Wait(mpi_request_list+(i-1),&mpi_status);
        }
        InfoDecomposition.Ntotal = NSampleRecvAll;
    }else{
        MPI_Request mpi_request_pos,mpi_request_list;
        MPI_Isend(InfoDecomposition.Pos,3*NSampleThisTime,MPI_DOUBLE,MPI_ROOT_RANK,
                TAG_PREDECOMPOSITION_POS,MPI_COMM_WORLD,&mpi_request_pos);
        MPI_Isend(InfoDecomposition.InteractionList,NSampleThisTime,MPI_INT,MPI_ROOT_RANK,
                TAG_PREDECOMPOSITION_LIST,MPI_COMM_WORLD,&mpi_request_list);
        MPI_Wait(&mpi_request_pos,&mpi_status);
        MPI_Wait(&mpi_request_list,&mpi_status);
    }

#if (InteractionListBase)
    if(MyID == MPI_ROOT_RANK)
        BiSectionQsortInteractionList();
#else
    if(MyID == MPI_ROOT_RANK)
        BiSectionQsort();
#endif

    int SizeBisection = (1<<MPIGetNumProcsPower())-1;
    MPI_Bcast(InfoBiSection,sizeof(struct StructInfoBiSection)*SizeBisection,MPI_BYTE,MPI_ROOT_RANK,MPI_COMM_WORLD);

    /*
    sprintf(fname,"Domain.%02d",MyID);
    DomainOutPut(fname);
    sprintf(fname,"DomainSample.%02d.%02d",NProcs,MyID);
    DomainSamplesOutPut(fname);
    */

    return;
}
#endif


void PreDomainDecompositionNew(const int mode){

    double TimingResultThisRoutine = GetElapsedTime();

    if(MPIGetNumProcs() == 1) 
        return;

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    MPI_Status mpi_status;

    const double rate = ((double)NProcs*PREDECOMPOSITION_SAMPLING_NUMBER)/(double)Pall.Ntotal_t;
    unsigned long int LocalSumIntList,GlobalSumIntList;
    LocalSumIntList = 0;
    for(int i=0;i<Pall.Ntotal;i++)
        LocalSumIntList += Pbody[i]->InteractionList+1;
    MPI_Allreduce(&LocalSumIntList,&GlobalSumIntList,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);

#define PREDECOMPOSITION_MINIMUM_SAMPLING_CONTROL_FACTOR (0.1)

    int NSample = MAX(Pall.Ntotal_t*rate*(((double)LocalSumIntList)/((double)GlobalSumIntList)),
                    Pall.Ntotal*rate/(1+PREDECOMPOSITION_MINIMUM_SAMPLING_CONTROL_FACTOR));
    NSample = MIN(NSample,Pall.Ntotal*rate/(1-PREDECOMPOSITION_MINIMUM_SAMPLING_CONTROL_FACTOR));
    int NElements[NProcs]; 
    memset(NElements,0,sizeof(int)*NProcs);
    //NElements[NProcs] = NSample;

    MPI_Gather(&NSample,1,MPI_INT,NElements,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);

    int NElementsAll = 0;
    for(int i=0;i<NProcs;i++){
        NElementsAll += NElements[i];
    }

    InfoDecomposition.Ntotal = NElementsAll;

    int RecvPosition[NProcs];
    int RecvDataSize[NProcs];
    int RecvPositionInteraction[NProcs];
    int RecvDataSizeInteraction[NProcs];

    if(MyID == MPI_ROOT_RANK){
        if(InfoDecomposition.Allocated < NElementsAll){
            InfoDecomposition.Allocated = (int)(ForAngelsShare*NElementsAll);
            free(InfoDecomposition.Pos);
            free(InfoDecomposition.InteractionList);
            InfoDecomposition.Pos = malloc(sizeof(double)*3*InfoDecomposition.Allocated);
            InfoDecomposition.InteractionList = malloc(sizeof(int)*InfoDecomposition.Allocated);
        }
        // count displacement
        RecvPosition[0] = 0;
        RecvDataSize[0] = 3*NElements[0];
        RecvPositionInteraction[0] = 0;
        RecvDataSizeInteraction[0] = NElements[0];
        for(int i=1;i<NProcs;i++){
            RecvPosition[i] = RecvPosition[i-1]+RecvDataSize[i-1];
            RecvDataSize[i] = 3*NElements[i];
            RecvPositionInteraction[i] = RecvPositionInteraction[i-1]+RecvDataSizeInteraction[i-1];
            RecvDataSizeInteraction[i] = NElements[i];
        }
    } else {
        if(InfoDecomposition.Allocated < NSample){
            InfoDecomposition.Allocated = (int)(ForAngelsShare*NSample);
            free(InfoDecomposition.Pos);
            free(InfoDecomposition.InteractionList);
            InfoDecomposition.Pos = malloc(sizeof(double)*3*InfoDecomposition.Allocated);
            InfoDecomposition.InteractionList = malloc(sizeof(int)*InfoDecomposition.Allocated);
        }
    }


///


    double PreDecompositionSendBufferPos[NSample][3];
    int PreDecompositionSendBufferInteractionList[NSample];

#ifdef USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //{
    int NSampleThisTime = ParticleSamplingWithInteractionList(NSample,PreDecompositionSendBufferPos,PreDecompositionSendBufferInteractionList);
#else //USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //}//{
    int NSampleThisTime = ParticleSamplingPositionOnly(NSample,PreDecompositionSendBufferPos);
#endif //USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //}
    assert(NSample == NSampleThisTime);

    MPI_Gatherv(PreDecompositionSendBufferPos,NSample*3,MPI_DOUBLE,
                InfoDecomposition.Pos,RecvDataSize,RecvPosition,MPI_DOUBLE,MPI_ROOT_RANK,MPI_COMM_WORLD);
#ifdef USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //{
    MPI_Gatherv(PreDecompositionSendBufferInteractionList,NSample,MPI_INT,
                InfoDecomposition.InteractionList,RecvDataSizeInteraction,RecvPositionInteraction,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
#endif //USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //}

#if 0
    {
    if(MPIGetMyID() == MPI_ROOT_RANK){
    FILE *fp;
    FileOpen(fp,"pred.dat","w");

    for(int i=0;i<InfoDecomposition.Ntotal;i++){
        fprintf(fp,"%g %g %g %d\n",InfoDecomposition.Pos[i][0],
        InfoDecomposition.Pos[i][1],InfoDecomposition.Pos[i][2],
        InfoDecomposition.InteractionList[i]);
    }
    fclose(fp);
    fflush(NULL);
    }
    }
#endif

    // Calc boundaries
    if(MyID == MPI_ROOT_RANK){
        if(mode == 0){
            BiSectionQsortFreeDirection();
            //BiSectionQsort();
        } else {
            BiSectionQsortFreeDirectionNoUpdate();
        }
    }

    int NDomain = (1<<MPIGetNumProcsPower())-1;
    MPI_Bcast(InfoBiSection,sizeof(struct StructInfoBiSection)*NDomain,MPI_BYTE,0,MPI_COMM_WORLD);

    TimingResults.DecompositionThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
}

static int ParticleSampling(void){

    int step = MAX(1,Pall.Ntotal/DecompNSample);

    int nsample = 0;  
    for(int i=0;i<Pall.Ntotal;i+=step){
        InfoDecomposition.Pos[nsample][0] = Pbody[i]->PosP[0];
        InfoDecomposition.Pos[nsample][1] = Pbody[i]->PosP[1];
        InfoDecomposition.Pos[nsample][2] = Pbody[i]->PosP[2];
        InfoDecomposition.InteractionList[nsample] = Pbody[i]->InteractionList;
        nsample ++;
        if(nsample == DecompNSample)
            break;
    }
    return (nsample); 
}

static int ParticleSamplingEnsemble(const int Element){

    int nsample = 0;
    while(nsample < Element){
        int index = (int)(Pall.Ntotal*gsl_rng_uniform(RandomGenerator));
        PreDecompositionSend[nsample].Pos[0] = Pbody[index]->PosP[0];
        PreDecompositionSend[nsample].Pos[1] = Pbody[index]->PosP[1];
        PreDecompositionSend[nsample].Pos[2] = Pbody[index]->PosP[2];
        //PreDecompositionSend[nsample].InteractionList = 1;
        PreDecompositionSend[nsample].InteractionList = Pbody[index]->InteractionList;
        nsample ++;
    }
    return (nsample); 
}

static int ParticleSamplingPositionOnly(const int Element, double Pos[restrict][3]){

    int nsample = 0;
    while(nsample < Element){
        int index = (int)(Pall.Ntotal*gsl_rng_uniform(RandomGenerator));
        Pos[nsample][0] = Pbody[index]->PosP[0];
        Pos[nsample][1] = Pbody[index]->PosP[1];
        Pos[nsample][2] = Pbody[index]->PosP[2];
        nsample ++;
    }
    return (nsample); 
}

static int ParticleSamplingWithInteractionList(const int Element, double Pos[restrict][3], int InteractionList[]){

    int nsample = 0;
    while(nsample < Element){
        int index = (int)(Pall.Ntotal*gsl_rng_uniform(RandomGenerator));
        Pos[nsample][0] = Pbody[index]->PosP[0];
        Pos[nsample][1] = Pbody[index]->PosP[1];
        Pos[nsample][2] = Pbody[index]->PosP[2];
        InteractionList[nsample] = sqrt(Pbody[index]->InteractionList+1);
        nsample ++;
    }
    return (nsample); 
}

static int ParticleSampling_(const int NSample){

    int nsample = 0;
    while(nsample < NSample){
        int index = (int)(Pall.Ntotal*gsl_rng_uniform(RandomGenerator));
        InfoDecomposition.Pos[nsample][0] = Pbody[index]->PosP[0];
        InfoDecomposition.Pos[nsample][1] = Pbody[index]->PosP[1];
        InfoDecomposition.Pos[nsample][2] = Pbody[index]->PosP[2];
        InfoDecomposition.InteractionList[nsample] = Pbody[index]->InteractionList;
        nsample ++;
    }
    return (nsample); 
}

static int PosCmpX(const void *x, const void *y){
    const struct StructPreDecompositionQsort *pointer1 = x;
    const struct StructPreDecompositionQsort *pointer2 = y;
    if( pointer1->Pos[0] > pointer2->Pos[0])
        return 1;
    else if( pointer1->Pos[0] < pointer2->Pos[0])
        return -1;
    else
        return 0;
    //return (pointer1->Pos[0] > pointer2->Pos[0] ? 1 : pointer1->Pos[0] > pointer2->Pos[0] ?  -1 : 0);
}
static int PosCmpY(const void *x, const void *y){
    const struct StructPreDecompositionQsort *pointer1 = x;
    const struct StructPreDecompositionQsort *pointer2 = y;
    if( pointer1->Pos[1] > pointer2->Pos[1])
        return 1;
    else if( pointer1->Pos[1] < pointer2->Pos[1])
        return -1;
    else
        return 0;
    //return (pointer1->Pos[1] > pointer2->Pos[1] ? 1 : pointer1->Pos[1] > pointer2->Pos[1] ?  -1 : 0);
}
static int PosCmpZ(const void *x, const void *y){
    const struct StructPreDecompositionQsort *pointer1 = x;
    const struct StructPreDecompositionQsort *pointer2 = y;
    if( pointer1->Pos[2] > pointer2->Pos[2])
        return 1;
    else if( pointer1->Pos[2] < pointer2->Pos[2])
        return -1;
    else
        return 0;
    //return (x->Pos[2] > y->Pos[2] ? 1 : x->Pos[2] > y->Pos[2] ?  -1 : 0);
}

static void BiSectionQsort(void){

    int Power = MPIGetNumProcsPower();

    struct StructPreDecompositionQsort *PreDecompositionQsort;
    PreDecompositionQsort = malloc(sizeof(struct StructPreDecompositionQsort)*InfoDecomposition.Ntotal);

    for(int i=0;i<InfoDecomposition.Ntotal;i++){
        PreDecompositionQsort[i].Pos[0] = InfoDecomposition.Pos[i][0];
        PreDecompositionQsort[i].Pos[1] = InfoDecomposition.Pos[i][1];
        PreDecompositionQsort[i].Pos[2] = InfoDecomposition.Pos[i][2];
    }

    int HeaderDelimiter = 0;

    int Nmax = 1<<Power;
    int Header[Power][Nmax];
    int Size[Power][Nmax];

    for(int i=0;i<Power;i++){
        int Axis = i%BisectionDimension;
        int Division = (1<<i);

        int HeaderDelimiterCopy = HeaderDelimiter;
        if(i==0){
            Header[i][0] = 0;
            Size[i][0] = InfoDecomposition.Ntotal;
        } else {
            for(int k=0;k<(Division>>1);k++){
                Header[i][2*k] = Header[i-1][k];
                Size[i][2*k] = Size[i-1][k]/2;
                Header[i][2*k+1] = Header[i-1][k]+Size[i][2*k];
                Size[i][2*k+1] = Size[i-1][k]-Size[i][2*k];

                InfoBiSection[HeaderDelimiter-(1<<(i-1))+k].Left = HeaderDelimiterCopy;
                InfoBiSection[HeaderDelimiter-(1<<(i-1))+k].Right = HeaderDelimiterCopy+1;

                HeaderDelimiterCopy += 2;
            }
        }

        for(int k=0;k<Division;k++){
            switch(Axis){
                case 0:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpX);
                    break;
                case 1:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpY);
                    break;
                case 2:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpZ);
                    break;
            }
            InfoBiSection[HeaderDelimiter].Pos
                = PreDecompositionQsort[Header[i][k]+Size[i][k]/2].Pos[Axis];
            InfoBiSection[HeaderDelimiter].Axis = Axis;
            HeaderDelimiter ++;
        }
    }

    free(PreDecompositionQsort);

    return;
}

//static int ReturnMaximumExtendedDirection(const int NParticles, const double Pos[][3]){
static int ReturnMaximumExtendedDirection(const int NParticles, struct StructPreDecompositionQsort PreDecomp[]){

    double Edges[2][3] = {{PreDecomp[0].Pos[0],PreDecomp[0].Pos[1],PreDecomp[0].Pos[2]},
                          {PreDecomp[0].Pos[0],PreDecomp[0].Pos[1],PreDecomp[0].Pos[2]}};

    for(int i=0;i<NParticles;i++){
        for(int k=0;k<3;k++){
            Edges[0][k] = fmin(PreDecomp[i].Pos[k],Edges[0][k]);
            Edges[1][k] = fmax(PreDecomp[i].Pos[k],Edges[1][k]);
        }
    }
    double Width[3] = {Edges[1][0]-Edges[0][0],
                       Edges[1][1]-Edges[0][1],
                       Edges[1][2]-Edges[0][2]};

    // int first = Width[0]>Width[1]?0:1;
    // int second = Width[2]>Widht[first]?2:1;
    // fprintf(stderr,"Width %g %g %g\n",Width[0],Width[1],Width[2]);

    //int second = Width[2]>Widht[Width[0]>Width[1]?0:1]?2:first;
    return Width[2]>Width[Width[0]>Width[1]?0:1]?2:Width[0]>Width[1]?0:1;
}


static void BiSectionQsortFreeDirection(void){

    int Power = MPIGetNumProcsPower();

    struct StructPreDecompositionQsort *PreDecompositionQsort;
    PreDecompositionQsort = malloc(sizeof(struct StructPreDecompositionQsort)*InfoDecomposition.Ntotal);

    for(int i=0;i<InfoDecomposition.Ntotal;i++){
        PreDecompositionQsort[i].Pos[0] = InfoDecomposition.Pos[i][0];
        PreDecompositionQsort[i].Pos[1] = InfoDecomposition.Pos[i][1];
        PreDecompositionQsort[i].Pos[2] = InfoDecomposition.Pos[i][2];
#ifdef USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //{
        PreDecompositionQsort[i].InteractionList = InfoDecomposition.InteractionList[i];
#endif // USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //}
    }

    int HeaderDelimiter = 0;

    int Nmax = 1<<Power;
    int Header[Power][Nmax];
    int Size[Power][Nmax];

    for(int i=0;i<Power;i++){
        int Division = (1<<i);

        int HeaderDelimiterCopy = HeaderDelimiter;
        if(i==0){
            Header[i][0] = 0;
            Size[i][0] = InfoDecomposition.Ntotal;
        } else {
            for(int k=0;k<(Division>>1);k++){
                Header[i][2*k] = Header[i-1][k];
                Size[i][2*k] = Size[i-1][k]/2;
                Header[i][2*k+1] = Header[i-1][k]+Size[i][2*k];
                Size[i][2*k+1] = Size[i-1][k]-Size[i][2*k];

                InfoBiSection[HeaderDelimiter-(1<<(i-1))+k].Left = HeaderDelimiterCopy;
                InfoBiSection[HeaderDelimiter-(1<<(i-1))+k].Right = HeaderDelimiterCopy+1;

                HeaderDelimiterCopy += 2;
            }
        }

        for(int k=0;k<Division;k++){
            int Axis = ReturnMaximumExtendedDirection(Size[i][k],PreDecompositionQsort+Header[i][k]);
            switch(Axis){
                case 0:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpX);
                    break;
                case 1:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpY);
                    break;
                case 2:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpZ);
                    break;
            }

            // dprintlmpi(Axis);
#ifdef USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //{
            long long int TotalInteractionList = 0;
            for(int l=0;l<Size[i][k];l++){
                TotalInteractionList += (long long int)(PreDecompositionQsort[Header[i][k]+l].InteractionList);
            }
            int CountInteractionList = 0;
            long long int SumInteractionList = 0;
            long long int HalfInteractionList = TotalInteractionList/2;
            // dprintlmpi(Size[i][k]);
            // dprintlmpi(HalfInteractionList);
            // dprintlmpi(TotalInteractionList);
            while(SumInteractionList<HalfInteractionList){
                SumInteractionList += 
                    (long long int)(PreDecompositionQsort[Header[i][k]+CountInteractionList].InteractionList);
                CountInteractionList ++;
            }
            InfoBiSection[HeaderDelimiter].Pos
                = PreDecompositionQsort[Header[i][k]+CountInteractionList].Pos[Axis];
#else //USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //}//{
            InfoBiSection[HeaderDelimiter].Pos
                = PreDecompositionQsort[Header[i][k]+Size[i][k]/2].Pos[Axis];
#endif //USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //}
            InfoBiSection[HeaderDelimiter].Axis = Axis;
            HeaderDelimiter ++;
        }
    }

    free(PreDecompositionQsort);

    return;
}

static void BiSectionQsortFreeDirectionNoUpdate(void){

    int Power = MPIGetNumProcsPower();

    struct StructPreDecompositionQsort *PreDecompositionQsort;
    PreDecompositionQsort = malloc(sizeof(struct StructPreDecompositionQsort)*InfoDecomposition.Ntotal);

    for(int i=0;i<InfoDecomposition.Ntotal;i++){
        PreDecompositionQsort[i].Pos[0] = InfoDecomposition.Pos[i][0];
        PreDecompositionQsort[i].Pos[1] = InfoDecomposition.Pos[i][1];
        PreDecompositionQsort[i].Pos[2] = InfoDecomposition.Pos[i][2];
#ifdef USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //{
        PreDecompositionQsort[i].InteractionList = InfoDecomposition.InteractionList[i];
#endif // USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //}
    }

    int HeaderDelimiter = 0;

    int Nmax = 1<<Power;
    int Header[Power][Nmax];
    int Size[Power][Nmax];

    for(int i=0;i<Power;i++){
        int Division = (1<<i);

        int HeaderDelimiterCopy = HeaderDelimiter;
        if(i==0){
            Header[i][0] = 0;
            Size[i][0] = InfoDecomposition.Ntotal;
        } else {
            for(int k=0;k<(Division>>1);k++){
                Header[i][2*k] = Header[i-1][k];
                Size[i][2*k] = Size[i-1][k]/2;
                Header[i][2*k+1] = Header[i-1][k]+Size[i][2*k];
                Size[i][2*k+1] = Size[i-1][k]-Size[i][2*k];

                InfoBiSection[HeaderDelimiter-(1<<(i-1))+k].Left = HeaderDelimiterCopy;
                InfoBiSection[HeaderDelimiter-(1<<(i-1))+k].Right = HeaderDelimiterCopy+1;

                HeaderDelimiterCopy += 2;
            }
        }

        for(int k=0;k<Division;k++){
            int Axis = InfoBiSection[HeaderDelimiter].Axis;
            switch(Axis){
                case 0:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpX);
                    break;
                case 1:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpY);
                    break;
                case 2:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpZ);
                    break;
            }
#ifdef USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //{
            long long int TotalInteractionList = 0;
            for(int l=0;l<Size[i][k];l++){
                TotalInteractionList += (long long int)(PreDecompositionQsort[Header[i][k]+l].InteractionList);
            }
            int CountInteractionList = 0;
            long long int SumInteractionList = 0;
            long long int HalfInteractionList = TotalInteractionList/2;
            // dprintlmpi(Size[i][k]);
            // dprintlmpi(HalfInteractionList);
            // dprintlmpi(TotalInteractionList);
            while(SumInteractionList<HalfInteractionList){
                SumInteractionList += 
                    (long long int)(PreDecompositionQsort[Header[i][k]+CountInteractionList].InteractionList);
                CountInteractionList ++;
            }
            InfoBiSection[HeaderDelimiter].Pos
                = PreDecompositionQsort[Header[i][k]+CountInteractionList].Pos[Axis];
#else //USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //}//{
            //InfoBiSection[HeaderDelimiter].Pos
                //= PreDecompositionQsort[Header[i][k]+Size[i][k]/2].Pos[Axis];
            InfoBiSection[HeaderDelimiter].Pos
                = PreDecompositionQsort[Header[i][k]+Size[i][k]/2].Pos[Axis];
#endif //USE_INTERACTIONLIST_FOR_PREDECOMPOSITION_SAMPLING //}
            HeaderDelimiter ++;
        }
    }

    free(PreDecompositionQsort);

    return;
}

static void BiSectionQsortInteractionList(void){

    static int Flag = 0;
    int Power = MPIGetNumProcsPower();

    struct StructPreDecompositionQsort *PreDecompositionQsort;
    PreDecompositionQsort = malloc(sizeof(struct StructPreDecompositionQsort)*InfoDecomposition.Ntotal);

    for(int i=0;i<InfoDecomposition.Ntotal;i++){
        PreDecompositionQsort[i].Pos[0] = InfoDecomposition.Pos[i][0];
        PreDecompositionQsort[i].Pos[1] = InfoDecomposition.Pos[i][1];
        PreDecompositionQsort[i].Pos[2] = InfoDecomposition.Pos[i][2];
        PreDecompositionQsort[i].InteractionList = InfoDecomposition.InteractionList[i];
    }

#if 0
    {
    FILE *fp;
    FileOpen(fp,"PreDecomposition.data","w");
    for(int i=0;i<InfoDecomposition.Ntotal;i++){
        fprintf(fp,"%g %g %g %d\n",PreDecompositionQsort[i].Pos[0],
            PreDecompositionQsort[i].Pos[1],PreDecompositionQsort[i].Pos[2],
                PreDecompositionQsort[i].InteractionList);
    }
    fclose(fp);
    }
#endif

    int HeaderDelimiter = 0;

    int Nmax = 1<<Power;
    int Header[Power][Nmax];
    int Size[Power][Nmax];

    for(int i=0;i<Power;i++){
        int Axis = i%BisectionDimension;
        int Division = (1<<i);

        int HeaderDelimiterCopy = HeaderDelimiter;
        if(i==0){
            Header[i][0] = 0;
            Size[i][0] = InfoDecomposition.Ntotal;
        } else {
            for(int k=0;k<(Division>>1);k++){
                Header[i][2*k] = Header[i-1][k];
                Size[i][2*k] = Size[i-1][k]/2;
                Header[i][2*k+1] = Header[i-1][k]+Size[i][2*k];
                Size[i][2*k+1] = Size[i-1][k]-Size[i][2*k];

                InfoBiSection[HeaderDelimiter-(1<<(i-1))+k].Left = HeaderDelimiterCopy;
                InfoBiSection[HeaderDelimiter-(1<<(i-1))+k].Right = HeaderDelimiterCopy+1;

                HeaderDelimiterCopy += 2;
            }
        }

        for(int k=0;k<Division;k++){
            switch(Axis){
                case 0:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpX);
                    break;
                case 1:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpY);
                    break;
                case 2:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpZ);
                    break;
            }


            long long int TotalInteractionList = 0;
            for(int l=0;l<Size[i][k];l++)
                TotalInteractionList += (long long int)(PreDecompositionQsort[Header[i][k]+l].InteractionList);
            //dlprintlmpi(TotalInteractionList);

            int CountInteractionList = 0;
            long long int SumInteractionList = 0;
            long long int HalfInteractionList = TotalInteractionList/2;
            while(SumInteractionList<HalfInteractionList){
                SumInteractionList += 
                    (long long int)(PreDecompositionQsort[Header[i][k]+CountInteractionList].InteractionList);
                CountInteractionList ++;
            }
            //dprintlmpi(CountInteractionList);
            //dlprintlmpi(SumInteractionList);

#if 0
            if(Flag == 0){
                InfoBiSection[HeaderDelimiter].Pos
                    = PreDecompositionQsort[Header[i][k]+CountInteractionList].Pos[Axis];
            }else{
                InfoBiSection[HeaderDelimiter].Pos
                    = 0.5*(InfoBiSection[HeaderDelimiter].Pos
                        +PreDecompositionQsort[Header[i][k]+CountInteractionList].Pos[Axis]);
            }
#else
            InfoBiSection[HeaderDelimiter].Pos
                = PreDecompositionQsort[Header[i][k]+CountInteractionList].Pos[Axis];
#endif
            fprintf(stderr,"[%02d] PreDecomp[%02d] %g\n",
                    MPIGetMyID(),HeaderDelimiter,InfoBiSection[HeaderDelimiter].Pos);
            HeaderDelimiter ++;
        }
    }
    assert(HeaderDelimiter==((1<<MPIGetNumProcsPower())-1));

    free(PreDecompositionQsort);
    Flag ++;

    return;
}

static void BiSectionQsortHybrid(void){

    //static int Flag = 0;
    int Power = MPIGetNumProcsPower();

    struct StructPreDecompositionQsort *PreDecompositionQsort;
    PreDecompositionQsort = malloc(sizeof(struct StructPreDecompositionQsort)*InfoDecomposition.Ntotal);

    for(int i=0;i<InfoDecomposition.Ntotal;i++){
        PreDecompositionQsort[i].Pos[0] = InfoDecomposition.Pos[i][0];
        PreDecompositionQsort[i].Pos[1] = InfoDecomposition.Pos[i][1];
        PreDecompositionQsort[i].Pos[2] = InfoDecomposition.Pos[i][2];
        PreDecompositionQsort[i].InteractionList = InfoDecomposition.InteractionList[i];
    }

    int HeaderDelimiter = 0;

    int Nmax = 1<<Power;
    int Header[Power+1][Nmax*2];
    int Size[Power+1][Nmax*2];

    Header[0][0] = 0;
    Size[0][0] = InfoDecomposition.Ntotal;
    for(int i=0;i<Power;i++){
        int Axis = i%BisectionDimension;
        int Division = (1<<i);

        int HeaderDelimiterCopy = HeaderDelimiter;
        /*
        if(i==0){
            Header[i][0] = 0;
            Size[i][0] = InfoDecomposition.Ntotal;
        } else {
            for(int k=0;k<(Division>>1);k++){
                Header[i][2*k] = Header[i-1][k];
                Size[i][2*k] = Size[i-1][k]/2;
                Header[i][2*k+1] = Header[i-1][k]+Size[i][2*k];
                Size[i][2*k+1] = Size[i-1][k]-Size[i][2*k];

                InfoBiSection[HeaderDelimiter-(1<<(i-1))+k].Left = HeaderDelimiterCopy;
                InfoBiSection[HeaderDelimiter-(1<<(i-1))+k].Right = HeaderDelimiterCopy+1;

                HeaderDelimiterCopy += 2;
            }
        }
        */
        if(i != 0){
            for(int k=0;k<(Division>>1);k++){
                InfoBiSection[HeaderDelimiter-(1<<(i-1))+k].Left = HeaderDelimiterCopy;
                InfoBiSection[HeaderDelimiter-(1<<(i-1))+k].Right = HeaderDelimiterCopy+1;
                HeaderDelimiterCopy += 2;
            }
        }

        for(int k=0;k<Division;k++){
            //fprintf(stderr,"[%02d] Power = %d, Division,Header,Size = %d %d %d, Axis  %d\n",
                    //MPIGetMyID(),i,k,Header[i][k],Size[i][k],Axis);
            switch(Axis){
                case 0:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpX);
                    break;
                case 1:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpY);
                    break;
                case 2:
                    qsort(PreDecompositionQsort+Header[i][k], Size[i][k],
                        sizeof(struct StructPreDecompositionQsort),
                            (int(*)(const void*, const void*))PosCmpZ);
                    break;
            }


            long long int TotalInteractionList = 0;
            for(int l=0;l<Size[i][k];l++)
                TotalInteractionList += 
                    (long long int)(1+PreDecompositionQsort[Header[i][k]+l].InteractionList/BaseInteractionLength);

            int CountInteractionList = 0;
            long long int SumInteractionList = 0;
            long long int HalfInteractionList = TotalInteractionList/2;
            while(SumInteractionList<HalfInteractionList){
                SumInteractionList += 
                    (long long int)(1+PreDecompositionQsort[Header[i][k]+CountInteractionList].InteractionList/BaseInteractionLength);
                CountInteractionList ++;
            }

            Size[i+1][2*k] = CountInteractionList;
            Size[i+1][2*k+1] = Size[i][k]-CountInteractionList;
            Header[i+1][2*k] = Header[i][k];
            Header[i+1][2*k+1] = Header[i][k]+CountInteractionList;

#if 0
            if(Flag == 0){
                InfoBiSection[HeaderDelimiter].Pos
                    = PreDecompositionQsort[Header[i][k]+CountInteractionList].Pos[Axis];
            }else{
                InfoBiSection[HeaderDelimiter].Pos
                    = 0.5*(InfoBiSection[HeaderDelimiter].Pos
                        +PreDecompositionQsort[Header[i][k]+CountInteractionList].Pos[Axis]);
            }
#else
            InfoBiSection[HeaderDelimiter].Pos
                = PreDecompositionQsort[Header[i][k]+CountInteractionList].Pos[Axis];
#endif
            //fprintf(stderr,"[%02d] PreDecomp[%02d] %g\n",
                    //MPIGetMyID(),HeaderDelimiter,InfoBiSection[HeaderDelimiter].Pos);
            HeaderDelimiter ++;
        }
    }
    assert(HeaderDelimiter == ((1<<MPIGetNumProcsPower())-1));

    free(PreDecompositionQsort);
    //Flag ++;

    return;
}

static void WriteDomainBisectionInfo(char *name){

    int SizeBisection = (1<<MPIGetNumProcsPower())-1;
    FILE *fp;
    FileOpen(fp,name,"w");
    for(int i=0;i<SizeBisection;i++)
        fprintf(fp,"%d %g %d %d\n",MPIGetMyID(),
                InfoBiSection[i].Pos,InfoBiSection[i].Right,InfoBiSection[i].Left);
    fclose(fp);
    return;
}

static void DomainOutPut(char *name){

    FILE *fp;

    FileOpen(fp,name,"w");
    for(int i=0;i<Pall.Ntotal;i++)
        fprintf(fp,"%ld %g %g %g\n",Pbody[i]->GlobalID,Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2]);
    fclose(fp);
    return;
}

static void DomainSamplesOutPut(char *name){

    FILE *fp;

    FileOpen(fp,name,"w");
    for(int i=0;i<InfoDecomposition.Ntotal;i++)
        fprintf(fp,"%g %g %g %d\n",
            InfoDecomposition.Pos[i][0],
                InfoDecomposition.Pos[i][1],
                    InfoDecomposition.Pos[i][2],
                        InfoDecomposition.InteractionList[i]);
    fclose(fp);

    return;
}


void InspectParticleFluctuation(void){

    int Nlocal = Pall.Ntotal;
    int NArray[MPIGetNumProcs()];

    MPI_Gather(&Nlocal,1,MPI_INT,NArray,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
 
    if(MPIGetMyID() == MPI_ROOT_RANK){
        int NProcs = MPIGetNumProcs();
        int Total = 0;
        int Max = NArray[0]; 
        int Min = NArray[0]; 
        for(int i=0;i<NProcs;i++){
            Total += NArray[i];
            Max = MAX(Max,NArray[i]);
            Min = MIN(Min,NArray[i]);
        }
        double Mean = (double)Total/(double)NProcs;

        double ArrayMu[MPIGetNumProcs()];
        double ArraySigma2[MPIGetNumProcs()];
        for(int i=0;i<NProcs;i++){
            ArrayMu[i] = NArray[i]-Mean;
            ArraySigma2[i] = SQ(ArrayMu[i]);
        }
        double Sigma2 = 0.e0;
        for(int i=0;i<NProcs;i++){
            Sigma2 += ArraySigma2[i];
        }
        Sigma2 /= (double)NProcs;

        fprintf(stderr,"Particle distribution infomation\n");
        fprintf(stderr," Mean = %g, Min = %d, Max = %d\n",Mean,Min,Max);
        fprintf(stderr," Sigma2 = %g, Sigma = %g\n",Sigma2,sqrt(Sigma2));
        fprintf(stderr," 1 Sigma range : %g < %g < %g\n",
                Mean-sqrt(Sigma2),Mean,Mean+sqrt(Sigma2));

        fflush(NULL);
    }

    return ;
}


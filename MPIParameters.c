#include "config.h"
#include "MPIParameters.h"

static int numprocs,myid,namelen,numgrapes,numprocspower;
static char processor_name[MPI_MAX_PROCESSOR_NAME];

void InitializeMPIEnv(int *argc, char **argv){

    MPI_Init(argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);

    MPISetMyID(myid);
    MPISetNumProcs(numprocs);
    MPISetNumProcsPower(numprocs);
    MPISetNumGrapes(MIN(MPIGetNumProcs(),4));
    MPISetNameLen(namelen);
    MPISetProcessorName(processor_name);

    return ;
}

void FinalizeMPIEnv(void){

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return ;
}


void MPISetMyID(int current_id){
    myid = current_id;
    return;
}

void MPISetNumProcs(int current_numprocs){
    numprocs = current_numprocs;
    return;
}

void MPISetNumGrapes(int current_numgrapes){
    numgrapes = current_numgrapes;
    return;
}

void MPISetNameLen(int current_namelen){
    namelen = current_namelen;
    return;
}

void MPISetProcessorName(char *current_processor_name){
    strcpy(processor_name,current_processor_name);
    return;
}

void MPISetNumProcsPower(int current_numprocs){

    int ndummy = 1;

    numprocspower = 0;
    if(current_numprocs%2 != 0){
        numprocspower = 0;
        return;
    }

    while(ndummy != current_numprocs){
        ndummy *= 2;
        numprocspower ++;
    }

    return;
}

int MPIGetMyID(void){
    return myid;
}

int MPIGetNumProcs(void){
    return numprocs;
}

int MPIGetNumGrapes(void){
    return numgrapes;
}

int MPIGetNameLen(void){
    return namelen;
}

char *MPIGetProcessorName(void){
    return processor_name;
}

int MPIGetNumProcsPower(void){
    return numprocspower;
}

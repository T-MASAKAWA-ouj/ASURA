#include "config.h"

void InitializeCommunicationOrder(void);

void InitializeCommunicationTable(void){

    int Nprocs = MPIGetNumProcs();
    CommunicationTable = malloc(sizeof(struct StructCommunicationTable)*Nprocs);

    InitializeCommunicationOrder();

    return;
}

void InitializeCommunicationOrder(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    int SendRecvRank[NProcs][NProcs];
    for(int i=0;i<NProcs;i++)
        SendRecvRank[0][i] = i;

    int levelmax = (int)log2(NProcs);
    for(int level=0;level<levelmax;level++){
        int remainder,division,header;
        remainder = MyID%(1<<level);
        division = MyID/(1<<level);

        header = (1<<level)*division + Parity(division)*(1<<level);
        for(int k=0;k<NProcs;k++){
            header = k+Parity(k/(1<<level))*(1<<level);
            for(int i=0;i<(1<<level);i++){
                SendRecvRank[(1<<level)+i][k] = SendRecvRank[i][header];
            }
        }
    }

    for(int i=0;i<NProcs-1;i++){
        CommunicationTable[i].SendRank = SendRecvRank[i+1][MyID];
        CommunicationTable[i].RecvRank = SendRecvRank[i+1][MyID];
    }

    return;
}

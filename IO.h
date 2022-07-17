/* IO.c */
#ifndef __IO_H_INCLUDED__
#define __IO_H_INCLUDED__
void GetRunStatus(const int argc, char **argv);
int GetRunStatusForSingleNodeAnalysis(const int argc, char **argv);
void CheckIOFileName(void);
void WriteAllData(void);
void ReadAllData(void);
void DataDump(void);
void DataFullDump(void);
void FileOutPutConstantInterval(void);
void BinaryDump(void);
void ParallelWriteAllData(char fname[]);
void ParallelDumpAllData(char fname[]);
void ParallelWriteAllDataASCIIFormat(char fname[]);
void ParallelReadAllData(void);
void ReadParallelDataOnSingleNode(void);
void OutPutASCIIDATA(void);
void WriteHydroDataASCIIFormat(const int mode);
void ShowHydroDataASCIIFormat(const int mode);
void WriteNeighborInformation(const int index, const int mode);
void ShowNeighborInformation(const int index, const int mode);
void WriteTimeSteps(void);
void WriteCurrentActiveParticles(const long int FileID, char suffix[]);

#endif // __IO_H_INCLUDED__

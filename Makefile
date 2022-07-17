#
# Makefile for ASURA.
#

UNAME = $(USER)
SYSTEM = $(HOSTNAME)

ifeq ($(SYSTEM),localhost.localdomain)
include Makefile.yutaka
endif
ifeq ($(SYSTEM),dhcp-219-191.mtk.nao.ac.jp)
include Makefile.yutaka
endif
ifeq ($(SYSTEM),cotponic)
include Makefile.cotponic
endif
ifeq ($(SYSTEM),sarrych)
include Makefile.sarrych
endif
ifeq ($(SYSTEM),calyc)
include Makefile.cotponic
endif
ifeq ($(SYSTEM),xc00)
include Makefile.xc
endif
ifeq ($(SYSTEM),xc01)
include Makefile.xc
endif
ifeq ($(SYSTEM),xc02)
include Makefile.xc
endif
ifeq ($(SYSTEM),xc03)
include Makefile.xc
endif
ifeq ($(SYSTEM),xc04)
include Makefile.xc
endif
ifeq ($(SYSTEM),an00.cfca.nao.ac.jp)
include Makefile.an00
endif
ifeq ($(SYSTEM),an01.cfca.nao.ac.jp)
include Makefile.an00
endif
ifeq ($(SYSTEM),an02.cfca.nao.ac.jp)
include Makefile.an00
endif
ifeq ($(SYSTEM),an03.cfca.nao.ac.jp)
include Makefile.an03
endif
ifeq ($(SYSTEM),an04.cfca.nao.ac.jp)
include Makefile.an03
endif
ifeq ($(SYSTEM),an05.cfca.nao.ac.jp)
include Makefile.an03
endif
ifeq ($(SYSTEM),an06.cfca.nao.ac.jp)
include Makefile.an03
endif
ifeq ($(SYSTEM),spaurh.mtk.nao.ac.jp)
include Makefile.spaurh
endif
ifeq ($(SYSTEM),pfs00)
include Makefile.pfs
endif
ifeq ($(SYSTEM),pfs01)
include Makefile.pfs
endif
ifeq ($(SYSTEM),spaurh00)
include Makefile.spaurh
endif
ifeq ($(SYSTEM),spaurh01)
include Makefile.spaurh
endif
include makeflags


CFLAGS = $(OPTIMIZE) $(OPTIONS) $(INCLUDEPATH) 

LDFLAGS += -lm


MOBJS	=	main.o 

OBJS =	DataStructures.o \
    StructureOperation.o \
    Run.o \
    Decomposition.o \
    PreDecomposition.o \
    PeriodicWrapping.o \
    ForceMisc.o \
    ForceDirect.o \
    ForceGRAPE.o \
    ForceTree.o \
    ForceParallelTreeGRAPE.o \
    ForceFromExternalPotentials.o \
    GRAPEEmulator.o \
    GRAPEManager.o \
    NeighborSearch.o \
    PlantGravityTree.o \
    PlantHydroTree.o \
    PlantHydroTreeImported.o \
    PlantStellarTree.o \
    TimeStep.o \
    Integral.o \
    CosmologicalTime.o \
    Cosmological.o \
    SizeDetermination.o \
    HydroMisc.o \
    HydroDensity.o \
    HydroAcc.o \
    HydroKernel.o \
    HydroExtraOperations.o \
    Cooling.o \
    Heating.o \
    H2.o \
    StarFormation.o \
    IMF.o \
    SNII.o \
    Delayed.o \
	StellarFeedback.o \
	ThermalConductivity.o \
    HIIregion.o \
    SinkParticle.o \
    SetUpTestRun.o \
    TestRuns.o \
    Read.o \
    IO.o \
    FOF.o \
    Unit.o \
    Logs.o \
    CommunicationTable.o \
    BufferOperation.o \
    OrderingKey.o \
	MPIParameters.o \
	ParallelOperation.o \
	RunLogs.o \
	StructuresForIO.o \
	Utilities.o \
	ElapsedTime.o \
	Allocators.o \
	ParticleAddRemoveOperation.o \
	RandomNumberGenerator.o \
	CheckStructures.o \
	BoundaryCondition.o \
	Exit.o \


MSRCS = main.c

SRCS =	DataStructures.c \
    StructureOperation.c \
    Run.c \
    Decomposition.c \
    PreDecomposition.c \
    PeriodicWrapping.c \
    ForceMisc.c \
    ForceDirect.c \
    ForceGRAPE.c \
    ForceTree.c \
    ForceParallelTreeGRAPE.c \
    ForceFromExternalPotentials.c \
    GRAPEEmulator.c \
    GRAPEManager.c \
    NeighborSearch.c \
    PlantGravityTree.c \
    PlantHydroTree.c \
    PlantHydroTreeImported.c \
    PlantStellarTree.c \
    TimeStep.c \
    Integral.c \
    CosmologicalTime.c \
    Cosmological.c \
    SizeDetermination.c \
    HydroMisc.c \
    HydroDensity.c \
    HydroAcc.c \
    HydroKernel.c \
    HydroExtraOperations.c \
    Cooling.c \
    Heating.c \
    H2.c \
    StarFormation.c \
    IMF.c \
    SNII.c \
    Delayed.c \
	StellarFeedback.c \
	ThermalConductivity.c \
    HIIregion.c \
    SinkParticle.c \
    SetUpTestRun.c \
    TestRuns.c \
    IO.c \
    Read.c \
    FOF.c \
    Unit.c \
    Logs.c \
    CommunicationTable.c \
    BufferOperation.c \
    OrderingKey.c \
	MPIParameters.c \
	ParallelOperation.c \
	RunLogs.c \
	StructuresForIO.c \
	Utilities.c \
	ElapsedTime.c \
	Allocators.c  \
	ParticleAddRemoveOperation.c \
	RandomNumberGenerator.c \
	CheckStructures.c \
	BoundaryCondition.c \
	Exit.c \

CHEADS = config.h \
	reconfig.h \
	Astro.h \
	DataStructures.h \
	Constants.h \
    Unit.h \
    Logs.h \
	MPIParameters.h \
	MPITags.h \
    BufferOperation.h \
	ParallelOperation.h \
	RunLogs.h \
	Utilities.h \
	ElapsedTime.h \
    HydroKernel.h \
    CommunicationTable.h \
	Allocators.h 

HEADS = $(CHEADS) \
    StructureOperation.h \
    Run.h \
    Decomposition.h \
    PreDecomposition.h \
    PeriodicWrapping.h \
    ForceMisc.h \
    ForceDirect.h \
    ForceGRAPE.h \
    ForceTree.h \
    ForceParallelTreeGRAPE.h \
    ForceFromExternalPotentials.h \
    GRAPEEmulator.h \
    GRAPEManager.h \
    NeighborSearch.h \
    PlantGravityTree.h \
    PlantHydroTree.h \
    PlantHydroTreeImported.h \
    PlantStellarTree.h \
    TimeStep.h \
    Integral.h \
    CosmologicalTime.h \
    Cosmological.h \
    SizeDetermination.h \
    HydroMisc.h \
    HydroDensity.h \
    HydroAcc.h \
    HydroExtraOperations.h \
    Cooling.h \
    CoolingTable.h \
    Heating.h \
    H2.h \
    StarFormation.h \
    IMF.h \
    SNII.h \
    Delayed.h \
	StellarFeedback.h \
	ThermalConductivity.h \
    HIIregion.h \
    SinkParticle.h \
    SetUpTestRun.h \
    TestRuns.h \
    IO.h \
    Read.h \
    FOF.h \
    OrderingKey.h \
	ParticleAddRemoveOperation.h \
	RandomNumberGenerator.h \
	CheckStructures.h \
	BoundaryCondition.h \
	Exit.h \

ALL:   start

starte:$(MOBJS) $(OBJS) 
	$(CC) $(CFLAGS) $(MOBJS) $(OBJS) $(LDFLAGS) -o asura.out

start:
	@echo '';\
	./AutoLog.rb makeflags;\
	./AutoMakeIO.rb ;\
	echo '=============================================================================';\
	echo ' make ASURA';\
	echo ' OPTIMIZATION PARAMETERS = ';\
	echo '   $(OPTIMIZE)';\
	echo ' OPTIONS = ';\
	echo '   $(OPTIONS)';\
	echo ' INCLUDEPATH = ';\
	echo '   $(INCLUDEPATH)';\
	echo ' LIBRARIES =';\
	echo '   $(LDFLAGS)';\
	echo ' Compile mode = $(GRAPE)';\
	echo '=============================================================================';\
	echo '';\
	echo '';\
	make starte ;\
	echo '';\
	echo '=============================================================================';\
	echo ' make ASURA';\
	echo ' OPTIMIZATION PARAMETERS = ';\
	echo '   $(OPTIMIZE)';\
	echo ' OPTIONS = ';\
	echo '   $(OPTIONS)';\
	echo ' INCLUDEPATH = ';\
	echo '   $(INCLUDEPATH)';\
	echo ' LIBRARIES =';\
	echo '   $(LDFLAGS)';\
	echo ' Compile mode = $(GRAPE)';\
	echo ' System = $(SYSTEM)';\
	echo '=============================================================================';\
	echo '';\
	echo ''


tag:
	ctags *.[c,h] AutoMakeIO.rb makeflags

clean:
	rm $(OBJS) main.o asura.out
#	rm $(OBJS) main.o asura.out asurapkg.tar.gz

#DATE = $(shell date '+%Y%m%d')
DATE = $(shell date '+%Y%m%d%H%M')
pkg:
	echo '$(DATE)' 
	tar cvf asurapkg.$(DATE).tar ./*.[ch] ./*.rb ./CoolingCurveSpaans2008.dat Makefile Makefile.* makeflags 
	gzip asurapkg.$(DATE).tar

strip:
	strip asura.out

main.o: main.c Makefile makeflags $(HEADS) 
Run.o: Run.c Makefile makeflags $(HEADS)

DataStructures.o: DataStructures.c Makefile makeflags $(HEADS)
StructureOperation.o: StructureOperation.c Makefile makeflags $(HEADS)

Decomposition.o: Decomposition.c Makefile makeflags $(HEADS)
PreDecomposition.o: PreDecomposition.c Makefile makeflags $(HEADS)
PeriodicWrapping.o: PeriodicWrapping.c Makefile makeflags $(HEADS)

ForceMisc.o: Makefile makeflags $(HEADS)
ForceDirect.o: Makefile makeflags $(HEADS)
ForceGRAPE.o: Makefile makeflags $(HEADS)
ForceTree.o: Makefile makeflags $(HEADS)
ForceParallelTreeGRAPE.o: Makefile makeflags $(HEADS) 
ForceFromExternalPotentials.o: Makefile makeflags $(HEADS)
GRAPEEmulator.o: Makefile makeflags $(HEADS)
GRAPEManager.o: Makefile makeflags $(HEADS)
NeighborSearch.o: Makefile makeflags $(HEADS)

PlantGravityTree.o: Makefile makeflags $(HEADS)
PlantHydroTree.o: Makefile makeflags $(HEADS)
PlantHydroTreeImported.o: Makefile makeflags $(HEADS)
PlantStellarTree.o: Makefile makeflags $(HEADS)

TimeStep.o: Makefile makeflags $(HEADS)
Integral.o: Makefile makeflags $(HEADS)
CosmologicalTime.o: Makefile makeflags $(HEADS)
Cosmological.o: Cosmological.c Makefile makeflags $(HEADS)

SizeDetermination.o: SizeDetermination.c Makefile makeflags $(HEADS)
HydroMisc.o: HydroMisc.c Makefile makeflags $(HEADS)
HydroDensity.o: HydroDensity.c Makefile makeflags $(HEADS)
HydroAcc.o: HydroAcc.c Makefile makeflags $(HEADS)
HydroKernel.o: HydroKernel.c Makefile makeflags $(HEADS)
HydroExtraOperations.o: HydroExtraOperations.c Makefile makeflags $(HEADS)
Cooling.o: Cooling.c CoolingTable.h CoolingTableSpaans2008.h CoolingCurveSpaans2008.dat Makefile makeflags $(HEADS)
Heating.o: Heating.c CoolingTableSpaans2008.h CoolingCurveSpaans2008.dat Makefile makeflags $(HEADS)
H2.o: H2.c CoolingTableSpaans2008.h CoolingCurveSpaans2008.dat Makefile makeflags $(HEADS)
StarFormation.o: StarFormation.c Makefile makeflags $(HEADS)
IMF.o: IMF.c Makefile makeflags $(HEADS)
SNII.o: SNII.c Makefile makeflags $(HEADS)
Delayed.o: Delayed.c Makefile makeflags $(HEADS)
StellarFeedback.o: StellarFeedback.c Makefile makeflags $(HEADS)
ThermalConductivity.o: ThermalConductivity.c Makefile makeflags $(HEADS)
HIIregion.o: HIIregion.c Makefile makeflags $(HEADS)
SinkParticle.o: SinkParticle.c Makefile makeflags $(HEADS)

SetUpTestRun.o: SetUpTestRun.c Makefile makeflags $(HEADS)
IO.o: IO.c StructuresForIO.c StructuresForIO.h AutoMakeIO.rb Makefile makeflags $(HEADS)
Read.o: Read.c Makefile makeflags $(HEADS)
BoundaryCondition.o: BoundaryCondition.c Makefile makeflags $(HEADS)
TestRuns.o: TestRuns.c Makefile makeflags $(HEADS)

CommunicationTable.o: CommunicationTable.c Makefile makeflags $(CHEADS)
BufferOperation.o: BufferOperation.c Makefile makeflags $(CHEADS) 
Allocators.o: Allocators.c Makefile makeflags $(CHEADS) 
MPIParameters.o: MPIParameters.c Makefile makeflags $(CHEADS) 
ParallelOperation.o: ParallelOperation.c Makefile makeflags $(CHEADS) 
ElapsedTime.o: ElapsedTime.c Makefile makeflags $(CHEADS) 
Unit.o: Makefile makeflags $(HEADS) 
Logs.o: Makefile makeflags $(HEADS) 
ParticleAddRemoveOperation.o: ParticleAddRemoveOperation.c Makefile makeflags $(CHEADS)
RandomNumberGenerator.o: RandomNumberGenerator.c Makefile makeflags $(CHEADS)
RunLogs.o: RunLogs.c AutoLog.rb Makefile makeflags $(CHEADS)
StructuresForIO.o: StructuresForIO.c AutoMakeIO.rb Makefile makeflags $(CHEADS)
Utilities.o: Utilities.c Makefile makeflags $(CHEADS)
CheckStructures.o: CheckStructures.c Makefile makeflags $(CHEADS)
Exit.o: Exit.c Makefile makeflags $(CHEADS)

FOF.o: FOF.c Makefile makeflags $(CHEADS)

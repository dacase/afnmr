SHELL		= /bin/sh

MEAD_OBJECTS =					\
	AnalyMaker.o				\
	AnalySlab.o				\
	AnalySphere.o				\
	Angle.o					\
	Atom.o					\
	AtomChargeSet.o				\
	AtomID.o				\
	AtomSet.o				\
	AvsSF_read_header.o			\
	AvsScalarField.o			\
	AvsScalarField_input.o			\
	Bigmem.o				\
	CLShell.o				\
	ChargeCubeRep.o				\
	ChargeDist.o				\
	CubeLatSpec.o				\
	Debye.o					\
	DielByAtoms.o				\
	DielCubeRep.o				\
	DielMembAtoms.o				\
	DielectricEnvironment.o			\
	DielectricSlab.o			\
	DielectricSphere.o			\
	ElectrolyteByAtoms.o			\
	ElectrolyteEnvironment.o		\
	ElstatMaker.o				\
	ElstatPot.o				\
	ElstatPotCombination.o			\
	Elstat_list_init.o			\
	ElySphere.o				\
	FDChargeIterator.o			\
	FDElstatMaker.o				\
	FDGridLevel.o				\
	FinDiffElstatPot.o			\
	FinDiffMethod.o				\
	ManyPointCharge.o			\
	MomentAnalysis.o                        \
	OnePointCharge.o			\
	Pair.o					\
	PairGenBorn.o				\
	PhysCond.o				\
	Polynomial.o				\
	Potat.o					\
	SAVanal_calc.o				\
	SAVcalc_cuberep.o			\
	SAVtag_points.o				\
	Sausage.o				\
	SphericalHarmonic.o                     \
	Shell.o					\
	SolvAccVol.o				\
	TPShell.o				\
	UniformDielectric.o			\
	UniformElectrolyte.o			\
	globals.o				\
	multPotxCharge.o			\
	sor.o

include ../config.h

install: solinprot

libmead.a: $(MEAD_OBJECTS)
	rm -f t$@
	$(AR) t$@ $(MEAD_OBJECTS) 
	$(RANLIB) t$@
	rm -f $@ && mv t$@ $@

solinprot: libmead.a
	$(CXX) -I.. $(COPTFLAGS) -DUSE_EXCEPTIONS -o $(BINDIR)/solinprot solinprot.cc libmead.a

clean:
	-$(RM) -f libmead.a tlibmead.a *.o

.SUFFIXES:
.SUFFIXES: .cc .o

# The rule for making .o files from C++ source files:
.cc.o:
	$(CXX) -c -I.. $(COPTFLAGS) -DUSE_EXCEPTIONS $<


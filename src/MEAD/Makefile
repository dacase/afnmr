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
	@echo "[MEAD]  AR $@"
	$(VB)rm -f t$@
	$(VB)$(AR) t$@ $(MEAD_OBJECTS) > /dev/null
	$(VB)$(RANLIB) t$@
	$(VB)rm -f $@ && mv t$@ $@

solinprot: libmead.a
	@echo "[MEAD]  CXX $@"
	$(VB)$(CXX) -I.. $(COPTFLAGS) $(LDFLAGS) -DUSE_EXCEPTIONS -o $(BINDIR)/solinprot solinprot.cc libmead.a

clean:
	-$(RM) -f libmead.a tlibmead.a *.o depend

.SUFFIXES:
.SUFFIXES: .cc .o

# The rule for making .o files from C++ source files:
.cc.o:
	@echo "[MEAD]  CXX $<"
	$(VB)$(CXX) -c -I.. $(COPTFLAGS) -DUSE_EXCEPTIONS $<

depend: $(MEAD_OBJECTS:.o=.cc)
	@echo "[MEAD]  make depend"
	$(VB)$(CXX) -I.. -MM $^ > $@

include depend

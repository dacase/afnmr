#include "MEAD/FinDiffElstatPot.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/PhysCond.h"
#include "MEAD/globals.h"
#include "MEAD/AnalyMaker.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/FDGridLevel.h"


FinDiffElstatPot::FinDiffElstatPot (DielectricEnvironment_lett *e,
				    ChargeDist_lett *r,
				    ElectrolyteEnvironment_lett* ely)
 : ElstatPot_lett(e, r, ely), eps(e), rho(r), electrolyte(ely),
   analytic_approx(AnalyMaker::maker(e,r,ely)), solved(0)
{
  this_owns_analytic_approx = true;  // Ugh!
  fine_grid_pot = 0;

  // Implement a default method..
  Coord origin;
  method.add_level (65, 1.0, origin);

}

FinDiffElstatPot::FinDiffElstatPot (FinDiffMethod fdm,
				    DielectricEnvironment_lett *e,
				    ChargeDist_lett *r,
				    ElectrolyteEnvironment_lett* ely)
 : ElstatPot_lett(e, r, ely), eps(e), rho(r), electrolyte(ely), method(fdm),
   analytic_approx(AnalyMaker::maker(e,r,ely)), solved(0)
{
  this_owns_analytic_approx = true;  // Ugh!
  fine_grid_pot = 0;
}

FinDiffElstatPot::FinDiffElstatPot
(FinDiffMethod fdm, DielectricEnvironment_lett *e, ChargeDist_lett *r,
 ElectrolyteEnvironment_lett* ely, AnalyticEP* aep)
 : ElstatPot_lett(e, r, ely), eps(e), rho(r), electrolyte(ely), method(fdm),
   analytic_approx(aep), solved(0)
{
  this_owns_analytic_approx = false;  // Ugh!
  fine_grid_pot = 0;
}

FinDiffElstatPot::~FinDiffElstatPot()
{
  blab2 << "FinDiffElstatPot destructor invoked" << endl;
  FDGridLevel *next_victim = fine_grid_pot;
  while (next_victim) {
    FDGridLevel *gl = next_victim;
    next_victim = gl->coarser;
    delete gl;
  }
  if (this_owns_analytic_approx)  // Ugh!
    delete analytic_approx;
}

void FinDiffElstatPot::solve_using_coarse_init(string fieldname)
{
  blab1 << "FinDiffElstatPot::solve entered" << endl;
  const CubeLatSpec *fdm = method.get_coarsest();
  if (!fdm) {
    error ("ERROR: FinDiffElstatPot::solve: method is empty, can't solve");
  }
  fine_grid_pot = new FDGridLevel(*fdm, rho, eps, electrolyte,
				  analytic_approx);
  if (fieldname.length() == 0) // No initial field file
    fine_grid_pot->solve();
  else
    fine_grid_pot->solve_using_initial(fieldname);

  while ((fdm = method.get_finer())) {
    FDGridLevel *glp = new FDGridLevel(*fdm, rho, eps, electrolyte,
				       fine_grid_pot);
    glp->solve();
    fine_grid_pot = glp;
  }

  // HMMM! should rho retain a "gridded" variant of itself for
  // altermate rho*phi methods to be explored later???

  solved = 1; // raise the solved flag to legalize certain ops.
  blab1 << "FinDiffElstatPot::solve exiting" << endl;
}


void FinDiffElstatPot::solve()
{
  blab1 << "FinDiffElstatPot::solve entered" << endl;
  const CubeLatSpec *fdm = method.get_coarsest();
  if (!fdm) {
    error ("ERROR: FinDiffElstatPot::solve: method is empty, can't solve");
  }
  fine_grid_pot = new FDGridLevel(*fdm, rho, eps, electrolyte,
				  analytic_approx);
  fine_grid_pot->solve();
  while ((fdm = method.get_finer())) {
    FDGridLevel *glp = new FDGridLevel(*fdm, rho, eps, electrolyte,
				       fine_grid_pot);
    glp->solve();
    fine_grid_pot = glp;
  }

  // HMMM! should rho retain a "gridded" variant of itself for
  // altermate rho*phi methods to be explored later???

  solved = 1; // raise the solved flag to legalize certain ops.
  blab1 << "FinDiffElstatPot::solve exiting" << endl;
}

// FIXME  This should be inlined?
float FinDiffElstatPot::value(Coord c) const
{
  if (solved) {
    return fine_grid_pot->potint(c);
  }
  else {
    cerr << "ERROR, FinDiffElstatPot::value called for unsolved object" << endl;
  }
  return 1.0;
}

// FIXME  This should be inlined?
Coord FinDiffElstatPot::field(Coord c) const
{
  if (solved) {
    return fine_grid_pot->fieldint(c);
  }
  else {
    cerr << "ERROR, FinDiffElstatPot::field called for unsolved object" << endl;
  }
  return Coord(0,0,0);
}

// FIXME  This should be inlined?
Coord FinDiffElstatPot::displacement(Coord c) const
{
  if (solved) {
    return fine_grid_pot->displacement_int(c);
  }
  else {
    cerr << "ERROR, FinDiffElstatPot::field called for unsolved object" << endl;
  }
  return Coord(0,0,0);
}

void
FinDiffElstatPot::write_coarse_field(const string& fieldname)
{
  FDGridLevel* coarsept = fine_grid_pot;
  if (coarsept == 0)
    ::error("ERROR: FinDiffElstatPot::write_coarse_field \n",
	    "called but no field lattices created.");
  while (coarsept->coarser) coarsept = coarsept->coarser;
  coarsept->write(fieldname);
}

float* FinDiffElstatPot::coarse_field_into_array (int* grid_dim) const
{
  FDGridLevel* coarsept = fine_grid_pot;
  if (coarsept == 0)
    ::error("ERROR: FinDiffElstatPot::coarse_field_into_array\n"
	    "called but no field lattices created.");
  while (coarsept->coarser) coarsept = coarsept->coarser;
  FinDiffMethod themethod(method);
  const CubeLatSpec cls(*themethod.get_coarsest());
  const int ngr = cls.get_grid_dim();
  const int ngr_cubed = ngr*ngr*ngr;
  float* retptr;
  retptr = new float [ngr_cubed];
  if (retptr==0) {
    cerr << "Array allocation failed in coarse_field_into_array";
    *grid_dim = -1;
  }
  else {
    coarsept->get_val_array(retptr);
    *grid_dim = ngr;
  }
  return retptr;
}

CubeLatSpec FinDiffElstatPot::coarse_lattice_spec () const
{
  FinDiffMethod themethod(method);
  return CubeLatSpec (*themethod.get_coarsest());
}

// FinDiffElstatPot.cc ends here

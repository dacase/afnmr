#ifndef RISM3D_SNGLPNT_H
#define RISM3D_SNGLPNT_H

// C program to provide single point 3D-RISM calculations from the
// commandline.

// Structures

struct mdOptions {
    char *pdb, *prmtop, *rst, *traj;
};

struct rismOptions {
    char *xvv, *guv, *cuv, *huv, *uuv, *asymp, *quv, *chgdist,
        *exchem, *solvene, *entropy,
        *exchemGF, *solveneGF, *entropyGF,
        *exchemISC, *solveneISC, *entropyISC,
        *exchemUC, *solveneUC, *entropyUC,
        *potUV, *electronMap, *volfmt, *periodic, *closure;
    int closureOrder, asympcorr;
    double buffer, solvcut;
    double grdspcX, grdspcY, grdspcZ, solvboxX, solvboxY, solvboxZ;
    int ngX, ngY, ngZ;
    double mdiis_del, mdiis_restart;
    int mdiis_nvec, mdiis_method;
    int maxstep, npropagate;
    int centering, zerofrc, apply_rism_force, polarDecomp, entropicDecomp;
    int gfCorrection, pcplusCorrection;
    int ntwrism, verbose, progress, saveprogress;
    int molReconstruct;
    double uccoeff1, uccoeff2, uccoeff3, uccoeff4;
    double biasPotential;
    int treeDCF, treeTCF, treeCoulomb;
    double treeDCFMAC, treeTCFMAC, treeCoulombMAC;
    int treeDCFOrder, treeTCFOrder, treeCoulombOrder;
    int treeDCFN0, treeTCFN0, treeCoulombN0;
    double asympKSpaceTolerance, chargeSmear;
    double ljTolerance;
    int selftest;
    int request_help;
};


// Global variables to reduce the number of parameters.

// Program name.
extern char *name;

extern int ntolerance;
extern double tolerance[5];

extern int nclosure;
extern char *closure[5];

/////////////////
//MD options
/////////////////
extern struct mdOptions mdOpt;

/////////////////
//3D-RISM options
/////////////////
extern struct rismOptions rismOpt;

/////////////////
//Helper values
/////////////////
extern int centeringIsSet;

// Parse function defined in rism_options.l
void parse_options( int argc, char **argv);

#endif
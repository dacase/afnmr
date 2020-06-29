// C program to provide single point 3D-RISM calculations from the
// commandline.

// All 3D-RISM options are supported except saveprogress.
// saveprogress will be enabled when a better format is introduced.

// Tyler Luchko 2010/12
// Jesse Johnson 2014 (periodic boundary conditions, code cleanup)
// Dave Case 2018 (hand converted NAB code to C)
// Jamie Smith 2019 (split into header and c file)

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "rism3d.snglpnt.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#ifdef BINTRAJ
#  include "../sff/AmberNetcdf.h"
#endif
#include "../sff/sff.h"
#include "../sff/timer.h"

#ifdef USE_INTERNAL_GETLINE
#include <getline.h>
#endif

#ifdef MPI
int mpierror(int);
int mpifinalize(void);
int mpiinit(int *, char **, int *, int *);
int mytaskid;
int numtasks;
#else
int mytaskid = 0;
int numtasks = 1;
#endif

FILE *nabout;


// Definitions of global variables

// Program name.
char *name;

int ntolerance=1;
double tolerance[5];

int nclosure=1;
char *closure[5];

/////////////////
//MD options
/////////////////
struct mdOptions mdOpt;

/////////////////
//3D-RISM options
/////////////////
struct rismOptions rismOpt;

/////////////////
//Helper values
/////////////////
int centeringIsSet;


///////////////////////////////////////////////////////////////////////////////
//Prints an error message to stderr
//IN:
//   mytaskid : MPI process id
//   message  : message to output to user
///////////////////////////////////////////////////////////////////////////////
int error(char *message)
{
    if (mytaskid == 0) {
        return fprintf(stderr, "ERROR:  %s\n", message);
    }
    return -1;
}

///////////////////////////////////////////////////////////////////////////////
//sets default options
//IN:
//   mdOpt : mdOptions structure
//   rismOpt : rismOptions structure
///////////////////////////////////////////////////////////////////////////////
int setDefaults(struct mdOptions *mdOpt, struct rismOptions *rismOpt)
{
    mdOpt->pdb = NULL;
    mdOpt->prmtop = NULL;
    mdOpt->rst = NULL;
    mdOpt->traj = NULL;
    rismOpt->xvv = NULL;
    rismOpt->guv = NULL;
    rismOpt->huv = NULL;
    rismOpt->cuv = NULL;
    rismOpt->uuv = NULL;
    rismOpt->asymp = NULL;
    rismOpt->quv = NULL;
    rismOpt->chgdist = NULL;
    rismOpt->exchem = NULL;
    rismOpt->solvene = NULL;
    rismOpt->entropy = NULL;
    rismOpt->exchemGF = NULL;
    rismOpt->solveneGF = NULL;
    rismOpt->entropyGF = NULL;
    rismOpt->exchemISC = NULL;
    rismOpt->solveneISC = NULL;
    rismOpt->entropyISC = NULL;
    rismOpt->exchemUC = NULL;
    rismOpt->solveneUC = NULL;
    rismOpt->entropyUC = NULL;
    rismOpt->potUV = NULL;
    rismOpt->electronMap = NULL;
    rismOpt->volfmt = "dx";
    closure[0] = "kh";
    rismOpt->closureOrder = 1;
    rismOpt->asympcorr = 1;
    rismOpt->buffer = 14;
    rismOpt->solvcut = -1;
    rismOpt->grdspcX = 0.5;
    rismOpt->grdspcY = 0.5;
    rismOpt->grdspcZ = 0.5;
    rismOpt->ngX = 0;
    rismOpt->ngY = 0;
    rismOpt->ngZ = 0;
    rismOpt->solvboxX = 0;
    rismOpt->solvboxY = 0;
    rismOpt->solvboxZ = 0;
    tolerance[0] = 1.e-5;
    rismOpt->mdiis_del = 0.7;
    rismOpt->mdiis_restart = 10;
    rismOpt->mdiis_nvec = 5;
    rismOpt->mdiis_method = 2;
    rismOpt->maxstep = 10000;
    rismOpt->npropagate = 5;
    rismOpt->centering = 1;
    rismOpt->zerofrc = 1;
    rismOpt->apply_rism_force = 0;
    rismOpt->polarDecomp = 0;
    rismOpt->entropicDecomp = 0;
    rismOpt->gfCorrection = 0;
    rismOpt->pcplusCorrection = 0;
    rismOpt->periodic = NULL;
    rismOpt->ntwrism = 0;
    rismOpt->verbose = 0;
    rismOpt->progress = 1;
    rismOpt->saveprogress = 0;
    rismOpt->molReconstruct = 0;
    rismOpt->uccoeff1 = 0.0;
    rismOpt->uccoeff2 = 0.0;
    rismOpt->uccoeff3 = 0.0;
    rismOpt->uccoeff4 = 0.0;
    rismOpt->biasPotential = 0.0;
    rismOpt->treeDCF = 1;
    rismOpt->treeTCF = 1;
    rismOpt->treeCoulomb = 0;
    rismOpt->treeDCFMAC = 0.1;
    rismOpt->treeTCFMAC = 0.1;
    rismOpt->treeCoulombMAC = 0.1;
    rismOpt->treeDCFOrder = 2;
    rismOpt->treeTCFOrder = 2;
    rismOpt->treeCoulombOrder = 2;
    rismOpt->treeDCFN0 = 500;
    rismOpt->treeTCFN0 = 500;
    rismOpt->treeCoulombN0 = 500;
    rismOpt->asympKSpaceTolerance = -1.;
    rismOpt->ljTolerance = -1.;
    rismOpt->chargeSmear = 1.;
    rismOpt->selftest = 0;
    rismOpt->request_help = 0;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
//Prints option values
//IN:
//   mdOpt : mdOptions structure
//   rismOpt : rismOptions structure
///////////////////////////////////////////////////////////////////////////////
int printOptions(struct mdOptions mdOpt, struct rismOptions rismOpt)
{
    setDefaults(&mdOpt, &rismOpt);

    fprintf(stderr, "\n");
    fprintf(stderr, "%-23s %s\n", "Key", "Default");
    fprintf(stderr, "%-23s %s\n", "---", "-------");

    fprintf(stderr, "%-23s %s\n", "--pdb", mdOpt.pdb);
    fprintf(stderr, "%-23s %s\n", "--prmtop", mdOpt.prmtop);
    fprintf(stderr, "%-23s %s\n", "--rst", mdOpt.rst);
    fprintf(stderr, "%-23s %s\n", "-y|--traj", mdOpt.traj);
    fprintf(stderr, "%-23s %s\n", "--xvv", rismOpt.xvv);
    fprintf(stderr, "%-23s %s\n", "--guv", rismOpt.guv);
    fprintf(stderr, "%-23s %s\n", "--huv", rismOpt.huv);
    fprintf(stderr, "%-23s %s\n", "--cuv", rismOpt.cuv);
    fprintf(stderr, "%-23s %s\n", "--uuv", rismOpt.uuv);
    fprintf(stderr, "%-23s %s\n", "--asymp", rismOpt.asymp);
    fprintf(stderr, "%-23s %s\n", "--quv", rismOpt.quv);
    fprintf(stderr, "%-23s %s\n", "--chgdist", rismOpt.chgdist);
    fprintf(stderr, "%-23s %s\n", "--exchem", rismOpt.exchem);
    fprintf(stderr, "%-23s %s\n", "--solvene", rismOpt.solvene);
    fprintf(stderr, "%-23s %s\n", "--entropy", rismOpt.entropy);
    fprintf(stderr, "%-23s %s\n", "--exchemGF", rismOpt.exchemGF);
    fprintf(stderr, "%-23s %s\n", "--solveneGF", rismOpt.solveneGF);
    fprintf(stderr, "%-23s %s\n", "--entropyGF", rismOpt.entropyGF);
    fprintf(stderr, "%-23s %s\n", "--exchemISC", rismOpt.exchemISC);
    fprintf(stderr, "%-23s %s\n", "--solveneISC", rismOpt.solveneISC);
    fprintf(stderr, "%-23s %s\n", "--entropyISC", rismOpt.entropyISC);
    fprintf(stderr, "%-23s %s\n", "--exchemUC", rismOpt.exchemUC);
    fprintf(stderr, "%-23s %s\n", "--solveneUC", rismOpt.solveneUC);
    fprintf(stderr, "%-23s %s\n", "--entropyUC", rismOpt.entropyUC);
    fprintf(stderr, "%-23s %s\n", "--potUV", rismOpt.potUV);
    fprintf(stderr, "%-23s %d\n", "--molReconstruct", rismOpt.molReconstruct);
    fprintf(stderr, "%-23s %s\n", "--volfmt", rismOpt.volfmt);
    fprintf(stderr, "%-23s %s\n", "--periodic", rismOpt.periodic);
    fprintf(stderr, "%-23s %s\n", "--closure", rismOpt.closure);
    fprintf(stderr, "%-23s %d\n", "--closureorder", rismOpt.closureOrder);
    fprintf(stderr, "%-23s %i\n", "--asympcorr", rismOpt.asympcorr);
    fprintf(stderr, "%-23s %f\n", "--buffer", rismOpt.buffer);
    fprintf(stderr, "%-23s %f\n", "--solvcut", rismOpt.solvcut);
    fprintf(stderr, "%-23s %f,%f,%f\n", "--grdspc", rismOpt.grdspcX,
            rismOpt.grdspcY, rismOpt.grdspcZ);
    fprintf(stderr, "%-23s %d,%d,%d\n", "--ng", rismOpt.ngX, rismOpt.ngY,
            rismOpt.ngZ);
    fprintf(stderr, "%-23s %f,%f,%f\n", "--solvbox", rismOpt.solvboxX,
            rismOpt.solvboxY, rismOpt.solvboxZ);
    fprintf(stderr,"%-23s %f\n", "--tolerance", tolerance[0]);
    fprintf(stderr, "%-23s %f\n", "--mdiis_del", rismOpt.mdiis_del);
    fprintf(stderr, "%-23s %f\n", "--mdiis_restart",
            rismOpt.mdiis_restart);
    fprintf(stderr, "%-23s %d\n", "--mdiis_nvec", rismOpt.mdiis_nvec);
    //  fprintf(stderr,"%-23s = %d\n", "--mdiis_method", rismOpt.mdiis_method);
    fprintf(stderr, "%-23s %d\n", "--maxstep", rismOpt.maxstep);
    fprintf(stderr, "%-23s %d\n", "--npropagate", rismOpt.npropagate);
    fprintf(stderr, "%-23s %d\n", "--centering", rismOpt.centering);
    //  fprintf(stderr,"%-23s = %d\n", "--zerofrc", rismOpt.zerofrc);
    //  fprintf(stderr,"%-23s = %d\n", "--apply_rism_force", rismOpt.apply_rism_force);
    fprintf(stderr, "%-23s %d\n", "--polarDecomp", rismOpt.polarDecomp);
    fprintf(stderr, "%-23s %d\n", "--entropicDecomp",
            rismOpt.entropicDecomp);
    fprintf(stderr, "%-23s %d\n", "--gf", rismOpt.gfCorrection);
    fprintf(stderr, "%-23s %d\n", "--pc+", rismOpt.pcplusCorrection);
    fprintf(stderr, "%-23s %f,%f,%f,%f\n", "--uccoeff", rismOpt.uccoeff1,
            rismOpt.uccoeff2, rismOpt.uccoeff3, rismOpt.uccoeff4);
    fprintf(stderr, "%-23s %d\n", "--verbose", rismOpt.verbose);
    //  fprintf(stderr,"%-23s = %d\n", "--progress", rismOpt.progress);
    //  fprintf(stderr,"%-23s = %d\n", "--saveprogress", rismOpt.saveprogress);
    fprintf(stderr,"%-23s %d\n", "--treeDCF", rismOpt.treeDCF);
    fprintf(stderr,"%-23s %d\n", "--treeTCF", rismOpt.treeTCF);
    fprintf(stderr,"%-23s %d\n", "--treeCoulomb", rismOpt.treeCoulomb);
    fprintf(stderr,"%-23s %d\n", "--treeDCFOrder", rismOpt.treeDCFOrder);
    fprintf(stderr,"%-23s %d\n", "--treeTCFOrder", rismOpt.treeTCFOrder);
    fprintf(stderr,"%-23s %d\n", "--treeCoulombOrder", rismOpt.treeCoulombOrder);
    fprintf(stderr,"%-23s %d\n", "--treeDCFN0", rismOpt.treeDCFN0);
    fprintf(stderr,"%-23s %d\n", "--treeTCFN0", rismOpt.treeTCFN0);
    fprintf(stderr,"%-23s %d\n", "--treeCoulombN0", rismOpt.treeCoulombN0);
    fprintf(stderr,"%-23s %f\n", "--treeDCFMAC", rismOpt.treeDCFMAC);
    fprintf(stderr,"%-23s %f\n", "--treeTCFMAC", rismOpt.treeTCFMAC);
    fprintf(stderr,"%-23s %f\n", "--treeCoulombMAC", rismOpt.treeCoulombMAC);
    fprintf(stderr,"%-23s %f\n", "--asympKSpaceTolerance", rismOpt.asympKSpaceTolerance);
    fprintf(stderr,"%-23s %f\n", "--ljTolerance", rismOpt.ljTolerance);
    fprintf(stderr,"%-23s %f\n", "--chargeSmear", rismOpt.chargeSmear);
    fprintf(stderr, "%-23s %i\n", "--selftest", rismOpt.selftest);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
//Presents usage information for the program to stderr and exits with status 1.
//IN:
//   name : name of the program
//   mytaskid : MPI process id
///////////////////////////////////////////////////////////////////////////////
int usage()
{
    char empty[2] = " ";
    int i;
    if (mytaskid == 0) {
        fprintf(stderr,
                "USAGE: %s --pdb pdbfile --prmtop prmtopfile [--rst rstfile] [-y|--traj netCDFfile]\n",
                name);
        fprintf(stderr,
                "       %s --xvv Xvv_filename [--guv Guv_rootname]\n",
                empty);
        fprintf(stderr,
                "       %s [--cuv Cuv_rootname] [--huv Huv_rootname]\n",
                empty);
        fprintf(stderr,
                "       %s [--uuv Uuv_rootname] [--asymp asymp_rootname]\n",
                empty);
        fprintf(stderr,
                "       %s [--quv Quv_rootname] [--chgdist chgdist_rootname]\n",
                empty);
        fprintf(stderr,
                "       %s [--exchem exchem_rootname] [--solvene solvene_rootname]\n",
                empty);
        fprintf(stderr, "       %s [--entropy entropy_rootname]\n", empty);
        fprintf(stderr,
                "       %s [--exchemGF exchemGF_rootname] [--solveneGF solveneGF_rootname]\n",
                empty);
        fprintf(stderr, "       %s [--entropyGF entropyGF_rootname]\n",
                empty);
        fprintf(stderr,
                "       %s [--exchemISC exchemISC_rootname] [--solveneISC solveneISC_rootname]\n",
                empty);
        fprintf(stderr, "       %s [--entropyISC entropyISC_rootname]\n",
                empty);
        fprintf(stderr,
                "       %s [--exchemUC exchemUC_rootname] [--solveneUC solveneUC_rootname]\n",
                empty);
        fprintf(stderr, "       %s [--entropyUC entropyUC_rootname]\n",
                empty);
        fprintf(stderr, "       %s [--potUV potUV_rootname]\n", empty);
        fprintf(stderr, "       %s [--molReconstruct]\n", empty);
        fprintf(stderr, "       %s [--volfmt volume_format]\n", empty);
        fprintf(stderr, "       %s [--periodic pme|ewald]\n", empty);
        fprintf(stderr,
                "       %s [--closure kh|hnc|pse(1|2|...)[,closure2[ ...]]]\n",
                empty);
        fprintf(stderr,
                "       %s [--[no]asympcorr] [--buffer distance] [--solvcut distance]\n",
                empty);
        fprintf(stderr,
                "       %s [--grdspc dx,dy,dz] [--ng nx,ny,nz] [-solvbox lx,ly,lz]\n",
                empty);
        fprintf(stderr,
                "       %s [--tolerance tol1[ tol2 [...]] [--mdiis_del step_size]\n",
                empty);
        fprintf(stderr,
                "       %s [--mdiis_restart threshold] [--mdiis_nvec vectors]\n",
                empty);
        fprintf(stderr,
                "       %s [--maxstep 1000] [--npropagate #_old_solutions]\n",
                empty);
        fprintf(stderr,
                "       %s [--[no]polarDecomp] [--[no]entropicDecomp]\n",
                empty);
        fprintf(stderr, "       %s [--gf] [--pc+]\n", empty);
        fprintf(stderr, "       %s [--centering -4..4]\n", empty);
        fprintf(stderr,"       %s [--periodic potential_abbrev]\n",empty);
        fprintf(stderr, "       %s [--uccoeff a0,b0,a1,b1]\n", empty);
        // fprintf(stderr,"       %s [--biasPotential biasValue]\n",empty);
	fprintf(stderr,"       %s [--[no]treeDCF] [--[no]treeTCF] [--[no]treeCoulomb]\n",
		empty);
	fprintf(stderr,"       %s [--treeDCFMAC multipole-accept] [--treeTCFMAC multipole-accept] [--treeCoulombMAC multipole-accept]\n",
		empty);
	fprintf(stderr,"       %s [--treeDCFOrder Taylor-order] [--treeTCFOrder Taylor-order] [--treeCoulombOrder Taylor-order]\n",
		empty);
	fprintf(stderr,"       %s [--treeDCFN0 max-leaf-size] [--treeTCFN0 max-leaf-size] [--treeCoulombN0 max-leaf-size]\n",
		empty);
	fprintf(stderr,"       %s [--asympKSpaceTolerance tol] [--chargeSmear smear]\n",empty);
	fprintf(stderr,"       %s [--ljTolerance tol]\n",empty);
        fprintf(stderr, "       %s [--verbose 0|1|2]\n", empty);
        fprintf(stderr, "       %s [--selftest\n", empty);
        printOptions(mdOpt, rismOpt);
    }
    exit(1);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
//counts the number of Comma Separated Values in a string
//IN:
//   str : string
//OUT:
//   Number of comma separated values.  Includes blank values so this can be
//   thought of as number of commas plus one.
///////////////////////////////////////////////////////////////////////////////
int numCSV(char *str)
{
    int i, count;
    for (i = 0, count = 0; str[i]; i++)
        count += (str[i] == ',');
    return count + 1;
}

///////////////////////////////////////////////////////////////////////////////
//Final check to ensure options are correctly set.
//IN:
//   mdOpt : mdOptions structure
//   rismOpt : rismOptions structure
///////////////////////////////////////////////////////////////////////////////
int checkOptions(struct mdOptions mdOpt, struct rismOptions rismOpt)
{
    if (rismOpt.request_help) usage();
    if (mdOpt.rst == NULL && mdOpt.traj == NULL) {
        error("a RESTART or TRAJ file is required\n");
        usage();
    }
    if (mdOpt.prmtop == NULL) {
        error("a PRMTOP file is required\n");
        usage();
    }
    if (rismOpt.xvv == NULL) {
        error("an XVV file is required\n");
        usage();
    }
    return 0;
}

int main(int argc, char **argv)
{

    nabout = stdout; //  default; otherwise set to a file pointer of your choice
#ifdef MPI
	if ( mpiinit(&argc, argv, &mytaskid, &numtasks) != 0 ) {
		printf("Error in mpiinit!");
		fflush(stdout);
		exit(1);
	}
#endif

#ifdef BINTRAJ
    struct AmberNetcdf nc;
#endif
    double time = 0., dgrad, fret;
    int i, j;
    double ftmp;
    FILE *asciiTraj;
    int iframe, ipos, trajDone;

    //set global variables
    asprintf(&name, "%s", argv[0]);
    centeringIsSet = 0;

    setDefaults(&mdOpt, &rismOpt);
    parse_options(argc, argv);
    checkOptions(mdOpt, rismOpt);

    ////////////////////
    //Run 3D-RISM
    ////////////////////

    PARMSTRUCT_T *prm = rdparm(mdOpt.prmtop);   // reads the prmtop file
    int natm = prm->Natom;
    int natm3 = 3 * natm;
    double *p_xyz = malloc(natm3 * (sizeof(double)));
    double *f_xyz = malloc(natm3 * (sizeof(double)));
    double *v_xyz = malloc(natm3 * (sizeof(double)));

    char *mmopt;

    if (rismOpt.periodic) {
        asprintf(&mmopt, "cut=%e", sqrt(400));
        mm_options(mmopt);
        free(mmopt);
    } else {
        asprintf(&mmopt, "cut=%e", sqrt(10e37));
        mm_options(mmopt);
        free(mmopt);
    }

    asprintf(&mmopt, "rism=1, ntpr_rism=1, xvvfile=%s", rismOpt.xvv);
    mm_options(mmopt);
    free(mmopt);

    if (rismOpt.guv) {
        asprintf(&mmopt, "guvfile=%s", rismOpt.guv);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.cuv) {
        asprintf(&mmopt, "cuvfile=%s", rismOpt.cuv);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.huv) {
        asprintf(&mmopt, "huvfile=%s", rismOpt.huv);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.uuv) {
        asprintf(&mmopt, "uuvfile=%s", rismOpt.uuv);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.asymp) {
        asprintf(&mmopt, "asympfile=%s", rismOpt.asymp);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.quv) {
        asprintf(&mmopt, "quvfile=%s", rismOpt.quv);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.chgdist) {
        asprintf(&mmopt, "chgdistfile=%s", rismOpt.chgdist);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.exchem) {
        asprintf(&mmopt, "exchemfile=%s", rismOpt.exchem);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.solvene) {
        asprintf(&mmopt, "solvenefile=%s", rismOpt.solvene);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.entropy) {
        asprintf(&mmopt, "entropyfile=%s", rismOpt.entropy);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.exchemGF) {
        asprintf(&mmopt, "exchemGFfile=%s", rismOpt.exchemGF);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.solveneGF) {
        asprintf(&mmopt, "solveneGFfile=%s", rismOpt.solveneGF);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.entropyGF) {
        asprintf(&mmopt, "entropyGFfile=%s", rismOpt.entropyGF);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.exchemISC) {
        asprintf(&mmopt, "exchemISCfile=%s", rismOpt.exchemISC);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.solveneISC) {
        asprintf(&mmopt, "solveneISCfile=%s", rismOpt.solveneISC);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.entropyISC) {
        asprintf(&mmopt, "entropyISCfile=%s", rismOpt.entropyISC);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.exchemUC) {
        asprintf(&mmopt, "exchemUCfile=%s", rismOpt.exchemUC);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.solveneUC) {
        asprintf(&mmopt, "solveneUCfile=%s", rismOpt.solveneUC);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.entropyUC) {
        asprintf(&mmopt, "entropyUCfile=%s", rismOpt.entropyUC);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.potUV) {
        asprintf(&mmopt, "potUVfile=%s", rismOpt.potUV);
        mm_options(mmopt);
        free(mmopt);
    }
    if (rismOpt.electronMap) {
        asprintf(&mmopt, "electronMapfile=%s", rismOpt.electronMap);
        mm_options(mmopt);
        free(mmopt);
    }

    asprintf(&mmopt, "molReconstruct=%d", rismOpt.molReconstruct);
    mm_options(mmopt);
    free(mmopt);
    
    asprintf(&mmopt, "volfmt=%s", rismOpt.volfmt);
    mm_options(mmopt);
    free(mmopt);

    asprintf(&mmopt, "rst=%s", mdOpt.rst);
    mm_options(mmopt);
    free(mmopt);

    char *tmpstr;
    if (nclosure >= 1) {
        asprintf(&tmpstr, "closure=%s", closure[0]);
    }
    for (i = 1; i < nclosure; i++) {
        asprintf(&tmpstr, "%s,%s", tmpstr, closure[i]);
    }

    asprintf(&mmopt, "%s, closureOrder=%i", tmpstr, rismOpt.closureOrder);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "buffer=%.15g, solvcut=%.15g", rismOpt.buffer,
             rismOpt.solvcut);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "grdspcx=%.15g, grdspcy=%.15g, grdspcz=%.15g",
             rismOpt.grdspcX, rismOpt.grdspcY, rismOpt.grdspcZ);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "solvboxx=%.15g, solvboxy=%.15g, solvboxz=%.15g",
             rismOpt.solvboxX, rismOpt.solvboxY, rismOpt.solvboxZ);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "ngx=%d, ngy=%d, ngz=%d", rismOpt.ngX, rismOpt.ngY,
             rismOpt.ngZ);
    mm_options(mmopt);
    free(mmopt);

    if (ntolerance >= 1) {
        asprintf(&tmpstr, "tolerance=%.15g", tolerance[0]);
    }
    for (i = 1; i < ntolerance; i++) {
        asprintf(&tmpstr, "%s,%.15g", tmpstr, tolerance[i]);
    }

    asprintf(&mmopt, "%s, mdiis_del=%.15g", tmpstr, rismOpt.mdiis_del);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "mdiis_restart=%.15g", rismOpt.mdiis_restart);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "mdiis_nvec=%d, mdiis_method=%d", rismOpt.mdiis_nvec,
             rismOpt.mdiis_method);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "maxstep=%d, npropagate=%d", rismOpt.maxstep,
             rismOpt.npropagate);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "centering=%d, zerofrc=%d, apply_rism_force=%d",
             rismOpt.centering, rismOpt.zerofrc, rismOpt.apply_rism_force);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "polarDecomp=%d", rismOpt.polarDecomp);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "entropicDecomp=%d", rismOpt.entropicDecomp);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "gfCorrection=%d", rismOpt.gfCorrection);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "pcplusCorrection=%d", rismOpt.pcplusCorrection);
    mm_options(mmopt);
    free(mmopt);
    if (rismOpt.periodic) {
        asprintf(&mmopt, "periodic=%s", rismOpt.periodic);
        mm_options(mmopt);
        free(mmopt);
    }
    asprintf(&mmopt, "ntwrism=%d, verbose=%d, progress=%d",
             rismOpt.ntwrism, rismOpt.verbose, rismOpt.progress);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "asympCorr=%d", rismOpt.asympcorr);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "selftest=%d", rismOpt.selftest);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "saveprogress=%d", rismOpt.saveprogress);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "uccoeff=%.15g,%.15g,%.15g,%.15g",
             rismOpt.uccoeff1, rismOpt.uccoeff2,
             rismOpt.uccoeff3, rismOpt.uccoeff4);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "treeDCF=%d, treeTCF=%d, treeCoulomb=%d",
	    rismOpt.treeDCF, rismOpt.treeTCF, rismOpt.treeCoulomb);
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "treeDCFMAC=%.15g, treeTCFMAC=%.15g, treeCoulombMAC=%.15g",
	    rismOpt.treeDCFMAC, rismOpt.treeTCFMAC, rismOpt.treeCoulombMAC); 
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "treeDCFOrder=%d, treeTCFOrder=%d, treeCoulombOrder=%d",
	    rismOpt.treeDCFOrder, rismOpt.treeTCFOrder, rismOpt.treeCoulombOrder); 
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "treeDCFN0=%d, treeTCFN0=%d, treeCoulombN0=%d",
	    rismOpt.treeDCFN0, rismOpt.treeTCFN0, rismOpt.treeCoulombN0); 
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "asympKSpaceTolerance=%.15g, chargeSmear=%.15g",
	    rismOpt.asympKSpaceTolerance, rismOpt.chargeSmear); 
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "ljTolerance=%.15g", rismOpt.ljTolerance); 
    mm_options(mmopt);
    free(mmopt);
    asprintf(&mmopt, "biasPotential=%.15g", rismOpt.biasPotential);
    mm_options(mmopt);
    free(mmopt);

    // nothing frozen or constrained for now:
    int frozen[natm], constrained[natm];
    memset(frozen, 0, natm * sizeof(int));
    memset(constrained, 0, natm * sizeof(int));

    mme_init_sff(prm, frozen, constrained, NULL, NULL);

    if (mdOpt.traj) {
#ifdef BINTRAJ
        // Try NetCDF first.
        if (netcdfLoad(&nc, mdOpt.traj) == 0) {
            netcdfLoad(&nc, mdOpt.traj);
            if (mytaskid == 0)
                printf("\nProcessing NetCDF trajectory: %s\n", mdOpt.traj);
            while (netcdfGetNextFrame(&nc, p_xyz, NULL)) {
                if (mytaskid == 0)
                    printf("\nFrame: %d of %d\n", nc.currentFrame,
                           nc.ncframe);
                mme(p_xyz, f_xyz, &nc.currentFrame);
            }
            netcdfClose(&nc);
        } else {
#endif
            // Assume ASCII.
            if (mytaskid == 0) {
                if (numtasks > 1) {
                    error
                        ("ASCII trajectories not supported for more than one process");
                    exit(1);
                }
                printf("\nProcessing ASCII trajectory: %s\n", mdOpt.traj);
            }
            trajDone = 0;
            asciiTraj = fopen(mdOpt.traj, "r");
            if (asciiTraj == NULL) {
                fprintf(stderr, "Failed to open '%s'", mdOpt.traj);
                exit(1);
            }
            char* line = NULL;
            size_t llength = 0;
            getline(&line, &llength, asciiTraj); // read header
            for (iframe = 1;; iframe++) {
                for (ipos = 0; ipos < 3 * natm; ipos++) {
                    if (fscanf(asciiTraj, "%lf", &p_xyz[ipos]) < 1) {
                        trajDone = 1;
                        break;
                    }
                }
                if (trajDone)
                    break;
                if (mytaskid == 0) {
                    printf("\nFrame: %d\n", iframe);
                }
                mme(p_xyz, f_xyz, &iframe);
            }
#ifdef BINTRAJ
        }
#endif

    } else if (mdOpt.rst) {

        if (mytaskid == 0)
             printf("\nProcessing restart file: %s\n", mdOpt.rst);
        getxv(mdOpt.rst, natm, time, p_xyz, v_xyz);
        if (mytaskid == 0)
            printf("\nRunning MME\n");
        int iter = 1;
        mme(p_xyz, f_xyz, &iter);
        if (mytaskid == 0)
            printf("\nFinished MME\n");

    } else {
        fprintf( stderr, "Must provide either a restart or trajectory file.");
        exit(1);
    }

    if (mytaskid == 0) {
        printf("\n3D-RISM processing complete.\n");
    }

    mme_rism_max_memory();
    mme_timer();
#ifdef MPI
    mpifinalize();
#endif

}

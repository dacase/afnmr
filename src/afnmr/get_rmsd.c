#include "../nabc/nabc.h"
FILE* nabout;

int main( int argc, char *argv[] )
{

    MOLECULE_T *mol, *molref;
    REAL_T backbone, nonH;    // computed rmsd between mol and molref
    char *heavy_exp = "::C*,N*,O*,P*";

	nabout = stdout; /*default*/

    mol = getpdb( argv[1], NULL );
    molref = getpdb( argv[3], NULL );
    superimpose( mol, argv[2], molref, argv[4] );
    rmsd( &mol, &argv[2], &molref, &argv[4], &backbone );
    superimpose( mol, "::C*,N*,O*,P*", molref, "::C*,N*,O*,P*" );
    rmsd( &mol, &heavy_exp, &molref, &heavy_exp, &nonH );
    printf( "backbone rmsd = %10.4f,  nonH rmsd = %10.4f\n", backbone, nonH );

	exit( 0 );
}

%option prefix="rism"
%option noyywrap

%{
extern int mytaskid;

#include <rism3d.snglpnt.h>
#include <assert.h>

#undef ECHO
#define ECHO yyout=stderr;if(mytaskid==0) fprintf(yyout, "\tcmdline:  %s\n",yytext)

char* setVarStr ( char *str);
int setVarInt ( char *str);
double setVarDbl ( char *str);
void setArrayStr ( char *str, char *a, int n, int len);
void setArrayInt ( char *str, int *a, int n);
void setArrayDbl ( char *str, double *a, int n);
char *token;

%}
Float [+-]?((([0-9]+\.?[0-9]*)|([0-9]*\.?[0-9]+))([eE][-+]?[0-9]+)?)
Int [+-]?[0-9]+
PosInt [+]?[1-9][0-9]*
NegInt [-]?[1-9][0-9]*
NonNegInt [+]?[0-9]+
Str [^,\ =]+
Assign (\ *=\ *|\ +)
ListSep (\ *,\ *)
Closure (kh|pse1|pse2|pse3|pse4|hnc)
%%

--pdb{Assign}{Str}               { mdOpt.pdb=setVarStr(yytext); }
--rst{Assign}{Str}               { mdOpt.rst=setVarStr(yytext); }
--traj{Assign}{Str}              { mdOpt.traj=setVarStr(yytext); }
--prmtop{Assign}{Str}            { mdOpt.prmtop=setVarStr(yytext); }

--xvv{Assign}{Str}               { rismOpt.xvv=setVarStr(yytext);}

--guv{Assign}{Str}               { rismOpt.guv=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--huv{Assign}{Str}               { rismOpt.huv=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--cuv{Assign}{Str}               { rismOpt.cuv=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--uuv{Assign}{Str}               { rismOpt.uuv=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--asymp{Assign}{Str}             { rismOpt.asymp=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--quv{Assign}{Str}               { rismOpt.quv=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--chgdist{Assign}{Str}           { rismOpt.chgdist=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--exchem{Assign}{Str}            { rismOpt.exchem=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--solvene{Assign}{Str}           { rismOpt.solvene=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--entropy{Assign}{Str}           { rismOpt.entropy=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--exchemGF{Assign}{Str}          { rismOpt.exchemGF=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--solveneGF{Assign}{Str}         { rismOpt.solveneGF=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--entropyGF{Assign}{Str}         { rismOpt.entropyGF=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--exchemISC{Assign}{Str}         { rismOpt.exchemISC=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--solveneISC{Assign}{Str}        { rismOpt.solveneISC=setVarStr(yytext);
                                   rismOpt.ntwrism=1; }
--entropyISC{Assign}{Str}        { rismOpt.entropyISC=setVarStr(yytext);
                                   rismOpt.ntwrism=1; }
--exchemUC{Assign}{Str}          { rismOpt.exchemUC=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--solveneUC{Assign}{Str}         { rismOpt.solveneUC=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--entropyUC{Assign}{Str}         { rismOpt.entropyUC=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--potUV{Assign}{Str}             { rismOpt.potUV=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--electronMap{Assign}{Str}       { rismOpt.electronMap=setVarStr(yytext); 
                                   rismOpt.ntwrism=1; }
--molReconstruct                 { rismOpt.molReconstruct=1; }
--volfmt{Assign}(ccp4|dx|xyzv)   { rismOpt.volfmt=setVarStr(yytext); }
--periodic{Assign}(pme|ewald|pmen|ewaldn)    { rismOpt.periodic=setVarStr(yytext); 
                                   rismOpt.centering=0; }

 /* read in a comma-separated array of closure names: */
--closure{Assign}({Closure}{ListSep}){0,5}{Closure} {
             token = strtok( yytext, "= " );
             token = strtok( NULL, ", " );
             int ic=0;
             while( token != NULL ){
                closure[ic] = malloc( strlen(token) + 1 );
                strncpy( closure[ic], token, 5 );
                ic++;
                token = strtok( NULL, ", " );
             }
             nclosure = ic;
}

 /* read in a comma-separated array of tolerances: */
--tolerance{Assign}({Float}{ListSep}){0,5}{Float} {
             token = strtok( yytext, "= " );
             token = strtok( NULL, ", " );
             int ic=0;
             while( token != NULL ){
                sscanf(token,"%lg",&tolerance[ic]);
                ic++;
                token = strtok( NULL, ", " );
             }
             ntolerance = ic;
}

--closureOrder{Assign}[123]       { rismOpt.closureOrder=setVarInt(yytext); }
--asympCorr{Assign}[01]           { rismOpt.asympcorr=setVarInt(yytext); }
--solvcut{Assign}{Float}          { rismOpt.solvcut=setVarDbl(yytext); }
--buffer{Assign}{Float}           { rismOpt.buffer=setVarDbl(yytext); }
--grdspc{Assign}{Float}           { rismOpt.grdspcX=setVarDbl(yytext); \
                                    rismOpt.grdspcY = rismOpt.grdspcX; \
                                    rismOpt.grdspcZ = rismOpt.grdspcX; }
--grdspc{Assign}{Float},{Float},{Float}  { \
       token=strtok(yytext,"= "); \
       token=strtok(NULL,",");    sscanf(token,"%lg",&rismOpt.grdspcX); 
       token=strtok(NULL,",");    sscanf(token,"%lg",&rismOpt.grdspcY); 
       token=strtok(NULL," ");    sscanf(token,"%lg",&rismOpt.grdspcZ); }
--grdspcx{Assign}{Float}          { rismOpt.grdspcX=setVarDbl(yytext); }
--grdspcy{Assign}{Float}          { rismOpt.grdspcY=setVarDbl(yytext); }
--grdspcz{Assign}{Float}          { rismOpt.grdspcZ=setVarDbl(yytext); }
--ng{Assign}{PosInt}              { rismOpt.ngX=setVarInt(yytext);  
                                    rismOpt.ngY = rismOpt.ngX; 
                                    rismOpt.ngZ = rismOpt.ngX; }
--ng{Assign}{PosInt},{PosInt},{PosInt}  { \
       token=strtok(yytext,"= "); \
       token=strtok(NULL,",");    sscanf(token,"%d",&rismOpt.ngX); 
       token=strtok(NULL,",");    sscanf(token,"%d",&rismOpt.ngY); 
       token=strtok(NULL," ");    sscanf(token,"%d",&rismOpt.ngZ); }
--ngx{Assign}{Int}                { rismOpt.ngX=setVarInt(yytext); }
--ngy{Assign}{Int}                { rismOpt.ngY=setVarInt(yytext); }
--ngz{Assign}{Int}                { rismOpt.ngZ=setVarInt(yytext); }
--solvbox{Assign}{Float}          { rismOpt.solvboxX=setVarDbl(yytext); 
                                    rismOpt.solvboxY = rismOpt.solvboxX; 
                                    rismOpt.solvboxZ = rismOpt.solvboxX; }
--solvbox{Assign}{Float},{Float},{Float}  { \
       token=strtok(yytext,"= "); \
       token=strtok(NULL,",");    sscanf(token,"%lg",&rismOpt.solvboxX); 
       token=strtok(NULL,",");    sscanf(token,"%lg",&rismOpt.solvboxY); 
       token=strtok(NULL," ");    sscanf(token,"%lg",&rismOpt.solvboxZ); }
--solvboxx{Assign}{Float}         { rismOpt.solvboxX=setVarDbl(yytext); }
--solvboxy{Assign}{Float}         { rismOpt.solvboxY=setVarDbl(yytext); }
--solvboxz{Assign}{Float}         { rismOpt.solvboxZ=setVarDbl(yytext); }
--mdiis_del{Assign}{Float}        { rismOpt.mdiis_del=setVarDbl(yytext); }
--mdiis_restart{Assign}{Float}    { rismOpt.mdiis_restart=setVarDbl(yytext); }
--mdiis_nvec{Assign}{PosInt}      { rismOpt.mdiis_nvec=setVarInt(yytext); }
--mdiis_method{Assign}[012]       { rismOpt.mdiis_method=setVarInt(yytext); }
--maxstep{Assign}{PosInt}         { rismOpt.maxstep=setVarInt(yytext); }
--npropagate{Assign}[0-5]         { rismOpt.npropagate=setVarInt(yytext); }
--centering{Assign}[0-4\-\+]+     { rismOpt.centering=setVarInt(yytext); }
--zerofrc{Assign}[01]             { rismOpt.zerofrc=setVarInt(yytext); }
--apply_rism_force{Assign}[01]    { rismOpt.apply_rism_force=setVarInt(yytext); }
--polarDecomp                     { rismOpt.polarDecomp=1; }
--entropicDecomp                  { rismOpt.entropicDecomp=1; }
--gf                              { rismOpt.gfCorrection=1; }
--pc\+                            { rismOpt.pcplusCorrection=1; }
--saveprogress{Assign}{NonNegInt} { rismOpt.saveprogress=setVarInt(yytext); }
--ntwrism{Assign}{Int}            { rismOpt.ntwrism=setVarInt(yytext); }
--verbose{Assign}[012]            { rismOpt.verbose=setVarInt(yytext); }
--progress                        { rismOpt.progress=1; }
--noprogress                      { rismOpt.progress=0; }
--uccoeff{Assign}{Float},{Float},{Float},{Float}  { \
       token=strtok(yytext,"= "); \
       token=strtok(NULL,","); sscanf(token,"%lg",&rismOpt.uccoeff1); 
       token=strtok(NULL,",");    sscanf(token,"%lg",&rismOpt.uccoeff2); 
       token=strtok(NULL,",");    sscanf(token,"%lg",&rismOpt.uccoeff3); 
       token=strtok(NULL," ");    sscanf(token,"%lg",&rismOpt.uccoeff4); }
--biasPotential{Assign}{Float}    { rismOpt.biasPotential=setVarDbl(yytext); }
--treeDCF                     { rismOpt.treeDCF=1; }
--treeTCF                     { rismOpt.treeTCF=1; }
--treeCoulomb                 { rismOpt.treeCoulomb=1; }
--notreeDCF                   { rismOpt.treeDCF=0; }
--notreeTCF                   { rismOpt.treeTCF=0; }
--notreeCoulomb               { rismOpt.treeCoulomb=0; }
--treeDCFOrder{Assign}{NonNegInt} { rismOpt.treeDCFOrder=setVarInt(yytext); }
--treeTCFOrder{Assign}{NonNegInt} { rismOpt.treeTCFOrder=setVarInt(yytext); }
--treeCoulombOrder{Assign}{NonNegInt} { rismOpt.treeCoulombOrder=setVarInt(yytext); }
--treeDCFN0{Assign}{NonNegInt}    { rismOpt.treeDCFN0=setVarInt(yytext); }
--treeTCFN0{Assign}{NonNegInt}    { rismOpt.treeTCFN0=setVarInt(yytext); }
--treeCoulombN0{Assign}{NonNegInt} { rismOpt.treeCoulombN0=setVarInt(yytext); }
--treeDCFMAC{Assign}{Float}  { rismOpt.treeDCFMAC=setVarDbl(yytext); }
--treeTCFMAC{Assign}{Float}  { rismOpt.treeTCFMAC=setVarDbl(yytext); }
--treeCoulombMAC{Assign}{Float} { rismOpt.treeCoulombMAC=setVarDbl(yytext); }
--asympKSpaceTolerance{Assign}{Float} { rismOpt.asympKSpaceTolerance=setVarDbl(yytext); }
--ljTolerance{Assign}{Float}      { rismOpt.ljTolerance=setVarDbl(yytext); }
--chargeSmear{Assign}{Float}       { rismOpt.chargeSmear=setVarDbl(yytext); }
--selftest                        { rismOpt.selftest=1; }
--help                            { rismOpt.request_help=1; }
-h                                { rismOpt.request_help=1; }

[\ \n\t,]+                                                ;

.                               {        
 char c;
 fprintf( stderr, "Error: unknown keyword or invalid value: %s", yytext );
 for(c=input();c != ' ' && c != '=' && c != EOF && c != '\0';c=input())
   fprintf( stderr, "%c",c);
 fprintf( stderr,"\n");
 fprintf( stderr, "Type 'rism3d.snglpnt --help' to see valid options\n" );
 exit( 1 );
}
%%

/* Parse the command-line:  */

void parse_options( int argc, char **argv)
{

   char cmdline[4096];
   int index=0;
   int argl;
   int i;

   /* assemble the command-line tokens back into a single string */
   for( i=1; i<argc; i++ ){
      argl = strlen( argv[i]);
      strcpy( &cmdline[index], argv[i] );
      cmdline[ index + argl ] = ' '; 
      index += argl + 1;
      assert( index < 4096 );
   }
   cmdline[ index - 1 ] = '\0';

   rism_scan_string( cmdline );
   rismlex();
   yy_delete_buffer(YY_CURRENT_BUFFER);
}

#define	PGM(name,purpose) /* name: purpose */
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

/*   char	*strcpy (), *malloc ();  */
#define	strdup(s)	strcpy(malloc((unsigned)strlen(s)+1),s)

char	*Argv0 = "";

#ifdef US_DEBUG
#undef US_DEBUG
#define	US_DEBUG(msg) fputs ("debug: msg\n", stderr);
#else
#define US_DEBUG(msg)
#endif

#define	WARNING(msg) fprintf (stderr, "\007%s: msg.\n", Argv0);
#define	USAGE(synopsis) {\
				fprintf (stderr, "\007Usage: %s synopsis\n", Argv0);\
				exit (1);\
				}
#define	ERRMSG0(msg) {\
				fprintf (stderr, "\007%s: msg.\n", Argv0);\
				exit (1);\
				}
#define	ERRMSG1(msg, arg1) {\
				fprintf (stderr, "\007%s: msg. %d\n", Argv0, arg1);\
				exit (1);\
				}
#define	ERRMSG2(msg, arg1, arg2)\
				{\
				fprintf (stderr, "\007%s: msg. %s %s\n", Argv0, arg1, arg2);\
				exit (1);\
				}
#define	ERRMSG3(msg, arg1, arg2, arg3)\
				{\
				fprintf (stderr, "\007%s: msg. %s %s %s\n", Argv0, arg1, arg2, arg3);\
				exit (1);\
				}
#define ERRDATA              ERRMSG0 (Not enough (or no) input data)
#define ERRMANY(stuff,many)  ERRMSG1 (Too much stuff; at most %d allowed, many)
#define	ERROPEN(file)        ERRMSG1 (Cannot open '%s', file)
#define	ERROPT(opt)          ERRMSG1 (Unknown option -- '%c', opt)
#define	ERRSPACE(whatever)   ERRMSG0 (No storage space left for whatever)
#define	ERRNUM(string)       ERRMSG1 ('%s' is not a number, string)
#define	ERRRAGGED            ERRMSG0 (Ragged input file)


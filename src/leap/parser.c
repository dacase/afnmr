/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     LVARIABLE = 258,
     LSTRING = 259,
     LNUMBER = 260,
     LASSIGN = 261,
     LENDOFCOMMAND = 262,
     LOPENLIST = 263,
     LCLOSELIST = 264,
     LOPENPAREN = 265,
     LCLOSEPAREN = 266,
     LQUIT = 267,
     LCOMMA = 268,
     LDOT = 269,
     LCOMMAND = 270,
     LDUMMY = 271,
     LNULL = 272,
     LNOTSINGLECHAR = 273
   };
#endif
/* Tokens.  */
#define LVARIABLE 258
#define LSTRING 259
#define LNUMBER 260
#define LASSIGN 261
#define LENDOFCOMMAND 262
#define LOPENLIST 263
#define LCLOSELIST 264
#define LOPENPAREN 265
#define LCLOSEPAREN 266
#define LQUIT 267
#define LCOMMA 268
#define LDOT 269
#define LCOMMAND 270
#define LDUMMY 271
#define LNULL 272
#define LNOTSINGLECHAR 273




/* Copy the first part of user declarations.  */
#line 74 "parser.y"

#include	<unistd.h>
#include        "basics.h"

#include        "classes.h"

#include        "dictionary.h"
#include        "parmLib.h"

#include        "commands.h"
#include	"block.h"
#include	"parser.h"

#include        "leap.h"
#include        "block.h"

#include        "help.h"

#define         MESSAGEFILTER   MESSPARSER




#define         NULLSTR         "null"

/*
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        GLOBAL VARIABLES

*/

                /* The global variable GfCurrentInput defines the file */
                /* from which input is currently read.  When the file */
                /* is empty, then switch back to stdin */

#define	MAXINPUT	1000
#define	MAXINPUTFILES	10		/* Maximum 10 input files */
					/* can be open at once */


char            GsInputLine[MAXINPUT] = "";
BOOL		GbLastLine = FALSE;
BOOL		bCmdDeleteObj;
int             GiInputPos = 0;
PARMLIB		GplAllParameters;
RESULTt		GrMainResult;
BLOCK		GbCommand = NULL;
BLOCK		GbExecute = NULL;
int		GiClipPrompts = 0;
BOOL		GbGraphicalEnvironment;
STRING		GsProgramName;

extern int	iMemDebug;

static	STRING	*SbFirstSourceFiles = NULL;
static	int	iFirstSource = 0;
static	BOOL	SbUseStartup = TRUE;



/*
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *
 *	The following is used by the parser to maintain a stack of
 *	files where input is received from.  If the file is NULL
 *	then input is received from the main program in the form
 *	of BLOCKS.  The main program is then responsible for reading
 *	the stdin ( in the command line interface ) or for gathering
 *	keypress events from X-Windows ( in the graphical interface ).
 *
 */


int             GiInputFileStackPos = 0;
FILE*           GfaInputFileStack[MAXINPUTFILES];



/*
 *----------------------------------------------------------------
 *
 *	Not quite GLOBAL variables used by the parser.
 */

                                /* Arguments to routines are passed through */
                                /* an array */
#define MAXARGS         10
#define MAXLISTNEXT     10


ATOM            aDummy;
ASSOC           aaArgs[MAXARGS];
int             iArgCount, i;

                                /* List stuff is used for input of nested */
                                /* lists */
#define MAXLISTNEST     10
ASSOC           aaLists[MAXLISTNEST];
int             iCurrentList = -1;
#define PUSHLIST()      iCurrentList++
#define POPLIST()       iCurrentList--
#define CURRENTLIST     aaLists[iCurrentList]


OBJEKT          o0;
double          dTemp;
ASSOC           aAssoc;
STRING          sTemp;
BOOL		bQuit = FALSE;
BOOL		bCommandFound = FALSE;

                /* There seems to be a problem with YACC not properly */
                /* declaring yylval and yyval */
typedef struct  {
	ASSOC		aVal;
	double		dVal;
	STRING		sVal;
	FUNCTION	fCallback;
} YYSTYPEt;

#define YYSTYPE YYSTYPEt


extern  OBJEKT oGetObject( char *sName );
extern  int     yyparse();

/*  avoid compiler warnings:  */
int yylex();
int yyerror( char * );



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 275 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  13
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   46

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  19
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  14
/* YYNRULES -- Number of rules.  */
#define YYNRULES  27
/* YYNRULES -- Number of states.  */
#define YYNSTATES  34

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   273

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     5,     7,     9,    12,    15,    18,    22,
      25,    27,    29,    30,    35,    37,    39,    41,    43,    45,
      46,    50,    51,    54,    56,    57,    59,    62
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      20,     0,    -1,    21,    -1,     7,    -1,    22,    -1,     1,
       7,    -1,    23,     7,    -1,    27,     7,    -1,     3,     6,
      24,    -1,     3,     6,    -1,    25,    -1,    27,    -1,    -1,
       8,    26,    29,     9,    -1,     5,    -1,     4,    -1,     3,
      -1,    16,    -1,    17,    -1,    -1,    30,    28,    31,    -1,
      -1,    29,    25,    -1,    15,    -1,    -1,    32,    -1,    31,
      32,    -1,    25,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   228,   228,   234,   235,   239,   246,   247,   250,   268,
     275,   276,   280,   279,   293,   302,   312,   320,   327,   339,
     338,   371,   372,   382,   385,   386,   387,   390
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "LVARIABLE", "LSTRING", "LNUMBER",
  "LASSIGN", "LENDOFCOMMAND", "LOPENLIST", "LCLOSELIST", "LOPENPAREN",
  "LCLOSEPAREN", "LQUIT", "LCOMMA", "LDOT", "LCOMMAND", "LDUMMY", "LNULL",
  "LNOTSINGLECHAR", "$accept", "input", "line", "instruct", "assign",
  "express", "rawexp", "@1", "function", "@2", "elements", "cmdname",
  "args", "arg", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    19,    20,    21,    21,    21,    22,    22,    23,    23,
      24,    24,    26,    25,    25,    25,    25,    25,    25,    28,
      27,    29,    29,    30,    31,    31,    31,    32
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     1,     2,     2,     2,     3,     2,
       1,     1,     0,     4,     1,     1,     1,     1,     1,     0,
       3,     0,     2,     1,     0,     1,     2,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     3,    23,     0,     2,     4,     0,     0,
      19,     5,     9,     1,     6,     7,    24,    16,    15,    14,
      12,    17,    18,     8,    10,    11,    27,    20,    25,    21,
      26,     0,    13,    22
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     5,     6,     7,     8,    23,    26,    29,     9,    16,
      31,    10,    27,    28
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -18
static const yytype_int8 yypact[] =
{
      31,    -3,    -1,   -18,   -18,     7,   -18,   -18,     1,     2,
     -18,   -18,    -2,   -18,   -18,   -18,    20,   -18,   -18,   -18,
     -18,   -18,   -18,   -18,   -18,   -18,   -18,    20,   -18,   -18,
     -18,    13,   -18,   -18
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -18,   -18,   -18,   -18,   -18,   -18,   -12,   -18,     0,   -18,
     -18,   -18,   -18,   -17
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      24,    17,    18,    19,    11,    12,    20,    13,    14,    15,
      30,     0,    25,     4,    21,    22,    17,    18,    19,    33,
       0,    20,    32,    17,    18,    19,     0,     0,    20,    21,
      22,     0,     1,     0,     2,     0,    21,    22,     3,     0,
       0,     0,     0,     0,     0,     0,     4
};

static const yytype_int8 yycheck[] =
{
      12,     3,     4,     5,     7,     6,     8,     0,     7,     7,
      27,    -1,    12,    15,    16,    17,     3,     4,     5,    31,
      -1,     8,     9,     3,     4,     5,    -1,    -1,     8,    16,
      17,    -1,     1,    -1,     3,    -1,    16,    17,     7,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    15
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     7,    15,    20,    21,    22,    23,    27,
      30,     7,     6,     0,     7,     7,    28,     3,     4,     5,
       8,    16,    17,    24,    25,    27,    25,    31,    32,    26,
      32,    29,     9,    25
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 229 "parser.y"
    {
			return 0;
		}
    break;

  case 4:
#line 236 "parser.y"
    {
                        bCommandFound = FALSE;
                        }
    break;

  case 5:
#line 240 "parser.y"
    {
                            VP0(( "\n" ));
                            yyerrok;
                        }
    break;

  case 8:
#line 251 "parser.y"
    {
                                /* Set the value of the variable */
                            aAssoc = (yyvsp[(3) - (3)].aVal);
			    if ( aAssoc != NULL ) {
                                MESSAGE(( "Assigning a value to %s\n", 
                                                        (yyvsp[(1) - (3)].sVal) ));
                                VariableSet( (yyvsp[(1) - (3)].sVal), oAssocObject(aAssoc) );
				MESSAGE(( "DEREF (assign) - %s\n",
							sAssocName(aAssoc) ));
                                DEREF( aAssoc );
			    } else {
				MESSAGE(("Not assigning value to %s - rmving\n",
							(yyvsp[(1) - (3)].sVal) ));
				VariableRemove( (yyvsp[(1) - (3)].sVal) );
			    }

                        }
    break;

  case 9:
#line 269 "parser.y"
    {
                            MESSAGE(( "Removing variable %s\n", (yyvsp[(1) - (2)].sVal) ));
                            VariableRemove( (yyvsp[(1) - (2)].sVal) );
                        }
    break;

  case 12:
#line 280 "parser.y"
    {
                            PUSHLIST();
                                /* Create an ASSOC for the list */
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            AssocSetObject( aAssoc, oCreate(LISTid) );
                            CURRENTLIST = aAssoc;
                        }
    break;

  case 13:
#line 289 "parser.y"
    {
                            (yyval.aVal) = (yyvsp[(3) - (4)].aVal);
                            POPLIST();
                        }
    break;

  case 14:
#line 294 "parser.y"
    {
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            o0 = oCreate(ODOUBLEid);
                            ODoubleSet( o0, (yyvsp[(1) - (1)].dVal) );
                            AssocSetObject( aAssoc, o0 );
                            (yyval.aVal) = aAssoc;
                        }
    break;

  case 15:
#line 303 "parser.y"
    {
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            o0 = oCreate(OSTRINGid);
                            OStringDefine( (OSTRING) o0, (yyvsp[(1) - (1)].sVal) );
                            AssocSetObject( aAssoc, o0 );
			    DEREF( o0 );	/* keeps count = 1 */
                            (yyval.aVal) = aAssoc;
                        }
    break;

  case 16:
#line 313 "parser.y"
    {
			    OBJEKT	o = oGetObject((yyvsp[(1) - (1)].sVal));
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, (yyvsp[(1) - (1)].sVal) );
                            AssocSetObject( aAssoc, o ); /* REF's o */
                            (yyval.aVal) = aAssoc;
                        }
    break;

  case 17:
#line 321 "parser.y"
    {
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            AssocSetObject( aAssoc, aDummy );
                            (yyval.aVal) = aAssoc;
                        }
    break;

  case 18:
#line 328 "parser.y"
    {
                            MESSAGE(( "Parsed a null\n" ));
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            AssocSetObject( aAssoc, NULL );
                            (yyval.aVal) = aAssoc;
                        }
    break;

  case 19:
#line 339 "parser.y"
    { iArgCount = 0;
                        }
    break;

  case 20:
#line 342 "parser.y"
    {
                                /* Execute the command */
			    MESSAGE(( "executing function\n"));
			    bCmdDeleteObj = FALSE;
                            o0 = (yyvsp[(1) - (3)].fCallback)( iArgCount, aaArgs );
			    if ( o0 != NULL ) {
                            	aAssoc = (ASSOC)oCreate(ASSOCid);
                            	AssocSetObject( aAssoc, o0 );
                            	(yyval.aVal) = aAssoc;
			    } else {
				MESSAGE(( "func == NULL---\n"));
				(yyval.aVal) = NULL;
			    }

                                /* DEREF each of the arguments */

                            for ( i=0; i<iArgCount; i++ ) {
				if ( bCmdDeleteObj ) {
					MESSAGE(( "bCmdDeleteObj---\n" ));
					DEREF( oAssocObject( aaArgs[i] ) );
				}
				MESSAGE(( "DEREF (function) - %s\n",
						sAssocName(aaArgs[i])));
                                DEREF( aaArgs[i] );
                            }
                        }
    break;

  case 22:
#line 373 "parser.y"
    {
                                /* Get the element and add it to the list */
                            MESSAGE(( "Adding to list:\n" ));
                            ListAddToEnd( (LIST)oAssocObject(CURRENTLIST), 
							(OBJEKT) (yyvsp[(2) - (2)].aVal) );
                            (yyval.aVal) = CURRENTLIST;
                        }
    break;

  case 27:
#line 391 "parser.y"
    {
                            aaArgs[iArgCount++] = (yyvsp[(1) - (1)].aVal);
                        }
    break;


/* Line 1267 of yacc.c.  */
#line 1676 "y.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 407 "parser.y"

/*------------------------------------------------------------

        ROUTINES

*/

static  BOOL    SbGotUngetc = FALSE;
static  char    ScUngetc;


/*
 *      yyerror
 *
 *      Respond to errors.
 */
int
yyerror( char *sStr )
{
    VPFATALEXIT(( "Error from the parser: %s\n", sStr ));
    return 1;
}

FILE *
fINPUTFILE()            
{
	if ( GiInputFileStackPos < 0 )
		return(NULL);
	return( GfaInputFileStack[GiInputFileStackPos] );
}



/*
 *	zbGetLine
 *
 *	Get the next line of input from the current input source.
 *	This may be GbCommand or GbExecute depending if
 *	the input is coming from the user or from an execute file.
 *	Return FALSE if there is no more lines to be received from 
 *	the BLOCK.
 */
BOOL
zbGetLine( char *sLine, BOOL *bPFromExecute )
{
BOOL		bGotBlock;
char		c;

    /* VPTRACEENTER(( "zbGetLine" )); */

    if ( fINPUTFILE() == NULL ) {
	*bPFromExecute = FALSE;
	/* VPTRACEEXIT (( "zbGetLine" )); */
	return(bBlockReadLine( GbCommand, sLine ));
    }

    *bPFromExecute = TRUE;

		/* If there is no line in the execute BLOCK */
		/* then read in another block from the current file */

    if ( bBlockEndOfRead(GbExecute) ) {
	    BlockEmpty( GbExecute );
	    bGotBlock = FALSE;
	    while ( !feof(fINPUTFILE()) && !bGotBlock ) {
		c = fgetc(fINPUTFILE());
		if ( feof(fINPUTFILE()) ) break;
		bGotBlock = bBlockAddChar( GbExecute, c );
	    }

		/* If a complete BLOCK was not read then */
		/* append a final '\n' character, if that doesn't */
		/* make this a complete BLOCK then there is an error */
		/* in the input, also we are at the end of the file */
		/* so pop the file off the stack and say that we have */
		/* a complete BLOCK */

	    if ( !bGotBlock ) {
	    	bGotBlock = bBlockAddChar( GbExecute, '\n' );
	    	if ( fINPUTFILE() != NULL )
	    		fclose( fINPUTFILE() );
	    	INPUTPOPFILE();
	    	VPTRACE(( "After pop in zbGetLine GiInputFileStackPos = %d\n",
	    	          GiInputFileStackPos ));
	    }
    }
    /* VPTRACEEXIT (( "zbGetLine" )); */
    return(bBlockReadLine( GbExecute, sLine ));
}


		
	


/*
 *      zcGetChar
 *
 *      Get the next character from the line buffer.
 *	If there are no more characters in the line buffer and
 *	GbLastLine is TRUE then return '\0', otherwise if the
 *	line buffer is empty, fill it and return the next character
 *	in it.
 */
char
zcGetChar()
{
char            c;
BOOL		bFromExecute;

                /* If there is a pushed character then return it */

    if ( SbGotUngetc ) {
        SbGotUngetc = FALSE;
        c = ScUngetc;
        goto DONE;
    }

		/* Now if the input line is empty then fill it */
    if ( GsInputLine[GiInputPos] == '\0' ) {
	if ( GbLastLine ) {
	    c = '\0';
	    GbLastLine = FALSE;
	    goto DONE;
	}
	GbLastLine = zbGetLine( GsInputLine, &bFromExecute );
	GiInputPos = 0;
	if ( bFromExecute ) {
	    for ( i=GiClipPrompts; i<=iINPUTSTACKDEPTH(); i++ ) VP2(( ">" ));
	    VP2(( " " ));
	    VP2(( "%s", GsInputLine ));
	} else {
	    VPLOG(( "> %s", GsInputLine ));
	}
    }

    c = GsInputLine[GiInputPos++];

DONE:
    return(c);

}



/*
 *      zUngetc
 *
 *      Push one character back to the file.
 */
void
zUngetc( char c )
{
    SbGotUngetc = TRUE;
    ScUngetc = c;
}






/*
 *      cGetChar
 *
 *      Return the next character, skipping over comments
 *      which start with '#' and end with '\n'.
 */
char
cGetChar()
{
char    c;

	while (1) {
        	c = zcGetChar();
        	if ( c == '#' )
            		while ( ( c=zcGetChar() ) != '\n' )  /* Nothing */ ;
		else
			break;
        } 
	return(c);
}

 
int
iOneCharToken( int c )
{
	switch(c) {
		case ';':
        		return(LENDOFCOMMAND);
		case '=':
        		return(LASSIGN);
		case '*':
        		return(LDUMMY);
		case '(':
        		return(LOPENPAREN);
		case ')':
        		return(LCLOSEPAREN);
		case '{':
        		return(LOPENLIST);
		case '}':
        		return(LCLOSELIST);
		default:
        		return(LNOTSINGLECHAR);
	}
}

/*
 *      yylex
 *
 *      Lexical analyzer for the parser.
 *      Read characters from stdin and return the token types
 *      read and place the value read into the global UNION
 *      yylval.
 *
 *      The things that it recognizes are:
 *              LDOUBLE         [-]###.###E## or ###
 *              LSTRING         "xx xx xxx" or '$'everything up to ' ' ',' ';'
 *              commands        xxxxxxx
 *              LVARIABLE       xxxxxxx which are not commands
 *              LTERMINATOR     ;
 *              LASSIGN         =
 *
 *	Modified 17 November 1992 - David A. Rivkin
 *		Added checking the alias table for command matches.
 *	Total rewrite October 1993 - Bill Ross
 *
 */
int
yylex()
{
STRING          sStr;
int             j, iMax, tok;
BOOL            bGotExp, bGotDot;
char            c;
STRING		sCmd;
STRING		sPossibleCmd;

                /* Skip over blanks, tabs, end of lines etc */

    while ( (c=cGetChar())==' '  ||  c=='\t'  ||  c=='\n'  ||  c == ',' );

    if ( c == '\0' ) 
	return(LENDOFCOMMAND);

    /*
     *  Check the 1-character possibilities: , ; = * ( ) { }
     */
    tok = iOneCharToken( c );
    if ( tok != LNOTSINGLECHAR ) {
        MESSAGE(( "Parsed /%c/\n", c ));
	return(tok);
    }

    /*
     *  it isn't a 1-char thing; read in the rest 
     *	and push back the terminating char
     */
    sStr[0] = c;
    for (j=1;;j++) {

	if ( j >= sizeof(STRING) )
	    DFATAL(( "string too long" ));

	c = cGetChar();
	/*
	 *  NULL terminates anything (?)
	 */
	if ( c == '\0' ) {
	    	sStr[j] = '\0';
	    	break;
	}
	/*
	 *  allow anything inside quotes; chop closing quote
	 */
	if ( sStr[0] == '"' ) {
		if ( c == '"' ) {
			sStr[j] = '\0';
			break;
		}
		sStr[j] = c;
		continue;
	}
	/*
	 *  whitespace is a delimiter outside of quotes
	 */
	if ( c == ' '  ||  c == '\t'  ||  c == '\n'  ||  c == ',' ) {
	    	sStr[j] = '\0';
	    	break;
	}
	/*
	 *  special case for $-type strings: allow embedded single-char
	 *	tokens, except ';'
	 */
	if ( sStr[0] == '$' ) {
		if ( c == ';' ) {
			zUngetc( c );
	    		sStr[j] = '\0';
	    		break;
		}
		sStr[j] = c;
		continue;
	}
	tok = iOneCharToken( c );
	if ( tok != LNOTSINGLECHAR ) {
		zUngetc( c );
	    	sStr[j] = '\0';
	    	break;
	}
	sStr[j] = c;
    }

    /*
     *  see if it's a number
     */
    bGotExp = FALSE;
    bGotDot = FALSE;
    if ( isdigit(sStr[0]) || sStr[0] == '-' || sStr[0] == '+' || 
				( sStr[0] == '.' && isdigit(sStr[1]) ) ) {

        for ( j=0; j<sizeof(STRING); j++ ) {
            MESSAGE(( "Thinking NUMBER got: %c\n", sStr[j] ));
	    switch ( sStr[j] ) {
		case '\0':
        	    if ( sscanf( sStr, "%lf", &yylval.dVal ) != 1 ) {
			VPWARN(( " Couldn't scan NUMBER from (%s)\n", sStr ));
			return(LNULL);
		    }
        	    MESSAGE(( "Parsed a number: %lf\n", yylval.dVal ));
        	    return(LNUMBER);
		case '.':
		    if ( bGotDot ) {
			VPWARN(( "(Multiple '.' in NUMBER-like thing (%s))\n", 
								sStr ));
			goto notnum;
		    }
		    if ( bGotExp ) {
			VPWARN(( 
			 "('.' follows exponent in NUMBER-like thing (%s))\n",
								sStr ));
			goto notnum;
		    }
        	    bGotDot = TRUE;
		    break;
		case 'e':
		case 'E':
            	    if ( bGotExp ) {
			VPWARN(( "(Multiple 'e' in NUMBER-like thing (%s))\n", 
								sStr ));
			goto notnum;
		    }
                    bGotExp = TRUE;
                    break;
		case '+':
		case '-':
                    break;
		default:
            	    if ( !isdigit(sStr[j]) ) {
			goto notnum;
		    }
            	    break;
	    }
	}
    }
notnum:
    /* 
     *  see if it's a string in quotes
     */
    if ( sStr[0] == '"' ) {
        strcpy( yylval.sVal, &sStr[1] );
        MESSAGE(( "Parsed a STRING: %s\n", sStr ));
        return(LSTRING);
    }

    /* 
     *  see if it's a string prefixed w/ '$'
     */
    if ( sStr[0] == '$' ) {
        strcpy( yylval.sVal, &sStr[1] );
        MESSAGE(( "Parsed a STRING: %s\n", sStr ));
        return(LSTRING);
    }

                /* LASTLY!!!!!!!! */
    /* 
     *  see if it's a variable/command 
     */
    strcpy( yylval.sVal, sStr );
    strcpy( sPossibleCmd, sStr );
    StringLower( sPossibleCmd );

    		/* Check if there is an alias that is an exact match */
    if ( (iMax = iVarArrayElementCount( GvaAlias )) ) {
	ALIAS		aAlias;
	aAlias = PVAI( GvaAlias, ALIASt, 0 );
	for ( i=0; i<iMax; i++, aAlias++ ) {
	    if ( strcmp( aAlias->sName, sPossibleCmd ) == 0 ) {
	    	strcpy( sPossibleCmd, aAlias->sCommand );
	    }
        }
    }
                /* Check if there is an exact match of the command */
                /* If a command has already been found for this input
                	line, then do not consider the string a command
                	but rather as a STRING variable */
                	
    if ( !bCommandFound ) {
	for ( j=0; strlen(cCommands[j].sName) != 0; j++ ) {
	    strcpy( sCmd, cCommands[j].sName );
	    StringLower( sCmd );
            if ( strcmp( sCmd, sPossibleCmd ) == 0 ) {
		yylval.fCallback = cCommands[j].fCallback;
		MESSAGE(( "Parsed a command: %s\n", sStr ));
		bCommandFound = TRUE;
		return(LCOMMAND);
	    }
        }
    }


                /* If the variable name is null then return LNULL */

    if ( strcmp( sStr, NULLSTR ) == 0 ) 
	return(LNULL);
    
                /* Return the variable name */

    strcpy( yylval.sVal, sStr );
    MESSAGE(( "Parsed a variable: %s\n", sStr ));
    return(LVARIABLE);

}




/*
 *      oGetObject
 *
 *      If the string is a variable then return the OBJEKT that
 *      it is attached to, otherwise if the string is
 *      a string like: 'unit.mol.res.atom' then parse the
 *      individual names and search the CONTAINERS for
 *      the subcontainers.
 */
OBJEKT
oGetObject( char *sName )
{
CONTAINER       cCont[5];
int             j, k, iSeq;
OBJEKT          oObj;
STRING          sLine, sHead;
BOOL		bDot, bAt, bPdbSeq;
STRING		sGroup;
LIST		lGroup;
OSTRING		osString;
LOOP		lRes;
RESIDUE		rRes;

    oObj = oVariable( sName );
    if ( oObj != NULL ) return(oObj);

        /* Now try to parse the name */
    strcpy( sLine, sName );
    bDot = FALSE;
    bAt = FALSE;
    bPdbSeq = FALSE;
    for ( k=0; k<strlen(sName); k++ ) {
	if ( sLine[k] == '.' ) {
	    sLine[k]=' ';
	    bDot = TRUE;
/*fprintf(stderr, "GOTDOT\n"); */
	} else if ( sLine[k] == '@' ) {
	    sLine[k] = ' ';
	    bAt = TRUE;
	} else if ( sLine[k] == '%' ) {
	    if ( !bDot ) {
	        sLine[k] = ' ';
	        bPdbSeq = TRUE;
	    }
	}
    }


    sRemoveFirstString( sLine, sHead );
    cCont[0] = (CONTAINER)oVariable(sHead);
    if ( cCont[0] == NULL ) {
/* fprintf(stderr, "STRING %s\n", sName); */
    	/* It is not an object variable so...
    	   return the whole string as a OSTRING */
	goto String;
    }
    
/*fprintf(stderr, "VARIABLE %c %c\n", bAt, bPdbSeq);*/
    if ( bPdbSeq ) {
	sRemoveLeadingSpaces( sLine );
	sRemoveFirstString( sLine, sHead );
	if ( strlen(sHead) == 0 ) {
	    goto String;
	}
	if ( isdigit(sHead[0]) ) {

			/* Make sure the rest are digits */
			/* If not return NULL */
	    for ( j=1; j<strlen(sHead); j++ ) {
		if ( !isdigit(sHead[j]) ) {
			    /* It is not an object variable so...
			       return the whole string as a OSTRING */
		    goto String;
		}
	    }
	
			/* Find the PDB sequence number */

	    cCont[1] = NULL;	
	    iSeq = atoi(sHead);
	    lRes = lLoop( (OBJEKT)cCont[0], RESIDUES );
	    while ( (rRes = (RESIDUE)oNext(&lRes)) ) {
		if ( iResiduePdbSequence(rRes)==iSeq ) {
		    cCont[1] = (CONTAINER) rRes;
		}
	    }

	    if ( !cCont[1] ) {
		goto String;
	    }

	    if ( bDot ) {
		sRemoveLeadingSpaces( sLine );
		sRemoveFirstString( sLine, sHead );
		if ( strlen(sHead) == 0 ) {
		    goto String;
		}
		if ( isdigit(sHead[0]) ) {

				/* Make sure the rest are digits */
				/* If not return NULL */
		    for ( j=1; j<strlen(sHead); j++ ) {
			if ( !isdigit(sHead[j]) ) {
				    /* It is not an object variable so...
				       return the whole string as a OSTRING */
			    goto String;
			}
		    }
			
		    iSeq = atoi(sHead);
		    cCont[2] = cContainerFindSequence( cCont[1], 
						DIRECTCONTENTS, iSeq );
		} else {
		    cCont[2] = 
			cContainerFindName( cCont[1], DIRECTCONTENTS, sHead );
		}
		return((OBJEKT)cCont[2]);
	    } else {
		return((OBJEKT)cCont[1]);
	    }
	} else {
	    goto String;
	}
    }


    if ( bDot ) {
	k = 0;
	do {
	    sRemoveLeadingSpaces( sLine );
	    sRemoveFirstString( sLine, sHead );
	    if ( strlen(sHead) == 0 ) break;
	    k++;
	    if ( cCont[k-1]->oHeader.cObjType == PARMSETid ) {
		/*  semi-HACK - parmsets don't have contents in this sense */
		cCont[k] = NULL;
	    } else if ( isdigit(sHead[0]) ) {

			    /* Make sure the rest are digits */
			    /* If not return NULL */
		for ( j=1; j<strlen(sHead); j++ ) {
		    if ( !isdigit(sHead[j]) ) {
   				/* It is not an object variable so...
			    	   return the whole string as a OSTRING */
			goto String;
    		    }
		}
		iSeq = atoi(sHead);
		cCont[k] = cContainerFindSequence( cCont[k-1], 
					    DIRECTCONTENTS, iSeq );
	    } else {
		cCont[k] = 
			cContainerFindName( cCont[k-1], DIRECTCONTENTS, sHead );
	    }
	    if ( cCont[k] == NULL ) break;
	} while ( strlen(sHead) != 0 ) ;
        if ( cCont[k] == NULL ) {
    		/* It is not an object variable so...
    		   return the whole string as a OSTRING */
		goto String;
    	}
	return((OBJEKT)cCont[k]);
    }

			/* If group notation then return the group */

    if ( bAt ) {
	sRemoveLeadingSpaces( sLine );
	sRemoveFirstString( sLine, sGroup );
	if ( iObjectType(cCont[0]) != UNITid ) {
    		/* It is not an object variable so...
    		   return the whole string as a OSTRING */
		goto String;
	}

	lGroup = lUnitGroup( (UNIT)cCont[0], sGroup );
	return((OBJEKT)lGroup);
    }

String:

   /* 
    *  It is not an object variable so...
    *	   return the whole string as a OSTRING
    *	   -- need to set refs to 0 so that
    *	      it will be freed later.. HACK
    */

    osString = (OSTRING)oCreate(OSTRINGid);
    OStringDefine( osString, sName );
    ((OBJEKT)(osString))->iReferences = 0;
    return((OBJEKT)osString );
}





    
/*
 *================================================================
 *
 *	Public routines
 */




/*
 *	ParseArguments
 *
 *	Parse the arguments that the user passes to LEaP from
 *	the command line arguments.
 */
void
ParseArguments( int argc, char *argv[] )
{
char		c;
extern	char	*optarg;

    while ( (c = getopt( argc, argv, "hsI:f:" )) != (char)(EOF) ) {
	switch (c) {
	    case 'h':
		printf( "Usage: %s [options]\n", argv[0] );
		printf( "Options:\n" );
		printf( " -h         Generate this message.\n" );
		printf( " -s         Ignore %s startup file.\n", LEAPRC );
		printf( " -I {dir}   Add {dir} to search path.\n" );
		printf( " -f {file}  Source {file}.\n" );
		exit(0);
	    case 's':
		printf( "-s: Ignoring startup file: %s\n", LEAPRC );
		SbUseStartup = FALSE;
		break;
	    case 'I':
		printf( "-I: Adding %s to search path.\n", optarg );
		BasicsAddDirectory( optarg, 1 );
		break;
	    case 'f':
		printf( "-f: Source %s.\n", optarg );
		if ( iFirstSource == 0 ) {
			MALLOC( SbFirstSourceFiles, STRING *, sizeof(STRING) );
			iFirstSource = 1;
		} else {
			iFirstSource++;
			REALLOC( SbFirstSourceFiles, STRING *, SbFirstSourceFiles,
					iFirstSource * sizeof(STRING));
		}
		strcpy( SbFirstSourceFiles[iFirstSource-1], optarg );
		break;
	}
    }
}






/*
 *	ParseInit
 *
 *	Initialize the parser.
 *	If SbStartup is TRUE the execute the LEAPRC script.
 */
void
ParseInit( RESULTt *rPResult )
{
FILE	*fStartup;
int	iFile;

    VPTRACEENTER(( "ParseInit" ));
    VP0(( "\nWelcome to LEaP!\n" ));

#ifdef  DEBUG
    VP0(( "LEaP is running in DEBUG mode!\n" ));
#endif

		/* Initialize memory manager debugging */
    INITMEMORYDEBUG();

    HelpInitialize();

                /* Initialize the first file in the stack to be stdin */
                
    GfaInputFileStack[0] = NULL;
    GbExecute = bBlockCreate();

                /* Create a few OBJEKTs that will be used by the parser */

    aDummy = (ATOM)oCreate(ATOMid);
    ContainerSetName( aDummy, "DUMMY" );
    GplAllParameters = plParmLibCreate();

    VariablesInit();
    GrMainResult.iCommand = CNONE;
    rPResult->iCommand = CNONE;
    
                /* Parse the LEAPRC file if bUseStartup is TRUE */

    if ( SbUseStartup ) {
        fStartup = FOPENNOCOMPLAIN( LEAPRC, "r" );
        if ( fStartup == NULL ) {
	    VP0(( "(no %s in search path)\n", LEAPRC ));
	} else {
	    /*
	     *  source the leaprc
	     */
	    VP0(( "Sourcing %s: %s\n", LEAPRC, GsBasicsFullName ));
	    INPUTPUSHFILE( fStartup );
	    GiClipPrompts = 1;
	    while ( fINPUTFILE() != NULL ) {
	        yyparse();
		if ( GrMainResult.iCommand == CQUIT ) {
	    	    if ( fINPUTFILE() != NULL )
	    	    	fclose( fINPUTFILE() );
	    	    INPUTPOPFILE();
		}
	    }
	    *rPResult = GrMainResult;
	    GiClipPrompts = 0;
	}
    }

		/* Parse the source files specified */
		/* on the command line using the -f options */

    for (iFile=0; iFile<iFirstSource; iFile++) {
	fStartup = FOPENCOMPLAIN( SbFirstSourceFiles[iFile], "r" );
	if ( fStartup != NULL ) {
	    VP0(( "Sourcing: %s\n", GsBasicsFullName ));
	    INPUTPUSHFILE( fStartup );
	    VPTRACE(( "After push GiInputFileStackPos = %d\n",
	               GiInputFileStackPos ));
	    GiClipPrompts = 1;
	    while ( fINPUTFILE() != NULL ) {
		yyparse();
		if ( GrMainResult.iCommand == CQUIT ) {
	    	    if ( fINPUTFILE() != NULL )
	    	    	fclose( fINPUTFILE() );
	    	    INPUTPOPFILE();
	    	    VPTRACE(( "After pop GiInputFileStackPos = %d\n",
	    	              GiInputFileStackPos ));
		}
	    }
	    *rPResult = GrMainResult;
	    GiClipPrompts = 0;
	} else {
	    exit(41);
	}
    }
    if ( SbFirstSourceFiles )
	FREE( SbFirstSourceFiles );
    VPTRACEEXIT (( "ParseInit" ));
}




/*
 *	ParseBlock
 *
 *	Parse a BLOCK containing one complete command.
 *	Return in rResult the result of the command.
 */
void
ParseBlock( BLOCK bBlock, RESULTt *rPResult )
{
		/* Set up the BLOCK from which to read the command */

    VPTRACEENTER(( "ParseBlock" ));
    MESSAGE(( "Parsing block: %s\n", sBlockText(bBlock) ));
    GbCommand = bBlock;
    BlockResetRead( GbCommand );

    GrMainResult.iCommand = CNONE;

		/* Parse the BLOCK */
		/* Keep parsing as long as the 'execute' command */
		/* keeps setting the GLOBAL variable GbContinueParsing */

    do {
	yyparse();
	if ( GrMainResult.iCommand == CQUIT ) {
	    if ( fINPUTFILE() != NULL )
	    	fclose( fINPUTFILE() );
	    INPUTPOPFILE();
	    VPTRACE(( "After pop GiInputFileStackPos = %d\n",
	              GiInputFileStackPos ));
	}
    } while ( fINPUTFILE() != NULL );

	/* Reset the bCommandFound variable as a new command may be
	   available in a new input */
    bCommandFound = FALSE;
    
	/* Copy the RESULT from the global result variable */

    *rPResult = GrMainResult;

    VPTRACEEXIT (( "ParseBlock" ));
}





/*
 *	ParseShutdown
 *
 *	Shutdown the parser, release all variables setup in ParseInit
 *	Only needed if debugging memory mgt, since prog mem is all
 *	freed when process exits anyway.
 */
void
ParseShutdown()
{
	if ( !iMemDebug )
		return;

	Destroy( (OBJEKT *)&aDummy );
	VariablesDestroy();
	ParmLibDestroy( &GplAllParameters );

	BlockDestroy( &GbExecute );

	HelpShutdown();
}


/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

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
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 2 "parser.y" /* yacc.c:339  */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../math/units.h"
#include "../structures/list.h"
#include "../structures/hashmap.h"
#include "../system/color.h"
#include "../system/error.h"
#include "../sensors/sensor.h"

extern FILE * yyin;

int yylex();
void yyerror(char * message);

typedef struct Trigger Trigger;
typedef struct Numeric Numeric;
typedef struct Specification Specification;

Charge        * make_charge(Numeric * wire, Specification * tag);
Trigger       * make_trigger(List * charges, List * options);
Specification * specify_trigger(char * id, bool less, Numeric * threshold, Trigger * trigger);
Specification * make_tag(char * id, List * options, List * args);

void build_sensor(char * id, Numeric * frequency, List * specifications);
void print_config();

#line 96 "y.tab.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "y.tab.h".  */
#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif
/* "%code requires" blocks.  */
#line 32 "parser.y" /* yacc.c:355  */

  
  #include <stdbool.h>
  #include "../math/units.h"
  #include "../structures/list.h"
  #include "../sensors/sensor.h"
  
  typedef struct Specification {
      
      char * id;
      List * options;
      List * args;
      
      bool should_destroy_options;
      bool should_destroy_args;
      
  } Specification;

#line 145 "y.tab.c" /* yacc.c:355  */

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    TRIGGER = 258,
    ID = 259,
    NUMERIC = 260
  };
#endif
/* Tokens.  */
#define TRIGGER 258
#define ID 259
#define NUMERIC 260

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 51 "parser.y" /* yacc.c:355  */

  char          * string;
  int           * integer;
  float         * decimal;
  Numeric       * numeric;
  
  List          * list;
  Charge        * charge;
  Trigger       * trigger;
  Specification * specification;

#line 179 "y.tab.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 196 "y.tab.c" /* yacc.c:358  */

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
#else
typedef signed char yytype_int8;
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
# elif ! defined YYSIZE_T
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
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
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
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
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
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
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
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  6
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   58

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  15
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  12
/* YYNRULES -- Number of rules.  */
#define YYNRULES  27
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  59

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   260

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    11,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    14,     8,
       9,     2,    10,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    12,     2,    13,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     6,     2,     7,     2,     2,     2,     2,
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
       5
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,    90,    90,    91,    94,    95,    98,    99,   102,   103,
     106,   107,   108,   111,   112,   113,   114,   117,   118,   121,
     122,   125,   126,   127,   130,   131,   134,   135
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "TRIGGER", "ID", "NUMERIC", "'{'", "'}'",
  "';'", "'<'", "'>'", "','", "'['", "']'", "':'", "$accept", "Config",
  "Sensors", "Sensor", "Specs", "Spec", "Actuator", "Charges", "Charge",
  "Tag", "Args", "Options", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   123,   125,    59,    60,
      62,    44,    91,    93,    58
};
# endif

#define YYPACT_NINF -40

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-40)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int8 yypact[] =
{
      -1,     4,    23,    -1,   -40,    10,   -40,   -40,     0,    24,
      29,    27,     1,   -40,    30,    34,    35,   -40,    22,    33,
     -40,   -40,    39,    40,   -40,    -4,   -40,    38,    38,   -40,
      41,     6,    21,   -40,   -40,   -40,    16,    43,    41,    36,
      -2,     7,   -40,    44,   -40,   -40,    19,   -40,   -40,    13,
      11,    45,   -40,   -40,   -40,   -40,    14,   -40,   -40
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       3,     0,     0,     2,     4,     0,     1,     5,     0,     0,
       0,     0,     0,     8,     0,     0,     0,     6,     0,     0,
       9,    12,     0,     0,    21,     0,     7,     0,     0,    26,
       0,     0,     0,    10,    11,    25,     0,     0,     0,    20,
       0,     0,    17,     0,    22,    27,     0,    19,    13,     0,
       0,     0,    24,    23,    14,    15,     0,    18,    16
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -40,   -40,   -40,    42,   -40,    46,    25,   -40,     3,    12,
      17,   -39
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     2,     3,     4,    12,    13,    33,    41,    42,    14,
      36,    31
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      29,    49,    29,     1,     9,     9,    48,    10,    19,     5,
      30,    56,    11,    11,    50,    29,     8,    37,    51,    55,
      38,    54,    58,     6,    37,    37,    39,    43,    40,    44,
      43,    18,    53,    15,    16,    24,    25,    17,    21,    22,
      23,    26,    27,    28,    32,     7,    35,    45,    11,    52,
      39,    47,     0,    34,    57,    46,     0,     0,    20
};

static const yytype_int8 yycheck[] =
{
       4,    40,     4,     4,     4,     4,     8,     7,     7,     5,
      14,    50,    12,    12,     7,     4,     6,    11,    11,     8,
      14,     8,     8,     0,    11,    11,     5,    11,     7,    13,
      11,     4,    13,     9,    10,    13,    14,     8,     8,     5,
       5,     8,     3,     3,     6,     3,     5,     4,    12,     5,
       5,    39,    -1,    28,    51,    38,    -1,    -1,    12
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     4,    16,    17,    18,     5,     0,    18,     6,     4,
       7,    12,    19,    20,    24,     9,    10,     8,     4,     7,
      20,     8,     5,     5,    13,    14,     8,     3,     3,     4,
      14,    26,     6,    21,    21,     5,    25,    11,    14,     5,
       7,    22,    23,    11,    13,     4,    25,    24,     8,    26,
       7,    11,     5,    13,     8,     8,    26,    23,     8
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    15,    16,    16,    17,    17,    18,    18,    19,    19,
      20,    20,    20,    21,    21,    21,    21,    22,    22,    23,
      23,    24,    24,    24,    25,    25,    26,    26
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     0,     1,     2,     5,     6,     1,     2,
       5,     5,     2,     3,     4,     4,     5,     1,     3,     2,
       1,     3,     6,     7,     3,     1,     1,     3
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

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
#ifndef YYINITDEPTH
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
static YYSIZE_T
yystrlen (const char *yystr)
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
static char *
yystpcpy (char *yydest, const char *yysrc)
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

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
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
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
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

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
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
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

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
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 3:
#line 91 "parser.y" /* yacc.c:1646  */
    { printf("\nExperiment does not involve sensors"); }
#line 1298 "y.tab.c" /* yacc.c:1646  */
    break;

  case 6:
#line 98 "parser.y" /* yacc.c:1646  */
    { build_sensor((yyvsp[-4].string), (yyvsp[-3].numeric), NULL);                      }
#line 1304 "y.tab.c" /* yacc.c:1646  */
    break;

  case 7:
#line 99 "parser.y" /* yacc.c:1646  */
    { build_sensor((yyvsp[-5].string), (yyvsp[-4].numeric),   (yyvsp[-2].list));                      }
#line 1310 "y.tab.c" /* yacc.c:1646  */
    break;

  case 8:
#line 102 "parser.y" /* yacc.c:1646  */
    { (yyval.list) = list_from(1, (yyvsp[0].specification));                           }
#line 1316 "y.tab.c" /* yacc.c:1646  */
    break;

  case 9:
#line 103 "parser.y" /* yacc.c:1646  */
    { list_insert((yyvsp[-1].list), (yyvsp[0].specification)); (yyval.list) = (yyvsp[-1].list);                    }
#line 1322 "y.tab.c" /* yacc.c:1646  */
    break;

  case 10:
#line 106 "parser.y" /* yacc.c:1646  */
    { (yyval.specification) = specify_trigger((yyvsp[-4].string),  true, (yyvsp[-2].numeric), (yyvsp[0].trigger));         }
#line 1328 "y.tab.c" /* yacc.c:1646  */
    break;

  case 11:
#line 107 "parser.y" /* yacc.c:1646  */
    { (yyval.specification) = specify_trigger((yyvsp[-4].string), false, (yyvsp[-2].numeric), (yyvsp[0].trigger));         }
#line 1334 "y.tab.c" /* yacc.c:1646  */
    break;

  case 12:
#line 108 "parser.y" /* yacc.c:1646  */
    { (yyval.specification) = (yyvsp[-1].specification);                                         }
#line 1340 "y.tab.c" /* yacc.c:1646  */
    break;

  case 13:
#line 111 "parser.y" /* yacc.c:1646  */
    { printf("Triggers must list pins\n"); exit(3);    }
#line 1346 "y.tab.c" /* yacc.c:1646  */
    break;

  case 14:
#line 112 "parser.y" /* yacc.c:1646  */
    { printf("Triggers must list pins\n"); exit(3);    }
#line 1352 "y.tab.c" /* yacc.c:1646  */
    break;

  case 15:
#line 113 "parser.y" /* yacc.c:1646  */
    { (yyval.trigger) = make_trigger((yyvsp[-2].list), NULL);                     }
#line 1358 "y.tab.c" /* yacc.c:1646  */
    break;

  case 16:
#line 114 "parser.y" /* yacc.c:1646  */
    { (yyval.trigger) = make_trigger((yyvsp[-3].list),   (yyvsp[-1].list));                     }
#line 1364 "y.tab.c" /* yacc.c:1646  */
    break;

  case 17:
#line 117 "parser.y" /* yacc.c:1646  */
    { (yyval.list) = list_from(1, (yyvsp[0].charge));                           }
#line 1370 "y.tab.c" /* yacc.c:1646  */
    break;

  case 18:
#line 118 "parser.y" /* yacc.c:1646  */
    { list_insert((yyvsp[-2].list), (yyvsp[0].charge)); (yyval.list) = (yyvsp[-2].list);                    }
#line 1376 "y.tab.c" /* yacc.c:1646  */
    break;

  case 19:
#line 121 "parser.y" /* yacc.c:1646  */
    { (yyval.charge) = make_charge((yyvsp[-1].numeric),   (yyvsp[0].specification));                      }
#line 1382 "y.tab.c" /* yacc.c:1646  */
    break;

  case 20:
#line 122 "parser.y" /* yacc.c:1646  */
    { (yyval.charge) = make_charge((yyvsp[0].numeric), NULL);                      }
#line 1388 "y.tab.c" /* yacc.c:1646  */
    break;

  case 21:
#line 125 "parser.y" /* yacc.c:1646  */
    { (yyval.specification) = make_tag((yyvsp[-1].string), NULL, NULL);                   }
#line 1394 "y.tab.c" /* yacc.c:1646  */
    break;

  case 22:
#line 126 "parser.y" /* yacc.c:1646  */
    { (yyval.specification) = make_tag((yyvsp[-4].string), NULL,   (yyvsp[-1].list));                   }
#line 1400 "y.tab.c" /* yacc.c:1646  */
    break;

  case 23:
#line 127 "parser.y" /* yacc.c:1646  */
    { (yyval.specification) = make_tag((yyvsp[-5].string),   (yyvsp[-3].list),   (yyvsp[-1].list));                   }
#line 1406 "y.tab.c" /* yacc.c:1646  */
    break;

  case 24:
#line 130 "parser.y" /* yacc.c:1646  */
    { list_insert((yyvsp[-2].list), (yyvsp[0].numeric)); (yyval.list) = (yyvsp[-2].list);                    }
#line 1412 "y.tab.c" /* yacc.c:1646  */
    break;

  case 25:
#line 131 "parser.y" /* yacc.c:1646  */
    { (yyval.list) = list_from(1, (yyvsp[0].numeric));                           }
#line 1418 "y.tab.c" /* yacc.c:1646  */
    break;

  case 26:
#line 134 "parser.y" /* yacc.c:1646  */
    { (yyval.list) = list_from(1, (yyvsp[0].string));                           }
#line 1424 "y.tab.c" /* yacc.c:1646  */
    break;

  case 27:
#line 135 "parser.y" /* yacc.c:1646  */
    { list_insert((yyvsp[-2].list), (yyvsp[0].string)); (yyval.list) = (yyvsp[-2].list);                    }
#line 1430 "y.tab.c" /* yacc.c:1646  */
    break;


#line 1434 "y.tab.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
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

  /* Else will try to reuse lookahead token after shifting the error
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

  /* Do not reclaim the symbols of the rule whose action triggered
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
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
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

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


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

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
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
  return yyresult;
}
#line 138 "parser.y" /* yacc.c:1906  */


static void represent_as_decimal(Numeric * numeric) {
    
    if (numeric -> is_decimal) return;
    
    numeric -> decimal    = (float) numeric -> integer;
    numeric -> is_decimal = true;
    
    if (!strcmp(numeric -> units, "i"))    // change generic unit if needed
      numeric -> units[0] = 'f';           // -----------------------------
}

static void specification_destroy(void * vspecification) {
    
    Specification * specification = vspecification;
    
    /*if (specification -> options && specification -> should_destroy_options)
        list_destroy(specification -> options);
    
    if (specification -> args && specification -> should_destroy_args)
    list_destroy(specification -> args);*/
    
    free(specification -> id);
    free(specification);
}

Charge * make_charge(Numeric * wire, Specification * tag) {
    /* note, pin 0 not allowed. Shouldn't matter since that's SDA.0, *
     * but it's important for future bug potential.                  */
  
    if (strcmp(wire -> units, "i")) {
        printf(RED "Charge wires must be integers" RESET);
	exit(ERROR_EXPERIMENTER);        
    }
    
    Charge * charge = calloc(1, sizeof(*charge));
    
    charge -> gpio = abs(wire -> integer);
    charge -> hot  = (wire -> integer > 0);
    
    if (tag) {
        if (strcmp(tag -> id, "pulse")) {
            printf(RED "Unknown tag %s for charge\n" RESET, tag -> id);
            exit(ERROR_EXPERIMENTER);
        }
        
        if (tag -> args -> size != 1) {
            printf(RED "Exactly 1 argument shoud be used for 'pulse' tag\n" RESET);
            exit(ERROR_EXPERIMENTER);
        }
        
        charge -> duration = (int) tag -> args -> head -> value;
        
        specification_destroy(tag);    // no longer needed
    }
    
    return charge;
}

Trigger * make_trigger(List * charges, List * options) {
    // creates a trigger, which becomes a specification later
    
    Trigger * trigger = calloc(1, sizeof(*trigger));
    
    trigger -> fired    = false;
    trigger -> singular = true;     // defaults
    trigger -> reverses = false;    // --------
    
    trigger -> charges = charges;
    
    if (!options) return trigger;
    
    for (iterate(options, char *, option)) {
        if      (!strcmp(option, "singular")) trigger -> singular = true;
        else if (!strcmp(option, "forever" )) trigger -> singular = false;
        else if (!strcmp(option, "reverses")) trigger -> reverses = true;
        else {
            printf("Unknown option " RED "%s\n" RESET, option);
            exit(ERROR_EXPERIMENTER);
        }
    }
    
    return trigger;
}


Specification * specify_trigger(char * id, bool less, Numeric * threshold, Trigger * trigger) {
    /* Wraps and modifies the trigger into a specification so that sensor construction *
     * may have a linked list consisting of nodes with the same content                */
    
    represent_as_decimal(threshold);
    
    trigger -> id        = id;
    trigger -> less      = less;
    trigger -> threshold = threshold;
    
    Specification * specification = malloc(sizeof(*specification));
    
    specification -> id      = strdup("trigger");
    specification -> options = list_from(1, trigger);
    
    return specification;
}

Specification * make_tag(char * id, List * options, List * args) {
  
    Specification * tag = malloc(sizeof(*tag));
    
    tag -> id      = id;
    tag -> options = options;
    tag -> args    = args;
    
    
    // [ print, color : 5hz ] - no color or freq required
    // [ smooth : ] -
    // [ calibrate target : 1.0, 3.4, ... ]
    
    
    
    /* Perform error checking on this tag now */
    
    if (!strcmp(id, "print")) {
      
        if (options) {
	  
            char * color      = (char *) options -> head -> value;
            char * color_code = get_color_by_name(color);
            
            if (!color_code) {
                printf(RED "Color name %s not recognized\n" RESET, color);
                exit(ERROR_EXPERIMENTER);
            }
        }
    }
    
    else if (!strcmp(id, "pulse")) {
        
        if (options)
            exit_printing("Pulsing does not support options at this time\n", ERROR_EXPERIMENTER);
	
        if (!args || args -> size != 1)
            exit_printing("Pulsing requires exactly 1 argument\n", ERROR_EXPERIMENTER);
    }
    
    else if (!strcmp(id, "smooth")) {
      
        if (options)
            exit_printing("Smoothing does not support options at this time\n", ERROR_EXPERIMENTER);
	
        if (args) {
	    Numeric * numeric = (Numeric *) list_get(args, 0);
            float value = numeric -> decimal;
	    
            if (strcmp(numeric -> units, "f") || value < 0.0f || value > 1.0f) {
                printf(RED "Autoregressive smoothing requires a # in [0.0, 1.0]\n" RESET);
                exit(ERROR_EXPERIMENTER);
            }
        }
    }
    
    else if (!strcmp(id, "calibrate")) {
      
        if (!args   ) exit_printing("Calibration curves require at least one constant\n", ERROR_EXPERIMENTER);
        if (!options) exit_printing("Calibration requires a target\n", ERROR_EXPERIMENTER);
	if (options -> size < 2) exit_printing("Calibration requires a unit\n", ERROR_EXPERIMENTER);
	
        char * curve = (char *) list_get(options, 1);
        
        if      (!strcmp(curve, "poly"));
        else if (!strcmp(curve, "hart")) {
            if (args -> size != 3)
	        exit_printing("The Steinhart and Hart Equation requires 3 constants\n", ERROR_EXPERIMENTER);
        }
        else {
            printf(RED "Unknown calibration curve " CYAN "%s\n" RESET, curve);
            exit(ERROR_EXPERIMENTER);
        }
    }
    
    else {
        printf(RED "Invalid tag id " CYAN "%s\n" RESET, id);
        exit(ERROR_EXPERIMENTER);
    }
    
    return tag;
}

void build_sensor(char * id, Numeric * frequency, List * specifications) {
  
  ProtoSensor * proto = hashmap_get(proto_sensors, id);
  
  if (!proto) {
    printf(RED "%s is not a sensor\n" RESET, id);
    exit(ERROR_EXPERIMENTER);
  }
  
  proto -> hertz = frequency -> integer;    // only support Hz for now
  proto -> requested = true;
  
  if (!specifications) return;
  
  List * triggers = list_create();
  
  for (iterate(specifications, Specification *, specification)) {
      
      if (!strcmp(specification -> id, "trigger")) {
	  list_insert(triggers, list_get(specification -> options, 0));
      }
      
      else if (!strcmp(specification -> id, "print")) {
          proto -> print = true;
      }
      
      else if (!strcmp(specification -> id, "smooth")) {
	  
          Numeric * numeric = list_get(specification -> args, 0);
	  
          float auto_regressive = numeric -> decimal;
          
          if (proto -> auto_regressive) {
              printf(RED "Duplicate autoregressive constant for %s\n" RESET, id);
              exit(ERROR_EXPERIMENTER);
          }
          
          proto -> auto_regressive = auto_regressive;
      }
      
      else if (!strcmp(specification -> id, "calibrate")) {
          
          char * target = list_get(specification -> options, 0);
	  char * curve  = list_get(specification -> options, 1);
	  char * unit   = list_get(specification -> options, 2);
	  
          if (!proto -> targets || !hashmap_exists(proto -> targets, target)) {
              printf(RED "Calibration target %s for %s is not known\n" RESET, target, id);
              exit(ERROR_EXPERIMENTER);
          }
	  
          if (hashmap_exists(proto -> calibrations, target)) {
              printf(RED "Duplicate target calibration for %s of %s found\n" RESET, target, id);
              exit(ERROR_EXPERIMENTER);
          }
	  
	  // force all calibration constants to be decimals
	  for (iterate(specification -> args, Numeric *, constant))
	    represent_as_decimal(constant);
	  
	  // make first node the calibration type
	  list_insert_first(specification -> args, curve);
	  
	  hashmap_add(proto -> calibrations, target, specification -> args);
	  hashmap_add(proto -> output_units, target,                  unit);
      }
  }
  
  
  
  for (iterate(triggers, Trigger *, trigger)) {
      
      // check trigger targets
      if (!proto -> targets || !hashmap_exists(proto -> targets, trigger -> id)) {
          printf(RED "Trigger target %s unknown\n" RESET, trigger -> id);
          exit(ERROR_EXPERIMENTER);
      }
      
      // duplicate reversing triggers
      
      if (trigger -> reverses) {
	
	  Trigger * opposite = calloc(1, sizeof(*opposite));
	  
	  opposite -> id       =  strdup(trigger -> id);
	  opposite -> less     = !trigger -> less;
	  opposite -> fired    =  trigger -> fired;
	  opposite -> singular =  trigger -> singular;
	  opposite -> reverses =  false;
	  	  
	  opposite -> threshold = trigger -> threshold;
	  opposite -> charges   = list_create();
	  
	  for (iterate(trigger -> charges, Charge *, charge)) {
	      
              Charge * anti_charge = malloc(sizeof(*charge));
	      
	      anti_charge -> gpio     =  charge -> gpio;
	      anti_charge -> hot      = !charge -> hot;
	      anti_charge -> duration =  charge -> duration;
	      
	      list_insert(opposite -> charges, anti_charge);
	  }
	  
	  list_insert(triggers, opposite);
      }
  }
  
  proto -> triggers = triggers;
  
  specifications -> free = specification_destroy;
  list_destroy(specifications);
  
  n_triggers += triggers -> size;
}

void print_config() {
  
  if (!n_triggers) {
    printf("\n\nExperiment does not involve triggers\n\n");
    return;
  }
  
  printf("\n\n" GRAY "%d" RESET " Triggers\n\n", n_triggers);
  
  for (iterate(proto_sensors -> all, ProtoSensor *, proto)) {
    
    if (!proto -> requested) continue;
    if (!proto -> triggers ) continue;
    
    printf(GREEN "%s\n" RESET, proto -> code_name);
    
    for (iterate(proto -> triggers, Trigger *, trigger)) {
      
      Numeric * numeric = trigger -> threshold;
      
      Numeric standard;
      to_standard_units(&standard, numeric);
      
      printf(CYAN "    %s" RESET, trigger -> id);
      
      if (trigger -> less)
	printf(GRAY " <" MAGENTA " %.3f%s " GRAY "(= " MAGENTA "%.3f%s" GRAY ")",
	       numeric -> decimal, numeric -> units, standard.decimal, standard.units);
      else
	printf(GRAY " >" MAGENTA " %.3f%s " GRAY "(= " MAGENTA "%.3f%s" GRAY ")",
	       numeric -> decimal, numeric -> units, standard.decimal, standard.units);
      
      
      printf(" { " RESET);
      
      for (iterate(trigger -> charges, Charge *, charge)) {
        
	if (charge -> hot) printf(YELLOW "+" RESET);
	else               printf(YELLOW "-" RESET);
	
	printf("%d ", charge -> gpio);
      }
      
      if (!trigger -> singular) printf(GRAY "}" YELLOW " *\n" RESET);
      else                      printf(GRAY "}\n");
    }
    
    printf("\n");
  }
}

/* Easel's foundation.
 *
 * Contents:
 *    1. Exception and fatal error handling.
 *    2. Memory allocation/deallocation conventions.
 *    3. Standard banner for Easel miniapplications.
 *    4. Improved replacements for some C library functions.
 *    5. Portable drop-in replacements for nonstandard C functions.
 *    6. Additional string functions, esl_str*()
 *    7. File path/name manipulation, including tmpfiles.
 *    8. Typed comparison functions.
 *    9. Unit tests.
 *   10. Test driver.
 *   11. Examples.
 *   12. Copyright and license.
 */

/* Easel */

#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>		/* for FILE */
#include <stdarg.h>		/* for va_list */
#include <math.h>               /* for HUGE_VAL */
#ifdef HAVE_STDINT_H
#include <stdint.h>		/* for uint32_t and the like (C99) */
#endif
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>		/* some systems allegedly put uints here */
#endif

/*** Start of inlined file: esl_stopwatch.h ***/
#ifndef eslSTOPWATCH_INCLUDED
#define eslSTOPWATCH_INCLUDED

#include <time.h>
#ifdef HAVE_TIMES
#include <sys/times.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* need for sysconf() */
#endif

typedef struct {
#ifdef eslSTOPWATCH_HIGHRES
  double     t0;                /* baseline wall time from Nadeau routine */
#elif  HAVE_TIMES
  clock_t    t0;		/* baseline wall time, POSIX times()      */
#else
  time_t     t0;                /* baseline wall time from ANSI time()    */
#endif

#ifdef HAVE_TIMES
  struct tms cpu0;		/* baseline CPU/system time, POSIX times()      */
#else
  clock_t cpu0;			/* baseline CPU time, fallback to ANSI clock()  */
#endif

  /* elapsed/user/sys are t-t0 results for the last time the
   * watch was Stop()'ed.
   */
  double elapsed;               /* elapsed wall time, seconds */
  double user;                  /* CPU time, seconds          */
  double sys;                   /* system time, seconds       */
} ESL_STOPWATCH;

ESL_STOPWATCH *esl_stopwatch_Create(void);
void           esl_stopwatch_Destroy(ESL_STOPWATCH *w);

int esl_stopwatch_Start(ESL_STOPWATCH *w);
int esl_stopwatch_Stop(ESL_STOPWATCH *w);
int esl_stopwatch_Display(FILE *fp, ESL_STOPWATCH *w, char *prefix);

double esl_stopwatch_GetElapsed(ESL_STOPWATCH *w);

int esl_stopwatch_Include(ESL_STOPWATCH *master, ESL_STOPWATCH *w);

#endif /*eslSTOPWATCH_INCLUDED*/

/*** Start of inlined file: easel.h ***/
#ifndef eslEASEL_INCLUDED
#define eslEASEL_INCLUDED

/*** End of inlined file: esl_config.h ***/



/*****************************************************************
 * 1. Macros implementing Easel's error handling conventions
 *****************************************************************/
/* Many objects contain a fixed length "errbuf" for failure
 * diagnostics: ESL_FAIL() and ESL_XFAIL() fill this buffer.
 */
#define eslERRBUFSIZE 128

/* ESL_FAIL()       - return an error message, without cleanup.
 * ESL_XFAIL()      - return an error message, with cleanup.
 * ESL_EXCEPTION()  - throwing an exception, without cleanup.
 * ESL_XEXCEPTION() - throwing an exception, with cleanup.
 *
 * The X versions (with cleanup) require the caller to have an
 * <int status> variable and a <ERROR:> goto target in scope,
 * which, yes, is a little hacky.
 *
 * Wrapping these macros in <while(0)> loops allows a statement:
 *       if (something) ESL_XEXCEPTION(code,mesg);
 * without the trailing semicolon becoming a null statement after
 * macro expansion.
 *
 * All esl_fail() does is vsnprintf() to the <errbuf>; the reason to
 * have ESL_FAIL and ESL_XFAIL call esl_fail() is to enable us to set
 * a debugging breakpoint in esl_fail(), so we can break execution at
 * a normal failure.
 */
/*::cexcerpt::error_macros::begin::*/
#define ESL_FAIL(code, errbuf, ...) do {				\
	esl_fail(errbuf, __VA_ARGS__);                                      \
	return code; }							\
  while (0)

#define ESL_XFAIL(code, errbuf, ...) do {				\
	status = code;							\
	esl_fail(errbuf, __VA_ARGS__);                                      \
	goto ERROR; }							\
  while (0)

#define ESL_EXCEPTION(code, ...)  do {					\
	esl_exception(code, FALSE, __FILE__, __LINE__, __VA_ARGS__);	\
	return code; }							\
  while (0)

#define ESL_XEXCEPTION(code, ...)  do {					\
	status = code;							\
	esl_exception(code, FALSE, __FILE__, __LINE__, __VA_ARGS__);	\
	goto ERROR; }							\
  while (0)

#define ESL_EXCEPTION_SYS(code, ...) do {				\
	esl_exception(code, TRUE, __FILE__, __LINE__, __VA_ARGS__);		\
	return code; }							\
  while (0)

#define ESL_XEXCEPTION_SYS(code, ...)  do {				\
	status = code;							\
	esl_exception(code, TRUE, __FILE__, __LINE__, __VA_ARGS__);	\
	goto ERROR; }							\
  while (0)
/*::cexcerpt::error_macros::end::*/

/* Return codes for error handler
 */
/*::cexcerpt::statuscodes::begin::*/
#define eslOK              0    /* no error/success             */
#define eslFAIL            1    /* failure                      */
#define eslEOL             2    /* end-of-line (often normal)   */
#define eslEOF             3    /* end-of-file (often normal)   */
#define eslEOD             4    /* end-of-data (often normal)   */
#define eslEMEM            5    /* malloc or realloc failed     */
#define eslENOTFOUND       6    /* file or key not found        */
#define eslEFORMAT         7    /* file format not correct      */
#define eslEAMBIGUOUS      8    /* an ambiguity of some sort    */
#define eslEDIVZERO        9    /* attempted div by zero        */
#define eslEINCOMPAT      10    /* incompatible parameters      */
#define eslEINVAL         11    /* invalid argument/parameter   */
#define eslESYS           12    /* generic system call failure  */
#define eslECORRUPT       13    /* unexpected data corruption   */
#define eslEINCONCEIVABLE 14    /* "can't happen" error         */
#define eslESYNTAX        15    /* invalid user input syntax    */
#define eslERANGE         16    /* value out of allowed range   */
#define eslEDUP           17    /* saw a duplicate of something */
#define eslENOHALT        18    /* a failure to converge        */
#define eslENORESULT      19    /* no result was obtained       */
#define eslENODATA        20    /* no data provided, file empty */
#define eslETYPE          21    /* invalid type of argument     */
#define eslEOVERWRITE     22    /* attempted to overwrite data  */
#define eslENOSPACE       23    /* ran out of some resource     */
#define eslEUNIMPLEMENTED 24    /* feature is unimplemented     */
#define eslENOFORMAT      25	/* couldn't guess file format   */
#define eslENOALPHABET    26	/* couldn't guess seq alphabet  */
#define eslEWRITE         27   	/* write failed (fprintf, etc)  */
#define eslEINACCURATE    28    /* return val may be inaccurate */
/*::cexcerpt::statuscodes::end::*/

/*****************************************************************
 * 2. Macros implementing Easel's memory allocation conventions
 *****************************************************************/
/* ESL_ALLOC(), ESL_RALLOC():
 *
 * Allocation and reallocation wrappers.
 * Both require <int status> in scope, and <ERROR:> goto target.
 * ESL_RALLOC() also requires <void *> ptr to be provided as <tmp>.
 *
 * ESL_REALLOC() is a newer version of ESL_RALLOC() which doesn't
 * need a tmp ptr. All ESL_RALLOC() calls can be safely converted
 * to ESL_REALLOC() calls.
 *
 * The result of malloc(0) is implementation-defined (either NULL or
 * a ptr that may not be dereferenced), a bit of a hole in the C
 * standard. In Easel, we want to avoid having NULL as a valid
 * non-error result of malloc(), because it confuses static analysis
 * tools when they see dereferences of possibly NULL pointers. We
 * therefore treat malloc(0) as an eslEMEM error.
 */
/*::cexcerpt::alloc_macros::begin::*/
#define ESL_ALLOC(p, size) do {\
	if ( size <= 0 ) { \
	   p = NULL; \
	   status = eslEMEM; \
	   esl_exception(status, FALSE, __FILE__, __LINE__, "zero malloc disallowed"); \
	   goto ERROR;\
	}\
	if ( ((p) = malloc(size)) == NULL)  { \
	   status = eslEMEM;\
	   esl_exception(status, FALSE, __FILE__, __LINE__, "malloc of size %d failed", size); \
	   goto ERROR;\
	 }} while (0)

#define ESL_RALLOC(p, tmp, newsize) do {\
	 if ((p) == NULL) { (tmp) = malloc(newsize);         }\
	 else             { (tmp) = realloc((p), (newsize)); }\
	 if ((tmp) != NULL) (p) = (tmp);\
	 else {\
	   status = eslEMEM;\
	   esl_exception(status, FALSE, __FILE__, __LINE__, "realloc for size %d failed", newsize);	\
	   goto ERROR;\
	 }} while (0)

#define ESL_REALLOC(p, newsize) do {\
	 void *esltmpp;\
	 if ((p) == NULL) { (esltmpp) = malloc(newsize);         }\
	 else             { (esltmpp) = realloc((p), (newsize)); }\
	 if ((esltmpp) != NULL) (p) = (esltmpp);\
	 else {\
	   status = eslEMEM;\
	   esl_exception(status, FALSE, __FILE__, __LINE__, "realloc for size %d failed", newsize); \
	   goto ERROR;\
	 }} while (0)
/*::cexcerpt::alloc_macros::end::*/

/* Convert MiB,GiB,TiB to bytes, using binary definitions (2^20, 2^30, 2^40):
 * Pedantically speaking, that's: mebibytes (MiB), gibibytes (GiB), tebibytes (TiB).
 * 1 TB = 10^12 bytes; 1 TiB = 2^40 bytes.
 */
#define ESL_MBYTES(x) ((x) * 1048576)
#define ESL_GBYTES(x) ((x) * 1024 * 1048576)
#define ESL_TBYTES(x) ((x) * 1024 * 1024 * 1048576)

/* Round integer <n> up to the nearest multiple of <m>.
 * Particularly useful when dealing w/ memory alignment issues.
 */
#define ESL_UPROUND(n, m)  ( ((n) + (m)-1) / (m) * (m))

/*****************************************************************
 * 3. Macros implementing Easel's function argument conventions
 *****************************************************************/

#define esl_byp_IsInternal(p) ((p) == NULL)
#define esl_byp_IsReturned(p) ((p) != NULL && (*p) == NULL)
#define esl_byp_IsProvided(p) ((p) != NULL && (*p) != NULL)

/* Sometimes a shared function API dictates arguments that a function
 * doesn't use, and we want to silence compiler warnings about this.
 * Putting ESL_UNUSED(x) in the function, for an unused argument <x>,
 * should silence the compiler, and should generate a no-op.
 */
#define ESL_UNUSED(x) (void)(sizeof((x)))

/*****************************************************************
 * 4. Macros implementing Easel's debugging output conventions
 *****************************************************************/
/* Debugging hooks, w/ three levels (1-3).
 */
#if eslDEBUGLEVEL >= 1		/* for ESL_DASSERT() macros */
#include <assert.h>
#endif

#if (eslDEBUGLEVEL >= 1)
#define ESL_DPRINTF1(x)  printf x
#define ESL_DASSERT1(x)  assert x
#else
#define ESL_DPRINTF1(x)
#define ESL_DASSERT1(x)
#endif
#if (eslDEBUGLEVEL >= 2)
#define ESL_DPRINTF2(x)  printf x
#define ESL_DASSERT2(x)  assert x
#else
#define ESL_DPRINTF2(x)
#define ESL_DASSERT2(x)
#endif
#if (eslDEBUGLEVEL >= 3)
#define ESL_DPRINTF3(x)  printf x
#define ESL_DASSERT3(x)  assert x
#else
#define ESL_DPRINTF3(x)
#define ESL_DASSERT3(x)
#endif

/*****************************************************************
 * 5. Defined constants
 *****************************************************************/

/* Making sure TRUE/FALSE are defined, for convenience */
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* Some basic mathematical constants.
 * Assuming IEEE754 math with 64-bit doubles (53-bit mantissas), we
 * want 17 significant decimal digits in our constants. More is
 * a waste (but we do it for some anyway).
 */
#define eslCONST_E     2.71828182845904523536028747135
#define eslCONST_PI    3.14159265358979323846264338328
#define eslCONST_EULER 0.57721566490153286060651209008
#define eslCONST_GOLD  1.61803398874989484820458683437
#define eslCONST_LOG2  0.69314718055994529
#define eslCONST_LOG2R 1.44269504088896341

/* Define <eslINFINITY>, <eslNaN> portably.
 * ANSI C99 makes it easy (at last!). The failovers to other pre-C99
 * methods are legacy now that we require a C99 compiler - but no harm
 * leaving them in Just In Case.
 */
#if defined (INFINITY)
#define eslINFINITY    INFINITY  /* C99 */
#elif defined (HUGE_VAL)
#define eslINFINITY    HUGE_VAL	 /* assume IEEE754 HUGE_VAL = infinity. ok? */
#else
#define eslINFINITY    (1.0/0.0) /* portable? */
#endif

#if defined (NAN)
#define eslNaN         NAN	/* C99 */
#else
#define eslNaN         (eslINFINITY/eslINFINITY) /* portably make a IEEE754 NaN */
#endif

/* Define crossovers for numerical approximations.
 */
/* log(1+x) ~ x and  1-e^x = -x approximation.
 * Same threshold appears to be optimal for float or double x. xref STL9/138.
 */
#define eslSMALLX1    5e-9

/*****************************************************************
 * 6. Basic support for Easel's digitized biosequences.
 *****************************************************************/

/* Most of this support is in the alphabet module, but we externalize
 * some into the easel foundation because ESL_INMAP is used in unaugmented
 * sqio, msa modules.
 *
 * A digital sequence residue (ESL_DSQ) is an unsigned 8-bit type
 * (0..255).  A valid digital residue has a value in the range 0..127
 * (Easel can represent alphabets of up to 128 different characters).
 * Values 128..255 are reserved for flags.
 *
 * An "inmap" is ESL_DSQ[128], or *ESL_DSQ allocated for 128 values;
 * it is a many-to-one construct for mapping 7-bit ASCII chars (in
 * range 0..127) either to new ASCII chars (in the case of raw
 * sequence input in sqio, msa) or to digital codes (in the alphabet
 * module).  Valid mapped values are 0..127; any value in range
 * 128..255 is some kind of flag.
 */
typedef uint8_t ESL_DSQ;
#define eslDSQ_SENTINEL 255	/* sentinel bytes 0,L+1 in a dsq */
#define eslDSQ_ILLEGAL  254	/* input symbol is unmapped and unexpected */
#define eslDSQ_IGNORED  253     /* input symbol is unmapped and ignored */
#define eslDSQ_EOL      252	/* input symbol marks end of a line */
#define eslDSQ_EOD      251     /* input symbol marks end of a seq record */

/* If you try to test sym > 0 && sym <= 127 below, instead of isascii(sym),
 * you'll get a compiler warning for an always-successful test regardless
 * of whether a char is signed or unsigned. So we trust that isascii() is
 * doing the Right Thing.
 */
#define esl_inmap_IsValid(inmap, sym)  (isascii(sym) && (inmap)[(int)sym] <= 127)

/*****************************************************************
 * 7. Miscellaneous.
 *****************************************************************/
/* A placeholder for helping w/ portability of filenames/paths.
 * I think, but have not tested, that:
 *   VMS:            #define DIRSLASH ']'
 *   ancient MacOS:  #define DIRSLASH ':'
 *   DOS:            #define DIRSLASH '\\'
 * Setting DIRSLASH correctly is probably not the only thing
 * that would need to be done to port to other OS's, but it's
 * probably a start.
 *
 * The code assumes that '.' is used for file name extensions,
 * such as "foo.bar".
 *
 * This gets used in easel.c's *_File*() functions.
 */
#define eslDIRSLASH '/'           /* UNIX directory paths have /foo/bar */

/* Some generic macros for swapping, min, and max.
 */
#define ESL_SWAP(x, y, type)  do { type esltmpxyz = (x); (x) = (y); (y) = esltmpxyz; } while (0)
#define ESL_MIN(a,b)          (((a)<(b))?(a):(b))
#define ESL_MAX(a,b)          (((a)>(b))?(a):(b))

static inline float esl_log (double x) { return (x == 0.0 ? -eslINFINITY : log(x));  } /* avoid fp exceptions; log(0) = -inf is fine */
static inline float esl_logf(float x)  { return (x == 0.0 ? -eslINFINITY : logf(x)); }
static inline float esl_log2f(float x) { return (x == 0.0 ? -eslINFINITY : eslCONST_LOG2R * logf(x)); }

/* Typedef: <esl_pos_t>
 *
 * <esl_pos_t> is a signed integer type suitable for safe casting
 * to EITHER an <off_t> or <size_t> on this system (i.e. as a position
 * in memory or in a file on disk), where we may use a -1 as a flag
 * (or even other negative numbers).
 *
 * <esl_pos_t> is for use for anything having to do with positions in
 * large buffers, strings, or files, where we want strict control
 * of integer range limits.
 *
 * Note that POSIX requires size_t to be unsigned, and off_t to be
 * signed.
 */
typedef int64_t esl_pos_t;

/* ESL_ANALYZER_NORETURN
 *    adds some optional support for clang static analysis.
 *    The static analyzer sometimes needs to be clued in when a
 *    function cannot return: fatal error handlers, for example.
 *    clang, gcc, and other gcc-like compilers support the __attribute__
 *    extension on function declarations. We detect this support
 *    at compile-time in the configure script. Functions that
 *    don't return are declared like:
 *       void fatal(char *msg, ...) ESL_ANALYZER_NORETURN;
 */
#ifndef ESL_ANALYZER_NORETURN
#ifdef  HAVE_FUNC_ATTRIBUTE_NORETURN
#define ESL_ANALYZER_NORETURN __attribute__((__noreturn__))
#else
#define ESL_ANALYZER_NORETURN
#endif
#endif

/*****************************************************************
 * 8. Void declarations of missing augmentations
 *****************************************************************/
#ifndef eslAUGMENT_ALPHABET
typedef void ESL_ALPHABET;
#endif
#ifndef eslAUGMENT_KEYHASH
typedef void ESL_KEYHASH;
#endif

/*****************************************************************
 * 9. The API declarations for easel.c
 *****************************************************************/

/* 1. Error handling. */
typedef void (*esl_exception_handler_f)(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, va_list argp);
void esl_fail(char *errbuf, const char *format, ...);
void esl_exception(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, ...);
void esl_exception_SetHandler(esl_exception_handler_f);
void esl_exception_ResetDefaultHandler(void);
void esl_nonfatal_handler(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, va_list argp);
void esl_fatal(const char *format, ...) ESL_ANALYZER_NORETURN;

/* 2. Memory allocation/deallocation conventions. */
void esl_Free2D(void  **p, int dim1);
void esl_Free3D(void ***p, int dim1, int dim2);

/* 3. Standard banner for Easel miniapplications. */
int  esl_banner    (FILE *fp, char *progname, char *banner);
int  esl_usage     (FILE *fp, char *progname, char *usage);
int  esl_dataheader(FILE *fp, ...);

/* 4. Improved replacements for some C library functions */
int  esl_fgets(char **buf, int *n, FILE *fp);
int  esl_strdup(const char *s, int64_t n, char **ret_dup);
int  esl_strcat(char **dest, int64_t ldest, const char *src, int64_t lsrc);
int  esl_strmapcat        (const ESL_DSQ *inmap, char **dest, int64_t *ldest, const char *src, esl_pos_t lsrc);
int  esl_strmapcat_noalloc(const ESL_DSQ *inmap,  char *dest, int64_t *ldest, const char *src, esl_pos_t lsrc);
int  esl_strtok    (char **s, char *delim, char **ret_tok);
int  esl_strtok_adv(char **s, char *delim, char **ret_tok, int *opt_toklen, char *opt_endchar);
int  esl_sprintf (char **ret_s, const char *format, ...);
int  esl_vsprintf(char **ret_s, const char *format, va_list *ap);
int  esl_strcmp(const char *s1, const char *s2);

/* 5.  Portable drop-in replacements for non-standard C functions */
#ifndef HAVE_STRCASECMP
#ifdef _MSC_VER
#define strcasecmp stricmp
#else
int  esl_strcasecmp(const char *s1, const char *s2);
#define strcasecmp esl_strcasecmp
#endif
#endif

/* 6. Additional string functions, esl_str*() */
int     esl_strchop(char *s, int64_t n);
int     esl_strdealign(char *s, const char *aseq, const char *gapchars, int64_t *opt_rlen);
int     esl_str_IsBlank(char *s);
int     esl_str_IsInteger(char *s);
int     esl_str_IsReal(char *s);
int64_t esl_str_GetMaxWidth(char **s, int n);

/* 7. File path/name manipulation functions, including tmpfiles */
int  esl_FileExists(const char *filename);
int  esl_FileTail(const char *path, int nosuffix, char **ret_file);
int  esl_file_Extension(char *filename, esl_pos_t n_ignore, char **ret_sfx, esl_pos_t *ret_n);
int  esl_FileConcat(const char *dir, const char *file, char **ret_path);
int  esl_FileNewSuffix(const char *filename, const char *sfx, char **ret_newpath);
int  esl_FileEnvOpen(const char *fname, const char *env,
			    FILE **ret_fp, char **ret_path);
int  esl_tmpfile(char *basename6X, FILE **ret_fp);
int  esl_tmpfile_named(char *basename6X, FILE **ret_fp);
int  esl_getcwd(char **ret_cwd);

/* 8. Typed comparison routines. */
int  esl_DCompare   (double a, double b, double tol);
int  esl_FCompare   (float  a, float  b, float  tol);
int  esl_DCompareAbs(double a, double b, double tol);
int  esl_FCompareAbs(float  a, float  b, float  tol);
int  esl_CCompare(char *s1, char *s2);

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/*** End of inlined file: easel.h ***/


/*** Start of inlined file: esl_alphabet.h ***/
#ifndef eslALPHABET_INCLUDED
#define eslALPHABET_INCLUDED

#include <ctype.h>		/* isascii() */

/* Flags for alphabet types.
 * Do not change, only add, because these codes are used in file formats.
 */
#define eslUNKNOWN     0        /* 0=unknown is easel-wide convention; don't change */
#define eslRNA         1
#define eslDNA         2
#define eslAMINO       3
#define eslCOINS       4	/* for toy examples      */
#define eslDICE        5	/* also for toy examples */
#define eslNONSTANDARD 6
/* ... if you add here, change esl_abc_ValidateType() too. */

/* Structure: ESL_ALPHABET
 */
typedef struct {
  int      type;	     /* eslDNA, eslRNA, eslAMINO, eslNONSTANDARD, etc.                 */
  int      K;		     /* uniq alphabet size: 4 or 20                                    */
  int      Kp;		     /* total size: alphabet + degen + gap + missing                   */
  char    *sym;              /* "ACGT-RYMKSWHBVDN*~", for instance    [0..Kp-1]                */
  ESL_DSQ  inmap[128];       /* inmap['A'] = 0, etc: dsq[] index for a symbol                  */
  char   **degen;            /* 1/0, which syms inc which res [0..Kp-1][0..K-1]                */
  int     *ndegen;	     /* # of degenerate residues per code  [0..Kp-1]                   */
  ESL_DSQ *complement;       /* maps sym to complements, [0..Kp-1]; NULL if <type> not DNA/RNA */
} ESL_ALPHABET;

/* 1. An ESL_ALPHABET object.
 */
ESL_ALPHABET *esl_alphabet_Create(int type);
ESL_ALPHABET *esl_alphabet_CreateCustom(const char *alphabet, int K, int Kp);
int           esl_alphabet_SetEquiv(ESL_ALPHABET *a, char sym, char c);
int           esl_alphabet_SetCaseInsensitive(ESL_ALPHABET *a);
int           esl_alphabet_SetDegeneracy(ESL_ALPHABET *a, char c, char *ds);
int           esl_alphabet_SetIgnored(ESL_ALPHABET *a, const char *ignoredchars);
size_t        esl_alphabet_Sizeof(ESL_ALPHABET *a);
void          esl_alphabet_Destroy(ESL_ALPHABET *a);

/* 2. Digitized sequences.
 */
int     esl_abc_CreateDsq(const ESL_ALPHABET *a, const char    *seq,        ESL_DSQ **ret_dsq);
int     esl_abc_Digitize (const ESL_ALPHABET *a, const char    *seq,        ESL_DSQ *dsq);
int     esl_abc_Textize  (const ESL_ALPHABET *a, const ESL_DSQ *dsq,  int64_t L, char   *seq);
int     esl_abc_TextizeN (const ESL_ALPHABET *a, const ESL_DSQ *dptr, int64_t L, char   *buf);
int     esl_abc_dsqcpy(const ESL_DSQ *dsq, int64_t L, ESL_DSQ *dcopy);
int     esl_abc_dsqdup(const ESL_DSQ *dsq, int64_t L, ESL_DSQ **ret_dup);
int     esl_abc_dsqcat        (const ESL_DSQ *inmap, ESL_DSQ **dsq, int64_t *L, const char *s, esl_pos_t n);
int     esl_abc_dsqcat_noalloc(const ESL_DSQ *inmap, ESL_DSQ  *dsq, int64_t *L, const char *s, esl_pos_t n);
int64_t esl_abc_dsqlen(const ESL_DSQ *dsq);
int64_t esl_abc_dsqrlen(const ESL_ALPHABET *a, const ESL_DSQ *dsq);
int     esl_abc_CDealign(const ESL_ALPHABET *abc, char    *s, const ESL_DSQ *ref_ax, int64_t *opt_rlen);
int     esl_abc_XDealign(const ESL_ALPHABET *abc, ESL_DSQ *x, const ESL_DSQ *ref_ax, int64_t *opt_rlen);
int     esl_abc_ConvertDegen2X(const ESL_ALPHABET *abc, ESL_DSQ *dsq);
int     esl_abc_revcomp(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int n);

/* 3. Other routines in the API.
 */
int    esl_abc_ValidateType(int type);
int    esl_abc_GuessAlphabet(const int64_t *ct, int *ret_type);
double esl_abc_Match       (const ESL_ALPHABET *a, ESL_DSQ x, ESL_DSQ y, double *p);
int    esl_abc_IAvgScore   (const ESL_ALPHABET *a, ESL_DSQ x, const int    *sc);
float  esl_abc_FAvgScore   (const ESL_ALPHABET *a, ESL_DSQ x, const float  *sc);
double esl_abc_DAvgScore   (const ESL_ALPHABET *a, ESL_DSQ x, const double *sc);
int    esl_abc_IExpectScore(const ESL_ALPHABET *a, ESL_DSQ x, const int    *sc, const float  *p);
float  esl_abc_FExpectScore(const ESL_ALPHABET *a, ESL_DSQ x, const float  *sc, const float  *p);
double esl_abc_DExpectScore(const ESL_ALPHABET *a, ESL_DSQ x, const double *sc, const double *p);

int    esl_abc_IAvgScVec   (const ESL_ALPHABET *a, int    *sc);
int    esl_abc_FAvgScVec   (const ESL_ALPHABET *a, float  *sc);
int    esl_abc_DAvgScVec   (const ESL_ALPHABET *a, double *sc);
int    esl_abc_IExpectScVec(const ESL_ALPHABET *a, int    *sc, const float  *p);
int    esl_abc_FExpectScVec(const ESL_ALPHABET *a, float  *sc, const float  *p);
int    esl_abc_DExpectScVec(const ESL_ALPHABET *a, double *sc, const double *p);
int    esl_abc_FCount      (const ESL_ALPHABET *a, float  *ct, ESL_DSQ x, float  wt);
int    esl_abc_DCount      (const ESL_ALPHABET *a, double *ct, ESL_DSQ x, double wt);
int    esl_abc_EncodeType  (char *typestring);
char  *esl_abc_DecodeType  (int type);
int    esl_abc_ValidateSeq(const ESL_ALPHABET *a, const char *seq, int64_t L, char *errbuf);

/* In the tests below, remember the rules of order in internal alphabets:
 *   Canonical alphabet   Gap   Degeneracies   Any    None    Missing
 *        0..K-1           K      K+1..Kp-4   (Kp-3)  (Kp-2)   (Kp-1)
 *         ACGT            -     RYMKSWHBVD     N       *        ~           DNA: K=4  Kp=18
 *  ACDEFGHIKLMNPQRSTVWY   -        BJZOU       X       *        ~       protein: K=20 Kp=29
 *
 * ESL_DSQ is an unsigned 8-bit type, so don't test for >= 0 or compilers will complain.
 */
#define esl_abc_DigitizeSymbol(a, c) ((a)->inmap[(int)c])
#define esl_abc_XIsValid(a, x)       ((x) < (a)->Kp)
#define esl_abc_XIsResidue(a, x)     ((x) < (a)->K || ((x) > (a)->K && (x) < (a)->Kp-2))
#define esl_abc_XIsCanonical(a, x)   ((x) < (a)->K)
#define esl_abc_XIsGap(a, x)         ((x) == (a)->K)
#define esl_abc_XIsDegenerate(a, x)  ((x) >  (a)->K && (x) < (a)->Kp-2)
#define esl_abc_XIsUnknown(a, x)     ((x) == (a)->Kp-3)
#define esl_abc_XIsNonresidue(a, x)  ((x) == (a)->Kp-2)
#define esl_abc_XIsMissing(a, x)     ((x) == (a)->Kp-1)
#define esl_abc_XGetGap(a)           ((a)->K)
#define esl_abc_XGetUnknown(a)       ((a)->Kp-3)
#define esl_abc_XGetNonresidue(a)    ((a)->Kp-2)
#define esl_abc_XGetMissing(a)       ((a)->Kp-1)

#define esl_abc_CIsValid(a, c)       (isascii(c) && (a)->inmap[(int)c] < (a)->Kp)
#define esl_abc_CIsResidue(a, c)     ((a)->inmap[(int)c] < (a)->K || ((a)->inmap[(int)c] > (a)->K && (a)->inmap[(int)c] < (a)->Kp-2))
#define esl_abc_CIsCanonical(a, c)   ((a)->inmap[(int)c] < (a)->K)
#define esl_abc_CIsGap(a, c)         ((a)->inmap[(int)c] == (a)->K)
#define esl_abc_CIsDegenerate(a, c)  ((a)->inmap[(int)c] > (a)->K  && (a)->inmap[(int)c] < (a)->Kp-2)
#define esl_abc_CIsUnknown(a, c)     ((a)->inmap[(int)c] == (a)->Kp-3)
#define esl_abc_CIsNonresidue(a, c)  ((a)->inmap[(int)c] == (a)->Kp-2)
#define esl_abc_CIsMissing(a, c)     ((a)->inmap[(int)c] == (a)->Kp-1)
#define esl_abc_CGetGap(a)           ((a)->sym[(int)(a)->K])
#define esl_abc_CGetUnknown(a)       ((a)->sym[(int)(a)->Kp-3])
#define esl_abc_CGetNonresidue(a)    ((a)->sym[(int)(a)->Kp-2])
#define esl_abc_CGetMissing(a)       ((a)->sym[(int)(a)->Kp-1])

#endif /*eslALPHABET_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_alphabet.h ***/


/*** Start of inlined file: esl_avx.h ***/
/* Vectorized utility routines for Intel AVX instructions and
 * compatible processors.
 *
 * This header file, unusually, provides many complete function
 * implementations; this is so that they can be inlined by the
 * compiler, for maximum efficiency.
 *
 * Contents:
 *    1. Inlined horizontal functions for 8 and 16-bit quantities
 *       in 256-bit vectors (__m256i)
 */

#ifndef eslAVX_INCLUDED
#define eslAVX_INCLUDED

#include <stdio.h>
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#ifdef HAVE_AVX2 // don't include on architectures that can't compile avx2
/* Function:  esl_avx_hmax_epu8()
 * Synopsis:  Return the unsigned max of the 32 elements in epu8 vector.
 *
 * Purpose:   Returns the maximum value of the 32 elements in
 *            an <epu8> vector.
 */
static inline uint8_t
esl_avx_hmax_epu8(__m256i a)
{

	  __m256i temp1_AVX = _mm256_permute2x128_si256(a, a, 0x01);
	  // Swap the 128-bit halves from a into temp1

	  __m256i temp2_AVX = _mm256_max_epu8(temp1_AVX, a); // each 8-bit field in temp2_AVX now has the max of the
	  //corresponding fields in the high and low halves of a

	  temp1_AVX = _mm256_shuffle_epi32(temp2_AVX, 0x4e);  // Swap the 64-bit halves of each 128-bit half of a
	  temp2_AVX = _mm256_max_epu8(temp1_AVX, temp2_AVX);  // Each 64-bit quantity in temp2 now has the max of the corresponding
	  // 8-bit fields from the 64-bit quarters of a

	  temp1_AVX = _mm256_shuffle_epi32(temp2_AVX, 0xb1);  // Swap the 32-bit halves of each 64-bit quarter of temp2_AVX
	  temp2_AVX = _mm256_max_epu8(temp1_AVX, temp2_AVX);  // Each 32-bit quantity in temp2 now has the max of the corresponding
	  // 8 bit fields from the 32-bit eighths of a

	  temp1_AVX = _mm256_shufflelo_epi16(temp2_AVX, 0xb1); // bottom 32-bits of temp1_AVX now contain the swapped 16-bit halves
	  // of the low 32 bits of temp2_AVX
	  temp2_AVX = _mm256_max_epu8(temp1_AVX, temp2_AVX);  //bottom 16 bits of temp2_AVX now contain the max of the 16-bit fields of a

	  uint8_t temp_stash = _mm256_extract_epi8(temp2_AVX, 1);
	  temp1_AVX = _mm256_insert_epi8(temp2_AVX, temp_stash, 0);  // low byte of temp1_AVX now has byte 2 of temp2_AVX
	  temp2_AVX = _mm256_max_epu8(temp1_AVX, temp2_AVX);  //bottom 16 bits of temp2_AVX now contain the max of the 16-bit fields of Dmaxv_AVX
	  return(_mm256_extract_epi8(temp2_AVX, 0));  // get low byte of temp2_AVX
}

/* Function:  esl_avx_hmax_epi16()
 * Synopsis:  Return the signed max of the 16 elements in epu8 vector.
 *
 * Purpose:   Returns the maximum value of the 32 elements in
 *            an <epu8> vector.
 */
static inline uint16_t
esl_avx_hmax_epi16(__m256i a)
{

	  __m256i temp1_AVX = _mm256_permute2x128_si256(a, a, 0x01);
	  // Swap the 128-bit halves from a into temp1

	  __m256i temp2_AVX = _mm256_max_epi16(temp1_AVX, a); // each 8-bit field in temp2_AVX now has the max of the
	  //corresponding fields in the high and low halves of a

	  temp1_AVX = _mm256_shuffle_epi32(temp2_AVX, 0x4e);  // Swap the 64-bit halves of each 128-bit half of a
	  temp2_AVX = _mm256_max_epi16(temp1_AVX, temp2_AVX);  // Each 64-bit quantity in temp2 now has the max of the corresponding
	  // 8-bit fields from the 64-bit quarters of a

	  temp1_AVX = _mm256_shuffle_epi32(temp2_AVX, 0xb1);  // Swap the 32-bit halves of each 64-bit quarter of temp2_AVX
	  temp2_AVX = _mm256_max_epi16(temp1_AVX, temp2_AVX);  // Each 32-bit quantity in temp2 now has the max of the corresponding
	  // 8 bit fields from the 32-bit eighths of a

	  temp1_AVX = _mm256_shufflelo_epi16(temp2_AVX, 0xb1); // bottom 32-bits of temp1_AVX now contain the swapped 16-bit halves
	  // of the low 32 bits of temp2_AVX
	  temp2_AVX = _mm256_max_epi16(temp1_AVX, temp2_AVX);  //bottom 16 bits of temp2_AVX now contain the max of the 16-bit fields of a

	  return(_mm256_extract_epi16(temp2_AVX, 0));  // get low 16 bits of temp2_AVX
}

// shifts vector left by num_bytes bytes.  Assumes that num_bytes < 16, and will fail horribly if not.
static inline __m256i esl_avx_leftshift(__m256i vector, int num_bytes){
   register __m256i temp_mask_AVX = _mm256_permute2x128_si256(vector, vector, _MM_SHUFFLE(0,0,3,0) );
   return(_mm256_alignr_epi8(vector, temp_mask_AVX,(16-num_bytes)));
}

#endif

#endif /*eslAVX_INCLUDED*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/*** End of inlined file: esl_avx.h ***/


/*** Start of inlined file: esl_buffer.h ***/
#ifndef eslBUFFER_INCLUDED
#define eslBUFFER_INCLUDED

#include <stdio.h>

#define eslBUFFER_PAGESIZE      4096    /* default for b->pagesize                       */
#define eslBUFFER_SLURPSIZE  4194304	/* switchover from slurping whole file to mmap() */

enum esl_buffer_mode_e {
  eslBUFFER_UNSET   = 0,
  eslBUFFER_STREAM  = 1,  /* chunk in mem[0..n-1] = input[baseoffset..baseoffset-n-1];  balloc>0; offset>=0; fp open  */
  eslBUFFER_CMDPIPE = 2,  /* chunk in mem[0..n-1] = input[baseoffset..baseoffset-n-1];  balloc>0; offset>=0; fp open  */
  eslBUFFER_FILE    = 3,  /* chunk in mem[0..n-1] = input[baseoffset..baseoffset-n-1];  balloc>0; offset>=0; fp open  */
  eslBUFFER_ALLFILE = 4,  /* whole file in mem[0..n-1];  balloc=0; offset=0;  fp=NULL  */
  eslBUFFER_MMAP    = 5,  /* whole file in mem[0..n-1];  balloc=0; offset=0;  fp=NULL  */
  eslBUFFER_STRING  = 6   /* whole str in mem[0..n-1];   balloc=0; offset=0;  fp=NULL  */
};

typedef struct {
  char      *mem;	          /* the buffer                                            */
  esl_pos_t  n;		          /* curr buf length; mem[0..n-1] contains valid bytes     */
  esl_pos_t  balloc;              /* curr buf alloc;  mem[0..balloc-1] may be used         */
  esl_pos_t  pos;	          /* curr buf parse position; n-pos = # of parseable bytes */
  esl_pos_t  baseoffset;          /* offset of byte mem[0] in the stream                   */

  esl_pos_t  anchor;	          /* buf[anchor..n-1] safe from overwrite [-1=unset]       */
  int        nanchor;		  /* number of anchors set at <anchor>                     */

  FILE      *fp;	          /* open stream; NULL if already entirely in memory       */
  char      *filename;	          /* for diagnostics. filename; or NULL (stdin, string)    */
  char      *cmdline;		  /* for diagnostics. NULL, or cmd for CMDPIPE             */

  esl_pos_t  pagesize;	          /* size of new <fp> reads. Guarantee: n-pos >= pagesize  */

  char     errmsg[eslERRBUFSIZE]; /* error message storage                                 */
  enum esl_buffer_mode_e mode_is; /* mode (stdin, cmdpipe, file, allfile, mmap, string)    */
} ESL_BUFFER;

/* 1. The ESL_BUFFER object: opening/closing.  */
int esl_buffer_Open      (const char *filename, const char *envvar, ESL_BUFFER **ret_bf);
int esl_buffer_OpenFile  (const char *filename,                     ESL_BUFFER **ret_bf);
int esl_buffer_OpenPipe  (const char *filename, const char *cmdfmt, ESL_BUFFER **ret_bf);
int esl_buffer_OpenMem   (const char *p,         esl_pos_t  n,      ESL_BUFFER **ret_bf);
int esl_buffer_OpenStream(FILE *fp,                                 ESL_BUFFER **ret_bf);
int esl_buffer_Close(ESL_BUFFER *bf);

/* 2. Positioning and anchoring an ESL_BUFFER. */
esl_pos_t esl_buffer_GetOffset      (ESL_BUFFER *bf);
int       esl_buffer_SetOffset      (ESL_BUFFER *bf, esl_pos_t offset);
int       esl_buffer_SetAnchor      (ESL_BUFFER *bf, esl_pos_t offset);
int       esl_buffer_SetStableAnchor(ESL_BUFFER *bf, esl_pos_t offset);
int       esl_buffer_RaiseAnchor    (ESL_BUFFER *bf, esl_pos_t offset);

/* 3. Raw access to the buffer */
int esl_buffer_Get(ESL_BUFFER *bf, char **ret_p, esl_pos_t *ret_n);
int esl_buffer_Set(ESL_BUFFER *bf, char  *p,     esl_pos_t  nused);

/* 4. Line-based parsing */
int esl_buffer_GetLine       (ESL_BUFFER *bf, char **opt_p, esl_pos_t *opt_n);
int esl_buffer_FetchLine     (ESL_BUFFER *bf, char **opt_p, esl_pos_t *opt_n);
int esl_buffer_FetchLineAsStr(ESL_BUFFER *bf, char **opt_s, esl_pos_t *opt_n);

/* 5. Token-based parsing */
int esl_buffer_GetToken       (ESL_BUFFER *bf, const char *sep, char **opt_p, esl_pos_t *opt_n);
int esl_buffer_FetchToken     (ESL_BUFFER *bf, const char *sep, char **opt_p, esl_pos_t *opt_n);
int esl_buffer_FetchTokenAsStr(ESL_BUFFER *bf, const char *sep, char **opt_p, esl_pos_t *opt_n);

/* 6. Binary (fread-like) parsing */
int esl_buffer_Read(ESL_BUFFER *bf, size_t nbytes, void *p);

/* Replaces functionality of esl_fileparser module as follows
 *
 *  esl_fileparser_Open()             -> esl_buffer_Open()
 *  esl_fileparser_Create()           -> esl_buffer_OpenStream()
 *  esl_fileparser_CreateMapped()     -> esl_buffer_OpenMem()
 *  esl_fileparser_SetCommentChar()   -> esl_buffer_SetCommentChar()
 *  esl_fileparser_GetToken()         -> esl_buffer_GetToken()
 *  esl_fileparser_NextLine()         -> do { esl_buffer_GetLine() } while esl_line_IsBlank();
 *  esl_fileparser_NextLinePeeked()   -> unneeded. esl_buffer_Peek*() functionality, syntax different
 *  esl_fileparser_GetTokenOnLine()   -> unneeded. esl_buffer_GetToken() has an idiom.
 *  esl_fileparser_GetRemainingLine() -> esl_buffer_GetLine()
 *  esl_fileparser_Destroy()          -> esl_buffer_Close()
 *  esl_fileparser_Close()            -> esl_buffer_Close()
 */

/* Replaces functionality of esl_recorder module as follows:
 *
 *  esl_recorder_Create()             -> esl_buffer_OpenStream()
 *  esl_recorder_ResizeTo()
 *  esl_recorder_GetFirst()           ->
 *  esl_recorder_GetLast()            ->
 *  esl_recorder_GetCurrent()         ->
 *  esl_recorder_GetNext()            ->
 *  esl_recorder_Destroy()            -> esl_buffer_Close()
 *  esl_recorder_Read()               -> esl_buffer_GetLine()
 *  esl_recorder_Position()
 *  esl_recorder_MarkBlock()
 *  esl_recorder_UnmarkBlock()
 *  esl_recorder_GetBlock()
 *
 */

#endif	/*eslBUFFER_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_buffer.h ***/


/*** Start of inlined file: esl_cluster.h ***/
#ifndef eslCLUSTER_INCLUDED
#define eslCLUSTER_INCLUDED

int esl_cluster_SingleLinkage(void *base, size_t n, size_t size,
				     int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
				     int *workspace, int *assignments, int *ret_C);
#endif /*eslCLUSTER_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_cluster.h ***/


/*** Start of inlined file: esl_composition.h ***/
int  esl_composition_BL62(double *f);
int  esl_composition_WAG (double *f);
int  esl_composition_SW34(double *f);
int  esl_composition_SW50(double *f);

/*** End of inlined file: esl_composition.h ***/


/*** Start of inlined file: esl_dirichlet.h ***/
#ifndef eslDIRICHLET_INCLUDED
#define eslDIRICHLET_INCLUDED

/* Structure: MIXDCHLET
 *
 * A mixture Dirichlet density, usually used as a prior
 * for a multinomial model (turning count vectors into probability
 * parameters).
 */
typedef struct {
 /*::cexcerpt::dirichlet_mixdchlet::begin::*/
  double  *pq;			/* mixture coefficients pq[0..N-1]          */
  double **alpha;               /* Dirichlet params alpha[0..N-1][0..K-1]   */
  int      N;			/* number of mixtures, e.g. 9 for Sjolander */
  int      K;			/* alphabet size, e.g. 20                   */
 /*::cexcerpt::dirichlet_mixdchlet::end::*/
} ESL_MIXDCHLET;

ESL_MIXDCHLET *esl_mixdchlet_Create(int N, int K);
int            esl_mixdchlet_Compare(ESL_MIXDCHLET *d1, ESL_MIXDCHLET *d2, double tol);
int            esl_mixdchlet_Copy(ESL_MIXDCHLET *d, ESL_MIXDCHLET *d_dst);
int            esl_mixdchlet_Dump(FILE *fp, ESL_MIXDCHLET *d);
void           esl_mixdchlet_Destroy(ESL_MIXDCHLET *pri);
int            esl_mixdchlet_MPParameters(double *c, int K,
 						 ESL_MIXDCHLET *pri, double *mix, double *p);
int            esl_mixdchlet_BILD_score(double *c, int K, int N, ESL_MIXDCHLET *pri,
												 double *mix, double *bg, double *q);

int esl_dirichlet_LogProbData(double *c, double *alpha, int K,
				     double *ret_answer);
int esl_dirichlet_LogProbData_Mixture(double *c, ESL_MIXDCHLET *d,
					     double *ret_answer);
int esl_dirichlet_LogProbProbs(double *p, double *alpha, int K,
				      double *ret_answer);

/* Optional fitting code, when augmented by minimizing module.
 */
#ifdef eslAUGMENT_MINIMIZER

/*** Start of inlined file: esl_minimizer.h ***/
#ifndef eslMINIMIZER_INCLUDED
#define eslMINIMIZER_INCLUDED

#define MAXITERATIONS 100

int esl_min_Bracket(double *a, double *d, double *u, int n,
			   double (*func)(double *, int, void *), void *prm,
			   double *ret_fa,
			   double *b, double *ret_bx, double *ret_fb,
			   double *c, double *ret_cx, double *ret_fc);
int esl_min_LineSearch(double *ori, double *d, double *u, int n,
			      double (*func)(double *, int, void *), void *prm,
			      double tol, double *b,
			      double *x, double *ret_xx, double *ret_fx);
int esl_min_ConjugateGradientDescent(double *x, double *u, int n,
					    double (*func)(double *, int, void *),
					    void (*dfunc)(double *, int, void *, double *),
					    void *prm, double tol, double *wrk, double *ret_fx);

#endif /*eslMINIMIZER_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_minimizer.h ***/


int esl_mixdchlet_Fit(double **c, int nc, ESL_MIXDCHLET *d, int be_verbose);
#ifdef eslAUGMENT_RANDOM

/*** Start of inlined file: esl_random.h ***/
#ifndef eslRANDOM_INCLUDED
#define eslRANDOM_INCLUDED

#include <stdint.h>

#define eslRND_FAST     0
#define eslRND_MERSENNE 1

typedef struct {
  int      type;		/* eslRND_FAST | eslRND_MERSENNE               */
  int      mti;			/* current position in mt[] table              */
  uint32_t mt[624];		/* state of the Mersenne Twister               */
  uint32_t x;			/* state of the Knuth generator                */
  uint32_t seed;		/* seed used to init the RNG                   */
} ESL_RANDOMNESS;

/* esl_rnd_Roll(a) chooses a uniformly distributed integer
 * in the range 0..a-1, given an initialized ESL_RANDOMNESS r.
 */
#define esl_rnd_Roll(r, a)    ((int) (esl_random(r) * (a)))

/* 1. The ESL_RANDOMNESS object.
 */
ESL_RANDOMNESS *esl_randomness_Create    (uint32_t seed);
ESL_RANDOMNESS *esl_randomness_CreateFast(uint32_t seed);
ESL_RANDOMNESS *esl_randomness_CreateTimeseeded(void); /* DEPRECATED */
void            esl_randomness_Destroy(ESL_RANDOMNESS *r);
int             esl_randomness_Init(ESL_RANDOMNESS *r, uint32_t seed);
uint32_t        esl_randomness_GetSeed(const ESL_RANDOMNESS *r);

/* 2. The generator, esl_random().
 */
double   esl_random       (ESL_RANDOMNESS *r);
uint32_t esl_random_uint32(ESL_RANDOMNESS *r);

/* 3. Debugging/development tools.
 */
int esl_randomness_Dump(FILE *fp, ESL_RANDOMNESS *r);

/* 4. Other fundamental sampling (including Gaussian, gamma).
 */
double esl_rnd_UniformPositive(ESL_RANDOMNESS *r);
double esl_rnd_Gaussian (ESL_RANDOMNESS *rng, double mean, double stddev);
double esl_rnd_Gamma    (ESL_RANDOMNESS *rng, double a);
void   esl_rnd_Dirichlet(ESL_RANDOMNESS *rng, const double *alpha, int K, double *p);
void   esl_rnd_mem      (ESL_RANDOMNESS *rng, void *buf, int n);

/* 5. Multinomial sampling from discrete probability n-vectors.
 */
int    esl_rnd_DChoose   (ESL_RANDOMNESS *r, const double *p,   int N);
int    esl_rnd_FChoose   (ESL_RANDOMNESS *r, const float  *p,   int N);
int    esl_rnd_DChooseCDF(ESL_RANDOMNESS *r, const double *cdf, int N);
int    esl_rnd_FChooseCDF(ESL_RANDOMNESS *r, const float  *cdf, int N);

#endif /*eslRANDOM_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_random.h ***/


int esl_mixdchlet_Fit_Multipass(ESL_RANDOMNESS *rng, double **c, int nc, int reps, ESL_MIXDCHLET *best_md, int verbose);
#endif /*eslAUGMENT_RANDOM*/
#endif /*eslAUGMENT_MINIMIZER*/

/* Optional sampling code, when augmented by random module.
 */
#ifdef eslAUGMENT_RANDOM
int esl_dirichlet_DSample(ESL_RANDOMNESS *r, double *alpha, int K, double *p);
int esl_dirichlet_FSample(ESL_RANDOMNESS *r, float  *alpha, int K, float  *p);
int esl_dirichlet_DSampleUniform(ESL_RANDOMNESS *r, int K, double *p);
int esl_dirichlet_FSampleUniform(ESL_RANDOMNESS *r, int K, float  *p);
int esl_dirichlet_SampleBeta(ESL_RANDOMNESS *r, double theta1,
				    double theta2, double *ret_answer);
#endif /*eslAUGMENT_RANDOM*/

/* Optional file input code, when augmented by fileparser module
 */
#ifdef eslAUGMENT_FILEPARSER

/*** Start of inlined file: esl_fileparser.h ***/
#ifndef eslFILEPARSER_INCLUDED
#define eslFILEPARSER_INCLUDED

#include <stdio.h>

typedef struct {
  FILE *fp;			/* open file pointer, for reading                  */
  char *buf;			/* current line; will be modified by esl_strtok(). */
  int   buflen;			/* current allocated length of buf                 */
  char *s;			/* used by esl_strtok(); current position in buf.  */
  char  commentchar;		/* often '#'                                       */

  char *filename;		/* name of opened file; or NULL (if just a stream) */
  int   linenumber;		/* what line is loaded into buf; 1..nlines         */
  char  errbuf[eslERRBUFSIZE];  /* for holding error diagnostics                   */

  int   is_buffer;              /* the file has been buffered into memory          */
  char *mem_buffer;             /* pointer to the buffered file                    */
  int   mem_size;               /* size of the buffered file                       */
  int   mem_pos;                /* current position in the buffer                  */
} ESL_FILEPARSER;

int  esl_fileparser_Open(const char *filename, const char *envvar, ESL_FILEPARSER **ret_efp);
ESL_FILEPARSER *esl_fileparser_Create(FILE *fp);
ESL_FILEPARSER *esl_fileparser_CreateMapped(void *buffer, int size);
int  esl_fileparser_SetCommentChar  (ESL_FILEPARSER *efp, char c);
int  esl_fileparser_GetToken        (ESL_FILEPARSER *efp, char **opt_tok, int *opt_toklen);
int  esl_fileparser_NextLine        (ESL_FILEPARSER *efp);
int  esl_fileparser_NextLinePeeked  (ESL_FILEPARSER *efp, char *prefix, int plen);
int  esl_fileparser_GetTokenOnLine  (ESL_FILEPARSER *efp, char **opt_tok, int *opt_toklen);
int  esl_fileparser_GetRemainingLine(ESL_FILEPARSER *efp, char **ret_s);
void esl_fileparser_Destroy         (ESL_FILEPARSER *efp);
void esl_fileparser_Close           (ESL_FILEPARSER *efp);

#endif /*eslFILEPARSER_INCLUDED */
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_fileparser.h ***/


int esl_mixdchlet_Read(ESL_FILEPARSER *efp,  ESL_MIXDCHLET **ret_pri);
int esl_mixdchlet_Write(FILE *fp,  ESL_MIXDCHLET *d);
#endif /*eslAUGMENT_FILEPARSER*/

#endif /*eslDIRICHLET_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_dirichlet.h ***/


/*** Start of inlined file: esl_distance.h ***/
#ifndef eslDISTANCE_INCLUDED
#define eslDISTANCE_INCLUDED

#ifdef eslAUGMENT_ALPHABET

#endif
#ifdef eslAUGMENT_DMATRIX

/*** Start of inlined file: esl_dmatrix.h ***/
#ifndef eslDMATRIX_INCLUDED
#define eslDMATRIX_INCLUDED

#include <stdio.h>

typedef struct {
  /*mx, mx[0] are allocated. */
/*::cexcerpt::dmatrix_obj::begin::*/
  double **mx;                  /* mx[i][j] is i'th row, j'th col */
  int      n;                   /* rows    */
  int      m;                   /* columns */
  enum { eslGENERAL, eslUPPER } type;
/*::cexcerpt::dmatrix_obj::end::*/
  int      ncells;		/* number of valid cells (nxm in standard matrix) */
} ESL_DMATRIX;

typedef struct {
  int     *pi;
  int      n;
} ESL_PERMUTATION;

/* 1. The ESL_DMATRIX object. */
ESL_DMATRIX *esl_dmatrix_Create(int n, int m);
ESL_DMATRIX *esl_dmatrix_CreateUpper(int n);
int          esl_dmatrix_Destroy(ESL_DMATRIX *A);
int          esl_dmatrix_Copy       (const ESL_DMATRIX *src, ESL_DMATRIX *dest);
ESL_DMATRIX *esl_dmatrix_Clone      (const ESL_DMATRIX *old);
int          esl_dmatrix_Compare    (const ESL_DMATRIX *A, const ESL_DMATRIX *B, double tol);
int          esl_dmatrix_CompareAbs (const ESL_DMATRIX *A, const ESL_DMATRIX *B, double tol);
int          esl_dmatrix_Set        (ESL_DMATRIX *A, double x);
int          esl_dmatrix_SetZero    (ESL_DMATRIX *A);
int          esl_dmatrix_SetIdentity(ESL_DMATRIX *A);

/* 2. Debugging/validation for ESL_DMATRIX. */
int          esl_dmatrix_Dump(FILE *ofp, const ESL_DMATRIX *A,
				     const char *rowlabel, const char *collabel);

/* 3. Visualization tools. */
int          esl_dmatrix_PlotHeatMap(FILE *fp, ESL_DMATRIX *D, double min, double max);

/* 4. The ESL_PERMUTATION object. */
ESL_PERMUTATION *esl_permutation_Create(int n);
int              esl_permutation_Destroy(ESL_PERMUTATION *P);
int              esl_permutation_Reuse(ESL_PERMUTATION *P);

/* 5. Debugging/validation for ESL_PERMUTATION. */
int              esl_permutation_Dump(FILE *ofp, const ESL_PERMUTATION *P,
					     const char *rowlabel, const char *collabel);

/* 6. The rest of the dmatrix API. */
double       esl_dmx_Max    (const ESL_DMATRIX *A);
double       esl_dmx_Min    (const ESL_DMATRIX *A);
double       esl_dmx_Sum    (const ESL_DMATRIX *A);
int          esl_dmx_MinMax(const ESL_DMATRIX *A, double *ret_min, double *ret_max);
int          esl_dmx_FrobeniusNorm(const ESL_DMATRIX *A, double *ret_fnorm);
int          esl_dmx_Multiply(const ESL_DMATRIX *A, const ESL_DMATRIX *B, ESL_DMATRIX *C);
int          esl_dmx_Exp(const ESL_DMATRIX *Q, double t, ESL_DMATRIX *P);
int          esl_dmx_Transpose(ESL_DMATRIX *A);
int          esl_dmx_Add(ESL_DMATRIX *A, const ESL_DMATRIX *B);
int          esl_dmx_Scale(ESL_DMATRIX *A, double k);
int          esl_dmx_AddScale(ESL_DMATRIX *A, double k, const ESL_DMATRIX *B);
int          esl_dmx_Permute_PA(const ESL_PERMUTATION *P, const ESL_DMATRIX *A, ESL_DMATRIX *B);
int          esl_dmx_LUP_decompose(ESL_DMATRIX *A, ESL_PERMUTATION *P);
int          esl_dmx_LU_separate(const ESL_DMATRIX *LU, ESL_DMATRIX *L, ESL_DMATRIX *U);
int          esl_dmx_Invert(const ESL_DMATRIX *A, ESL_DMATRIX *Ai);

/* 7. Optional: interoperability with GSL */
#ifdef HAVE_LIBGSL
#include <gsl/gsl_matrix.h>
int          esl_dmx_MorphGSL(const ESL_DMATRIX *E, gsl_matrix **ret_G);
int          esl_dmx_UnmorphGSL(const gsl_matrix *G, ESL_DMATRIX **ret_E);
#endif

/* 8. Optional: interfaces to LAPACK  */
#ifdef HAVE_LIBLAPACK
int esl_dmx_Diagonalize(const ESL_DMATRIX *A, double **ret_Er, double **ret_Ei, ESL_DMATRIX **ret_UL, ESL_DMATRIX **ret_UR);
#endif

#endif /*eslDMATRIX_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_dmatrix.h ***/


#endif
#ifdef eslAUGMENT_RANDOM

#endif

/* 1. Pairwise distances for aligned text sequences.
 */
int esl_dst_CPairId(const char *asq1, const char *asq2,
			   double *opt_pid, int *opt_nid, int *opt_n);
int esl_dst_CPairmatch(const char *asq1, const char *asq2,
			      double *opt_pmatch, int *opt_nmatch, int *opt_n);
int esl_dst_CJukesCantor(int K, const char *as1, const char *as2,
				double *opt_distance, double *opt_variance);

/* 2. Pairwise distances for aligned digital seqs.
 */
#ifdef eslAUGMENT_ALPHABET
int esl_dst_XPairId(const ESL_ALPHABET *abc, const ESL_DSQ *ax1, const ESL_DSQ *ax2,
			   double *opt_pid, int *opt_nid, int *opt_n);
int esl_dst_XPairMatch(const ESL_ALPHABET *abc, const ESL_DSQ *ax1, const ESL_DSQ *ax2,
			      double *opt_distance, int *opt_nmatch, int *opt_n);
int esl_dst_XJukesCantor(const ESL_ALPHABET *abc, const ESL_DSQ *ax, const ESL_DSQ *ay,
				double *opt_distance, double *opt_variance);
#endif

/* 3. Distance matrices for aligned text sequences.
 */
#ifdef eslAUGMENT_DMATRIX
int esl_dst_CPairIdMx     (char **as, int N, ESL_DMATRIX **ret_S);
int esl_dst_CDiffMx       (char **as, int N, ESL_DMATRIX **ret_D);
int esl_dst_CJukesCantorMx(int K, char **as, int N, ESL_DMATRIX **opt_D, ESL_DMATRIX **opt_V);
#endif

/* 4. Distance matrices for aligned digital sequences.
 */
#if defined(eslAUGMENT_DMATRIX) && defined(eslAUGMENT_ALPHABET)
int esl_dst_XPairIdMx(const ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_S);
int esl_dst_XDiffMx  (const ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_D);

int esl_dst_XJukesCantorMx(const ESL_ALPHABET *abc, ESL_DSQ **ax, int nseq,
				  ESL_DMATRIX **opt_D, ESL_DMATRIX **opt_V);
#endif

/*  5. Average pairwise identity for multiple alignments.
 */
#ifdef eslAUGMENT_RANDOM
int esl_dst_CAverageId   (char **as, int nseq, int max_comparisons, double *ret_id);
int esl_dst_CAverageMatch(char **as, int N, int max_comparisons, double *ret_match);
#endif
#if defined(eslAUGMENT_RANDOM) && defined(eslAUGMENT_ALPHABET)
int esl_dst_XAverageId   (const ESL_ALPHABET *abc, ESL_DSQ **ax, int N, int max_comparisons, double *ret_id);
int esl_dst_XAverageMatch(const ESL_ALPHABET *abc, ESL_DSQ **ax, int N, int max_comparisons, double *ret_match);

#endif

#endif /*eslDISTANCE_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_distance.h ***/


/*** Start of inlined file: esl_dsqdata.h ***/
#ifndef eslDSQDATA_INCLUDED
#define eslDSQDATA_INCLUDED


/*** Start of inlined file: esl_sqio.h ***/
#ifndef eslSQIO_INCLUDED
#define eslSQIO_INCLUDED

#include <stdio.h>


/*** Start of inlined file: esl_sqio_ascii.h ***/
#ifndef eslSQIO_ASCII_INCLUDED
#define eslSQIO_ASCII_INCLUDED

#include <stdio.h>

/*** Start of inlined file: esl_sq.h ***/
#ifndef eslSQ_INCLUDED
#define eslSQ_INCLUDED

#ifdef eslAUGMENT_ALPHABET

#endif
#ifdef eslAUGMENT_MSA

/*** Start of inlined file: esl_msa.h ***/
#ifndef eslMSA_INCLUDED
#define eslMSA_INCLUDED

#include <stdio.h>


/*** Start of inlined file: esl_keyhash.h ***/
#ifndef eslKEYHASH_INCLUDED
#define eslKEYHASH_INCLUDED

#include <stdio.h>		/* for FILE */

/* ESL_KEYHASH:
 *    a dynamically resized hash structure;
 *    contains a hash table and associated data
 *
 * Each key string is associated with an index i = (0..nkeys-1).
 * Key strings are stored in one array, in smem.
 * Each key has an offset in this array, key_offset[i].
 * Thus key number <i> is at: smem + key_offset[i].
 *
 * The keys are hashed, and stored in linked lists in
 * a hashtable by their index i = (0..nkeys-1), with -1
 * as a sentinel for end-of-list.
 *
 * hashtable[0..hashsize-1] = head of linked list;
 *                            index of first elem in list (0..nkeys-1),
 *                            or -1 if empty.
 * nxt[0..nkeys-1] = next elem in list (0..nkeys-1), or -1 if none.
 *
 * Thus a typical loop, looking for a <key>:
 *    uint32_t val = jenkins_hash(key, kh->hashsize);
 *    for (i = kh->hashtable[val]; i != -1; i = kh->nxt[i])
 *      if (strcmp(key, kh->smem + kh->key_offset[i]) == 0) found_it;
 *
 */
typedef struct {
  int      *hashtable;          /* hashtable[0..hashsize-1] = index of first elem, or -1 */
  uint32_t  hashsize;	        /* size of the hash table                                */

  int      *key_offset;		/* key [idx=0..nkeys-1] starts at smem + key_offset[idx] */
  int      *nxt;		/* nxt [idx=0..nkeys-1], next "pointers" in hash table   */
  int       nkeys;		/* number of keys stored                                 */
  int       kalloc;		/* number of keys allocated for                          */

  char *smem;	 	        /* Array of memory for storing key strings (w/ \0's)     */
  int   salloc;			/* current allocated size of <key_mem>                   */
  int   sn; 			/* current used size of key strings, inclusive \0's      */
} ESL_KEYHASH;

ESL_KEYHASH *esl_keyhash_Create(void);
ESL_KEYHASH *esl_keyhash_CreateCustom(uint32_t hashsize, int kalloc, int salloc);
ESL_KEYHASH *esl_keyhash_Clone(const ESL_KEYHASH *kh);
char *       esl_keyhash_Get(const ESL_KEYHASH *kh, int idx);
int          esl_keyhash_GetNumber(const ESL_KEYHASH *kh);
size_t       esl_keyhash_Sizeof(const ESL_KEYHASH *kh);
int          esl_keyhash_Reuse(ESL_KEYHASH *kh);
void         esl_keyhash_Destroy(ESL_KEYHASH *kh);
void         esl_keyhash_Dump(FILE *fp, const ESL_KEYHASH *kh);

int  esl_keyhash_Store (      ESL_KEYHASH *kh, const char *key, esl_pos_t n, int *ret_index);
int  esl_keyhash_Lookup(const ESL_KEYHASH *kh, const char *key, esl_pos_t n, int *ret_index);

#endif /* eslKEYHASH_INCLUDED */
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_keyhash.h ***/


/*** Start of inlined file: esl_ssi.h ***/
#ifndef eslSSI_INCLUDED
#define eslSSI_INCLUDED
#ifdef eslAUGMENT_SSI

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#define eslSSI_MAXFILES 32767	     /* 2^15-1 */
#define eslSSI_MAXKEYS  2147483647L  /* 2^31-1 */
#define eslSSI_MAXRAM   256	     /* >256MB indices trigger external sort */

#ifndef HAVE_FSEEKO
#define fseeko fseek
#define ftello ftell
#endif

/* ESL_SSI
 * Using an existing SSI index file.
 */
typedef struct {
  FILE      *fp;              /* open SSI index file                 */
  uint32_t   flags;	      /* optional behavior flags             */
  uint32_t   offsz;	      /* sizeof(off_t)'s in the SSI file     */
  uint16_t   nfiles;          /* number of files = 16 bit int        */
  uint64_t   nprimary;        /* number of primary keys              */
  uint64_t   nsecondary;      /* number of secondary keys            */
  uint32_t   flen;            /* length of filenames (inc '\0')      */
  uint32_t   plen;            /* length of primary keys (inc '\0')   */
  uint32_t   slen;            /* length of secondary keys (inc '\0') */
  uint32_t   frecsize;        /* # bytes in a file record            */
  uint32_t   precsize;        /* # bytes in a primary key record     */
  uint32_t   srecsize;        /* # bytes in a secondary key record   */
  off_t      foffset;         /* disk offset, start of file records  */
  off_t      poffset;         /* disk offset, start of pri key recs  */
  off_t      soffset;         /* disk offset, start of sec key recs  */

  /* File information:  */
  char     **filename;        /* list of file names [0..nfiles-1]    */
  uint32_t  *fileformat;      /* file formats                        */
  uint32_t  *fileflags;	      /* optional per-file behavior flags    */
  uint32_t  *bpl;             /* bytes per line in file              */
  uint32_t  *rpl;             /* residues per line in file           */
} ESL_SSI;

/* Flags for the <ssi->fileflags> bit vectors. */
#define eslSSI_FASTSUBSEQ   (1<<0)    /* we can do fast subseq lookup calculations on this file */

/* ESL_NEWSSI
 * Used to create a new SSI index.
 */
typedef struct {		/* Primary key data: */
  char      *key;               /* key name          */
  uint16_t   fnum;		/* file number       */
  off_t      r_off;		/* record offset     */
  off_t      d_off;		/* data offset       */
  int64_t    len;		/* sequence length   */
} ESL_PKEY;

typedef struct {		/* Secondary key data: */
  char        *key;             /* secondary key name  */
  char        *pkey;            /* primary key name    */
} ESL_SKEY;

typedef struct {
  char       *ssifile;		/* name of the SSI file we're creating    */
  FILE       *ssifp;		/* open SSI file being created            */
  int         external;	        /* TRUE if pkeys and skeys are on disk    */
  int         max_ram;	        /* threshold in MB to trigger sort */

  char      **filenames;
  uint32_t   *fileformat;
  uint32_t   *bpl;
  uint32_t   *rpl;
  uint32_t    flen;		/* length of longest filename, inc '\0' */
  uint16_t    nfiles;		/* can store up to 2^15-1 (32767) files */

  ESL_PKEY   *pkeys;
  uint32_t    plen;	        /* length of longest pkey, including '\0'    */
  uint64_t    nprimary;		/* can store up to 2^63-1 = 9.2e18 keys      */
  char       *ptmpfile;		/* primary key tmpfile name, for sort */
  FILE       *ptmp;	        /* handle on open ptmpfile */

  ESL_SKEY   *skeys;
  uint32_t    slen;        	/* length of longest skey, including '\0' */
  uint64_t    nsecondary;
  char       *stmpfile;		/* secondary key tmpfile name, for sort */
  FILE       *stmp;	        /* handle on open ptmpfile */

  char        errbuf[eslERRBUFSIZE];
} ESL_NEWSSI;

#define eslSSI_FCHUNK  16	/* chunk size for file name reallocation */
#define eslSSI_KCHUNK  128	/* and for key reallocation              */

/* 1. Using (reading) SSI indices */
int  esl_ssi_Open(const char *filename, ESL_SSI **ret_ssi);
void esl_ssi_Close(ESL_SSI *ssi);
int  esl_ssi_FindName(ESL_SSI *ssi, const char *key,
			     uint16_t *ret_fh, off_t *ret_roff, off_t *opt_doff, int64_t *opt_L);
int  esl_ssi_FindNumber(ESL_SSI *ssi, int64_t nkey,
			       uint16_t *opt_fh, off_t *opt_roff, off_t *opt_doff, int64_t *opt_L, char **opt_pkey);
int  esl_ssi_FindSubseq(ESL_SSI *ssi, const char *key, int64_t requested_start,
			       uint16_t *ret_fh, off_t *ret_roff, off_t *ret_doff, int64_t *ret_L, int64_t *ret_actual_start);
int  esl_ssi_FileInfo(ESL_SSI *ssi, uint16_t fh, char **ret_filename, int *ret_format);

/* 2. Creating (writing) SSI indices. */
int  esl_newssi_Open(const char *ssifile, int allow_overwrite, ESL_NEWSSI **ret_newssi);
int  esl_newssi_AddFile  (ESL_NEWSSI *ns, const char *filename, int fmt, uint16_t *ret_fh);
int  esl_newssi_SetSubseq(ESL_NEWSSI *ns, uint16_t fh, uint32_t bpl, uint32_t rpl);
int  esl_newssi_AddKey   (ESL_NEWSSI *ns, const char *key, uint16_t fh, off_t r_off, off_t d_off, int64_t L);
int  esl_newssi_AddAlias (ESL_NEWSSI *ns, const char *alias, const char *key);
int  esl_newssi_Write    (ESL_NEWSSI *ns);
void esl_newssi_Close    (ESL_NEWSSI *ns);

/* 3. Portable binary i/o. */
void     esl_byteswap(char *swap, int nbytes);
uint16_t esl_ntoh16(uint16_t netshort);
uint32_t esl_ntoh32(uint32_t netlong);
uint64_t esl_ntoh64(uint64_t net_int64);
uint16_t esl_hton16(uint16_t hostshort);
uint32_t esl_hton32(uint32_t hostlong);
uint64_t esl_hton64(uint64_t host_int64);
int      esl_fread_u16(FILE *fp, uint16_t *ret_result);
int      esl_fread_u32(FILE *fp, uint32_t *ret_result);
int      esl_fread_u64(FILE *fp, uint64_t *ret_result);
int      esl_fread_i16(FILE *fp, int16_t  *ret_result);
int      esl_fread_i32(FILE *fp, int32_t  *ret_result);
int      esl_fread_i64(FILE *fp, int64_t  *ret_result);
int      esl_fwrite_u16(FILE *fp, uint16_t n);
int      esl_fwrite_u32(FILE *fp, uint32_t n);
int      esl_fwrite_u64(FILE *fp, uint64_t n);
int      esl_fwrite_i16(FILE *fp, int16_t  n);
int      esl_fwrite_i32(FILE *fp, int32_t  n);
int      esl_fwrite_i64(FILE *fp, int64_t  n);
int	esl_fread_offset(FILE *fp, int mode, off_t *ret_offset);
int      esl_fwrite_offset(FILE *fp, off_t offset);
#endif /* eslAUGMENT_SSI*/
#ifndef eslAUGMENT_SSI
typedef void ESL_SSI;
typedef void ESL_NEWSSI;
#endif
#endif /* eslSSI_INCLUDED */
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_ssi.h ***/

/* The following constants define the Pfam/Rfam cutoff set we propagate
 * from Stockholm format msa's into HMMER and Infernal models.
 */
/*::cexcerpt::msa_cutoffs::begin::*/
#define eslMSA_TC1     0
#define eslMSA_TC2     1
#define eslMSA_GA1     2
#define eslMSA_GA2     3
#define eslMSA_NC1     4
#define eslMSA_NC2     5
#define eslMSA_NCUTS   6
/*::cexcerpt::msa_cutoffs::end::*/

/* Object: ESL_MSA
 *
 * A multiple sequence alignment.
 */
typedef struct {
  /* Mandatory information associated with the alignment.
   * (The important stuff.)
   */
  /*::cexcerpt::msa_mandatory::begin::*/
  char  **aseq;       /* alignment itself, [0..nseq-1][0..alen-1], \0-terminated */
  char  **sqname;     /* sequence names [0..nseq-1][], \0-terminated             */
  double *wgt;        /* sequence weights [0..nseq-1], default 1.0               */
  int64_t alen;       /* length of alignment (columns); or (if growable) -1      */
  int     nseq;       /* number of seqs in alignment; or (if growable) blocksize */
  int     flags;      /* flags for what info has been set                        */
  /*::cexcerpt::msa_mandatory::end::*/

#ifdef eslAUGMENT_ALPHABET
  /* When augmented w/ digital alphabets, we can store pre-digitized data in
   * ax[][], instead of the text info in aseq[][].
   */
  ESL_ALPHABET  *abc;    	/* reference ptr to alphabet            */
  ESL_DSQ      **ax;		/* digitized aseqs [0..nseq-1][1..alen] */
#endif

  /* Optional information that we understand, and that we might have.
   * (The occasionally useful stuff.)
   */
  /*::cexcerpt::msa_optional::begin::*/
  char  *name;      /* name of alignment, or NULL                                           */
  char  *desc;      /* description of alignment, or NULL                                    */
  char  *acc;       /* accession of alignment, or NULL                                      */
  char  *au;        /* "author" information, or NULL                                        */
  char  *ss_cons;   /* consensus sec structure, or NULL;  [0..alen-1], even in digital mode */
  char  *sa_cons;   /* consensus surface access, or NULL; [0..alen-1], even in digital mode */
  char  *pp_cons;   /* consensus posterior prob, or NULL; [0..alen-1], even in digital mode */
  char  *rf;        /* reference coord system, or NULL;   [0..alen-1], even in digital mode */
  char  *mm;        /* model mask, or NULL;   [0..alen-1], even in digital mode             */
  char **sqacc;     /* accession numbers for sequences i                                    */
  char **sqdesc;    /* description lines for sequences i                                    */
  char **ss;        /* per-seq secondary structures, or NULL                                */
  char **sa;        /* per-seq surface accessibilities, or NULL                             */
  char **pp;        /* posterior prob per residue, or NULL                                  */
  float  cutoff[eslMSA_NCUTS];  /* NC/TC/GA cutoffs propagated to Pfam/Rfam                 */
  int    cutset[eslMSA_NCUTS];  /* TRUE if a cutoff is set; else FALSE                      */
  /*::cexcerpt::msa_optional::end::*/

  /* Info needed for maintenance of the data structure
   * (internal stuff.)
   */
  int      sqalloc;		/* # seqs currently allocated for           */
  int64_t *sqlen;               /* individual seq lengths during parsing    */
  int64_t *sslen;               /* individual ss lengths during parsing     */
  int64_t *salen;               /* individual sa lengths during parsing     */
  int64_t *pplen;               /* individual pp lengths during parsing     */
  int      lastidx;		/* last index we saw; use for guessing next */

  /* Optional information, especially Stockholm markup.
   * (The stuff we don't understand, but we can regurgitate.)
   *
   * That is, we know what type of information it is, but it's
   * either (interpreted as) free-text comment, or it's Stockholm
   * markup with unfamiliar tags.
   */
  char  **comment;              /* free text comments, or NULL      */
  int     ncomment;		/* number of comment lines          */
  int     alloc_ncomment;	/* number of comment lines alloc'ed */

  char  **gf_tag;               /* markup tags for unparsed #=GF lines  */
  char  **gf;                   /* annotations for unparsed #=GF lines  */
  int     ngf;			/* number of unparsed #=GF lines        */
  int     alloc_ngf;		/* number of gf lines alloc'ed          */

  char  **gs_tag;               /* markup tags for unparsed #=GS lines     */
  char ***gs;                   /* [0..ngs-1][0..nseq-1][free text] markup */
  int     ngs;                  /* number of #=GS tag types                */

  char  **gc_tag;               /* markup tags for unparsed #=GC lines  */
  char  **gc;                   /* [0..ngc-1][0..alen-1] markup         */
  int     ngc;                  /* number of #=GC tag types             */

  char  **gr_tag;               /* markup tags for unparsed #=GR lines     */
  char ***gr;                   /* [0..ngr-1][0..nseq-1][0..alen-1] markup */
  int     ngr;			/* number of #=GR tag types                */

  /* Optional augmentation w/ keyhashes.
   * This can significantly speed up parsing of large alignments
   * with many (>1,000) sequences.
   */
#ifdef eslAUGMENT_KEYHASH
  ESL_KEYHASH  *index;	        /* name ->seqidx hash table */
  ESL_KEYHASH  *gs_idx;         /* hash of #=GS tag types   */
  ESL_KEYHASH  *gc_idx;         /* hash of #=GC tag types   */
  ESL_KEYHASH  *gr_idx;         /* hash of #=GR tag types   */
#endif /*eslAUGMENT_KEYHASH*/

#ifdef eslAUGMENT_SSI
  off_t         offset;		/* disk offset to start of 1st line of this MSA's record */
#endif
} ESL_MSA;

/* Flags for msa->flags */
#define eslMSA_HASWGTS (1 << 0)  /* 1 if wgts were set, 0 if default 1.0's */
#define eslMSA_DIGITAL (1 << 1)	 /* if ax[][] is used instead of aseq[][]  */

/* Declarations of the API */

/* 1. The ESL_MSA object */
ESL_MSA *esl_msa_Create(int nseq, int64_t alen);
int      esl_msa_Expand(ESL_MSA *msa);
int      esl_msa_Copy (const ESL_MSA *msa, ESL_MSA *new);
ESL_MSA *esl_msa_Clone(const ESL_MSA *msa);
void     esl_msa_Destroy(ESL_MSA *msa);

/* 2. Digital mode MSA's (augmentation: alphabet) */
#ifdef eslAUGMENT_ALPHABET
int      esl_msa_GuessAlphabet(const ESL_MSA *msa, int *ret_type);
ESL_MSA *esl_msa_CreateDigital(const ESL_ALPHABET *abc, int nseq, int64_t alen);
int      esl_msa_Digitize(const ESL_ALPHABET *abc, ESL_MSA *msa, char *errmsg);
int      esl_msa_Textize(ESL_MSA *msa);
int      esl_msa_ConvertDegen2X(ESL_MSA *msa);
#endif /*eslAUGMENT_ALPHABET*/

/* 3. Setting or checking data fields in an ESL_MSA */
int esl_msa_SetName          (ESL_MSA *msa, const char *s, esl_pos_t n);
int esl_msa_SetDesc          (ESL_MSA *msa, const char *s, esl_pos_t n);
int esl_msa_SetAccession     (ESL_MSA *msa, const char *s, esl_pos_t n);
int esl_msa_SetAuthor        (ESL_MSA *msa, const char *s, esl_pos_t n);
int esl_msa_SetSeqName       (ESL_MSA *msa, int idx, const char *s, esl_pos_t n);
int esl_msa_SetSeqAccession  (ESL_MSA *msa, int idx, const char *s, esl_pos_t n);
int esl_msa_SetSeqDescription(ESL_MSA *msa, int idx, const char *s, esl_pos_t n);
int esl_msa_SetDefaultWeights(ESL_MSA *msa);

int esl_msa_FormatName          (ESL_MSA *msa, const char *name,    ...);
int esl_msa_FormatDesc          (ESL_MSA *msa, const char *desc,    ...);
int esl_msa_FormatAccession     (ESL_MSA *msa, const char *acc,     ...);
int esl_msa_FormatAuthor        (ESL_MSA *msa, const char *author,  ...);
int esl_msa_FormatSeqName       (ESL_MSA *msa, int idx, const char *name, ...);
int esl_msa_FormatSeqAccession  (ESL_MSA *msa, int idx, const char *acc, ...);
int esl_msa_FormatSeqDescription(ESL_MSA *msa, int idx, const char *desc, ...);

int esl_msa_AddComment(ESL_MSA *msa, char *p,   esl_pos_t n);
int esl_msa_AddGF     (ESL_MSA *msa, char *tag, esl_pos_t taglen,            char *value, esl_pos_t vlen);
int esl_msa_AddGS     (ESL_MSA *msa, char *tag, esl_pos_t taglen, int sqidx, char *value, esl_pos_t vlen);
int esl_msa_AppendGC  (ESL_MSA *msa, char *tag, char *value);
int esl_msa_AppendGR  (ESL_MSA *msa, char *tag, int sqidx, char *value);

int esl_msa_CheckUniqueNames(const ESL_MSA *msa);

/* 4. Miscellaneous functions for manipulating MSAs */
int esl_msa_ReasonableRF(ESL_MSA *msa, double symfrac, int useconsseq, char *rfline);
int esl_msa_MarkFragments(ESL_MSA *msa, double fragthresh);
int esl_msa_SequenceSubset(const ESL_MSA *msa, const int *useme, ESL_MSA **ret_new);
int esl_msa_ColumnSubset (ESL_MSA *msa, char *errbuf, const int *useme);
int esl_msa_MinimGaps    (ESL_MSA *msa, char *errbuf, const char *gaps, int consider_rf);
int esl_msa_MinimGapsText(ESL_MSA *msa, char *errbuf, const char *gaps, int consider_rf, int fix_bps);
int esl_msa_NoGaps       (ESL_MSA *msa, char *errbuf, const char *gaps);
int esl_msa_NoGapsText   (ESL_MSA *msa, char *errbuf, const char *gaps, int fix_bps);
int esl_msa_SymConvert(ESL_MSA *msa, const char *oldsyms, const char *newsyms);
int esl_msa_Checksum(const ESL_MSA *msa, uint32_t *ret_checksum);

int esl_msa_RemoveBrokenBasepairsFromSS(char *ss, char *errbuf, int len, const int *useme);
int esl_msa_RemoveBrokenBasepairs(ESL_MSA *msa, char *errbuf, const int *useme);

int esl_msa_ReverseComplement(ESL_MSA *msa);
#ifdef eslAUGMENT_KEYHASH
int esl_msa_Hash(ESL_MSA *msa);
#endif

/* 5. Debugging, testing, development */
int      esl_msa_Validate(const ESL_MSA *msa, char *errmsg);
ESL_MSA *esl_msa_CreateFromString(const char *s, int fmt);
int      esl_msa_Compare         (ESL_MSA *a1, ESL_MSA *a2);
int      esl_msa_CompareMandatory(ESL_MSA *a1, ESL_MSA *a2);
int      esl_msa_CompareOptional (ESL_MSA *a1, ESL_MSA *a2);
#endif /*eslMSA_INCLUDED*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/

/*** End of inlined file: esl_msa.h ***/


#endif
#if defined eslAUGMENT_RANDOM && defined eslAUGMENT_RANDOMSEQ


/*** Start of inlined file: esl_randomseq.h ***/
#ifndef eslRANDOMSEQ_INCLUDED
#define eslRANDOMSEQ_INCLUDED

/* Control flag passed to esl_rsq_Sample():                          */
#define eslRSQ_SAMPLE_ALNUM  1	/* isalpha | isdigit                 */
#define eslRSQ_SAMPLE_ALPHA  2	/* islower | isupper                 */
#define eslRSQ_SAMPLE_LOWER  3	/* ASCII: a-z                        */
#define eslRSQ_SAMPLE_UPPER  4	/* ASCII: A-Z                        */
#define eslRSQ_SAMPLE_DIGIT  5	/* 0-9                               */
#define eslRSQ_SAMPLE_XDIGIT 6	/* 0-9, a-f, A-F                     */
#define eslRSQ_SAMPLE_CNTRL  7	/* ASCII: 0..0x1F, plus 0x7F (DEL)   */
#define eslRSQ_SAMPLE_GRAPH  8  /* isprint && ! ' ' (space)          */
#define eslRSQ_SAMPLE_SPACE  9	/* ' ', '\f', '\n', '\r', '\t', '\v' */
#define eslRSQ_SAMPLE_BLANK  10	/* ' ', '\t'                         */
#define eslRSQ_SAMPLE_PRINT  11 /* ASCII: 0x20 ' ' through 0x7E '~'  */
#define eslRSQ_SAMPLE_PUNCT  12	/* isprint && !(isspace || isalnum)  */

/* 1. Generating simple random character strings. */
int esl_rsq_Sample(ESL_RANDOMNESS *rng, int allowed_chars_flag, int L, char **ret_s);

/* 2. Generating iid sequences. */
int esl_rsq_IID  (ESL_RANDOMNESS *r, const char *alphabet, const double *p, int K, int L, char *s);
int esl_rsq_fIID (ESL_RANDOMNESS *r, const char *alphabet, const float  *p, int K, int L, char *s);

/* 3. Shuffling sequences. */
int esl_rsq_CShuffle       (ESL_RANDOMNESS *r, const char *s,        char *shuffled);
int esl_rsq_CShuffleDP     (ESL_RANDOMNESS *r, const char *s,        char *shuffled);
int esl_rsq_CShuffleKmers  (ESL_RANDOMNESS *r, const char *s, int K, char *shuffled);
int esl_rsq_CReverse       (const char *s, char *rev);
int esl_rsq_CShuffleWindows(ESL_RANDOMNESS *r, const char *s, int w, char *shuffled);

/* 4. Randomizing sequences */
int esl_rsq_CMarkov0  (ESL_RANDOMNESS *r, const char *s, char *markoved);
int esl_rsq_CMarkov1  (ESL_RANDOMNESS *r, const char *s, char *markoved);

/* 5. Generating iid sequences (digital mode). */
int esl_rsq_xIID       (ESL_RANDOMNESS *r, const double *p, int K, int L, ESL_DSQ *dsq);
int esl_rsq_xfIID      (ESL_RANDOMNESS *r, const float  *p, int K, int L, ESL_DSQ *dsq);
int esl_rsq_SampleDirty(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, double **byp_p, int L, ESL_DSQ *dsq);

/* 6. Shuffling sequences (digital mode). */
int esl_rsq_XShuffle       (ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L,        ESL_DSQ *shuffled);
int esl_rsq_XShuffleDP     (ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, int K, ESL_DSQ *shuffled);
int esl_rsq_XShuffleKmers  (ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, int K, ESL_DSQ *shuffled);
int esl_rsq_XReverse(const ESL_DSQ *dsq, int L, ESL_DSQ *rev);
int esl_rsq_XShuffleWindows(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, int w, ESL_DSQ *shuffled);

/* 7. Randomizing sequences (digital mode) */
int esl_rsq_XMarkov0  (ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, int K, ESL_DSQ *markoved);
int esl_rsq_XMarkov1  (ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, int K, ESL_DSQ *markoved);

#endif /*eslRANDOMSEQ_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/*** End of inlined file: esl_randomseq.h ***/

#endif

/* ESL_SQ - a biosequence
 *
 * Can be either in text mode <seq>, or digital mode <dsq>.
 * One of them has to be NULL, and the other contains the data.
 *
 * When in text mode, <ss> and <seq> can hold up to <n=salloc-1>
 * residues and a terminal '\0', and both are indexed <0..n-1>.
 *
 * When in digital mode, <ss> and <dsq> can hold up to <n=salloc-2>
 * residues; both are indexed <1..n>, and positions 0 and n+1 are
 * sentinel bytes. The digital sequence <dsq> uses <eslDSQ_SENTINEL>
 * as its sentinels; as a hack, <ss> uses '\0' as sentinels.  This
 * means that <sq->ss+1> is a valid NUL-terminated C string, but
 * <sq->ss> itself would be a string of length 0 because of the
 * leading NUL sentinel. Watch out for this.
 *
 * To save on allocation calls, the structure is designed to be reused
 * for subsequent sequences, rather than free'd and reallocated -
 * thus, we keep track of the allocated sizes of all the strings.
 *
 * Notes on when we need to reallocate:
 *    - In a text mode sequence (seq 0..n-1), byte salloc-1 is
 *      reserved for the NUL, so the sequence is full when
 *      n == salloc-1.
 *
 *    - In a digital mode sequence (dsq 1..n), bytes 0 and salloc-1
 *      are reserved for sentinel bytes, so the reallocation condition
 *      is when n == salloc-2.
 *
 * At least for now, the only way to set the <ss> structure annotation
 * field is by a CreateFrom(), by extraction from an MSA, or manually
 * (by directly allocating a sequence's <ss> field).
 *
 * A sequence object will usually be holding a complete (full length)
 * sequence. Three other cases arise less frequently:
 *
 * 1. We're a subsequence extracted from a source sequence.
 *    <sourcename> is the name of the source.
 *    <L> is the length of the source (and coords are 1..L).
 *    The subsequence is from <start>..<end> on the source.
 *    The length of the subsequence <n> is abs(<end>-<start>)+1.
 *    <start> can be greater than <end> for a nucleic acid sequence;
 *    in this case, the subsequence is reverse complemented.
 *
 * 2. We're a window on a source sequence.
 *    This is similar to being a subsequence, with the added
 *    wrinkle that we're scanning over a long source sequence
 *    in overlapping segments, defined by a "previous context"
 *    <C> and a "new window" <W> (the whole sequence is n=C+W
 *    residues long):
 *                       s  C          W      e
 *    current window:    |------||------------|
 *    next window read:                |------||------------|
 *                                     s  C           W     e
 *    Here, dsq[1..n] is source[s..e]; each newly read
 *    window starts at dsq[C+1], and is preceded by C
 *    residues of context.
 *
 * 3. We're just after information about the sequence, not the
 *    sequence itself; everything except the per-residue information
 *    (such as <dsq/seq> and <ss>). We do this when SSI indexing,
 *    for example, so we don't have to read entire huge seqs into
 *    memory just to calculate their lengths for the index.
 *
 * Note/TODO: use of "\0" empty string to indicate lack of optional
 * acc, desc info is now deprecated. Cannot distinguish empty string
 * from lack of annotation. Should use NULL ptr instead. Fix this in
 * future.  (21 Nov 09 xref J5/114)
 *
 */
typedef struct {
  /*::cexcerpt::sq_sq::begin::*/
  char    *name;           /* name; one word, no whitespace ("\0" if no name)  */
  char    *acc;            /* optional accession (1 word) ("\0" if none)       */
  char    *desc;           /* description line ("\0" if no description)        */
  int32_t  tax_id;         /* NCBI taxonomy id (-1 if none)                    */
  char    *seq;            /* sequence [0..n-1], or NULL if digital            */
  ESL_DSQ *dsq;            /* digitized sequence [1..n], or NULL if text       */
  char    *ss;             /* optional sec structure [0..n-1], [1..n], or NULL */
  int64_t  n;              /* length of seq (or dsq) and ss                    */
  /*::cexcerpt::sq_sq::end::*/

  /* Coordinate info for:                                       seq       subseq     window     info */
  /*                                                           ----       ------     ------    ----- */
  int64_t  start;  /* coord of seq[0],dsq[1] on source  [1..L]    1      1<=i<=L    1<=i<=L      0   */
  int64_t  end;    /* coord of seq[n-1],dsq[n] on source[1..L]    L      1<=j<=L    1<=j<=L      0   */
  int64_t  C;      /* # of context residues for a window          0            0        n-W      0   */
  int64_t  W;      /* window width                                L            n        n-C      0   */
  int64_t  L;      /* source sequence length in residues          L     L (or -1)   L (or -1)    L   */
  /* and   n: length of seq (or dsq) and ss actually stored:      L   abs(j-i)+1        C+W      0   */
  /* In all the above bookkeeping, a -1 means "unknown" */
  char    *source; /* name of the source of a subseq/window; or MSA name; or ""*/

  /* Memory allocation bookkeeping:  (all inclusive of \0;  >= strlen()+1)     */
  int      nalloc;         /* allocated length of name                         */
  int      aalloc;         /* allocated length of accession                    */
  int      dalloc;         /* allocated length of description                  */
  int64_t  salloc;         /* alloc for seq or dsq, and ss if present          */
  int      srcalloc;	   /* allocated length for source name                 */

  /* Disk offset bookkeeping:                                                  */
  int64_t  idx;           /* ctr for which # seq this is; -1 if not counting   */
  off_t    roff;          /* record offset (start of record); -1 if none       */
  off_t    hoff;          /* offset to last byte of header; -1 if unknown      */
  off_t    doff;          /* data offset (start of sequence data); -1 if none  */
  off_t    eoff;          /* offset to last byte of record; -1 if unknown      */

  /* Optional information for extra residue markups.
   * The number of them, and their tags are arbitrary
   */
  char  **xr_tag;          /* markup tags for extra residue markups [0..ntr-1][free-text], [0..ntr-1][free-text], or NULL */
  char  **xr;              /* annotations for extra residue markups [0..ntr-1][0..n-1],    [0..ntr-1][1..n],      or NULL */
  int     nxr;             /* number of extra residue markups                                                             */

  /* Copy of a pointer to the alphabet, if digital mode */
#if defined(eslAUGMENT_ALPHABET)
  const ESL_ALPHABET *abc; /* reference to the alphabet for <dsq>              */
#else
  const void         *abc; /* void reference, if we're not even augmented      */
#endif
} ESL_SQ;

typedef struct {
  int      count;       /* number of <ESL_SQ> objects in the block */
  int      listSize;    /* maximum number elements in the list     */
  int      complete;    /*TRUE if the the final ESL_SQ element on the block is complete, FALSE if it's only a partial winow of the full sequence*/
  int64_t  first_seqidx;/*unique identifier of the first ESL_SQ object on list;  the seqidx of the i'th entry on list is first_seqidx+i */
  ESL_SQ  *list;        /* array of <ESL_SQ> objects               */
} ESL_SQ_BLOCK;

/* These control default initial allocation sizes in an ESL_SQ.     */
#define eslSQ_NAMECHUNK   32	// allocation unit for name, source
#define eslSQ_ACCCHUNK    32	// allocation unit for accession
#define eslSQ_DESCCHUNK  128	// allocation unit for description
#define eslSQ_SEQCHUNK   256	// allocation unit for seqs
								//  .. dsqdata assumes _SEQCHUNK >= 4

ESL_SQ *esl_sq_Create(void);
ESL_SQ *esl_sq_CreateFrom(const char *name, const char *seq,
				 const char *desc, const char *acc, const char *ss);
int     esl_sq_Grow  (ESL_SQ *sq, int64_t *ret_nsafe);
int     esl_sq_GrowTo(ESL_SQ *sq, int64_t  n);
int     esl_sq_Copy(const ESL_SQ *src, ESL_SQ *dst);
int     esl_sq_Compare  (ESL_SQ *sq1, ESL_SQ *sq2);
int     esl_sq_Reuse    (ESL_SQ *sq);
int     esl_sq_IsDigital(const ESL_SQ *sq);
int     esl_sq_IsText   (const ESL_SQ *sq);
void    esl_sq_Destroy  (ESL_SQ *sq);

int     esl_sq_SetName        (ESL_SQ *sq, const char *name);
int     esl_sq_SetAccession   (ESL_SQ *sq, const char *acc);
int     esl_sq_SetDesc        (ESL_SQ *sq, const char *desc);
int     esl_sq_SetSource      (ESL_SQ *sq, const char *source);
int     esl_sq_FormatName     (ESL_SQ *sq, const char *name,   ...);
int     esl_sq_FormatAccession(ESL_SQ *sq, const char *acc,    ...);
int     esl_sq_FormatDesc     (ESL_SQ *sq, const char *desc,   ...);
int     esl_sq_FormatSource   (ESL_SQ *sq, const char *source, ...);
int     esl_sq_AppendDesc     (ESL_SQ *sq, const char *desc);
int     esl_sq_SetCoordComplete(ESL_SQ *sq, int64_t L);
int     esl_sq_CAddResidue (ESL_SQ *sq, char c);
int     esl_sq_ReverseComplement(ESL_SQ *sq);
int     esl_sq_Checksum(const ESL_SQ *sq, uint32_t *ret_checksum);
int     esl_sq_CountResidues(const ESL_SQ *sq, int start, int L, float *f);

#ifdef eslAUGMENT_ALPHABET
ESL_SQ *esl_sq_CreateDigital(const ESL_ALPHABET *abc);
ESL_SQ *esl_sq_CreateDigitalFrom(const ESL_ALPHABET *abc, const char *name, const ESL_DSQ *dsq,
					int64_t L, const char *desc, const char *acc,  const char *ss);
int     esl_sq_Digitize(const ESL_ALPHABET *abc, ESL_SQ *sq);
int     esl_sq_Textize(ESL_SQ *sq);
int     esl_sq_GuessAlphabet(ESL_SQ *sq, int *ret_type);
int     esl_sq_XAddResidue(ESL_SQ *sq, ESL_DSQ x);
int     esl_sq_ConvertDegen2X(ESL_SQ *sq);
#endif

#ifdef eslAUGMENT_MSA
int     esl_sq_GetFromMSA  (const ESL_MSA *msa, int which, ESL_SQ *sq);
int     esl_sq_FetchFromMSA(const ESL_MSA *msa, int which, ESL_SQ **ret_sq);
#endif

ESL_SQ_BLOCK *esl_sq_CreateBlock(int count);
int esl_sq_BlockGrowTo(ESL_SQ_BLOCK *sqblock, int newsize, int do_digital, const ESL_ALPHABET *abc);
#ifdef eslAUGMENT_ALPHABET
ESL_SQ_BLOCK *esl_sq_CreateDigitalBlock(int count, const ESL_ALPHABET *abc);
#endif
void          esl_sq_DestroyBlock(ESL_SQ_BLOCK *sqBlock);

#if defined eslAUGMENT_RANDOM && defined eslAUGMENT_RANDOMSEQ
int esl_sq_Sample(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, int maxL, ESL_SQ **ret_sq);
#endif /* eslAUGMENT_RANDOM && eslAUGMENT_RANDOMSEQ */

#endif /*eslSQ_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_sq.h ***/


#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef eslAUGMENT_MSA


/*** Start of inlined file: esl_msafile.h ***/
#ifndef eslMSAFILE_INCLUDED
#define eslMSAFILE_INCLUDED

#include <stdio.h>

/* Object: ESL_MSAFILE_FMTDATA
 *
 * Additional (often optional) information about variants of some file
 * formats. Not much in here right now - but figured this might need
 * to expand in the future, best to have the mechanism in place.
 *
 * Used in three ways:
 *   1. When opening an MSA file in a known format (as opposed to
 *      guessing an unknown format), caller may provide an <ESL_MSAFILE_FMTDATA>
 *      structure containing any additional constraints on the format.
 *      The new <afp> will copy this information into <afp->fmtd>.
 *   2. When opening an MSA file in an unknown format (calling GuessFileFormat()),
 *      format-specific autodetectors fill in <afp->fmtd> with any additional
 *      constraints.
 *   3. When writing an MSA file, caller may provide additional constraints on
 *      the format; notably <fmtd->rpl>, the number of residues per line,
 *      used for many formats.
 *
 * TODO: If this fills up with more information, we should eventually
 *       consolidate the format code too; create ESL_MSAFORMAT structure
 *       to hold both integer code and optional information; implement
 *       it in esl_msaformat.[ch]; put format guessing routines there;
 *       rename eslMSAFILE_* -> eslMSAFORMAT_*. For now, not worth the
 *       time, because it's really only a placeholder dealing with a small
 *       PHYLIP-specific format issue. <format>, <fmtd> are generally
 *       an ordered pair, to facilitate eventual replacement w/ single
 *       <fmt>. [SRE, 19 Jul 11]
 */
typedef struct {
  int namewidth;   /* PHYLIP only:     width of the name field (usually 10, but can vary) unset=0 */
  int rpl;	   /* several formats: residues per line                                  unset=0 */
} ESL_MSAFILE_FMTDATA;

/* Object: ESL_MSAFILE
 *
 * An alignment file open for parsing.
 */
typedef struct {
  ESL_BUFFER          *bf;            /* input file/data being parsed                          */

  int32_t              format;	      /* format of alignment file we're reading                */
  ESL_MSAFILE_FMTDATA  fmtd;          /* additional (often optional) format-specific details.  */

  char                *line;	      /* line read from <bf> by <esl_msafile_GetLine()>        */
  esl_pos_t            n;	      /* length of line in bytes (line is not NUL-terminated)  */
  int64_t              linenumber;    /* input linenumber for diagnostics; -1 if we lose track */
  esl_pos_t            lineoffset;    /* offset of start of <line> in <bf>; -1 if line unset   */

  ESL_DSQ              inmap[128];    /* input map, 0..127                                     */
  const ESL_ALPHABET  *abc;	      /* non-NULL if augmented and in digital mode             */
  ESL_SSI             *ssi;	      /* open SSI index; or NULL, if none or not augmented     */
  char                 errmsg[eslERRBUFSIZE];   /* user-directed message for normal errors     */
} ESL_MSAFILE;

/* Alignment file format codes.
 * Must coexist with sqio unaligned file format codes.
 * Rules:
 *     - 0 is an unknown/unassigned format
 *     - <=100 reserved for unaligned formats
 *     - >100 reserved for aligned formats
 */
#define eslMSAFILE_UNKNOWN     0    /* unknown format                              */
#define eslMSAFILE_STOCKHOLM   101  /* Stockholm format, interleaved               */
#define eslMSAFILE_PFAM        102  /* Pfam/Rfam one-line-per-seq Stockholm format */
#define eslMSAFILE_A2M         103  /* UCSC SAM's fasta-like a2m format            */
#define eslMSAFILE_PSIBLAST    104  /* NCBI PSI-BLAST alignment format             */
#define eslMSAFILE_SELEX       105  /* old SELEX format (largely obsolete)         */
#define eslMSAFILE_AFA         106  /* aligned FASTA format                        */
#define eslMSAFILE_CLUSTAL     107  /* CLUSTAL format                              */
#define eslMSAFILE_CLUSTALLIKE 108  /* CLUSTAL-like formats (MUSCLE, PROBCONS)     */
#define eslMSAFILE_PHYLIP      109  /* interleaved PHYLIP format                   */
#define eslMSAFILE_PHYLIPS     110  /* sequential PHYLIP format                    */

/* 1. Opening/closing an ESL_MSAFILE */
int   esl_msafile_Open      (ESL_ALPHABET **byp_abc, const char *msafile, const char *env, int format, ESL_MSAFILE_FMTDATA *fmtd, ESL_MSAFILE **ret_afp);
int   esl_msafile_OpenMem   (ESL_ALPHABET **byp_abc, const char *p, esl_pos_t n,           int format, ESL_MSAFILE_FMTDATA *fmtd, ESL_MSAFILE **ret_afp);
int   esl_msafile_OpenBuffer(ESL_ALPHABET **byp_abc, ESL_BUFFER *bf,                       int format, ESL_MSAFILE_FMTDATA *fmtd, ESL_MSAFILE **ret_afp);
void  esl_msafile_OpenFailure(ESL_MSAFILE *afp, int status);
int   esl_msafile_SetDigital (ESL_MSAFILE *afp, const ESL_ALPHABET *abc);
void  esl_msafile_Close(ESL_MSAFILE *afp);

/* 2. ESL_MSAFILE_FMTDATA: optional extra constraints on formats */
int   esl_msafile_fmtdata_Init(ESL_MSAFILE_FMTDATA *fmtd);
int   esl_msafile_fmtdata_Copy(ESL_MSAFILE_FMTDATA *src,  ESL_MSAFILE_FMTDATA *dst);

/* 3. Utilities for different file formats */
int   esl_msafile_GuessFileFormat(ESL_BUFFER *bf, int *ret_fmtcode, ESL_MSAFILE_FMTDATA *fmtd, char *errbuf);
int   esl_msafile_IsMultiRecord(int fmt);
int   esl_msafile_EncodeFormat(char *fmtstring);
char *esl_msafile_DecodeFormat(int fmt);

/* 4. Utilities for different alphabets */
#ifdef eslAUGMENT_ALPHABET
int esl_msafile_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type);
#endif

/* 5. Random access in a MSA flatfile database */
#ifdef eslAUGMENT_SSI
int esl_msafile_PositionByKey(ESL_MSAFILE *afp, const char *key);
#endif

/* 6. Reading an MSA from an ESL_MSAFILE */
int  esl_msafile_Read(ESL_MSAFILE *afp, ESL_MSA **ret_msa);
void esl_msafile_ReadFailure(ESL_MSAFILE *afp, int status);

/* 7. Writing an MSA to a stream */
int esl_msafile_Write(FILE *fp, ESL_MSA *msa, int fmt);

/* 8. Utilities for specific parsers */
int esl_msafile_GetLine(ESL_MSAFILE *afp, char **opt_p, esl_pos_t *opt_n);
int esl_msafile_PutLine(ESL_MSAFILE *afp);


/*** Start of inlined file: esl_msafile_a2m.h ***/
#ifndef eslMSAFILE_A2M_INCLUDED
#define eslMSAFILE_A2M_INCLUDED

int esl_msafile_a2m_SetInmap     (ESL_MSAFILE *afp);
int esl_msafile_a2m_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type);
int esl_msafile_a2m_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
int esl_msafile_a2m_Write        (FILE *fp,    const ESL_MSA *msa);

#endif /* eslMSAFILE_A2M_INCLUDED */

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_msafile_a2m.h ***/


/*** Start of inlined file: esl_msafile_afa.h ***/
#ifndef eslMSAFILE_AFA_INCLUDED
#define eslMSAFILE_AFA_INCLUDED

int esl_msafile_afa_SetInmap     (ESL_MSAFILE *afp);
int esl_msafile_afa_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type);
int esl_msafile_afa_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
int esl_msafile_afa_Write        (FILE *fp, const ESL_MSA *msa);

#endif /* eslMSAFILE_AFA_INCLUDED */

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_msafile_afa.h ***/


/*** Start of inlined file: esl_msafile_clustal.h ***/
#ifndef eslMSAFILE_CLUSTAL_INCLUDED
#define eslMSAFILE_CLUSTAL_INCLUDED

int esl_msafile_clustal_SetInmap     (ESL_MSAFILE *afp);
int esl_msafile_clustal_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type);
int esl_msafile_clustal_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
int esl_msafile_clustal_Write        (FILE *fp,    const ESL_MSA *msa, int fmt);

#endif /* eslMSAFILE_CLUSTAL_INCLUDED */

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_msafile_clustal.h ***/


/*** Start of inlined file: esl_msafile_phylip.h ***/
#ifndef eslMSAFILE_PHYLIP_INCLUDED
#define eslMSAFILE_PHYLIP_INCLUDED

int esl_msafile_phylip_SetInmap     (ESL_MSAFILE *afp);
int esl_msafile_phylip_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type);
int esl_msafile_phylip_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
int esl_msafile_phylip_Write        (FILE *fp, const ESL_MSA *msa, int format, ESL_MSAFILE_FMTDATA *opt_fmtd);

int esl_msafile_phylip_CheckFileFormat(ESL_BUFFER *bf, int *ret_format, int *ret_namewidth);

#endif /* eslMSAFILE_PHYLIP_INCLUDED */
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_msafile_phylip.h ***/


/*** Start of inlined file: esl_msafile_psiblast.h ***/
#ifndef eslMSAFILE_PSIBLAST_INCLUDED
#define eslMSAFILE_PSIBLAST_INCLUDED

int esl_msafile_psiblast_SetInmap     (ESL_MSAFILE *afp);
int esl_msafile_psiblast_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type);
int esl_msafile_psiblast_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
int esl_msafile_psiblast_Write        (FILE *fp, const ESL_MSA *msa);

#endif /* eslMSAFILE_PSIBLAST_INCLUDED */
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_msafile_psiblast.h ***/


/*** Start of inlined file: esl_msafile_selex.h ***/
#ifndef eslMSAFILE_SELEX_INCLUDED
#define eslMSAFILE_SELEX_INCLUDED

int esl_msafile_selex_SetInmap     (ESL_MSAFILE *afp);
int esl_msafile_selex_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type);
int esl_msafile_selex_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
int esl_msafile_selex_Write        (FILE *fp,    const ESL_MSA *msa);

#endif /* eslMSAFILE_SELEX_INCLUDED */

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_msafile_selex.h ***/


/*** Start of inlined file: esl_msafile_stockholm.h ***/
#ifndef eslMSAFILE_STOCKHOLM_INCLUDED
#define eslMSAFILE_STOCKHOLM_INCLUDED

int esl_msafile_stockholm_SetInmap     (ESL_MSAFILE *afp);
int esl_msafile_stockholm_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type);
int esl_msafile_stockholm_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
int esl_msafile_stockholm_Write        (FILE *fp, const ESL_MSA *msa, int fmt);

#endif /*eslMSAFILE_STOCKHOLM_INCLUDED*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_msafile_stockholm.h ***/

#endif /*eslMSAFILE_INCLUDED*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_msafile.h ***/

#endif

/* set the max residue count to 1 meg when reading a block */
#define MAX_RESIDUE_COUNT (1024 * 1024)

/* forward declaration */
struct esl_sqio_s;

/* ESL_SQASCII:
 * An open sequence file for reading.
 */
typedef struct esl_sqascii_s {

  FILE *fp;           	      /* Open file ptr                            */
  char  errbuf[eslERRBUFSIZE];/* parse error mesg.  Size must match msa.h */

  int   do_gzip;	      /* TRUE if we're reading from gzip -dc pipe */
  int   do_stdin;	      /* TRUE if we're reading from stdin         */
  int   do_buffer;            /* TRUE if we're reading from a buffer      */

  /* all input first gets buffered in memory; this gives us enough
   * recall to use Guess*() functions even in nonrewindable streams
   */
  char    *mem;		      /* buffered input                           */
  int      allocm;	      /* <mem> size, multiples of eslREADBUFSIZE  */
  int      mn;		      /* number of chars in <mem> (up to allocm)  */
  int      mpos;	      /* pos of next <buf> to load from <mem>     */
  off_t    moff;	      /* disk offset to start of <mem>            */
  int      is_recording;      /* TRUE if we need to keep buffering more   */

  /* input is either character-based [fread()] or line-based (esl_fgets())*/
  char    *buf;		      /* buffer for fread() or fgets() input      */
  off_t    boff;	      /* disk offset to start of buffer           */
  int      balloc;	      /* allocated size of buf                    */
  int      nc;		      /* #chars in buf (usually full, less at EOF)*/
  int      bpos;	      /* current position in the buffer (0..nc-1) */
  int64_t  L;		      /* #residues seen so far in current seq     */
  int64_t  linenumber;	      /* What line of the file  (1..N; -1=unknown)*/
  off_t    bookmark_offset;   /* bookmark fwd position before reversing...*/
  int64_t  bookmark_linenum;  /* in both linenumber and disk offset       */

  /* Format-specific configuration                                           */
  int   is_linebased;	      /* TRUE for fgets() parsers; FALSE for fread() */
  int   eof_is_ok;	      /* TRUE if record can end on EOF               */
  int  (*parse_header)(struct esl_sqio_s *, ESL_SQ *sq);
  int  (*skip_header) (struct esl_sqio_s *, ESL_SQ *sq);
  int  (*parse_end)   (struct esl_sqio_s *, ESL_SQ *sq);

  /* MSA augmentation confers reading MSA files as sequential seq files. */
#if defined(eslAUGMENT_MSA)
  ESL_MSAFILE  *afp;	      /* open ESL_MSAFILE for reading           */
  ESL_MSA      *msa;	      /* preloaded alignment to draw seqs from  */
  int           idx;	      /* index of next seq to return, 0..nseq-1 */
#else
  void        *afp; 	      /* NULL */
  void        *msa;           /* NULL */
  int          idx;           /* 0    */
#endif /*eslAUGMENT_MSA*/

  /* SSI augmentation confers random access of records in a seq file        */
  char    *ssifile;	      /* path to expected SSI index file            */
  int      rpl;		      /* residues per line in file; -1=unset 0=inval*/
  int      bpl;		      /* bytes per line in file; -1=unset, 0=inval  */
  int      currpl;	      /* residues on current line (-1=unknown)      */
  int      curbpl;	      /* bytes on current line    (-1=unknown)      */
  int      prvrpl;	      /* residues on previous line                  */
  int      prvbpl;	      /* bytes on previous line                     */
#if defined(eslAUGMENT_SSI)
  ESL_SSI *ssi;		/* open ESL_SSI index, or NULL if none     */
#else
  void    *ssi;		/* NULL */
#endif /*eslAUGMENT_SSI*/
} ESL_SQASCII_DATA;

int  esl_sqascii_Open(char *seqfile, int format, struct esl_sqio_s *sqfp);
int  esl_sqascii_WriteFasta(FILE *fp, ESL_SQ *s, int update);
int  esl_sqascii_Parse(char *buf, int size, ESL_SQ *s, int format);

#endif /*eslSQIO_ASCII_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_sqio_ascii.h ***/

#ifdef eslAUGMENT_NCBI

/*** Start of inlined file: esl_sqio_ncbi.h ***/
#ifndef eslSQIO_NCBI_INCLUDED
#define eslSQIO_NCBI_INCLUDED

#include <stdio.h>

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

/* forward declaration */
struct esl_sqio_s;

/* set the max residue count to 1 meg when reading a block */
#define MAX_RESIDUE_COUNT (1024 * 1024)

#define MAX_DB_VOLUMES   100

/* ESL_SQNCBI_VOLUME:
 * Information for the volume
 */
typedef struct esl_sqncbi_vol_s {
  char      *name;                 /* name of the volume                       */

  uint32_t   start_seq;            /* starting sequence number                 */
  uint32_t   end_seq;              /* ending sequence number                   */

  uint32_t   hdr_off;              /* disk offset in .pin to header index      */
  uint32_t   seq_off;              /* disk offset to .pin to sequence index    */
  uint32_t   amb_off;              /* disk offset to .pin to ambiguous index   */
} ESL_SQNCBI_VOLUME;

/* ESL_SQNCBI:
 * An open sequence file for reading.
 */
typedef struct esl_sqncbi_s {
  FILE      *fppin;                /* Open .pin file ptr                       */
  FILE      *fpphr;                /* Open .phr file ptr                       */
  FILE      *fppsq;                /* Open .psq file ptr                       */
  char       errbuf[eslERRBUFSIZE];/* parse error mesg.  Size must match msa.h */

  char      *title;                /* database title                           */
  int        version;              /* database version                         */
  char      *timestamp;            /* time stamp of database creation          */

  uint32_t   num_seq;              /* number of sequences in the database      */
  uint64_t   total_res;            /* total number of residues                 */
  uint32_t   max_seq;              /* longest sequence in the database         */

  uint32_t   hdr_off;              /* disk offset in .pin to header index      */
  uint32_t   seq_off;              /* disk offset to .pin to sequence index    */
  uint32_t   amb_off;              /* disk offset to .pin to ambiguous index   */

  int        index;                /* current sequence index in the database   */
  uint32_t   vol_index;            /* current volume index (-1 if no volumes)  */
  uint32_t   roff;                 /* record offset (start of header)          */
  uint32_t   hoff;                 /* offset to last byte of header            */
  uint32_t   doff;                 /* data offset (start of sequence data)     */
  uint32_t   eoff;                 /* offset to last byte of sequence          */

  uint32_t   index_start;          /* start of indexes currently loaded        */
  uint32_t   index_end;            /* end of indexes currently loaded          */
  uint32_t  *hdr_indexes;          /* block of header indexes from .pin        */
  uint32_t  *seq_indexes;          /* block of header indexes from .pin        */
  uint32_t  *amb_indexes;          /* block of header indexes from .pin        */

  /* volume information */
  uint32_t   volumes;              /* number of volumes                        */
  ESL_SQNCBI_VOLUME vols[MAX_DB_VOLUMES];

  /* information for the current header */
  unsigned char *hdr_buf;          /* buffer for holding unparsed header       */
  unsigned char *hdr_ptr;          /* current parser position                  */
  int            hdr_alloced;      /* size of the allocated buffer             */

  char          *name_ptr;         /* pointer to name NOT NULL TERMINATED      */
  int32_t        name_size;        /* length of the name                       */
  char          *acc_ptr;          /* pointer to accession NOT NULL TERMINATED */
  int32_t        acc_size;         /* length of the accession                  */
  int32_t        int_id;           /* integer sequence id                      */
  char          *str_id_ptr;       /* pointer to id NOT NULL TERMINATED        */
  int32_t        str_id_size;      /* length of the id                         */

  /* information on the current sequence */
  uint32_t       seq_apos;         /* position of ambiguity table              */
  uint32_t       seq_alen;         /* size of ambiguity table                  */
  uint32_t       seq_cpos;         /* current position in ambiguity table      */
  int32_t        seq_L;            /* true sequence length                     */

  /* alphabet used to convert ncbi to hmmer to ascii */
  int            alphatype;        /* amino or dna                             */
  char          *alphasym;         /* string of residues                       */

} ESL_SQNCBI_DATA;

int  esl_sqncbi_Open(char *seqfile, int format, struct esl_sqio_s *sqfp);

#endif /*eslSQIO_NCBI_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_sqio_ncbi.h ***/


#endif

#ifdef eslAUGMENT_ALPHABET

#endif
#ifdef eslAUGMENT_MSA

#endif

/*::cexcerpt::sq_sqio_data::begin::*/
/* ESL_SQDATA:
 * Data for different sequence formats.
 */
typedef union {
  ESL_SQASCII_DATA ascii;
#ifdef eslAUGMENT_NCBI
  ESL_SQNCBI_DATA  ncbi;
#endif
} ESL_SQDATA;
/*::cexcerpt::sq_sqio_data::end::*/

/* ESL_SQFILE:
 * An open sequence file for reading.
 */
typedef struct esl_sqio_s {
  char *filename;	      /* Name of file (for diagnostics)           */

  /* In digital mode, we have an alphabet ptr                             */
  int   do_digital;	      /* TRUE if we're reading in digital mode    */
#if defined(eslAUGMENT_ALPHABET)
  const ESL_ALPHABET *abc;
#else
  void               *abc;
#endif

  /* Format-specific configuration                                        */
  int     format;	      /* Format code of this file                 */
  ESL_DSQ inmap[128];	      /* an input map, 0..127                     */

  /* function pointers to format specific routines                        */
  int   (*position)        (struct esl_sqio_s *sqfp, off_t offset);
  void  (*close)           (struct esl_sqio_s *sqfp);

  int   (*set_digital)     (struct esl_sqio_s *sqfp, const ESL_ALPHABET *abc);
  int   (*guess_alphabet)  (struct esl_sqio_s *sqfp, int *ret_type);

  int   (*read)            (struct esl_sqio_s *sqfp, ESL_SQ *sq);
  int   (*read_info)       (struct esl_sqio_s *sqfp, ESL_SQ *sq);
  int   (*read_seq)        (struct esl_sqio_s *sqfp, ESL_SQ *sq);
  int   (*read_window)     (struct esl_sqio_s *sqfp, int C, int W, ESL_SQ *sq);
  int   (*echo)            (struct esl_sqio_s *sqfp, const ESL_SQ *sq, FILE *ofp);

  int   (*read_block)      (struct esl_sqio_s *sqfp, ESL_SQ_BLOCK *sqBlock, int max_residues, int max_sequences, int long_target);

#ifdef eslAUGMENT_SSI
  int   (*open_ssi)        (struct esl_sqio_s *sqfp, const char *ssifile_hint);
  int   (*pos_by_key)      (struct esl_sqio_s *sqfp, const char *key);
  int   (*pos_by_number)   (struct esl_sqio_s *sqfp, int which);

  int   (*fetch)           (struct esl_sqio_s *sqfp, const char *key, ESL_SQ *sq);
  int   (*fetch_info)      (struct esl_sqio_s *sqfp, const char *key, ESL_SQ *sq);
  int   (*fetch_subseq)    (struct esl_sqio_s *sqfp, const char *source, int64_t start, int64_t end, ESL_SQ *sq);
#endif

  int   (*is_rewindable)   (const struct esl_sqio_s *sqfp);
  const char *(*get_error) (const struct esl_sqio_s *sqfp);

  ESL_SQDATA data;            /* format specific data                     */
} ESL_SQFILE;

#if defined(eslAUGMENT_ALPHABET)
/* ESL_SQCACHE:
 * A entire database cached into memory.
 */
typedef struct esl_sqcache_s {
  char               *filename;    /* Name of file (for diagnostics)              */
  int                 format;      /* Format code of this file                    */

  const ESL_ALPHABET *abc;         /* alphabet for database                       */

  uint32_t            seq_count;   /* number of sequences                         */
  uint64_t            res_count;   /* number of residues                          */
  uint32_t            max_seq;     /* longest sequence                            */

  ESL_SQ             *sq_list;     /* list of cached sequences [0 .. seq_count-1] */

  void               *residue_mem; /* memory holding the residues                 */
  void               *header_mem;  /* memory holding the header strings           */

  uint64_t            res_size;    /* size of residue memory allocation           */
  uint64_t            hdr_size;    /* size of header memory allocation            */
} ESL_SQCACHE;
#endif

/*::cexcerpt::sq_sqio_format::begin::*/
/* Unaligned file format codes
 * These codes are coordinated with the msa module.
 *   - 0 is an unknown/unassigned format (eslSQFILE_UNKNOWN, eslMSAFILE_UNKNOWN)
 *   - <=100 is reserved for sqio, for unaligned formats
 *   - >100  is reserved for msa, for aligned formats
 */
#define eslSQFILE_UNKNOWN      0
#define eslSQFILE_FASTA        1   // FASTA format
#define eslSQFILE_EMBL         2   // EMBL DNA sequence
#define eslSQFILE_GENBANK      3   // Genbank DNA sequence
#define eslSQFILE_DDBJ         4   // DDBJ (currently identical to GenBank parser)
#define eslSQFILE_UNIPROT      5   // UniProt (currently identical to EMBL parser)
#define eslSQFILE_NCBI         6   // NCBI blast db, v4, single file
#define eslSQFILE_DAEMON       7   // Farrar format, hmmpgmd queries: fasta + // terminator
#define eslSQFILE_HMMPGMD      8   // Farrar hmmpgmd database format: fasta + # header
#define eslSQFILE_FMINDEX      9   // Pressed FM-index format used in HMMER
/*::cexcerpt::sq_sqio_format::end::*/

/* eslREADBUFSIZE is the fixed size of a block to bring in at one time,
 * in character-based (fread()) parsers (like the FASTA parser).
 */
#define eslREADBUFSIZE  4096

int  esl_sqfile_Open(const char *seqfile, int fmt, const char *env, ESL_SQFILE **ret_sqfp);
void esl_sqfile_Close(ESL_SQFILE *sqfp);

#ifdef eslAUGMENT_ALPHABET
int  esl_sqfile_OpenDigital(const ESL_ALPHABET *abc, const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp);
int  esl_sqfile_SetDigital(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc);
int  esl_sqfile_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type);
#endif

int   esl_sqio_Read        (ESL_SQFILE *sqfp, ESL_SQ *sq);
int   esl_sqio_ReadInfo    (ESL_SQFILE *sqfp, ESL_SQ *sq);
int   esl_sqio_ReadWindow  (ESL_SQFILE *sqfp, int C, int W, ESL_SQ *sq);
int   esl_sqio_ReadSequence(ESL_SQFILE *sqfp, ESL_SQ *sq);
int   esl_sqio_ReadBlock   (ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock, int max_residues, int max_sequences, int long_target);
int   esl_sqio_Parse       (char *buffer, int size, ESL_SQ *s, int format);

int   esl_sqio_Write       (FILE *fp, ESL_SQ *s, int format, int update);
int   esl_sqio_Echo        (ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp);

const char  *esl_sqfile_GetErrorBuf(const ESL_SQFILE *sqfp);
int   esl_sqfile_IsRewindable(const ESL_SQFILE *sqfp);
int   esl_sqio_IsAlignment(int fmt);
int   esl_sqio_EncodeFormat(char *fmtstring);
char *esl_sqio_DecodeFormat(int fmt);
int   esl_sqfile_Position(ESL_SQFILE *sqfp, off_t offset);
int   esl_sqio_Ignore(ESL_SQFILE *sqfp, const char *ignoredchars);
int   esl_sqio_AcceptAs(ESL_SQFILE *sqfp, char *xchars, char readas);

#ifdef eslAUGMENT_SSI
int   esl_sqfile_OpenSSI         (ESL_SQFILE *sqfp, const char *ssifile_hint);
int   esl_sqfile_PositionByKey   (ESL_SQFILE *sqfp, const char *key);
int   esl_sqfile_PositionByNumber(ESL_SQFILE *sqfp, int which);

int   esl_sqio_Fetch      (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq);
int   esl_sqio_FetchInfo  (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq);
int   esl_sqio_FetchSubseq(ESL_SQFILE *sqfp, const char *source, int64_t start, int64_t end, ESL_SQ *sq);
#endif

int   esl_sqfile_Cache(const ESL_ALPHABET *abc, const char *seqfile, int fmt, const char *env, ESL_SQCACHE **ret_sqcache);
void  esl_sqfile_Free(ESL_SQCACHE *sqcache);

#endif /*eslSQIO_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_sqio.h ***/

#ifdef HAVE_PTHREAD

#include <stdio.h>
#include <stdint.h>
#include <pthread.h>

/* ESL_DSQDATA_CHUNK
 * A data chunk returned by esl_dsqdata_Read().
 */
typedef struct esl_dsqdata_chunk_s {
  int64_t   i0;           // Chunk contains sequences i0..i0+N-1 from the database, 0-offset
  int       N;            // Chunk contains N sequences

  ESL_DSQ **dsq;          // Pointers to each of the N sequences
  char    **name;         // Names, \0 terminated.  Ptr into <metadata> buffer.
  char    **acc;          // Optional accessions, \0 terminated;   "\0" if none.
  char    **desc;         // Optional descriptions, \0 terminated; "\0" if none
  int32_t  *taxid;        // NCBI taxonomy identifiers. (>=1 is a taxid; -1 means none)
  int64_t  *L;            // Sequence lengths, in residues. The unpacker figures these out.

  /* Memory management */
  char     *smem;         // Unpacked (dsq[]) and packed (psq) data ptrs share this allocation.
  uint32_t *psq;          // Pointer into smem; packed data fread()'s go here.
  int       pn;           // how many uint32's are loaded in <psq>
  char     *metadata;     // Raw fread() buffer of all name/acc/desc/taxid data.
  int       mdalloc;      // Current allocation size for <metadata> in bytes
  struct esl_dsqdata_chunk_s *nxt; // Chunks can be put in linked lists
} ESL_DSQDATA_CHUNK;

/* ESL_DSQDATA_RECORD
 * The dsqi index file is composed of an array of these, aside from its header.
 */
typedef struct esl_dsqdata_record_s {
  int64_t  metadata_end;
  int64_t  psq_end;
} ESL_DSQDATA_RECORD;

/* ESL_DSQDATA
 * The object created by esl_dsqdata_Open() and used by esl_dsqdata_Read()
 * to read chunks of sequence data from the database.
 */
typedef struct esl_dsqdata_s {
  char         *basename;    // Basename of the four dsqdata data files
  FILE         *stubfp;      // Open <basename> stub file
  FILE         *ifp;         // Open basename.dsqi index file
  FILE         *sfp;         // Open basename.dsqs sequence file
  FILE         *mfp;         // Open basename.dsqm metadata file
  ESL_ALPHABET *abc_r;       // Copy of ptr to the alphabet the caller told us to read in.

  /* Header information from dsqi index file
   *  .. dsqm, dsqs have magic and uniquetag for integrity checking
   *  .. and stub file has uniquetag as text.
   */
  uint32_t     magic;       // Binary magic format code, for detecting byteswapping
  uint32_t     uniquetag;   // Random number tag that links the four files
  uint32_t     flags;       // Currently unused (0); reserved for future bitflags
  uint32_t     max_namelen; // Max name length in the dataset
  uint32_t     max_acclen;  //  .. and max accession length
  uint32_t     max_desclen; //  .. and max description length
  uint64_t     max_seqlen;  //  .. and max seq length. 64b = bring on Paris japonica.
  uint64_t     nseq;        // Total number of sequences in the dataset
  uint64_t     nres;        //  .. and total number of residues

  /* Control parameters. */
  int          chunk_maxseq;    // default = eslDSQDATA_CHUNK_MAXSEQ
  int          chunk_maxpacket; // default = eslDSQDATA_CHUNK_MAXPACKET
  int          do_byteswap;     // TRUE if we need to byteswap (bigendian <=> littleendian)
  int          pack5;           // TRUE if we're using all 5bit packing; FALSE for mixed 2+5bit

  /* Managing the reader's threaded producer/consumer pipeline:
   * consisting of 1 loader thread and 1 unpacker thread that we manage,
   * and <nconsumers> consumer threads that caller created to get
   * successive chunks with esl_dsqdata_Read().
   */
  int                nconsumers;               // Caller's reading with this # of readers

  pthread_t          loader_t;                 // Loader thread id
  ESL_DSQDATA_CHUNK *loader_outbox;            // A loaded chunk goes here, for unpacker
  pthread_mutex_t    loader_outbox_mutex;      // mutex protecting the outbox
  pthread_cond_t     loader_outbox_full_cv;    // signal to unpacker that next chunk is ready
  pthread_cond_t     loader_outbox_empty_cv;   // signal from unpacker that it's got the chunk

  pthread_t          unpacker_t;               // Unpacker thread id
  ESL_DSQDATA_CHUNK *unpacker_outbox;          // Unpacked chunk goes here, for _Read()
  pthread_mutex_t    unpacker_outbox_mutex;    // mutex protecting the outbox
  pthread_cond_t     unpacker_outbox_full_cv;  // signal to _Read() that chunk is ready (or at_eof)
  pthread_cond_t     unpacker_outbox_empty_cv; // signal from _Read() that it's got the chunk
  int                at_eof;                   // flag that goes up at end of the input file;
											   //  .. <at_eof> change is in unpacker's mutex
  ESL_DSQDATA_CHUNK *recycling;                // Linked list of chunk memory for reuse
  pthread_mutex_t    recycling_mutex;          // mutex protecting the recycling list
  pthread_cond_t     recycling_cv;             // signal to loader that a chunk is available

  /* Error handling.
   * Pthread variables don't define a value for "unset", so for pristine cleanup after
   * errors, we must use separate booleans to track which thread resources are  created.
   */
  int  lt_c;  int lom_c;  int lof_c;  int loe_c;
  int  ut_c;  int uom_c;  int uof_c;  int uoe_c;
  int  rm_c;  int r_c;
  char errbuf[eslERRBUFSIZE];   // User-directed error message in case of a failed open or read.
} ESL_DSQDATA;

/* Defaults for size of eslDSQDATA_CHUNK
 */
#define eslDSQDATA_CHUNK_MAXSEQ       4096      // max number of sequences
#define eslDSQDATA_CHUNK_MAXPACKET  262144      // max number of uint32 sequence packets

/* Reading the control bits on a packet v
 */
#define eslDSQDATA_EOD   (1 << 31)
#define eslDSQDATA_5BIT  (1 << 30)
#define ESL_DSQDATA_EOD(v)   ((v) & eslDSQDATA_EOD)
#define ESL_DSQDATA_5BIT(v)  ((v) & eslDSQDATA_5BIT)

/* Functions in the API
 */
int  esl_dsqdata_Open   (ESL_ALPHABET **byp_abc, char *basename, int nconsumers, ESL_DSQDATA **ret_dd);
int  esl_dsqdata_Read   (ESL_DSQDATA *dd, ESL_DSQDATA_CHUNK **ret_chu);
int  esl_dsqdata_Recycle(ESL_DSQDATA *dd, ESL_DSQDATA_CHUNK *chu);
int  esl_dsqdata_Close  (ESL_DSQDATA *dd);

int  esl_dsqdata_Write  (ESL_SQFILE *sqfp, char *basename, char *errbuf);

#endif /*eslDSQDATA_INCLUDED*/

#endif // HAVE_PTHREAD
/*** End of inlined file: esl_dsqdata.h ***/

/* Structure: ESL_HISTOGRAM
 *
 * Keeps a score histogram, in which scores are counted into bins of
 * size (width) w.
 *   histogram starts at bmin <  floor(xmin/w) * w
 *   histogram ends at   bmax >= ceil(xmax/w)*w
 *   nb = (bmax-bmin)/w
 *   each score x is counted into bin b = nb - (int) (bmax-x)/w
 *   each bin b contains scores bw+bmin < x <= (b+1)w + bmin
 *
 * Anything having to do with the counts themselves (obs, n, etc)
 * is a uint64_t, with range 0..2^64-1  (up to 2e19).
 */
typedef struct {
  /* The histogram is kept as counts in fixed-width bins.
   */
  uint64_t *obs;	/* observed counts in bin b, 0..nb-1 (dynamic)      */
  int       nb;         /* number of bins                                   */
  double    w;		/* fixed width of each bin                          */
  double    bmin, bmax;	/* histogram bounds: all x satisfy bmin < x <= bmax */
  int       imin, imax;	/* smallest, largest bin that contain obs[i] > 0    */

  /* Optionally, in a "full" h, we can also keep all the raw samples in x.
   */
  double    xmin, xmax;	/* smallest, largest sample value x observed        */
  uint64_t  n;          /* total number of raw data samples                 */
  double   *x;		/* optional: raw sample values x[0..n-1]            */
  uint64_t  nalloc;	/* current allocated size of x                      */

  /* The binned data might be censored (either truly, or virtually).
   * This information has to be made available to a binned/censored
   * parameter fitting function, and to goodness-of-fit tests.
   */
  double   phi;		/* censoring value; all x_i > phi                   */
  int      cmin;	/* smallest bin index that contains uncensored data */
  uint64_t z;		/* # of censored values <= phi                      */
  uint64_t Nc;	        /* # samples in complete data (including unobs)     */
  uint64_t No;		/* # of samples in observed data                    */

  /* Expected binned counts are set by SetExpect() or SetExpectedTail().
   */
  double *expect;	/* expected counts in bin b, 0..nb-1 (not resized)  */
  int     emin;		/* smallest bin index that contains expected counts */
  double  tailbase;	/* for tail fits: fitted x > tailbase               */
  double  tailmass;	/* for tail fits: fractional prob in the tail       */

  /* Some status flags
   */
  int is_full;		/* TRUE when we're keeping raw data in x           */
  int is_done;		/* TRUE if we prevent more Add()'s                 */
  int is_sorted;	/* TRUE if x is sorted smallest-to-largest         */
  int is_tailfit;	/* TRUE if expected dist only describes tail       */
  int is_rounded;	/* TRUE if values aren't more accurate than bins   */
  enum { COMPLETE, VIRTUAL_CENSORED, TRUE_CENSORED } dataset_is;

} ESL_HISTOGRAM;

/*** Start of inlined file: esl_exponential.h ***/
#ifndef eslEXPONENTIAL_INCLUDED
#define eslEXPONENTIAL_INCLUDED

double esl_exp_pdf    (double x, double mu, double lambda);
double esl_exp_logpdf (double x, double mu, double lambda);
double esl_exp_cdf    (double x, double mu, double lambda);
double esl_exp_logcdf (double x, double mu, double lambda);
double esl_exp_surv   (double x, double mu, double lambda);
double esl_exp_logsurv(double x, double mu, double lambda);
double esl_exp_invcdf (double p, double mu, double lambda);
double esl_exp_invsurv(double p, double mu, double lambda);

double esl_exp_generic_pdf   (double x, void *params);
double esl_exp_generic_cdf   (double x, void *params);
double esl_exp_generic_surv  (double x, void *params);
double esl_exp_generic_invcdf(double p, void *params);

int    esl_exp_Plot(FILE *fp, double mu, double lambda,
			   double (*func)(double x, double mu, double lambda),
			   double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
double esl_exp_Sample(ESL_RANDOMNESS *r, double mu, double lambda);
#endif

int esl_exp_FitComplete     (double *x, int n, double *ret_mu, double *ret_lambda);
int esl_exp_FitCompleteScale(double *x, int n, double      mu, double *ret_lambda);

#ifdef eslAUGMENT_HISTOGRAM
int esl_exp_FitCompleteBinned(ESL_HISTOGRAM *h,
				     double *ret_mu, double *ret_lambda);
#endif

#endif /*eslEXPONENTIAL_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_exponential.h ***/


/*** Start of inlined file: esl_gamma.h ***/
#ifndef eslGAMMA_INCLUDED
#define eslGAMMA_INCLUDED

#ifdef eslAUGMENT_HISTOGRAM

/*** Start of inlined file: esl_histogram.h ***/
#ifndef eslHISTOGRAM_INCLUDED
#define eslHISTOGRAM_INCLUDED

#include <math.h>   /* floor() is in one of the macros */



#define esl_histogram_Bin2LBound(h,b)  ((h)->w*(b) + (h)->bmin)
#define esl_histogram_Bin2UBound(h,b)  ((h)->w*((b)+1) + (h)->bmin)

/* Creating/destroying histograms and collecting data:
 */
ESL_HISTOGRAM *esl_histogram_Create    (double bmin, double bmax, double w);
ESL_HISTOGRAM *esl_histogram_CreateFull(double bmin, double bmax, double w);
void           esl_histogram_Destroy  (ESL_HISTOGRAM *h);
int            esl_histogram_Score2Bin(ESL_HISTOGRAM *h, double x, int *ret_b);
int            esl_histogram_Add      (ESL_HISTOGRAM *h, double x);

/* Declarations about the binned data before parameter fitting:
 */
int esl_histogram_DeclareCensoring(ESL_HISTOGRAM *h, int z, double phi);
int esl_histogram_DeclareRounding (ESL_HISTOGRAM *h);
int esl_histogram_SetTail         (ESL_HISTOGRAM *h, double phi,
					  double *ret_newmass);
int esl_histogram_SetTailByMass   (ESL_HISTOGRAM *h, double pmass,
					  double *ret_newmass);

/* Accessing data samples in a full histogram:
 */
int esl_histogram_GetRank(ESL_HISTOGRAM *h, int rank, double *ret_x);
int esl_histogram_GetData(ESL_HISTOGRAM *h, double **ret_x, int *ret_n);
int esl_histogram_GetTail(ESL_HISTOGRAM *h, double phi, double **ret_x,
				 int *ret_n, int *ret_z);
int esl_histogram_GetTailByMass(ESL_HISTOGRAM *h, double pmass,
				       double **ret_x, int *ret_n, int *ret_z);

/* Setting expected binned counts:
 */
int esl_histogram_SetExpect(ESL_HISTOGRAM *h,
				   double (*cdf)(double x, void *params),
				   void *params);
int esl_histogram_SetExpectedTail(ESL_HISTOGRAM *h, double base_val,
					 double pmass,
					 double (*cdf)(double x, void *params),
					 void *params);

/* Output/display of binned data:
 */
int esl_histogram_Write       (FILE *fp, ESL_HISTOGRAM *h);
int esl_histogram_Plot        (FILE *fp, ESL_HISTOGRAM *h);
int esl_histogram_PlotSurvival(FILE *fp, ESL_HISTOGRAM *h);
int esl_histogram_PlotQQ      (FILE *fp, ESL_HISTOGRAM *h,
				      double (*invcdf)(double, void *), void *params);

/* Goodness of fit testing
 */
int esl_histogram_Goodness(ESL_HISTOGRAM *h, int nfitted,
				  int *ret_nbins,
				  double *ret_G,  double *ret_Gp,
				  double *ret_X2, double *ret_X2p);

#endif /*eslHISTOGRAM_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_histogram.h ***/


#endif

double esl_gam_pdf    (double x, double mu, double lambda, double tau);
double esl_gam_logpdf (double x, double mu, double lambda, double tau);
double esl_gam_cdf    (double x, double mu, double lambda, double tau);
double esl_gam_logcdf (double x, double mu, double lambda, double tau);
double esl_gam_surv   (double x, double mu, double lambda, double tau);
double esl_gam_logsurv(double x, double mu, double lambda, double tau);
double esl_gam_invcdf (double p, double mu, double lambda, double tau);

double esl_gam_generic_pdf   (double x, void *params);
double esl_gam_generic_cdf   (double x, void *params);
double esl_gam_generic_surv  (double x, void *params);
double esl_gam_generic_invcdf(double x, void *params);

int esl_gam_Plot(FILE *fp, double mu, double lambda, double tau,
			double (*func)(double x, double mu, double lambda, double tau),
			double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
double esl_gam_Sample(ESL_RANDOMNESS *r, double mu, double lambda, double tau);
#endif

int esl_gam_FitComplete(double *x, int n, double mu, double *ret_lambda, double *ret_tau);

#ifdef eslAUGMENT_HISTOGRAM
int esl_gam_FitCompleteBinned(ESL_HISTOGRAM *h, double *ret_mu, double *ret_lambda, double *ret_tau);
#endif

#endif /*eslGAMMA_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_gamma.h ***/


/*** Start of inlined file: esl_gencode.h ***/
#ifndef eslGENCODE_INCLUDED
#define eslGENCODE_INCLUDED


/*** Start of inlined file: esl_getopts.h ***/
#ifndef eslGETOPTS_INCLUDED
#define eslGETOPTS_INCLUDED

/* Object: ESL_OPTIONS
 *
 * The application main.c defines an array of <ESL_OPTIONS> structures to
 * define what configuration options are used. The array is
 * terminated by a structure containing { NULL, NULL, NULL, 0, NULL,
 * NULL, NULL, NULL} (or more simply, just 0 in all 8 fields.)
 */
/*::cexcerpt::options_object::begin::*/
typedef struct {
  char *name;           /* either short "-a" or long "--foo" style               */
  int   type;           /* arg type, for type checking: (eslARG_INT, etc.)       */
  char *defval;         /* default setting, or NULL ("default" is a C keyword)   */
  char *envvar;         /* associated environ var ("BLASTDB"), or NULL           */
  char *range;          /* for range checking arg: ("0<=x<=1", etc.)             */
  char *toggle_opts;    /* comma-sep'd optlist: turn these off if this opt is on */
  char *required_opts;  /* comma-sep'd optlist: these must also be set           */
  char *incompat_opts;  /* comma-sep'd optlist: these must not be set            */
  char *help;           /* help/usage string                                     */
  int   docgrouptag;    /* integer tag for documentation groups                  */
} ESL_OPTIONS;
/*::cexcerpt::options_object::end::*/

/* Argument types: the "type" variable in <ESL_OPTIONS>.
 */
#define eslARG_NONE      0	/* option takes no argument (so, is boolean)   */
#define eslARG_INT       1	/* arg convertable by atoi()               <n> */
#define eslARG_REAL      2	/* arg convertable by atof()               <x> */
#define eslARG_CHAR      3	/* arg is a single character               <c> */
#define eslARG_STRING    4	/* unchecked arg type                      <s> */
#define eslARG_INFILE    5      /* input file - same as string, shown as   <f> */
#define eslARG_OUTFILE   6      /* output file - same as string, shown as  <f> */

/* Object: ESL_GETOPTS
 *
 * An <ESL_GETOPTS> object is created to parse configuration
 * from command line options, config file(s), and environment
 * variables.
 */
typedef struct {
  ESL_OPTIONS *opt;       /* array of app-defined options              */
  int          nopts;     /* number of options                         */

  int    argc;		  /* argc from command line                    */
  char **argv;		  /* argv from command line                    */
  int    optind;	  /* position in argc; eventually 1st arg idx  */
  int    nfiles;	  /* # of cfgfiles that have been processed    */

  char **val;		  /* config'ed val for each option (as string) */
  int   *setby;		  /* array [0..nopts-1] for who set option i   */
  int   *valloc;          /* 0, or length of alloc for val[i]          */

  char  *optstring;	  /* internal: ptr into string of 1-char opts in argv[]          */
  char  *spoof;	    	  /* internal allocation: ProcessSpoof() stores cmdline          */
  char **spoof_argv;	  /* internal allocation: ProcessSpoof()'s ptrs into its cmdline */

  char  errbuf[eslERRBUFSIZE];	/* buffer for reporting user error     */
} ESL_GETOPTS;

/* Possible values of the <setby> variable in ESL_GETOPTS.
 * Additionally, values of >3 also indicate a config file, in order
 * of _ProcessConfigFile() calls (that is, setby=3 is the first
 * config file, setby=4 is the second, etc.).
 */
#define eslARG_SETBY_DEFAULT  0
#define eslARG_SETBY_CMDLINE  1
#define eslARG_SETBY_ENV      2
#define eslARG_SETBY_CFGFILE  3

/* The visible API.
 */
ESL_GETOPTS *esl_getopts_Create(ESL_OPTIONS *opt);
ESL_GETOPTS *esl_getopts_CreateDefaultApp(ESL_OPTIONS *options, int nargs, int argc, char **argv, char *banner, char *usage);
int          esl_getopts_Reuse  (ESL_GETOPTS *g);
void         esl_getopts_Destroy(ESL_GETOPTS *g);
void         esl_getopts_Dump(FILE *ofp, ESL_GETOPTS *g);

int esl_opt_ProcessConfigfile (ESL_GETOPTS *g, char *filename, FILE *fp);
int esl_opt_ProcessEnvironment(ESL_GETOPTS *g);
int esl_opt_ProcessCmdline    (ESL_GETOPTS *g, int argc, char **argv);
int esl_opt_ProcessSpoof      (ESL_GETOPTS *g, const char *cmdline);
int esl_opt_VerifyConfig      (ESL_GETOPTS *g);
int esl_opt_ArgNumber   (const ESL_GETOPTS *g);
int esl_opt_SpoofCmdline(const ESL_GETOPTS *g, char **ret_cmdline);

int esl_opt_GetSetter(const ESL_GETOPTS *g, char *optname);

int    esl_opt_IsDefault (const ESL_GETOPTS *g, char *optname);
int    esl_opt_IsOn      (const ESL_GETOPTS *g, char *optname);
int    esl_opt_IsUsed    (const ESL_GETOPTS *g, char *optname);

int    esl_opt_GetBoolean(const ESL_GETOPTS *g, char *optname);
int    esl_opt_GetInteger(const ESL_GETOPTS *g, char *optname);
double esl_opt_GetReal   (const ESL_GETOPTS *g, char *optname);
char   esl_opt_GetChar   (const ESL_GETOPTS *g, char *optname);
char  *esl_opt_GetString (const ESL_GETOPTS *g, char *optname);
char  *esl_opt_GetArg    (const ESL_GETOPTS *g, int which);

int esl_opt_DisplayHelp(FILE *ofp, ESL_GETOPTS *go, int docgroup, int indent, int textwidth);

#endif /*eslGETOPTS_INCLUDED*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/*** End of inlined file: esl_getopts.h ***/

typedef struct {
  int     transl_table;      // NCBI transl_table number, or -1. Only set for a standard NCBI table, with _Set(); _Read() from file doesn't set this.
  char    desc[128];         // Description, or "".                ... ditto

  ESL_DSQ basic[64];         // Basic code table. aacode[0..63; pos1^16 + pos2^4 + pos3] = residue code for amino acid, 0..19 or the Nonresidue code. No degeneracies.
  int8_t  is_initiator[64];  // TRUE for allowed initiator codons; FALSE if not

  const ESL_ALPHABET *nt_abc;  // A reference to nucleic alphabet that caller is maintaining elsewhere
  const ESL_ALPHABET *aa_abc;  // A reference to amino alphabet that caller is maintaining
} ESL_GENCODE;

/* struct esl_gencode_workstate_s
 *   keeps state in DNA sequence <sq>, allowing us to process a sequence
 *   either in a single gulp (using ReadSeq) or in overlapping windows
 *   (using ReadWindow).
 *
 *   also contains one-time configuration information for translation
 */
typedef struct esl_gencode_workstate_s {
  /* stateful info (which may get updated with each new seq, strand, and/or window): */
  ESL_SQ *psq[3];     // Growing ORFs in each frame
  int8_t  in_orf[3];  // TRUE|FALSE: TRUE if we're growing an ORF in this frame
  int     apos;       // 1..L:  current nucleotide we're on (starting a codon) in <sq>
  int     frame;      // 0..2:  which frame <apos> is in
  int     codon;      // 0..63: Digitized codon for apos,apos+1,apos+2
  int     inval;      // 0..3:  how many apos increments we need to get past an ambiguous nucleotide
  int     is_revcomp; // TRUE|FALSE: TRUE if we're doing reverse complement strand
  int     orfcount;   // >=0:   How many ORFs we've processed so far

  ESL_SQ_BLOCK  *orf_block; // block of sequences to which to write ORFs

  /* one-time configuration information (from options) */
  int     do_watson;         // TRUE|FALSE:  TRUE if we translate the top strand
  int     do_crick;          // TRUE|FALSE:  TRUE if we translate the reverse complement strand
  int     using_initiators;  // TRUE|FALSE : TRUE if -m or -M, only valid initiators can start an ORF, and initiator codon always translates to Met
  int     minlen;            // >=0: minimum orf length that process_orf will deal with
  FILE   *outfp;             // default stdout: where to write output ORF data
  int     outformat;         // default eslSQFILE_FASTA: sqfile format to write ORFs in
} ESL_GENCODE_WORKSTATE;

/* Create/Destroy workstate */
void esl_gencode_WorkstateDestroy(ESL_GENCODE_WORKSTATE *wrk);
ESL_GENCODE_WORKSTATE * esl_gencode_WorkstateCreate(ESL_GETOPTS *go, ESL_GENCODE *gcode);

/* the ESL_GENCODE genetic code object */
ESL_GENCODE *esl_gencode_Create(const ESL_ALPHABET *nt_abc, const ESL_ALPHABET *aa_abc);
void         esl_gencode_Destroy            (ESL_GENCODE *gcode);
int          esl_gencode_Set                (ESL_GENCODE *gcode,  int ncbi_transl_table);
int          esl_gencode_SetInitiatorAny    (ESL_GENCODE *gcode);
int          esl_gencode_SetInitiatorOnlyAUG(ESL_GENCODE *gcode);

/* reading and writing genetic codes in NCBI format */
int          esl_gencode_Read(ESL_FILEPARSER *efp, const ESL_ALPHABET *nucleic_abc, const ESL_ALPHABET *amino_abc, ESL_GENCODE **ret_gcode);
int          esl_gencode_Write(FILE *ofp, const ESL_GENCODE *gcode, int add_comment);

/* DNA->protein digital translation, allowing ambiguity chars */
int   esl_gencode_GetTranslation(const ESL_GENCODE *gcode, ESL_DSQ *dsqp);
int   esl_gencode_IsInitiator   (const ESL_GENCODE *gcode, ESL_DSQ *dsqp);

/* Debugging/development utilities */
char *esl_gencode_DecodeDigicodon(const ESL_GENCODE *gcode, int digicodon, char *codon);
int   esl_gencode_DumpAltCodeTable(FILE *ofp);
int   esl_gencode_Compare(const ESL_GENCODE *gc1, const ESL_GENCODE *gc2, int metadata_too);

/* Functions for processing ORFs  */
int esl_gencode_ProcessOrf(ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq);
void esl_gencode_ProcessStart(ESL_GENCODE *gcode, ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq);
int esl_gencode_ProcessPiece(ESL_GENCODE *gcode, ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq);
int esl_gencode_ProcessEnd(ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq);

#endif	/*eslGENCODE_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_gencode.h ***/


/*** Start of inlined file: esl_gev.h ***/
#ifndef eslGEV_INCLUDED
#define eslGEV_INCLUDED

double esl_gev_pdf    (double x, double mu, double lambda, double alpha);
double esl_gev_logpdf (double x, double mu, double lambda, double alpha);
double esl_gev_cdf    (double x, double mu, double lambda, double alpha);
double esl_gev_logcdf (double x, double mu, double lambda, double alpha);
double esl_gev_surv   (double x, double mu, double lambda, double alpha);
double esl_gev_logsurv(double x, double mu, double lambda, double alpha);
double esl_gev_invcdf (double p, double mu, double lambda, double alpha);

double esl_gev_generic_pdf   (double x, void *params);
double esl_gev_generic_cdf   (double x, void *params);
double esl_gev_generic_surv  (double x, void *params);
double esl_gev_generic_invcdf(double p, void *params);

int    esl_gev_Plot(FILE *fp, double mu, double lambda, double alpha,
			   double (*func)(double x, double mu, double lambda, double alpha),
			   double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
double esl_gev_Sample(ESL_RANDOMNESS *r, double mu, double lambda, double alpha);
#endif

#ifdef eslAUGMENT_MINIMIZER
int esl_gev_FitComplete(double *x, int n,
			       double *ret_mu, double *ret_lambda,
			       double *ret_alpha);
int esl_gev_FitCensored(double *x, int n, int z, double phi,
			       double *ret_mu, double *ret_lambda,
			       double *ret_alpha);
#endif /*eslAUGMENT_MINIMIZER*/

#endif /*eslGEV_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_gev.h ***/


/*** Start of inlined file: esl_gumbel.h ***/
#ifndef eslGUMBEL_INCLUDED
#define eslGUMBEL_INCLUDED

double  esl_gumbel_pdf    (double x, double mu, double lambda);
double  esl_gumbel_logpdf (double x, double mu, double lambda);
double  esl_gumbel_cdf    (double x, double mu, double lambda);
double  esl_gumbel_logcdf (double x, double mu, double lambda);
double  esl_gumbel_surv   (double x, double mu, double lambda);
double  esl_gumbel_logsurv(double x, double mu, double lambda);
double  esl_gumbel_invcdf (double p, double mu, double lambda);
double  esl_gumbel_invsurv(double p, double mu, double lambda);

double  esl_gumbel_generic_pdf   (double x, void *params);
double  esl_gumbel_generic_cdf   (double x, void *params);
double  esl_gumbel_generic_surv  (double x, void *params);
double  esl_gumbel_generic_invcdf(double p, void *params);

int esl_gumbel_Plot(FILE *fp, double mu, double lambda,
			   double (*func)(double x, double mu, double lambda),
			   double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
double esl_gumbel_Sample(ESL_RANDOMNESS *r, double mu, double lambda);
#endif

int esl_gumbel_FitComplete   (double *x, int n,                    double *ret_mu, double *ret_lambda);
int esl_gumbel_FitCompleteLoc(double *x, int n,                    double lambda,  double *ret_mu);
int esl_gumbel_FitCensored   (double *x, int n, int z, double phi, double *ret_mu, double *ret_lambda);
int esl_gumbel_FitCensoredLoc(double *x, int n, int z, double phi, double lambda,  double *ret_mu);
#ifdef eslAUGMENT_MINIMIZER
int esl_gumbel_FitTruncated  (double *x, int n,        double phi, double *ret_mu, double *ret_lambda);
#endif

#endif /*eslGUMBEL_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_gumbel.h ***/


/*** Start of inlined file: esl_heap.h ***/
#ifndef eslHEAP_INCLUDED
#define eslHEAP_INCLUDED

#define eslHEAP_INITALLOC 128

#define eslHEAP_MIN   0
#define eslHEAP_MAX   1

#define ESL_HEAP_PARENT(i)  ( ((i)-1) / 2  )
#define ESL_HEAP_LEFT(i)    ( ((i)*2) + 1 )
#define ESL_HEAP_RIGHT(i)   ( ((i)+1) * 2  )

typedef struct esl_heap_s {
  int *idata;

  int  n;
  int  nalloc;
  int  maxormin;		/* eslHEAP_MAX | eslHEAP_MIN */
} ESL_HEAP;

ESL_HEAP *esl_heap_ICreate   (int maxormin);
int       esl_heap_GetCount  (ESL_HEAP *hp);
int       esl_heap_IGetTopVal(ESL_HEAP *hp);
int       esl_heap_Reuse     (ESL_HEAP *hp);
void      esl_heap_Destroy   (ESL_HEAP *hp);

int       esl_heap_IInsert(ESL_HEAP *hp, int val);

int       esl_heap_IExtractTop(ESL_HEAP *hp, int *ret_val);

int       esl_heap_IGetTop(ESL_HEAP *hp);

#endif /*eslHEAP_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/*** End of inlined file: esl_heap.h ***/


/*** Start of inlined file: esl_hmm.h ***/
#ifndef eslHMM_INCLUDED
#define eslHMM_INCLUDED

typedef struct {
  int     M;                    /* number of states in the model          */
  int     K;                    /* size of alphabet (redundant w/ abc->K) */
  float  *pi;                   /* initial (begin) distribution (0..M)    */
  float **t;                    /* Mx(M+1) state transition probabilities */
  float **e;                    /* MxK emission probabilities             */

  float **eo;			/* K'xM emission odds ratios              */
  const ESL_ALPHABET *abc;      /* ptr to alphabet                        */
} ESL_HMM;

typedef struct {
  float **dp;			/* [0..L][0..M-1] DP matrix                              */
  float  *sc;			/* [0..L+1] scale factors (log probs)                    */
  int     M;			/* actual model dimension (0..M-1)                       */
  int     L;			/* actual sequence dimension (1..L)                      */

  float    *dp_mem;		/* memory allocated for the resizable DP matrix          */
  int       allocR;		/* current allocated # of rows: L+1 <= validR <= allocR  */
  int       validR; 		/* # of dp rows actually pointing at DP memory           */
  int       allocM;		/* current set row width; M <= allocM                    */
  uint64_t  ncells;		/* total allocation of dp_mem; ncells >= (validR)(allocM)*/
} ESL_HMX;

ESL_HMM *esl_hmm_Create(const ESL_ALPHABET *abc, int M);
ESL_HMM *esl_hmm_Clone(const ESL_HMM *hmm);
int      esl_hmm_Configure(ESL_HMM *hmm, float *fq);
int      esl_hmm_SetDegeneracies(ESL_HMM *hmm);
void     esl_hmm_Destroy(ESL_HMM *hmm);

ESL_HMX *esl_hmx_Create(int allocL, int allocM);
int      esl_hmx_GrowTo (ESL_HMX *mx, int L, int M);
void     esl_hmx_Destroy(ESL_HMX *mx);

int      esl_hmm_Emit(ESL_RANDOMNESS *r, const ESL_HMM *hmm, ESL_DSQ **opt_dsq, int **opt_path, int *opt_L);
int      esl_hmm_Forward(const ESL_DSQ *dsq, int L, const ESL_HMM *hmm, ESL_HMX *fwd, float *opt_sc);
int      esl_hmm_Backward(const ESL_DSQ *dsq, int L, const ESL_HMM *hmm, ESL_HMX *bck, float *opt_sc);

#endif /*eslHMM_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_hmm.h ***/


/*** Start of inlined file: esl_hyperexp.h ***/

typedef struct {
  double *q;			/* mixture coefficients   [0..K-1]*/
  double *lambda;		/* scale params           [0..K-1]*/
  double *wrk;			/* tmp K-vector, for logpdf calc  */
  double  mu;			/* location (x offset) parameter  */
  int     K;			/* # of components                */
  char   *fixlambda;		/* TRUE to constrain a lambda val */
  int     fixmix;		/* TRUE to constrain the q's      */
} ESL_HYPEREXP;

ESL_HYPEREXP *esl_hyperexp_Create(int K);
void          esl_hyperexp_Destroy(ESL_HYPEREXP *h);
int           esl_hyperexp_Copy(ESL_HYPEREXP *src, ESL_HYPEREXP *dest);
int           esl_hyperexp_FixedUniformMixture(ESL_HYPEREXP *h);
int           esl_hyperexp_SortComponents(ESL_HYPEREXP *h);
int           esl_hyperexp_Write(FILE *fp, ESL_HYPEREXP *hxp);
int           esl_hyperexp_Dump(FILE *fp, ESL_HYPEREXP *hxp);
#ifdef eslAUGMENT_FILEPARSER
int           esl_hyperexp_Read(ESL_FILEPARSER *ef, ESL_HYPEREXP **ret_hxp);
int           esl_hyperexp_ReadFile(char *filename, ESL_HYPEREXP **ret_hxp);
#endif

double  esl_hxp_pdf    (double x, ESL_HYPEREXP *h);
double  esl_hxp_logpdf (double x, ESL_HYPEREXP *h);
double  esl_hxp_cdf    (double x, ESL_HYPEREXP *h);
double  esl_hxp_logcdf (double x, ESL_HYPEREXP *h);
double  esl_hxp_surv   (double x, ESL_HYPEREXP *h);
double  esl_hxp_logsurv(double x, ESL_HYPEREXP *h);
double  esl_hxp_invcdf (double p, ESL_HYPEREXP *h);

double  esl_hxp_generic_pdf   (double x, void *params);
double  esl_hxp_generic_cdf   (double x, void *params);
double  esl_hxp_generic_surv  (double x, void *params);
double  esl_hxp_generic_invcdf(double x, void *params);

int esl_hxp_Plot(FILE *fp, ESL_HYPEREXP *h,
			double (*func)(double x, ESL_HYPEREXP *h),
			double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
double esl_hxp_Sample(ESL_RANDOMNESS *r, ESL_HYPEREXP *h);
#endif
#ifdef eslAUGMENT_MINIMIZER
int esl_hxp_FitGuess   (double *x, int n, ESL_HYPEREXP *h);
int esl_hxp_FitComplete(double *x, int n, ESL_HYPEREXP *h);
#ifdef eslAUGMENT_HISTOGRAM
int esl_hxp_FitGuessBinned   (ESL_HISTOGRAM *g, ESL_HYPEREXP *h);
int esl_hxp_FitCompleteBinned(ESL_HISTOGRAM *g, ESL_HYPEREXP *h);
#endif
#endif

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/*** End of inlined file: esl_hyperexp.h ***/


/*** Start of inlined file: esl_mem.h ***/
#ifndef eslMEM_INCLUDED
#define eslMEM_INCLUDED

int       esl_mem_strtoi32(char *p, esl_pos_t n, int base, int *opt_nc, int32_t *opt_val);
int       esl_memnewline(const char *p, esl_pos_t n, esl_pos_t *ret_nline, int *ret_nterm);
int       esl_memtok(char **p, esl_pos_t *n, const char *delim, char **ret_tok, esl_pos_t *ret_toklen);
esl_pos_t esl_memspn (char *p, esl_pos_t n, const char *allow);
esl_pos_t esl_memcspn(char *p, esl_pos_t n, const char *disallow);
int       esl_memstrcmp     (const char *p, esl_pos_t n, const char *s);
int       esl_memstrpfx     (const char *p, esl_pos_t n, const char *s);
int       esl_memstrcontains(const char *p, esl_pos_t n, const char *s);
int       esl_memstrdup(const char *p, esl_pos_t n, char **ret_s);
int       esl_memstrcpy(const char *p, esl_pos_t n, char *dest);
int       esl_memtof(const char *p, esl_pos_t n, float  *ret_val);
int       esl_memtod(const char *p, esl_pos_t n, double *ret_val);
int       esl_mem_IsReal(const char *p, esl_pos_t n);

#endif /*eslMEM_INCLUDED*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_mem.h ***/


/*** Start of inlined file: esl_mixgev.h ***/
#ifndef eslMIXGEV_INCLUDED
#define eslMIXGEV_INCLUDED

#ifdef eslAUGMENT_RANDOM

#endif

typedef struct {
  double *q;			/* mixture coefficients      [0..K-1]*/
  double *mu;			/* location parameters       [0..K-1]*/
  double *lambda;		/* scale parameters          [0..K-1]*/
  double *alpha;		/* shape parameters          [0..K-1]*/
  double *wrk;			/* tmp vector needed for logpdf calc */
  int    *isgumbel;		/* flag:TRUE to constrain k to Gumbel*/
  int     K;			/* # of components                   */
} ESL_MIXGEV;

ESL_MIXGEV *esl_mixgev_Create(int K);
void        esl_mixgev_Destroy(ESL_MIXGEV *mg);
int         esl_mixgev_Copy(ESL_MIXGEV *dest, ESL_MIXGEV *src);
int         esl_mixgev_ForceGumbel(ESL_MIXGEV *mg, int which);

double      esl_mixgev_pdf    (double x, ESL_MIXGEV *mg);
double      esl_mixgev_logpdf (double x, ESL_MIXGEV *mg);
double      esl_mixgev_cdf    (double x, ESL_MIXGEV *mg);
double      esl_mixgev_logcdf (double x, ESL_MIXGEV *mg);
double      esl_mixgev_surv   (double x, ESL_MIXGEV *mg);
double      esl_mixgev_logsurv(double x, ESL_MIXGEV *mg);
double      esl_mixgev_invcdf (double p, ESL_MIXGEV *mg);

double      esl_mixgev_generic_pdf   (double x, void *params);
double      esl_mixgev_generic_cdf   (double x, void *params);
double      esl_mixgev_generic_surv  (double x, void *params);
double      esl_mixgev_generic_invcdf(double p, void *params);

int         esl_mixgev_Plot(FILE *fp, ESL_MIXGEV *mg,
				   double (*func)(double x, ESL_MIXGEV *mg),
				   double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
double      esl_mixgev_Sample(ESL_RANDOMNESS *r, ESL_MIXGEV *mg);
int         esl_mixgev_FitGuess(ESL_RANDOMNESS *r, double *x, int n,
				       ESL_MIXGEV *mg);
#endif /*eslAUGMENT_RANDOM*/

#ifdef eslAUGMENT_MINIMIZER
int         esl_mixgev_FitComplete(double *x, int n, ESL_MIXGEV *mg);
#endif /*eslAUGMENT_MINIMIZER*/

#endif /*eslMIXGEV_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_mixgev.h ***/


/*** Start of inlined file: esl_mpi.h ***/
#if defined(HAVE_MPI) && defined(eslLIBRARY)
#ifndef eslMPI_INCLUDED
#define eslMPI_INCLUDED
#include <mpi.h>



/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_stopwatch.h ***/

/* Many MPI implementations are not MPI2.2 compliant, and do not
 * support new MPI2.2 datatypes; work around that absence. [J10/152]
 * This configuration is better here than esl_config.h.in, because
 * we need to #include <mpi.h> first to see if the system MPI does
 * the right thing, and esl_config.h.in is intended to be included
 * BEFORE any system includes.
 */
#if MPI_VERSION < 2 || MPI_SUBVERSION < 2
#ifndef MPI_INT64_T
#define MPI_INT64_T  MPI_LONG_LONG_INT
#endif
#ifndef MPI_UINT64_T
#define MPI_UINT64_T MPI_UNSIGNED_LONG_LONG
#endif
#ifndef MPI_UINT32_T
#define MPI_UINT32_T MPI_UNSIGNED
#endif
#ifndef MPI_INT16_T
#define MPI_INT16_T  MPI_SHORT
#endif
#ifndef MPI_UINT8_T
#define MPI_UINT8_T  MPI_UNSIGNED_CHAR
#endif
#endif /*MPI_VERSION,MPI_SUBVERSION*/

/* 1. Communicating optional arrays */
int esl_mpi_PackOpt(void *inbuf, int incount, MPI_Datatype type, void *pack_buf,
			   int pack_buf_size, int *position, MPI_Comm comm);
int esl_mpi_PackOptSize(void *inbuf, int incount, MPI_Datatype type, MPI_Comm comm, int *ret_n);
int esl_mpi_UnpackOpt(void *pack_buf, int pack_buf_size, int *pos, void **outbuf,
			     int *opt_n, MPI_Datatype type, MPI_Comm comm);

/* 2. Communicating ESL_SQ (single sequences) */
int esl_sq_MPISend(ESL_SQ *sq, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
int esl_sq_MPIPackSize(ESL_SQ *sq, MPI_Comm comm, int *ret_n);
int esl_sq_MPIPack(ESL_SQ *sq, char *buf, int n, int *pos, MPI_Comm comm);
int esl_sq_MPIUnpack(const ESL_ALPHABET *abc, char *buf, int n, int *pos, MPI_Comm comm, ESL_SQ **ret_sq);
int esl_sq_MPIRecv(int source, int tag, MPI_Comm comm, const ESL_ALPHABET *abc,
			  char **buf, int *nalloc, ESL_SQ **ret_sq);

/* 3. Communicating ESL_MSA (multiple sequence alignments) */
int esl_msa_MPISend(const ESL_MSA *msa, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
int esl_msa_MPIPackSize(const ESL_MSA *msa, MPI_Comm comm, int *ret_n);
int esl_msa_MPIPack(const ESL_MSA *msa, char *buf, int n, int *position, MPI_Comm comm);
int esl_msa_MPIUnpack(const ESL_ALPHABET *abc, char *buf, int n, int *pos, MPI_Comm comm, ESL_MSA **ret_msa);
int esl_msa_MPIRecv(int source, int tag, MPI_Comm comm, const ESL_ALPHABET *abc, char **buf, int *nalloc, ESL_MSA **ret_msa);

/* 4. Communicating ESL_STOPWATCH (process timing) */
int esl_stopwatch_MPIReduce(ESL_STOPWATCH *w, int root, MPI_Comm comm);

#endif /*eslMPI_INCLUDED*/
#endif /*HAVE_MPI && eslLIBRARY*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *
 * SRE, Sat Jun  2 09:07:25 2007 [Janelia]
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/*** End of inlined file: esl_mpi.h ***/


/*** Start of inlined file: esl_msacluster.h ***/
#ifndef eslMSACLUSTER_INCLUDED
#define eslMSACLUSTER_INCLUDED

int esl_msacluster_SingleLinkage(const ESL_MSA *msa, double maxid,
					int **opt_c, int **opt_nin, int *opt_nc);

#endif /*eslMSACLUSTER_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/*** End of inlined file: esl_msacluster.h ***/


/*** Start of inlined file: esl_msafile2.h ***/
#ifndef eslMSAFILE2_INCLUDED
#define eslMSAFILE2_INCLUDED

#ifdef eslAUGMENT_ALPHABET

#endif
#ifdef eslAUGMENT_KEYHASH

#endif
#ifdef eslAUGMENT_SSI

#endif

/* Object: ESL_MSAFILE2
 *
 * Defines an alignment file that we open for reading,
 * in our legacy version. See ESL_MSAFILE (esl_msafile.c) for the
 * preferred version.
 */
typedef struct {
  FILE *f;                      /* open file pointer                         */
  char *fname;			/* name of file. used for diagnostic output  */
  int   linenumber;		/* what line are we on in the file           */
  char  errbuf[eslERRBUFSIZE];  /* buffer for holding parse error info       */

  char *buf;			/* buffer for line input w/ sre_fgets()      */
  int   buflen;			/* current allocated length for buf          */

  int   do_gzip;		/* TRUE if f is "gzip -dc |" (will pclose(f))*/
  int   do_stdin;		/* TRUE if f is stdin (won't close f)        */
  int   format;			/* format of alignment file we're reading    */

  int   do_digital;		/* TRUE to digitize seqs directly into ax    */
#if defined(eslAUGMENT_ALPHABET)
  const ESL_ALPHABET *abc;	/* AUGMENTATION (alphabet): digitized input  */
#else
  void               *abc;
#endif

#if defined(eslAUGMENT_SSI)		/* AUGMENTATION: SSI indexing of an MSA db   */
  ESL_SSI *ssi;		        /* open SSI index file; or NULL, if none.    */
#else
  void    *ssi;
#endif

  ESL_MSA *msa_cache;		/* occasional lookahead at next MSA; GuessAlphabet() */
} ESL_MSAFILE2;

/* 1. The ESL_MSAFILE2 object */
int  esl_msafile2_Open(const char *filename, const char *env, ESL_MSAFILE2 **ret_afp);
#ifdef eslAUGMENT_ALPHABET
int  esl_msafile2_OpenDigital(const ESL_ALPHABET *abc, const char *filename, const char *env, ESL_MSAFILE2 **ret_afp);
#endif
void esl_msafile2_Close(ESL_MSAFILE2 *afp);

/* 2. Memory efficient reading/writing in Pfam format (augmentation: keyhash, for regurgitating some but not all seqs) */
int   esl_msafile2_ReadInfoPfam(ESL_MSAFILE2 *afp, FILE *listfp, ESL_ALPHABET *abc, int64_t known_alen, char *known_rf, char *known_ss_cons, ESL_MSA **ret_msa,
				       int *opt_nseq, int64_t *opt_alen, int *opt_ngs, int *opt_maxname, int *opt_maxgf, int *opt_maxgc, int *opt_maxgr,
				       double ***opt_abc_ct, double ***opt_pp_ct, double ****opt_bp_ct, int **opt_spos_ct, int **opt_epos_ct);
#ifdef eslAUGMENT_KEYHASH
int   esl_msafile2_RegurgitatePfam(ESL_MSAFILE2 *afp, FILE *ofp, int maxname, int maxgf, int maxgc, int maxgr,
					  int do_header, int do_trailer, int do_blanks, int do_comments, int do_gf,
					  int do_gs, int do_gc, int do_gr, int do_aseq, ESL_KEYHASH *seqs2regurg, ESL_KEYHASH *seqs2skip,
					  int *useme, int *add2me, int exp_alen, char gapchar2add, int *opt_nseq_read, int *opt_nseq_written);
#endif

#endif //eslMSAFILE2_INCLUDED
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_msafile2.h ***/


/*** Start of inlined file: esl_msashuffle.h ***/
#ifndef eslMSASHUFFLE_INCLUDED
#define eslMSASHUFFLE_INCLUDED

#ifdef eslAUGMENT_ALPHABET

#endif

/* 1. Randomizing MSAs by column. */
int esl_msashuffle_Shuffle  (ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_MSA *shuf);
int esl_msashuffle_Bootstrap(ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_MSA *bootsample);

/* 2. Permuting the sequence order */
int esl_msashuffle_PermuteSequenceOrder(ESL_RANDOMNESS *r, ESL_MSA *msa);

/* 3. Shuffling pairwise (QRNA) alignments */
#ifdef eslAUGMENT_ALPHABET
int esl_msashuffle_CQRNA(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, char    *x, char    *y, char    *xs, char    *ys);
int esl_msashuffle_XQRNA(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, ESL_DSQ *x, ESL_DSQ *y, ESL_DSQ *xs, ESL_DSQ *ys);
#endif /*eslAUGMENT_ALPHABET*/

#endif /*eslMSASHUFFLE_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_msashuffle.h ***/


/*** Start of inlined file: esl_msaweight.h ***/
#ifndef eslMSAWEIGHT_INCLUDED
#define eslMSAWEIGHT_INCLUDED

int esl_msaweight_GSC(ESL_MSA *msa);
int esl_msaweight_PB(ESL_MSA *msa);
int esl_msaweight_BLOSUM(ESL_MSA *msa, double maxid);
int esl_msaweight_IDFilter(const ESL_MSA *msa, double maxid, ESL_MSA **ret_newmsa);

#endif /*eslMSAWEIGHT_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_msaweight.h ***/


/*** Start of inlined file: esl_normal.h ***/
#ifndef eslNORMAL_INCLUDED
#define eslNORMAL_INCLUDED

double esl_normal_pdf   (double x, double mu, double sigma);
double esl_normal_logpdf(double x, double mu, double sigma);
double esl_normal_cdf   (double x, double mu, double sigma);
double esl_normal_surv  (double x, double mu, double sigma);

double esl_normal_generic_pdf (double x, void *params);
double esl_normal_generic_cdf (double x, void *params);
double esl_normal_generic_surv(double x, void *params);

#endif /*eslNORMAL_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_normal.h ***/


/*** Start of inlined file: esl_paml.h ***/
#ifndef eslPAML_INCLUDED
#define eslPAML_INCLUDED

#include <stdio.h>

int esl_paml_ReadE(FILE *fp, ESL_DMATRIX *E, double *pi);

#endif /*eslPAML_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_paml.h ***/


/*** Start of inlined file: esl_ratematrix.h ***/
#ifndef eslRATEMATRIX_INCLUDED
#define eslRATEMATRIX_INCLUDED

/* 1. Setting standard rate matrix models. */
int esl_rmx_SetWAG(ESL_DMATRIX *Q, double *pi);
int esl_rmx_SetJukesCantor(ESL_DMATRIX *Q);
int esl_rmx_SetKimura(ESL_DMATRIX *Q, double alpha, double beta);
int esl_rmx_SetF81(ESL_DMATRIX *Q, double *pi);
int esl_rmx_SetHKY(ESL_DMATRIX *Q, double *pi, double alpha, double beta);

/* 2. Debugging routines for validating or dumping rate matrices. */
int esl_rmx_ValidateP(ESL_DMATRIX *P, double tol, char *errbuf);
int esl_rmx_ValidateQ(ESL_DMATRIX *Q, double tol, char *errbuf);

/* 3. Other routines in the exposed ratematrix API. */
int    esl_rmx_ScaleTo(ESL_DMATRIX *Q, double *pi, double unit);
int    esl_rmx_E2Q(ESL_DMATRIX *E, double *pi, ESL_DMATRIX *Q);
double esl_rmx_RelativeEntropy(ESL_DMATRIX *P, double *pi);
double esl_rmx_ExpectedScore  (ESL_DMATRIX *P, double *pi);

#endif /*eslRATEMATRIX_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_ratematrix.h ***/


/*** Start of inlined file: esl_recorder.h ***/
#ifndef eslRECORDER_INCLUDED
#define eslRECORDER_INCLUDED

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#include <stdio.h>

/* Object: ESL_RECORDER
 *
 * A history of a line-based input stream.
 * Allows (limited) rewinding in nonrewindable input streams;
 * also allows block-based parsing (as opposed to line-based).
 *
 * The history is kept in a rolling array of string ptrs. The
 * bookkeeping involved in indexing this array can be confusing.
 *
 * Lines in the file are numbered 0..N-1.
 * (N isn't known; we'll be reading sequentially.)
 * Lines in the recorder are 0..nalloc-1.
 *
 * The recorder keeps track of how many lines it has read so far, in
 * <nread>.
 *
 * The recorder can be backed up to any previous line. It sets <ncurr>
 * to be the number of lines it *appears* to have read so far;
 * that is, the next line it will return to the caller, upon a call
 * to esl_recorder_Read(), is line <ncurr>.
 *
 * A window of MIN(nread, nalloc) lines is stored;
 * consisting of line numbers MAX(baseline, nread-nalloc) .. nread-1).
 *
 * A line n in the file (0..n..N-1) corresponds to
 * an index i in the recorder by these transforms:
 *    i = (n-baseline) % nalloc
 *    n = i + MAX(baseline, nread-nalloc)
 *
 * Normally the baseline for the modulo calculation is just 0.
 *
 * The line array is circularly permuted (out of order) when
 * (nread-baseline) / nalloc != 0.
 */
typedef struct {
  FILE    *fp;		/* stream that we're reading line by line           */

  char   **line;	/* lines from input, line[0..nalloc-1]              */
  int      nalloc;	/* max number of lines remembered                   */
  int     *lalloc;	/* alloc for each line[0..nalloc-1][0..lalloc[i]-1] */
  off_t   *offset;	/* disk offsets to starts of each line              */

  int      nread;       /* max # of lines read from file in any pass [1..]  */
  int      ncurr;       /* # of lines into file in current pass      [1..]  */

  int      baseline;	/* line origin for n<->i transform [0..]            */
  int      markline;	/* line origin for start of current block [-1;0..]  */
} ESL_RECORDER;

ESL_RECORDER *esl_recorder_Create    (FILE *fp, int maxlines);
int           esl_recorder_ResizeTo  (ESL_RECORDER *rc, int new_maxlines);
int           esl_recorder_GetFirst  (ESL_RECORDER *rc);
int           esl_recorder_GetLast   (ESL_RECORDER *rc);
int           esl_recorder_GetCurrent(ESL_RECORDER *rc);
int           esl_recorder_GetNext   (ESL_RECORDER *rc);
void          esl_recorder_Destroy   (ESL_RECORDER *rc);

int           esl_recorder_Read(ESL_RECORDER *rc, char **opt_line);
int           esl_recorder_Position(ESL_RECORDER *rc, int linenumber);
int           esl_recorder_MarkBlock(ESL_RECORDER *rc, int markline);
int           esl_recorder_UnmarkBlock(ESL_RECORDER *rc);
int           esl_recorder_GetBlock(ESL_RECORDER *rc, char ***opt_lines, int **opt_lalloc, off_t **opt_offset, int *opt_nlines);

#endif /*eslRECORDER_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_recorder.h ***/


/*** Start of inlined file: esl_regexp.h ***/
#ifndef eslREGEXP_INCLUDED
#define eslREGEXP_INCLUDED

/* ESL_REGEXP_NSUB specifies the maximum number of () expressions
 * in a regexp. The whole regexp counts as one, so 16 allows for
 * parsing out up to 15 tokens from the match.
 */
#define ESL_REGEXP_NSUB 16

/* The esl__regexp structure is from the original Spencer code.
 * It's wrapped by the ESL_REGEXP structure, below.
 */
typedef struct {
  char *startp[ESL_REGEXP_NSUB]; /* ptrs to starts of submatches on target string */
  char *endp[ESL_REGEXP_NSUB];   /* ptrs to 1 char after ends of submatches */
  char regstart;		 /* Internal use only. */
  char reganch;		         /* Internal use only. */
  char *regmust;		 /* Internal use only. */
  int regmlen;		         /* Internal use only. */
  char program[1];	         /* Unwarranted chumminess with compiler. */
} esl__regexp;

/* This looks sort of stupid, wrapping a single ptr in a structure, but we
 * want the machine to be persistent even if different NDFAs are
 * compiled and used. Without this persistency, we would have to
 * create/destroy every time we used a different pattern, instead of
 * one create/destroy per block of code that uses regex matching
 * functionaility.
 *
 * Plus, if we ever need to keep other persistent info
 * beyond Spencer's compiled NDFA (which we'd rather not mess
 * with), we have a place to put it.
 */
typedef struct {
  esl__regexp *ndfa;	 /* a compiled regexp */
} ESL_REGEXP;

/* Declaration of functions in the API
 */

ESL_REGEXP *esl_regexp_Create(void);
void        esl_regexp_Destroy(ESL_REGEXP *machine);

int  esl_regexp_Match(ESL_REGEXP *machine, const char *pattern, const char *s);
int  esl_regexp_Compile(ESL_REGEXP *machine, const char *pattern);
int  esl_regexp_MultipleMatches(ESL_REGEXP *machine, char **sptr);

char *esl_regexp_SubmatchDup(ESL_REGEXP *machine, int elem);
int   esl_regexp_SubmatchCopy(ESL_REGEXP *machine, int elem, char *buffer, int nc);
int   esl_regexp_SubmatchCoords(ESL_REGEXP *machine, char *origin, int elem,
				       int *ret_start, int *ret_end);
int   esl_regexp_ParseCoordString(const char *cstring, uint32_t *ret_start, uint32_t *ret_end);

#endif /*eslREGEXP_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_regexp.h ***/


/*** Start of inlined file: esl_rootfinder.h ***/
#ifndef ESL_ROOTFINDER_INCLUDED
#define ESL_ROOTFINDER_INCLUDED

typedef struct {
  int   (*func)(double, void*, double*);
  int   (*fdf) (double, void*, double*, double*);
  void   *params;

  double xl;
  double fl;
  double xr;
  double fr;

  double x0;
  double f0;

  double x;
  double fx;
  double dfx;
  int    iter;

  double abs_tolerance;
  double rel_tolerance;
  double residual_tol;
  int    max_iter;
} ESL_ROOTFINDER;

ESL_ROOTFINDER *esl_rootfinder_Create   (int (*func)(double, void*, double*),          void *params);
ESL_ROOTFINDER *esl_rootfinder_CreateFDF(int (*fdf) (double, void*, double*, double*), void *params);

int esl_rootfinder_SetBrackets(ESL_ROOTFINDER *R, double xl, double xr);
int esl_rootfinder_SetAbsoluteTolerance(ESL_ROOTFINDER *R, double tol);
int esl_rootfinder_SetRelativeTolerance(ESL_ROOTFINDER *R, double tol);
int esl_rootfinder_SetResidualTolerance(ESL_ROOTFINDER *R, double tol);
int esl_rootfinder_SetMaxIterations(ESL_ROOTFINDER *R, int maxiter);
void esl_rootfinder_Destroy(ESL_ROOTFINDER *R);

int esl_root_Bisection(ESL_ROOTFINDER *R, double xl, double xr, double *ret_x);
int esl_root_NewtonRaphson(ESL_ROOTFINDER *R, double guess, double *ret_x);

#endif /*eslROOTFINDER_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_rootfinder.h ***/


/*** Start of inlined file: esl_scorematrix.h ***/
#ifndef eslSCOREMATRIX_INCLUDED
#define eslSCOREMATRIX_INCLUDED

/* ESL_SCOREMATRIX:
 * allocation is in one array in s[0].
 *
 * i,j can range from 0..Kp-1, including all characters valid in the alphabet.
 * Only values for 0..K-1 (canonical alphabet) are mandatory.
 */
typedef struct {
  int **s;			/* s[i][j] is the score of aligning residue i,j; i,j range 0..Kp-1 */
  int   K;			/* size of base alphabet (duplicate of S->abc_r->K) */
  int   Kp;			/* full size of s[][], including degeneracies (duplicate of S->abc_r->Kp) */

  /* bookkeeping for degenerate residues */
  char *isval;			/* array 0..Kp-1: which residues of alphabet have valid scores in S. */
  const ESL_ALPHABET *abc_r;	/* reference to the alphabet: includes K, Kp, and sym order */

  /* bookkeeping that lets us output exactly the residue order we read in a matrix file */
  int   nc;			/* number of residues with scores (inclusive of *, if present) */
  char *outorder;		/* NUL-terminated string 0..nc-1 giving order of residues in col/row labels   */

  char *name;			/* optional: name of score matrix; or NULL */
  char *path;			/* optional: full path to file that score matrix was read from; or NULL  */
} ESL_SCOREMATRIX;

/* 1. The ESL_SCOREMATRIX object. */
ESL_SCOREMATRIX *esl_scorematrix_Create(const ESL_ALPHABET *abc);
int              esl_scorematrix_Copy(const ESL_SCOREMATRIX *src, ESL_SCOREMATRIX *dest);
ESL_SCOREMATRIX *esl_scorematrix_Clone(const ESL_SCOREMATRIX *S);
int              esl_scorematrix_Compare(const ESL_SCOREMATRIX *S1, const ESL_SCOREMATRIX *S2);
int              esl_scorematrix_CompareCanon(const ESL_SCOREMATRIX *S1, const ESL_SCOREMATRIX *S2);
int              esl_scorematrix_Max(const ESL_SCOREMATRIX *S);
int              esl_scorematrix_Min(const ESL_SCOREMATRIX *S);
int              esl_scorematrix_IsSymmetric(const ESL_SCOREMATRIX *S);
int              esl_scorematrix_ExpectedScore(ESL_SCOREMATRIX *S, double *fi, double *fj, double *ret_E);
int              esl_scorematrix_RelEntropy(const ESL_SCOREMATRIX *S, const double *fi, const double *fj,
						   double lambda, double *ret_D);
int              esl_scorematrix_JointToConditionalOnQuery(const ESL_ALPHABET *abc, ESL_DMATRIX *P);
void             esl_scorematrix_Destroy(ESL_SCOREMATRIX *S);

/* 2. Some classic score matrices */
int              esl_scorematrix_Set(const char *name, ESL_SCOREMATRIX *S);
int              esl_scorematrix_SetIdentity(ESL_SCOREMATRIX *S);

/* 3. Deriving a score matrix probabilistically */
int              esl_scorematrix_SetFromProbs(ESL_SCOREMATRIX *S, double lambda, const ESL_DMATRIX *P,
						     const double *fi, const double *fj);
int              esl_scorematrix_SetWAG(ESL_SCOREMATRIX *S, double lambda, double t);

/* 4. Reading/writing score matrices. */
int  esl_scorematrix_Read(ESL_FILEPARSER *efp, const ESL_ALPHABET *abc, ESL_SCOREMATRIX **ret_S);
int  esl_scorematrix_Write(FILE *fp, const ESL_SCOREMATRIX *S);

/* 5. Implicit probabilistic basis, I: given bg. */
int esl_scorematrix_ProbifyGivenBG(const ESL_SCOREMATRIX *S, const double *fi, const double *fj,
					  double *opt_lambda, ESL_DMATRIX **opt_P);

/* 6. Implicit probabilistic basis, II: bg unknown. */
int esl_scorematrix_Probify(const ESL_SCOREMATRIX *S, ESL_DMATRIX **opt_P,
				   double **opt_fi, double **opt_fj, double *opt_lambda);

#endif /*eslSCOREMATRIX_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_scorematrix.h ***/


/*** Start of inlined file: esl_sse.h ***/
/* Vectorized routines for Intel/AMD, using Streaming SIMD Extensions (SSE).
 *
 * This header file, unusually, provides many complete function
 * implementations; this is so that they can be inlined by the
 * compiler, for maximum efficiency.
 *
 * Contents:
 *    1. Function declarations (from esl_sse.c)
 *    2. Inlined utilities for ps vectors (4 floats in __m128)
 *    3. Inlined utilities for epu8 vectors (16 uchars in __m128i)
 */
#ifdef HAVE_SSE2
#ifndef eslSSE_INCLUDED
#define eslSSE_INCLUDED

#include <stdio.h>
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

/* Some compilers (gcc 3.4) did not implement SSE2 cast functions
 * on the theory that they're unnecessary no-ops -- but then
 * code that has proper SSE cast calls doesn't compile. Provide
 * the no-ops.
 */
#ifndef HAVE_SSE2_CAST
#define _mm_castps_si128(x) (__m128i)(x)
#define _mm_castsi128_ps(x) (__m128)(x)
#endif

/*****************************************************************
 * 1. Function declarations (from esl_sse.c)
 *****************************************************************/
__m128  esl_sse_logf(__m128 x);
__m128  esl_sse_expf(__m128 x);
void    esl_sse_dump_ps(FILE *fp, __m128 v);

/*****************************************************************
 * 2. Inline utilities for ps vectors (4 floats in __m128)
 *****************************************************************/

/* Function:  esl_sse_select_ps()
 * Synopsis:  SSE equivalent of <vec_sel()>
 *
 * Purpose:   Vector select. Returns a vector <r[z] = a[z]> where <mask[z]>
 *            is all 0's; <r[z] = b[z]> where <mask[z]> is all 1's.
 *
 *            Useful for avoiding conditional branches. For example,
 *            to implement \ccode{if (a > 0) a += a;}:
 *
 *            \begin{cchunk}
 *              mask = _mm_cmpgt_ps(a, _mm_setzero_ps());
 *              twoa = _mm_add_ps(a, a);
 *              a    = esl_sse_select_ps(a, twoa, mask);
 *            \end{cchunk}
 *
 * Notes:     As recommended by the Altivec/SSE Migration Guide,
 *            Apple Computer, Inc.
 */
static inline __m128
esl_sse_select_ps(__m128 a, __m128 b, __m128 mask)
{
  b = _mm_and_ps(b, mask);
  a = _mm_andnot_ps(mask, a);
  return _mm_or_ps(a,b);
}

/* Function:  esl_sse_any_gt_ps()
 * Synopsis:  Returns TRUE if any a[z] > b[z]
 *
 * Purpose:   Returns TRUE if any a[z] > b[z] in two
 *            <ps> vectors of floats.
 *
 * Xref:      From Apple Altivec/SSE migration guide.
 */
static inline int
esl_sse_any_gt_ps(__m128 a, __m128 b)
{
  __m128 mask    = _mm_cmpgt_ps(a,b);
  int   maskbits = _mm_movemask_ps( mask );
  return maskbits != 0;
}

/* Function:  esl_sse_hmax_ps()
 * Synopsis:  Find the maximum of elements in a vector.
 *
 * Purpose:   Find the maximum valued element in the four float elements
 *            in <a>, and return that maximum value in <*ret_max>.
 *
 * Xref:      J3/90 for benchmarking of some alternative implementations.
 */
static inline void
esl_sse_hmax_ps(__m128 a, float *ret_max)
{
  a = _mm_max_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(0, 3, 2, 1)));
  a = _mm_max_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(1, 0, 3, 2)));
  _mm_store_ss(ret_max, a);
}

/* Function:  esl_sse_hmin_ps()
 * Synopsis:  Find the minimum of elements in a vector.
 *
 * Purpose:   Find the minimum valued element in the four float elements
 *            in <a> and return that minimum value in <*ret_min>.
 */
static inline void
esl_sse_hmin_ps(__m128 a, float *ret_min)
{
  a = _mm_min_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(0, 3, 2, 1)));
  a = _mm_min_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(1, 0, 3, 2)));
  _mm_store_ss(ret_min, a);
}

/* Function:  esl_sse_hsum_ps()
 * Synopsis:  Takes the horizontal sum of elements in a vector.
 *
 * Purpose:   Add the four float elements in vector <a>; return
 *            that sum in <*ret_sum>.
 */
static inline void
esl_sse_hsum_ps(__m128 a, float *ret_sum)
{
  a = _mm_add_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(0, 3, 2, 1)));
  a = _mm_add_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(1, 0, 3, 2)));
  _mm_store_ss(ret_sum, a);
}

/* Function:  esl_sse_rightshift_ps()
 * Synopsis:  Shift vector elements to the right.
 *
 * Purpose:   Returns a vector containing
 *            <{ b[0] a[0] a[1] a[2] }>:
 *            i.e. shift the values in <a> to the
 *            right, and load the first value of
 *            <b> into the first slot.
 */
static inline __m128
esl_sse_rightshift_ps(__m128 a, __m128 b)
{
  return _mm_move_ss(_mm_shuffle_ps(a, a, _MM_SHUFFLE(2, 1, 0, 0)), b);
}

/* Function:  esl_sse_leftshift_ps()
 * Synopsis:  Shift vector elements to the left.
 *
 * Purpose:   Returns a vector containing
 *            <{ a[1] a[2] a[3] b[0]}>:
 *            i.e. shift the values in <a> to the
 *            left and load the first value of
 *            <b> into the first slot.
 */
static inline __m128
esl_sse_leftshift_ps(__m128 a, __m128 b)
{
  register __m128 v = _mm_move_ss(a, b);                 /* now b[0] a[1] a[2] a[3] */
  return _mm_shuffle_ps(v, v, _MM_SHUFFLE(0, 3, 2, 1));  /* now a[1] a[2] a[3] b[0] */
}

/*****************************************************************
 * 3. Inlined utilities for epu8 vectors (16 uchars in __m128i)
 *****************************************************************/

/* Function:  esl_sse_any_gt_epu8()
 * Synopsis:  Returns TRUE if any a[z] > b[z].
 *
 * Purpose:   Return TRUE if any <a[z] > b[z]> for <z=0..15>
 *            in two <epu8> vectors of unsigned chars.
 *
 *            We need this incantation because SSE provides
 *            no <cmpgt_epu8> instruction.
 *
 *            For equality tests, note that <cmpeq_epi8> works fine
 *            for unsigned ints though there is no <cmpeq_epu8>
 *            instruction either).
 *
 *            See vec_any_gt
 */
static inline int
esl_sse_any_gt_epu8(__m128i a, __m128i b)
{
  __m128i mask    = _mm_cmpeq_epi8(_mm_max_epu8(a,b), b); /* anywhere a>b, mask[z] = 0x0; elsewhere 0xff */
  int   maskbits  = _mm_movemask_epi8(_mm_xor_si128(mask,  _mm_cmpeq_epi8(mask, mask))); /* the xor incantation is a bitwise inversion */
  return maskbits != 0;
}
static inline int
esl_sse_any_gt_epi16(__m128i a, __m128i b)
{
  return (_mm_movemask_epi8(_mm_cmpgt_epi16(a,b)) != 0);
}

/* Function:  esl_sse_hmax_epu8()
 * Synopsis:  Return the max of the 16 elements in epu8 vector.
 *
 * Purpose:   Returns the maximum value of the 16 elements in
 *            an <epu8> vector.
 */
static inline uint8_t
esl_sse_hmax_epu8(__m128i a)
{
  a = _mm_max_epu8(a, _mm_srli_si128(a, 8));
  a = _mm_max_epu8(a, _mm_srli_si128(a, 4));
  a = _mm_max_epu8(a, _mm_srli_si128(a, 2));
  a = _mm_max_epu8(a, _mm_srli_si128(a, 1));
  return (uint8_t) _mm_extract_epi16(a, 0);   /* only low-order 8 bits set; so _epi16 or _epi8 equiv; _epi8 is SSE4.1 */
}

/* Function:  esl_sse_hmax_epi16()
 * Synopsis:  Return the max of the 8 elements in epi16 vector.
 *
 * Purpose:   Returns the maximum value of the 16 elements in
 *            an <epu8> vector.
 */
static inline int16_t
esl_sse_hmax_epi16(__m128i a)
{
  a = _mm_max_epi16(a, _mm_srli_si128(a, 8));
  a = _mm_max_epi16(a, _mm_srli_si128(a, 4));
  a = _mm_max_epi16(a, _mm_srli_si128(a, 2));
  return (int16_t) _mm_extract_epi16(a, 0);   /* only low-order 8 bits set; so _epi16 or _epi8 equiv; _epi8 is SSE4.1 */
}

#endif /*eslSSE_INCLUDED*/
#endif /*HAVE_SSE2*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/*** End of inlined file: esl_sse.h ***/


/*** Start of inlined file: esl_stack.h ***/
#ifndef eslSTACK_INCLUDED
#define eslSTACK_INCLUDED

#define ESL_STACK_INITALLOC 128	/* initial allocation; realloc by doubling  */

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif
#ifdef eslAUGMENT_RANDOM

#endif /*eslAUGMENT_RANDOM*/

typedef struct esl_stack_s {
  int   *idata;			/* integer data stack                       */
  void **pdata;			/* pointer data stack                       */
  char  *cdata;			/* character data stack                     */

  int  n;			/* current (topmost) elem in data           */
  int  nalloc;			/* # of elems allocated right now           */

#ifdef HAVE_PTHREAD
  int              do_mutex;	/* TRUE if we need to mutex-protect this stack */
  int              do_cond;	/* TRUE if pushers want to notify poppers      */
  pthread_mutex_t *mutex;	/* protect while operating on stacks           */
  pthread_cond_t  *cond;	/* for pushers to notify poppers               */
#endif
} ESL_STACK;

ESL_STACK *esl_stack_ICreate(void);
ESL_STACK *esl_stack_CCreate(void);
ESL_STACK *esl_stack_PCreate(void);

int        esl_stack_Reuse(ESL_STACK *s);
void       esl_stack_Destroy(ESL_STACK *s);

int esl_stack_IPush(ESL_STACK *ns, int x);
int esl_stack_CPush(ESL_STACK *cs, char c);
int esl_stack_PPush(ESL_STACK *ps, void *p);

int esl_stack_IPop(ESL_STACK *ns, int   *ret_x);
int esl_stack_CPop(ESL_STACK *cs, char  *ret_c);
int esl_stack_PPop(ESL_STACK *ps, void **ret_p);

int esl_stack_ObjectCount(ESL_STACK *s);

char *esl_stack_Convert2String(ESL_STACK *cs);
int   esl_stack_DiscardTopN(ESL_STACK *s, int n);
int   esl_stack_DiscardSelected(ESL_STACK *s, int (*discard_func)(void *, void *), void *param);

#ifdef eslAUGMENT_RANDOM
int esl_stack_Shuffle(ESL_RANDOMNESS *r, ESL_STACK *s);
#endif /*eslAUGMENT_RANDOM*/

#ifdef HAVE_PTHREAD
int esl_stack_UseMutex   (ESL_STACK *s);
int esl_stack_UseCond    (ESL_STACK *s);
int esl_stack_ReleaseCond(ESL_STACK *s);
#endif
#endif /*eslSTACK_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_stack.h ***/


/*** Start of inlined file: esl_stats.h ***/
#ifndef eslSTATS_INCLUDED
#define eslSTATS_INCLUDED

/*****************************************************************
 * Splitting IEEE754 double-precision float into two uint32_t
 *****************************************************************
 *
 * Currently we only need these macros for one function,
 * esl_stats_erfc(). The Sun Microsystems erfc() code that we've
 * borrowed splits an IEEE754 double into two unsigned 32-bit
 * integers. It uses arcane trickery to deal with endianness at
 * runtime, using incantations like these:
 *    n0 = ((*(int*)&one)>>29)^1     0|1 = bigendian | littleendian
 *    hx = *(n0+(int*)&x);           get high word
 *    (1-n0+(int*)&z) = 0;           set low word
 *
 * Not only is this arcane and dubious, static code checking (using
 * the clang/llvm checker) doesn't like it. I found an improvement
 * in a library called zenilib, at:
 *  http://www-personal.umich.edu/~bazald/l/api/math__private_8h_source.html
 *
 * Here we do the same thing in an ANSI-respecting way using unions,
 * with endianness detected at compile time.
 *
 * The zenilib code also appears to derive from (C) Sun Microsystems
 * code.
 *
 * SRE TODO: insert license/copyright info here
 */

#ifdef WORDS_BIGENDIAN
typedef union
{
 double val;
 struct {
   uint32_t msw;
   uint32_t lsw;
 } parts;
} esl_double_split_t;
#else /* else we're littleendian, such as Intel */
typedef union
{
 double val;
 struct {
   uint32_t lsw;
   uint32_t msw;
 } parts;
} esl_double_split_t;
#endif /*WORDS_BIGENDIAN*/

#define ESL_GET_WORDS(ix0, ix1, d) \
  do { \
	esl_double_split_t esltmp_ds;  \
	esltmp_ds.val = (d);           \
	(ix0) = esltmp_ds.parts.msw;   \
	(ix1) = esltmp_ds.parts.lsw;   \
  } while (0)

#define ESL_GET_HIGHWORD(ix0, d)  \
  do { \
	esl_double_split_t esltmp_ds; \
	esltmp_ds.val = (d);          \
	(ix0) = esltmp_ds.parts.msw;  \
  } while (0)

#define ESL_GET_LOWWORD(ix0, d)   \
  do { \
	esl_double_split_t esltmp_ds; \
	esltmp_ds.val = (d);          \
	(ix0) = esltmp_ds.parts.lsw;  \
  } while (0)

#define ESL_SET_WORDS(d, ix0, ix1) \
   do { \
	esl_double_split_t esltmp_ds;  \
	esltmp_ds.parts.msw = (ix0);   \
	esltmp_ds.parts.lsw = (ix1);   \
	(d) = esltmp_ds.val;           \
   } while (0)

#define ESL_SET_HIGHWORD(d, ix0)  \
   do { \
	esl_double_split_t esltmp_ds; \
	esltmp_ds.val = (d);          \
	esltmp_ds.parts.msw = (ix0);  \
	(d) = esltmp_ds.val;          \
  } while (0)

#define ESL_SET_LOWWORD(d, ix1)   \
   do { \
	esl_double_split_t esltmp_ds; \
	esltmp_ds.val = (d);          \
	esltmp_ds.parts.lsw = (ix1);  \
	(d) = esltmp_ds.val;          \
  } while (0)

/*****************************************************************
 * Function declarations
 *****************************************************************/

/* 1. Summary statistics calculations */
int esl_stats_DMean(const double *x, int n, double *opt_mean, double *opt_var);
int esl_stats_FMean(const float  *x, int n, double *opt_mean, double *opt_var);
int esl_stats_IMean(const int    *x, int n, double *opt_mean, double *opt_var);

/* 2. Special functions */
int    esl_stats_LogGamma(double x, double *ret_answer);
int    esl_stats_Psi(double x, double *ret_answer);
int    esl_stats_IncompleteGamma(double a, double x, double *ret_pax, double *ret_qax);
double esl_stats_erfc(double x);

/* 3. Standard statistical tests */
int esl_stats_GTest(int ca, int na, int cb, int nb, double *ret_G, double *ret_P);
int esl_stats_ChiSquaredTest(int v, double x, double *ret_answer);

/* 4. Data fitting */
int esl_stats_LinearRegression(const double *x, const double *y, const double *sigma, int n,
				      double *opt_a,       double *opt_b,
				      double *opt_sigma_a, double *opt_sigma_b, double *opt_cov_ab,
				      double *opt_cc,      double *opt_Q);

/*****************************************************************
 * Portability
 *****************************************************************/

#ifndef HAVE_ERFC
#define erfc(x)  esl_stats_erfc(x)
#endif

#endif /*eslSTATS_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/*** End of inlined file: esl_stats.h ***/


/*** Start of inlined file: esl_stretchexp.h ***/
#ifndef eslSTRETCHEXP_INCLUDED
#define eslSTRETCHEXP_INCLUDED

double esl_sxp_pdf    (double x, double mu, double lambda, double tau);
double esl_sxp_logpdf (double x, double mu, double lambda, double tau);
double esl_sxp_cdf    (double x, double mu, double lambda, double tau);
double esl_sxp_logcdf (double x, double mu, double lambda, double tau);
double esl_sxp_surv   (double x, double mu, double lambda, double tau);
double esl_sxp_logsurv(double x, double mu, double lambda, double tau);
double esl_sxp_invcdf (double p, double mu, double lambda, double tau);

double esl_sxp_generic_pdf   (double x, void *params);
double esl_sxp_generic_cdf   (double x, void *params);
double esl_sxp_generic_surv  (double x, void *params);
double esl_sxp_generic_invcdf(double p, void *params);

int esl_sxp_Plot(FILE *fp, double mu, double lambda, double tau,
			double (*func)(double x, double mu, double lambda, double tau),
			double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
double esl_sxp_Sample(ESL_RANDOMNESS *r, double mu, double lambda, double tau);
#endif

#ifdef eslAUGMENT_MINIMIZER
int esl_sxp_FitComplete(double *x, int n,
			       double *ret_mu, double *ret_lambda, double *ret_tau);
#ifdef eslAUGMENT_HISTOGRAM
int esl_sxp_FitCompleteBinned(ESL_HISTOGRAM *g,
				     double *ret_mu, double *ret_lambda, double *ret_tau);
#endif /*eslAUGMENT_HISTOGRAM*/
#endif /*eslAUGMENT_MINIMIZER*/

#endif /*eslSTRETCHEXP_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_stretchexp.h ***/


/*** Start of inlined file: esl_threads.h ***/
#ifndef eslTHREADS_INCLUDED
#define eslTHREADS_INCLUDED

#include <pthread.h>

typedef struct {
  int             threadCount;      /* number of active worker threads                           */
  pthread_t      *threadId;	    /* threadId for each worker thread; [0..threadCount-1]       */
  void          **data;		    /* data pointer for each worker thread; [0..threadCount-1]   */

  int             startThread;      /* number of worker threads currently blocked at start mutex */
  pthread_mutex_t startMutex;	    /* the starting gate                                         */
  pthread_cond_t  startCond;	    /* the signal that workers are synchronized and may start    */

  void           (*func)(void *);   /* each worker thread runs this function; arg is to data[]   */
} ESL_THREADS;

ESL_THREADS *esl_threads_Create(void (*func)(void *));
void         esl_threads_Destroy(ESL_THREADS *obj);

int esl_threads_AddThread     (ESL_THREADS *obj, void *data);
int esl_threads_GetWorkerCount(ESL_THREADS *obj);
int esl_threads_WaitForStart  (ESL_THREADS *obj);
int esl_threads_WaitForFinish (ESL_THREADS *obj);

int   esl_threads_Started (ESL_THREADS *obj, int *ret_workeridx);
void *esl_threads_GetData (ESL_THREADS *obj, int workeridx);
int   esl_threads_Finished(ESL_THREADS *obj, int workeridx);

int esl_threads_CPUCount(int *ret_ncpu);

#endif /*eslTHREADS_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_threads.h ***/


/*** Start of inlined file: esl_tree.h ***/
#ifndef eslTREE_INCLUDED
#define eslTREE_INCLUDED

/* Object: ESL_TREE
 *
 * All trees are represented as rooted trees, starting from
 * node 0. For N taxa, there are N-1 internal nodes, numbered
 * 0..N-2. Taxa on leaves are numbered 0..N-1, and represented
 * in <left> and <right> as negative numbers.
 *
 */
typedef struct {
  int   N;		/* number of taxa */

  /* (Mandatory) information for the internal nodes of a rooted tree.
   * There are N-1 nodes, numbered 0..N-2, with the root at 0,
   * so each array below is indexed [0..N-2].
   * When an internal node has a left or right branch to a taxon,
   * <left>/<right> are <= 0, negative <taxon #>; if they're to
   * be used as array indices, flip their sign.
   * There is no ambiguity between taxon 0/root node 0, because
   * a taxon can't be a parent, and the root node can't be a child.
   * For an unrooted tree, by convention, taxon 0 is the outgroup: T->left[0] = 0,
   * and T->rd[0] = 0.0.
   */
  int    *parent;	/* index of parent of node: values are 0..N-2; parent of root 0 = 0 */
  int    *left;		/* index of left child:  values are -(N-1)..0=taxa; 1..N-2=nodes */
  int    *right;	/* index of right child: values are -(N-1)..0=taxa; 1..N-2=nodes */
  double *ld;	        /* left branch length under node: values are >= 0 */
  double *rd;	        /* right branch length under node: values are >= 0 */
						/* in linkage trees, ld[x]=rd[x]= "height" (linkage value) of node, not branch lengths */

  /* Derived (optional) information, that we can reconstruct if
   * we need to from the mandatory info above.
   */
  int    *taxaparent;   /* for taxa  [0..N-1]: index of its parent node, 0..N-2. [esl_tree_SetTaxaParents()] */
  int    *cladesize;	/* for nodes [0..N-2]: # taxa in this clade, 1..N        [esl_tree_SetCladesizes()]  */

  /* Optional information */
  char  **taxonlabel;	  /* labels for taxa: [0..N-1] array of char strings */
  char  **nodelabel;	  /* labels for nodes: [0..N-2] array of char strings */

  /* Tree mode options. */
  int   is_linkage_tree;	 /* TRUE if this is a linkage tree; if FALSE, it's an additive tree */

  /* Tree output options. */
  int   show_unrooted;	         /* TRUE to output 'root' as a trifurcation (a la PHYLIP) */
  int   show_node_labels;        /* TRUE to output labels for interior nodes */
  int   show_root_branchlength;  /* TRUE to show 0.0 branch length to root node (a la TreeAlign) */
  int   show_branchlengths;	 /* TRUE to output branch lengths */
  int   show_quoted_labels;	 /* TRUE to output ALL labels as quoted labels */
  int   show_numeric_taxonlabels;/* TRUE to output taxa labels as their 0..N-1 indices if no other taxonlabel is present */

  /* Memory allocation information, when growing a tree (on input, for example)
   */
  int     nalloc;	/* current allocated # of taxa */

} ESL_TREE;

/* UPGMA, average-link, minimum-link, and maximum-link clustering are
 * all implemented by one algorithm, cluster_engine(), in esl_tree.c.
 * We define some flags (within the scope of the tree module) to
 * control the behavior, as we call the algorithm engine from four
 * different API functions.
 */
#define eslUPGMA            0
#define eslWPGMA            1
#define eslSINGLE_LINKAGE   2
#define eslCOMPLETE_LINKAGE 3

/* 1. The ESL_TREE object.
 */
ESL_TREE *esl_tree_Create(int ntaxa);
ESL_TREE *esl_tree_CreateGrowable(int nalloc);
ESL_TREE *esl_tree_CreateFromString(char *s);
int       esl_tree_Grow(ESL_TREE *T);
int       esl_tree_SetTaxaParents(ESL_TREE *T);
int       esl_tree_SetCladesizes(ESL_TREE *T);
int       esl_tree_SetTaxonlabels(ESL_TREE *T, char **names);
int       esl_tree_RenumberNodes(ESL_TREE *T);
int       esl_tree_VerifyUltrametric(ESL_TREE *T);
int       esl_tree_Validate(ESL_TREE *T, char *errbuf);
void      esl_tree_Destroy(ESL_TREE *T);

/* 2. Newick format i/o
 */
int  esl_tree_WriteNewick(FILE *fp, ESL_TREE *T);
int  esl_tree_ReadNewick(FILE *fp, char *errbuf, ESL_TREE **ret_T);

/* 3. Tree comparison algorithms.
 */
int esl_tree_Compare(ESL_TREE *T1, ESL_TREE *T2);

/* 4. Clustering algorithms for distance-based tree construction.
 */
int esl_tree_UPGMA(ESL_DMATRIX *D, ESL_TREE **ret_T);
int esl_tree_WPGMA(ESL_DMATRIX *D, ESL_TREE **ret_T);
int esl_tree_SingleLinkage(ESL_DMATRIX *D, ESL_TREE **ret_T);
int esl_tree_CompleteLinkage(ESL_DMATRIX *D, ESL_TREE **ret_T);

/* 5. Generating simulated trees.
 */
int esl_tree_Simulate(ESL_RANDOMNESS *r, int N, ESL_TREE **ret_T);
int esl_tree_ToDistanceMatrix(ESL_TREE *T, ESL_DMATRIX **ret_D);

#endif /*eslTREE_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_tree.h ***/


/*** Start of inlined file: esl_vectorops.h ***/
#ifndef eslVECTOROPS_INCLUDED
#define eslVECTOROPS_INCLUDED

void   esl_vec_DSet(double *vec, int n, double value);
void   esl_vec_FSet(float  *vec, int n, float  value);
void   esl_vec_ISet(int    *vec, int n, int    value);

void   esl_vec_DScale(double *vec, int n, double scale);
void   esl_vec_FScale(float  *vec, int n, float  scale);
void   esl_vec_IScale(int    *vec, int n, int    scale);

void   esl_vec_DIncrement(double *v, int n, double x);
void   esl_vec_FIncrement(float  *v, int n, float  x);
void   esl_vec_IIncrement(int    *v, int n, int    x);

double esl_vec_DSum(double *vec, int n);
float  esl_vec_FSum(float  *vec, int n);
int    esl_vec_ISum(int    *vec, int n);

void   esl_vec_DAdd(double *vec1, const double *vec2, int n);
void   esl_vec_FAdd(float  *vec1, const float  *vec2, int n);
void   esl_vec_IAdd(int    *vec1, const int    *vec2, int n);

void   esl_vec_DAddScaled(double *vec1, double *vec2, double a, int n);
void   esl_vec_FAddScaled(float  *vec1, float  *vec2, float  a, int n);
void   esl_vec_IAddScaled(int    *vec1, int    *vec2, int    a, int n);

void   esl_vec_DCopy(const double *src, const int n, double *dest);
void   esl_vec_FCopy(const float  *src, const int n, float  *dest);
void   esl_vec_ICopy(const int    *src, const int n, int    *dest);

int    esl_vec_DCompare(const double *vec1, const double *vec2, int n, double tol);
int    esl_vec_FCompare(const float  *vec1, const float  *vec2, int n, float tol);
int    esl_vec_ICompare(const int    *vec1, const int    *vec2, int n);

void   esl_vec_DSwap(double *vec1, double *vec2, int n);
void   esl_vec_FSwap(float  *vec1, float  *vec2, int n);
void   esl_vec_ISwap(int    *vec1, int    *vec2, int n);

void   esl_vec_DReverse(double *vec, double *rev, int n);
void   esl_vec_FReverse(float  *vec, float  *rev, int n);
void   esl_vec_IReverse(int    *vec, int    *rev, int n);
void   esl_vec_CReverse(char   *vec, char   *rev, int n);

double esl_vec_DDot(double *vec1, double *vec2, int n);
float  esl_vec_FDot(float  *vec1, float  *vec2, int n);
int    esl_vec_IDot(int    *vec1, int    *vec2, int n);

double esl_vec_DMax(const double *vec, int n);
float  esl_vec_FMax(const float  *vec, int n);
int    esl_vec_IMax(const int    *vec, int n);

double esl_vec_DMin(const double *vec, int n);
float  esl_vec_FMin(const float  *vec, int n);
int    esl_vec_IMin(const int    *vec, int n);

int    esl_vec_DArgMax(const double *vec, int n);
int    esl_vec_FArgMax(const float  *vec, int n);
int    esl_vec_IArgMax(const int    *vec, int n);

int    esl_vec_DArgMin(const double *vec, int n);
int    esl_vec_FArgMin(const float  *vec, int n);
int    esl_vec_IArgMin(const int    *vec, int n);

void   esl_vec_DSortIncreasing(double *vec, int n);
void   esl_vec_FSortIncreasing(float  *vec, int n);
void   esl_vec_ISortIncreasing(int    *vec, int n);

void   esl_vec_DSortDecreasing(double *vec, int n);
void   esl_vec_FSortDecreasing(float  *vec, int n);
void   esl_vec_ISortDecreasing(int    *vec, int n);

int    esl_vec_DDump(FILE *ofp, double *v, int n, char *label);
int    esl_vec_FDump(FILE *ofp, float *v,  int n, char *label);
int    esl_vec_IDump(FILE *ofp, int *v,    int n, char *label);

void   esl_vec_D2F(double *src, int n, float  *dst);
void   esl_vec_F2D(float  *src, int n, double *dst);
void   esl_vec_I2F(int    *src, int n, float  *dst);
void   esl_vec_I2D(int    *src, int n, double *dst);

void   esl_vec_DNorm(double *vec, int n);
void   esl_vec_FNorm(float  *vec, int n);

void   esl_vec_DLog(double *vec, int n);
void   esl_vec_FLog(float  *vec, int n);

double esl_vec_DEntropy(const double *p, int n);
float  esl_vec_FEntropy(const float  *p, int n);

double esl_vec_DRelEntropy(const double *p, const double *f, int n);
float  esl_vec_FRelEntropy(const float  *p, const float  *f, int n);

void   esl_vec_DExp(double *vec, int n);
void   esl_vec_FExp(float  *vec, int n);

double esl_vec_DLogSum(double *vec, int n);
float  esl_vec_FLogSum(float  *vec, int n);

void   esl_vec_DLogNorm(double *vec, int n);
void   esl_vec_FLogNorm(float  *vec, int n);

void   esl_vec_DCDF(double *p, int n, double *cdf);
void   esl_vec_FCDF(float  *p, int n, float  *cdf);

int    esl_vec_DValidate(double *vec, int n, double tol, char *errbuf);
int    esl_vec_FValidate(float  *vec, int n, float  tol, char *errbuf);

int    esl_vec_DLogValidate(double *vec, int n, double tol, char *errbuf);
int    esl_vec_FLogValidate(float  *vec, int n, float  tol, char *errbuf);

#ifdef eslAUGMENT_RANDOM

int esl_vec_DShuffle(ESL_RANDOMNESS *r, double *v, int n);
int esl_vec_FShuffle(ESL_RANDOMNESS *r, float  *v, int n);
int esl_vec_IShuffle(ESL_RANDOMNESS *r, int    *v, int n);
#endif

#endif /* eslVECTOROPS_INCLUDED */

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/*** End of inlined file: esl_vectorops.h ***/


/*** Start of inlined file: esl_vmx.h ***/
#ifdef HAVE_VMX
/* Vectorized routines for PowerPC, using Altivec.
 *
 * SVN $Id$
 * SVN $URL$
 */
#ifndef eslVMX_INCLUDED
#define eslVMX_INCLUDED

#include <stdio.h>
#ifndef __APPLE_ALTIVEC__
#include <altivec.h>
#endif

vector float esl_vmx_logf(vector float x);
vector float esl_vmx_expf(vector float x);
void         esl_vmx_dump_vecfloat(FILE *fp, vector float v);

/* Function:  esl_vmx_set_float()
 * Synopsis:  Fills float vector with x.
 *
 * Purpose:   Sets all elements in the vector <float> to x.
 */
static inline vector float
esl_vmx_set_float(float x)
{
  vector float v;
  vector unsigned char p;

  v = vec_lde(0, &x);
  p = vec_lvsl(0, &x);
  v = vec_perm(v, v, p);
  v = vec_splat(v, 0);
  return v;
}

/* Function:  esl_vmx_set_s16()
 * Synopsis:  Fills short vector with x.
 *
 * Purpose:   Sets all elements in the vector <signed short> to x.
 */
static inline vector signed short
esl_vmx_set_s16(signed short x)
{
  vector signed short v;
  vector unsigned char p;

  v = vec_lde(0, &x);
  p = vec_lvsl(0, &x);
  v = vec_perm(v, v, p);
  v = vec_splat(v, 0);
  return v;
}

/* Function:  esl_vmx_set_u8()
 * Synopsis:  Fills byte vector with x.
 *
 * Purpose:   Sets all elements in the vector <unsigned char> to x.
 */
static inline vector unsigned char
esl_vmx_set_u8(unsigned char x)
{
  vector unsigned char v;
  vector unsigned char p;

  v = vec_lde(0, &x);
  p = vec_lvsl(0, &x);
  v = vec_perm(v, v, p);
  v = vec_splat(v, 0);
  return v;
}

/* Function:  esl_vmx_hsum_float()
 * Synopsis:  Returns sum of all floats.
 *
 * Purpose:   Resturns the sum of all elements in the vector <float>.
 */
static inline float
esl_vmx_hsum_float(vector float v)
{
  float f;

  v = vec_add(v, vec_sld(v, v, 4));
  v = vec_add(v, vec_sld(v, v, 8));
  vec_ste(v, 0, &f);

  return f;
}

/* Function:  esl_vmx_hsum_s16()
 * Synopsis:  Returns sum of all shorts.
 *
 * Purpose:   Resturns the sum of all elements in the vector <signed short>.
 */
static inline signed short
esl_vmx_hsum_s16(vector signed short v)
{
  signed short s;

  v = vec_add(v, vec_sld(v, v, 2));
  v = vec_add(v, vec_sld(v, v, 4));
  v = vec_add(v, vec_sld(v, v, 8));
  vec_ste(v, 0, &s);

  return s;
}

/* Function:  esl_vmx_hmax_float()
 * Synopsis:  Returns max of all floats.
 *
 * Purpose:   Resturns the maximum element in the vector <float>.
 */
static inline float
esl_vmx_hmax_float(vector float v)
{
  float f;

  v = vec_max(v, vec_sld(v, v, 4));
  v = vec_max(v, vec_sld(v, v, 8));
  vec_ste(v, 0, &f);

  return f;
}

/* Function:  esl_vmx_hmax_s16()
 * Synopsis:  Returns max of all shorts.
 *
 * Purpose:   Resturns the maximum element in the vector <signed short>.
 */
static inline signed short
esl_vmx_hmax_s16(vector signed short v)
{
  signed short s;

  v = vec_max(v, vec_sld(v, v, 2));
  v = vec_max(v, vec_sld(v, v, 4));
  v = vec_max(v, vec_sld(v, v, 8));
  vec_ste(v, 0, &s);

  return s;
}

/* Function:  esl_vmx_hmax_u8()
 * Synopsis:  Returns max of all bytes.
 *
 * Purpose:   Resturns the maximum element in the vector <unsigned char>.
 */
static inline unsigned char
esl_vmx_hmax_u8(vector unsigned char v)
{
  unsigned char s;

  v = vec_max(v, vec_sld(v, v, 1));
  v = vec_max(v, vec_sld(v, v, 2));
  v = vec_max(v, vec_sld(v, v, 4));
  v = vec_max(v, vec_sld(v, v, 8));
  vec_ste(v, 0, &s);

  return s;
}

#endif /*eslVMX_INCLUDED*/
#endif /*HAVE_VMX*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_vmx.h ***/


/*** Start of inlined file: esl_weibull.h ***/
#ifndef eslWEIBULL_INCLUDED
#define eslWEIBULL_INCLUDED

double esl_wei_pdf    (double x, double mu, double lambda, double tau);
double esl_wei_logpdf (double x, double mu, double lambda, double tau);
double esl_wei_cdf    (double x, double mu, double lambda, double tau);
double esl_wei_logcdf (double x, double mu, double lambda, double tau);
double esl_wei_surv   (double x, double mu, double lambda, double tau);
double esl_wei_logsurv(double x, double mu, double lambda, double tau);
double esl_wei_invcdf (double p, double mu, double lambda, double tau);

double esl_wei_generic_pdf   (double x, void *params);
double esl_wei_generic_cdf   (double x, void *params);
double esl_wei_generic_surv  (double x, void *params);
double esl_wei_generic_invcdf(double p, void *params);

int esl_wei_Plot(FILE *fp, double mu, double lambda, double tau,
			double (*func)(double x, double mu, double lambda, double tau),
			double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
double esl_wei_Sample(ESL_RANDOMNESS *r, double mu, double lambda, double tau);
#endif

#ifdef eslAUGMENT_MINIMIZER
int esl_wei_FitComplete(double *x, int n, double *ret_mu,
			       double *ret_lambda, double *ret_tau);
#ifdef eslAUGMENT_HISTOGRAM
int esl_wei_FitCompleteBinned(ESL_HISTOGRAM *h, double *ret_mu,
				     double *ret_lambda, double *ret_tau);
#endif /*eslAUGMENT_HISTOGRAM*/
#endif /*eslAUGMENT_MINIMIZER*/

#endif /*eslWEIBULL_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_weibull.h ***/


/*** Start of inlined file: esl_workqueue.h ***/
#ifndef eslWORKQUEUE_INCLUDED
#define eslWORKQUEUE_INCLUDED

typedef struct {
  pthread_mutex_t  queueMutex;          /* mutex for queue serialization                           */
  pthread_cond_t   readerQueueCond;     /* condition variable used to wake up the producer         */
  pthread_cond_t   workerQueueCond;     /* condition variable used to wake up the consumers        */

  void           **readerQueue;         /* list of objects the the workers have completed          */
  int              readerQueueCnt;      /* number of objects in the queue                          */
  int              readerQueueHead;     /* first object in the queue                               */

  void           **workerQueue;         /* list of objects ready to be processed by worker threads */
  int              workerQueueCnt;      /* number of objects in the queue                          */
  int              workerQueueHead;     /* first object in the queue                               */

  int              queueSize;           /* max number of items a queue will hold                   */
  int              pendingWorkers;      /* number of consumers waiting for work                    */
} ESL_WORK_QUEUE;

ESL_WORK_QUEUE *esl_workqueue_Create(int size);
void            esl_workqueue_Destroy(ESL_WORK_QUEUE *queue);

int esl_workqueue_Init    (ESL_WORK_QUEUE *queue, void *ptr);
int esl_workqueue_Complete(ESL_WORK_QUEUE *queue);
int esl_workqueue_Reset   (ESL_WORK_QUEUE *queue);

int esl_workqueue_Remove(ESL_WORK_QUEUE *queue, void **obj);

int esl_workqueue_ReaderUpdate(ESL_WORK_QUEUE *queue, void *in, void **out);
int esl_workqueue_WorkerUpdate(ESL_WORK_QUEUE *queue, void *in, void **out);

int esl_workqueue_Dump(ESL_WORK_QUEUE *queue);

#endif /*eslWORKQUEUE_INCLUDED*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_workqueue.h ***/


/*** Start of inlined file: esl_wuss.h ***/
#ifndef eslWUSS_INCLUDED
#define eslWUSS_INCLUDED

int esl_wuss2ct(char *ss, int len, int *ct);
int esl_ct2wuss(int *ct, int n, char *ss);
int esl_ct2simplewuss(int *ct, int n, char *ss);
int esl_wuss2kh(char *ss, char *kh);
int esl_kh2wuss(char *kh, char *ss);
int esl_wuss_full(char *oldss, char *newss);
int esl_wuss_nopseudo(char *ss1, char *ss2);
int esl_wuss_reverse(char *ss, char *new);

#endif /*eslWUSS_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version 0.43; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute
 * Other copyrights also apply. See the LICENSE file for a full list.
 *
 * Easel is open source software, distributed under the BSD license. See
 * the LICENSE file for more details.
 *****************************************************************/

/*** End of inlined file: esl_wuss.h ***/


/*** Start of inlined file: interface_gsl.h ***/
#ifdef HAVE_LIBGSL
/* interface_gsl.h
 * Easel's interfaces to the GNU Scientific Library
 *
 * SRE, Tue Jul 13 15:36:48 2004
 * SVN $Id$
 * SVN $URL$
 */
#ifndef eslINTERFACE_GSL_INCLUDED
#define eslINTERFACE_GSL_INCLUDED

#include <stdlib.h>
#include <easel/easel.h>
#include <easel/dmatrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_eigen.h>

int esl_GSL_MatrixInversion(ESL_DMATRIX *A, ESL_DMATRIX **ret_Ai);

#endif /*eslINTERFACE_GSL_INCLUDED*/
#endif /*HAVE_LIBGSL*/

/*** End of inlined file: interface_gsl.h ***/


/*** Start of inlined file: interface_lapack.h ***/
#ifdef HAVE_LIBLAPACK
/* interface_lapack.h
 *
 * SRE, Tue Jul 13 15:11:51 2004 [St. Louis]
 * SVN $Id$
 * SVN $URL$
 */
#ifndef eslINTERFACE_LAPACK_INCLUDED
#define eslINTERFACE_LAPACK_INCLUDED

/* This is the C interface to the Fortran77 dgeev routine,
 * provided by the LAPACK library:
 */
void  dgeev_(char *jobvl, char *jobvr, int *n, double *a,
					int *lda, double *wr, double *wi, double *vl,
					int *ldvl, double *vr, int *ldvr,
					double *work, int *lwork, int *info);

/* and this is our C interface to the lapack call:
 */
int esl_lapack_dgeev(ESL_DMATRIX *A, double **ret_Er, double **ret_Ei, ESL_DMATRIX **ret_VL, ESL_DMATRIX **ret_VR);

#endif /*eslINTERFACE_LAPACK_INCLUDED*/
#endif /*HAVE_LIBLAPACK*/

/*** End of inlined file: interface_lapack.h ***/

#endif /*eslEASEL_INCLUDED*/
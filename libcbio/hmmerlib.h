/* The all-encompassing include file for HMMER.
 * All-encompassing because there's a lot of crossdependency.
 * There's some opportunity for modularity, but not a lot.
 *
 *    1. P7_HMM:         a core model.
 *    2. P7_PROFILE:     a scoring profile, and its implicit model.
 *    3. P7_BG:          a null (background) model.
 *    4. P7_TRACE:       a traceback path (alignment of seq to profile).
 *    5. P7_HMMFILE:     an HMM save file or database, open for reading.
 *    6. P7_GMX:         a "generic" dynamic programming matrix
 *    7. P7_PRIOR:       mixture Dirichlet prior for profile HMMs
 *    8. P7_SPENSEMBLE:  segment pair ensembles for domain locations
 *    9. P7_ALIDISPLAY:  an alignment formatted for printing
 *   10. P7_DOMAINDEF:   reusably managing workflow in annotating domains
 *   11. P7_TOPHITS:     ranking lists of top-scoring hits
 *   12. P7_SCOREDATA:     data used in diagonal recovery and extension
 *   13. P7_HMM_WINDOW:  data used to track lists of sequence windows
 *   14. Inclusion of the architecture-specific optimized implementation.
 *   16. P7_PIPELINE:    H3's accelerated seq/profile comparison pipeline
 *   17. P7_BUILDER:     configuration options for new HMM construction.
 *   18. Declaration of functions in HMMER's exposed API.
 *   19. Copyright and license information.
 *
 * Also, see impl_{sse,vmx}/impl_{sse,vmx}.h for additional API
 * specific to the acceleration layer; in particular, the P7_OPROFILE
 * structure for an optimized profile.
 */
#ifndef P7_HMMERH_INCLUDED
#define P7_HMMERH_INCLUDED

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HMMER_THREADS
#include <pthread.h>
#endif

#include "easellib.h"
#include "p7_config.h"

/* Search modes. */
#define p7_NO_MODE   0
#define p7_LOCAL     1		/* multihit local:  "fs" mode   */
#define p7_GLOCAL    2		/* multihit glocal: "ls" mode   */
#define p7_UNILOCAL  3		/* unihit local: "sw" mode      */
#define p7_UNIGLOCAL 4		/* unihit glocal: "s" mode      */

#define p7_IsLocal(mode)  (mode == p7_LOCAL || mode == p7_UNILOCAL)
#define p7_IsMulti(mode)  (mode == p7_LOCAL || mode == p7_GLOCAL)

#define p7_NEVPARAM 6	/* number of statistical parameters stored in models                      */
#define p7_NCUTOFFS 6	/* number of Pfam score cutoffs stored in models                          */
#define p7_NOFFSETS 3	/* number of disk offsets stored in models for hmmscan's fast model input */
enum p7_evparams_e {    p7_MMU  = 0, p7_MLAMBDA = 1,     p7_VMU = 2,  p7_VLAMBDA = 3, p7_FTAU = 4, p7_FLAMBDA = 5 };
enum p7_cutoffs_e  {     p7_GA1 = 0,     p7_GA2 = 1,     p7_TC1 = 2,      p7_TC2 = 3,  p7_NC1 = 4,     p7_NC2 = 5 };
enum p7_offsets_e  { p7_MOFFSET = 0, p7_FOFFSET = 1, p7_POFFSET = 2 };

#define p7_EVPARAM_UNSET -99999.0f  /* if evparam[0] is unset, then all unset                         */
#define p7_CUTOFF_UNSET  -99999.0f  /* if cutoff[XX1] is unset, then cutoff[XX2] unset, XX={GA,TC,NC} */
#define p7_COMPO_UNSET   -1.0f      /* if compo[0] is unset, then all unset                           */

/* Option flags when creating multiple alignments with p7_tracealign_*() */
#define p7_DEFAULT             0
#define p7_DIGITIZE            (1<<0)
#define p7_ALL_CONSENSUS_COLS  (1<<1)
#define p7_TRIM                (1<<2)

/* Option flags when creating faux traces with p7_trace_FauxFromMSA() */
#define p7_MSA_COORDS	       (1<<0) /* default: i = unaligned seq residue coords     */

/* Which strand(s) should be searched */
enum p7_strands_e {    p7_STRAND_TOPONLY  = 0, p7_STRAND_BOTTOMONLY = 1,  p7_STRAND_BOTH = 2};

/*****************************************************************
 * 1. P7_HMM: a core model.
 *****************************************************************/

/* Bit flags used in <hmm->flags>: optional annotation in an HMM
 *
 * Flags marked with ! may not be changed nor used for other meanings,
 * because they're codes used by HMMER2 (and earlier) that must be
 * preserved for reverse compatibility with old HMMER files.
 *
 * Why use flags? (So I don't ask this question of myself again:)
 *   1. The way we allocate an HMM, we need to know if we're allocating
 *      M-width annotation fields (RF, CS, CA, MAP) before we read the
 *      annotation from the file.
 *   2. Historically, H2 used flags, so we still need to read H2 flags
 *      for backwards compatibility; so we may as well keep using them.
 */
#define p7H_HASBITS (1<<0)    /* obsolete (was: model has log-odds scores)       !*/
#define p7H_DESC    (1<<1)    /* description exists (legacy; xref SRE:J5/114)    !*/
#define p7H_RF      (1<<2)    /* #RF annotation available                        !*/
#define p7H_CS      (1<<3)    /* #CS annotation available                        !*/
#define p7H_XRAY    (1<<4)    /* obsolete (was: structural data available)       !*/
#define p7H_HASPROB (1<<5)    /* obsolete (was: model in probability form)       !*/
#define p7H_HASDNA  (1<<6)    /* obsolete (was: protein HMM->DNA seq params set) !*/
#define p7H_STATS   (1<<7)    /* model has E-value statistics calibrated         !*/
#define p7H_MAP     (1<<8)    /* alignment map is available                      !*/
#define p7H_ACC     (1<<9)    /* accession is available (legacy; xref SRE:J5/114)!*/
#define p7H_GA      (1<<10)   /* gathering thresholds available                  !*/
#define p7H_TC      (1<<11)   /* trusted cutoffs available                       !*/
#define p7H_NC      (1<<12)   /* noise cutoffs available                         !*/
#define p7H_CA      (1<<13)   /* surface accessibilities available               !*/
#define p7H_COMPO   (1<<14)   /* model-specific residue composition available     */
#define p7H_CHKSUM  (1<<15)   /* model has an alignment checksum                  */
#define p7H_CONS    (1<<16)   /* consensus residue line available                 */
#define p7H_MMASK   (1<<17)   /* #MM annotation available                        !*/

/* Indices of Plan7 main model state transitions, hmm->t[k][] */
enum p7h_transitions_e {
  p7H_MM = 0,
  p7H_MI = 1,
  p7H_MD = 2,
  p7H_IM = 3,
  p7H_II = 4,
  p7H_DM = 5,
  p7H_DD = 6
};
#define p7H_NTRANSITIONS 7

/* How the hmm->t[k] vector is interpreted as separate probability vectors. */
#define P7H_TMAT(hmm, k) ((hmm)->t[k])
#define P7H_TINS(hmm, k) ((hmm)->t[k]+3)
#define P7H_TDEL(hmm, k) ((hmm)->t[k]+5)
#define p7H_NTMAT 3
#define p7H_NTDEL 2
#define p7H_NTINS 2

/* Some notes:
 *   0. The model might be either in counts or probability form.
 *   1. t[0] is special: t[0][TMM,TMI,TMD] are the begin->M_1,I_0,D_1 entry probabilities,
 *      t[0][TIM,TII] are the I_0 transitions, and delete state 0 doesn't
 *      exist. Therefore D[0] transitions and mat[0] emissions are unused.
 *      To simplify some normalization code, we adopt a convention that these are set
 *      to valid probability distributions: 1.0 for t[0][TDM] and mat[0][0],
 *      and 0 for the rest.
 *   2. t[M] is also special: TMD and TDD are 0 because there is no next delete state;
 *      TDM is therefore 1.0 by definition. TMM and TDM are interpreted as the
 *      M->E and D->E end transitions. t[M][TDM] must be 1.0, therefore.
 */
typedef struct p7_hmm_s {
  /*::cexcerpt::plan7_core::begin::*/
  int     M;                    /* length of the model (# nodes)                           */
  float **t;                    /* transition prob's. t[(0),1..M][0..p7H_NTRANSITIONS-1]   */
  float **mat;                  /* match emissions.  mat[1..M][0..K-1]                     */
  float **ins;                  /* insert emissions. ins[1..M][0..K-1]                     */
  /*::cexcerpt::plan7_core::end::*/

  /* Annotation. Everything but <name> is optional. Flags are set when
   * optional values are set. All the char *'s are proper nul-terminated
   * strings, not just arrays. (hmm->map is an int array).
   */
  char    *name;                 /* name of the model                     (mandatory)      */ /* String, \0-terminated   */
  char    *acc;	                 /* accession number of model (Pfam)      (optional: NULL) */ /* String, \0-terminated   */
  char    *desc;                 /* brief (1-line) description of model   (optional: NULL) */ /* String, \0-terminated   */
  char    *rf;                   /* reference line from alignment 1..M    (p7H_RF)         */ /* String; 0=' ', M+1='\0' */
  char    *mm;                   /* model mask line from alignment 1..M   (p7H_MM)         */ /* String; 0=' ', M+1='\0' */
  char    *consensus;		         /* consensus residue line        1..M    (p7H_CONS)       */ /* String; 0=' ', M+1='\0' */
  char    *cs;                   /* consensus structure line      1..M    (p7H_CS)         */ /* String; 0=' ', M+1='\0' */
  char    *ca;	                 /* consensus accessibility line  1..M    (p7H_CA)         */ /* String; 0=' ', M+1='\0' */

  char    *comlog;               /* command line(s) that built model      (optional: NULL) */ /* String, \0-terminated   */
  int      nseq;	         /* number of training sequences          (optional: -1)   */
  float    eff_nseq;             /* effective number of seqs (<= nseq)    (optional: -1)   */
  int	   max_length;           /* upper bound length, all but 1e-7 prob (optional: -1)   */
  char    *ctime;	         /* creation date                         (optional: NULL) */
  int     *map;	                 /* map of alignment cols onto model 1..M (p7H_MAP)        */ /* Array; map[0]=0 */
  uint32_t checksum;             /* checksum of training sequences        (p7H_CHKSUM)     */
  float    evparam[p7_NEVPARAM]; /* E-value params                        (p7H_STATS)      */
  float    cutoff[p7_NCUTOFFS];  /* Pfam score cutoffs                    (p7H_{GA,TC,NC}) */
  float    compo[p7_MAXABET];    /* model bg residue comp                 (p7H_COMPO)      */

  off_t    offset;               /* HMM record offset on disk                              */
  const ESL_ALPHABET *abc;       /* ptr to alphabet info (hmm->abc->K is alphabet size)    */
  int      flags;                /* status flags                                           */
} P7_HMM;

/*****************************************************************
 * 2. P7_PROFILE: a scoring profile, and its implicit model.
 *****************************************************************/

/* Indices for special state types in the length model, gm->xsc[x][]
 */
enum p7p_xstates_e {
  p7P_E = 0,
  p7P_N = 1,
  p7P_J = 2,
  p7P_C = 3
};
#define p7P_NXSTATES 4

/* Indices for transitions from the length modeling scores gm->xsc[][x]
 */
enum p7p_xtransitions_e {
  p7P_LOOP = 0,
  p7P_MOVE = 1
};
#define p7P_NXTRANS 2

/* Indices for transition scores gm->tsc[k][] */
/* order is optimized for dynamic programming */
enum p7p_tsc_e {
  p7P_MM = 0,
  p7P_IM = 1,
  p7P_DM = 2,
  p7P_BM = 3,
  p7P_MD = 4,
  p7P_DD = 5,
  p7P_MI = 6,
  p7P_II = 7,
};
#define p7P_NTRANS 8

/* Indices for residue emission score vectors
 */
enum p7p_rsc_e {
  p7P_MSC = 0,
  p7P_ISC = 1
};
#define p7P_NR 2

/* Accessing transition, emission scores */
/* _BM is specially stored off-by-one: [k-1][p7P_BM] is score for entering at Mk */
#define p7P_TSC(gm, k, s) ((gm)->tsc[(k) * p7P_NTRANS + (s)])
#define p7P_MSC(gm, k, x) ((gm)->rsc[x][(k) * p7P_NR + p7P_MSC])
#define p7P_ISC(gm, k, x) ((gm)->rsc[x][(k) * p7P_NR + p7P_ISC])

typedef struct p7_profile_s {
  float  *tsc;          /* transitions  [0.1..M-1][0..p7P_NTRANS-1], hand-indexed  */
  float **rsc;          /* emissions [0..Kp-1][0.1..M][p7P_NR], hand-indexed       */
  float   xsc[p7P_NXSTATES][p7P_NXTRANS]; /* special transitions [NECJ][LOOP,MOVE] */

  int     mode;        	/* configured algorithm mode (e.g. p7_LOCAL)               */
  int     L;		/* current configured target seq length                    */
  int     allocM;	/* max # of nodes allocated in this structure              */
  int     M;		/* number of nodes in the model                            */
  int     max_length;	/* calculated upper bound on emitted seq length            */
  float   nj;		/* expected # of uses of J; precalculated from loop config */

  /* Info, most of which is a copy from parent HMM:                                       */
  char  *name;			/* unique name of model                                   */
  char  *acc;			/* unique accession of model, or NULL                     */
  char  *desc;                  /* brief (1-line) description of model, or NULL           */
  char  *rf;                    /* reference line from alignment 1..M; *rf=0 means unused */
  char  *mm;                    /* modelmask line           1..M; *ref=0: unused     */
  char  *cs;                    /* consensus structure line      1..M, *cs=0 means unused */
  char  *consensus;		/* consensus residues to display in alignments, 1..M      */
  float  evparam[p7_NEVPARAM]; 	/* parameters for determining E-values, or UNSET          */
  float  cutoff[p7_NCUTOFFS]; 	/* per-seq/per-domain bit score cutoffs, or UNSET         */
  float  compo[p7_MAXABET];	/* per-model HMM filter composition, or UNSET             */

  /* Disk offset information for hmmpfam's fast model retrieval                           */
  off_t  offs[p7_NOFFSETS];     /* p7_{MFP}OFFSET, or -1                                  */

  off_t  roff;                  /* record offset (start of record); -1 if none            */
  off_t  eoff;                  /* offset to last byte of record; -1 if unknown           */

  const ESL_ALPHABET *abc;	/* copy of pointer to appropriate alphabet                */
} P7_PROFILE;

/*****************************************************************
 * 3. P7_BG: a null (background) model.
 *****************************************************************/

/* This really contains three different things:
 *
 *   - the "null1" model, a one-state HMM consisting of background
 *     frequencies <f> and a parameter <p1> for a target-length
 *     dependent geometric;
 *
 *   - the "bias filter" <fhmm> a two-state HMM composed from null1's
 *     background <f> and the model's mean composition <compo>. This
 *     model is constructed dynamically, every time a new profile is
 *     considered;
 *
 *   - a single term <omega> that's needed by the "null2" model to set
 *     a balance between the null1 and null2 scoring terms.  The null2
 *     model is otherwise defined by construction, in p7_domaindef.c.
 *
 * Someday we might pull this apart into two or three separate
 * objects.
 */
typedef struct p7_bg_s {
  float   *f;		/* null1 background residue frequencies [0..K-1]: set at initialization    */
  float    p1;		/* null1's transition prob: p7_bg_SetLength() sets this from target seq L  */

  ESL_HMM *fhmm;	/* bias filter: p7_bg_SetFilter() sets this, from model's mean composition */

  float    omega;	/* the "prior" on null2/null3: set at initialization (one omega for both null types)  */

  const ESL_ALPHABET *abc;	/* reference to alphabet in use: set at initialization             */
} P7_BG;

/*****************************************************************
 * 4. P7_TRACE:  a traceback (alignment of seq to profile).
 *****************************************************************/

/* Traceback structure for alignment of a model to a sequence.
 *
 * A traceback only makes sense in a triplet (tr, gm, dsq), for a
 * given profile or HMM (with nodes 1..M) and a given digital sequence
 * (with positions 1..L).
 *
 * A traceback may be relative to a profile (usually) or to a core
 * model (as a special case in model construction; see build.c). You
 * can tell the difference by looking at the first statetype,
 * tr->st[0]; if it's a p7T_S, it's for a profile, and if it's p7T_B,
 * it's for a core model.
 *
 * A "profile" trace uniquely has S,N,C,T,J states and their
 * transitions; it also can have B->Mk and Mk->E internal entry/exit
 * transitions for local alignments. It may not contain X states.
 *
 * A "core" trace may contain I0, IM, and D1 states and their
 * transitions. A "core" trace can also have B->X->{MDI}k and
 * {MDI}k->X->E transitions as a special hack in a build procedure, to
 * deal with the case of a local alignment fragment implied by an
 * input alignment, which is "impossible" for a core model.
 * X "states" only appear in core traces, and only at these
 * entry/exit places; some code depends on this.
 *
 * A profile's N,C,J states emit on transition, not on state, so a
 * path of N emits 0 residues, NN emits 1 residue, NNN emits 2
 * residues, and so on. By convention, the trace always associates an
 * emission-on-transition with the trailing (destination) state, so
 * the first N, C, or J is stored in a trace as a nonemitter (i=0).
 *
 * A i coords in a traceback are usually 1..L with respect to an
 * unaligned digital target sequence, but in the special case of
 * traces faked from existing MSAs (as in hmmbuild), the coords may
 * be 1..alen relative to an MSA's columns.
 */

/* State types */
enum p7t_statetype_e {
  p7T_BOGUS =  0,
  p7T_M     =  1,
  p7T_D     =  2,
  p7T_I     =  3,
  p7T_S     =  4,
  p7T_N     =  5,
  p7T_B     =  6,
  p7T_E     =  7,
  p7T_C     =  8,
  p7T_T     =  9,
  p7T_J     = 10,
  p7T_X     = 11, 	/* missing data: used esp. for local entry/exits */
};
#define p7T_NSTATETYPES 12

typedef struct p7_trace_s {
  int    N;		/* length of traceback                       */
  int    nalloc;        /* allocated length of traceback             */
  char  *st;		/* state type code                   [0..N-1]*/
  int   *k;		/* node index; 1..M if M,D,I; else 0 [0..N-1]*/
  int   *i;		/* pos emitted in dsq, 1..L; else 0  [0..N-1]*/
  float *pp;		/* posterior prob of x_i; else 0     [0..N-1]*/
  int    M;		/* model length M (maximum k)                */
  int    L;		/* sequence length L (maximum i)             */

  /* The following section is data generated by "indexing" a trace's domains */
  int   ndom;		/* number of domains in trace (= # of B or E states) */
  int  *tfrom,   *tto;	/* locations of B/E states in trace (0..tr->N-1)     */
  int  *sqfrom,  *sqto;	/* first/last M-emitted residue on sequence (1..L)   */
  int  *hmmfrom, *hmmto;/* first/last M state on model (1..M)                */
  int   ndomalloc;	/* current allocated size of these stacks            */

} P7_TRACE;

/*****************************************************************
 * 5. P7_HMMFILE:  an HMM save file or database, open for reading.
 *****************************************************************/

/* These tags need to be in temporal order, so we can do tests
 * like "if (format >= p7_HMMFILE_3b) ..."
 */
enum p7_hmmfile_formats_e {
  p7_HMMFILE_20 = 0,
  p7_HMMFILE_3a = 1,
  p7_HMMFILE_3b = 2,
  p7_HMMFILE_3c = 3,
  p7_HMMFILE_3d = 4,
  p7_HMMFILE_3e = 5,
  p7_HMMFILE_3f = 6,
};

typedef struct p7_hmmfile_s {
  FILE         *f;		 /* pointer to stream for reading models                 */
  char         *fname;	         /* (fully qualified) name of the HMM file; [STDIN] if - */
  ESL_SSI      *ssi;		 /* open SSI index for model file <f>; NULL if none.     */

  int           do_gzip;	/* TRUE if f is "gzip -dc |" (will pclose(f))           */
  int           do_stdin;       /* TRUE if f is stdin (won't close f)                   */
  int           newly_opened;	/* TRUE if we just opened the stream (and parsed magic) */
  int           is_pressed;	/* TRUE if a pressed HMM database file (Pfam or equiv)  */

  int            format;	/* HMM file format code */
  int           (*parser)(struct p7_hmmfile_s *, ESL_ALPHABET **, P7_HMM **);
  ESL_FILEPARSER *efp;

  /* If <is_pressed>, we can read optimized profiles directly, via:  */
  FILE         *ffp;		/* MSV part of the optimized profile */
  FILE         *pfp;		/* rest of the optimized profile     */

#ifdef HMMER_THREADS
  int              syncRead;
  pthread_mutex_t  readMutex;
#endif

  char          errbuf[eslERRBUFSIZE];
} P7_HMMFILE;

/* note on <fname>, above:
 * this is the actual name of the HMM file being read.
 *
 * The way p7_hmmfile_Open() works, it will preferentially look for
 * hmmpress'ed binary files. If you open "foo", it will first try to
 * open "foo.h3m" and <fname> will be "foo.h3m". "foo" does not even
 * have to exist. If a parsing error occurs, you want <fname> to
 * be "foo.h3m", so error messages report blame correctly.
 * In the special case of reading from stdin, <fname> is "[STDIN]".
 */

/*****************************************************************
 * 6. P7_GMX: a "generic" dynamic programming matrix
 *****************************************************************/

enum p7g_scells_e {
  p7G_M = 0,
  p7G_I = 1,
  p7G_D = 2,
};
#define p7G_NSCELLS 3

enum p7g_xcells_e {
  p7G_E  = 0,
  p7G_N  = 1,
  p7G_J  = 2,
  p7G_B  = 3,
  p7G_C  = 4
};
#define p7G_NXCELLS 5

typedef struct p7_gmx_s {
  int  M;		/* actual model dimension (model 1..M)    */
  int  L;		/* actual sequence dimension (seq 1..L)   */

  int      allocR;      /* current allocated # of rows : L+1 <= validR <= allocR                */
  int      validR;	/* # of rows actually pointing at DP memory                             */
  int      allocW;	/* current set row width :  M+1 <= allocW                               */
  uint64_t ncells;	/* total # of allocated cells in 2D matrix : ncells >= (validR)(allocW) */

  float **dp;           /* logically [0.1..L][0.1..M][0..p7G_NSCELLS-1]; indexed [i][k*p7G_NSCELLS+s] */
  float  *xmx;          /* logically [0.1..L][0..p7G_NXCELLS-1]; indexed [i*p7G_NXCELLS+s]            */

  float  *dp_mem;
} P7_GMX;

/* Macros below implement indexing idioms for generic DP routines.
 * They require the following setup, for profile <gm> and matrix <gx>:
 *   float const *tsc = gm->tsc;
 *   float      **dp  = gx->dp;
 *   float       *xmx = gx->xmx;
 * and for each row i (target residue x_i in digital seq <dsq>):
 *   float const *rsc = gm->rsc[dsq[i]];
 */
#define MMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_M])
#define IMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_I])
#define DMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_D])
#define XMX(i,s) (xmx[(i) * p7G_NXCELLS + (s)])

#define TSC(s,k) (tsc[(k) * p7P_NTRANS + (s)])
#define MSC(k)   (rsc[(k) * p7P_NR     + p7P_MSC])
#define ISC(k)   (rsc[(k) * p7P_NR     + p7P_ISC])

/* Flags that control P7_GMX debugging dumps */
#define p7_HIDE_SPECIALS (1<<0)
#define p7_SHOW_LOG      (1<<1)

/*****************************************************************
 * 7. P7_PRIOR: mixture Dirichlet prior for profile HMMs
 *****************************************************************/

typedef struct p7_prior_s {
  ESL_MIXDCHLET *tm;		/*  match transitions */
  ESL_MIXDCHLET *ti;		/* insert transitions */
  ESL_MIXDCHLET *td;		/* delete transitions */
  ESL_MIXDCHLET *em;		/*  match emissions   */
  ESL_MIXDCHLET *ei;		/* insert emissions   */
} P7_PRIOR;

/*****************************************************************
 * 8. P7_SPENSEMBLE: segment pair ensembles for domain locations
 *****************************************************************/

/* struct p7_spcoord_s:
 *    a coord quad defining a segment pair.
 */
struct p7_spcoord_s {
  int idx; 	/* backreference index: which trace a seg came from, or which cluster a domain came from */
  int i, j;	/* start,end in a target sequence (1..L)  */
  int k, m;     /* start,end in a query model (1..M)      */
  float prob;	/* posterior probability of segment       */
};

/* Structure: P7_SPENSEMBLE
 *
 * Collection and clustering of an ensemble of sampled segment pairs,
 * in order to define domain locations using their posterior
 * probability distribution (as opposed to Viterbi MAP tracebacks).
 */
typedef struct p7_spensemble_s {
  /* Section 1: a collected ensemble of segment pairs                                       */
  int                  nsamples;    /* number of sampled traces                             */
  struct p7_spcoord_s *sp;	    /* array of sampled seg pairs; [0..n-1]                 */
  int                  nalloc;	    /* allocated size of <sp>                               */
  int                  n;	    /* number of seg pairs in <sp>                          */

  /* Section 2: then the ensemble is clustered by single-linkage clustering                 */
  int *workspace;                   /* temp space for Easel SLC algorithm: 2*n              */
  int *assignment;                  /* each seg pair's cluster index: [0..n-1] = (0..nc-1)  */
  int  nc;	                    /* number of different clusters                         */

  /* Section 3: then endpoint distribution is examined within each large cluster            */
  int *epc;	                    /* array counting frequency of each endpoint            */
  int  epc_alloc;	            /* allocated width of <epc>                             */

  /* Section 4: finally each large cluster is resolved into domain coords                   */
  struct p7_spcoord_s *sigc;	    /* array of coords for each domain, [0..nsigc-1]        */
  int                  nsigc;	    /* number of "significant" clusters, domains            */
  int                  nsigc_alloc; /* current allocated max for nsigc                      */
} P7_SPENSEMBLE;

/*****************************************************************
 * 9. P7_ALIDISPLAY: an alignment formatted for printing
 *****************************************************************/

/* Structure: P7_ALIDISPLAY
 *
 * Alignment of a sequence domain to an HMM, formatted for printing.
 *
 * For an alignment of L residues and names C chars long, requires
 * 6L + 2C + 30 bytes; for typical case of L=100,C=10, that's
 * <0.7 Kb.
 */
typedef struct p7_alidisplay_s {
  char *rfline;                 /* reference coord info; or NULL        */
  char *mmline;                 /* modelmask coord info; or NULL        */
  char *csline;                 /* consensus structure info; or NULL    */
  char *model;                  /* aligned query consensus sequence     */
  char *mline;                  /* "identities", conservation +'s, etc. */
  char *aseq;                   /* aligned target sequence              */
  char *ntseq;                  /* nucleotide target sequence if nhmmscant */
  char *ppline;			        /* posterior prob annotation; or NULL   */
  int   N;			            /* length of strings                    */

  char *hmmname;		/* name of HMM                          */
  char *hmmacc;			/* accession of HMM; or [0]='\0'        */
  char *hmmdesc;		/* description of HMM; or [0]='\0'      */
  int   hmmfrom;		/* start position on HMM (1..M, or -1)  */
  int   hmmto;			/* end position on HMM (1..M, or -1)    */
  int   M;			/* length of model                      */

  char *sqname;			/* name of target sequence              */
  char *sqacc;			/* accession of target seq; or [0]='\0' */
  char *sqdesc;			/* description of targ seq; or [0]='\0' */
  long  sqfrom;			/* start position on sequence (1..L)    */
  long  sqto;		    /* end position on sequence   (1..L)    */
  long  L;			/* length of sequence                   */

  int   memsize;                /* size of allocated block of memory    */
  char *mem;			/* memory used for the char data above  */
} P7_ALIDISPLAY;

/*****************************************************************
 * 10. P7_DOMAINDEF: reusably managing workflow in defining domains
 *****************************************************************/

typedef struct p7_dom_s {
  int            ienv, jenv;
  int            iali, jali;
  int            iorf, jorf; /*Used in translated search to capture the range in the DNA sequence of the ORF containing the match to a protein query */
  float          envsc;  	/* Forward score in envelope ienv..jenv; NATS; without null2 correction       */
  float          domcorrection;	/* null2 score when calculating a per-domain score; NATS                      */
  float          dombias;	/* FLogsum(0, log(bg->omega) + domcorrection): null2 score contribution; NATS */
  float          oasc;		/* optimal accuracy score (units: expected # residues correctly aligned)      */
  float          bitscore;	/* overall score in BITS, null corrected, if this were the only domain in seq */
  double         lnP;	        /* log(P-value) of the bitscore                                               */
  int            is_reported;	/* TRUE if domain meets reporting thresholds                                  */
  int            is_included;	/* TRUE if domain meets inclusion thresholds                                  */
  float         *scores_per_pos; /* score in BITS that each position in the alignment contributes to an overall viterbi score */
  P7_ALIDISPLAY *ad;
} P7_DOMAIN;

/* Structure: P7_DOMAINDEF
 *
 * This is a container for all the necessary information for domain
 * definition procedures in <p7_domaindef.c>, including a bunch of
 * heuristic thresholds. The structure is reusable to minimize the
 * number of allocation/free cycles that need to be done when
 * processing a large number of sequences. You create the structure
 * with <p7_domaindef_Create()>; after you're done with defining
 * domains on a sequence, you call <p7_domaindef_Reuse()> before using
 * it on the next sequence; and when you're completely done, you free
 * it with <p7_domaindef_Destroy()>. All memory management is handled
 * internally; you don't need to reallocate anything yourself.
 */
typedef struct p7_domaindef_s {
  /* for posteriors of being in a domain, B, E */
  float *mocc;			/* mocc[i=1..L] = prob that i is emitted by core model (is in a domain)       */
  float *btot; 			/* btot[i=1..L] = cumulative expected times that domain starts at or before i */
  float *etot;			/* etot[i=1..L] = cumulative expected times that domain ends at or before i   */
  int    L;
  int    Lalloc;

  /* the ad hoc null2 model: 1..L nat scores for each residue, log f'(x_i) / f(x_i) */
  float *n2sc;

  /* rng and reusable memory for stochastic tracebacks */
  ESL_RANDOMNESS *r;		/* random number generator                                 */
  int             do_reseeding;	/* TRUE to reset the RNG, make results reproducible        */
  P7_SPENSEMBLE  *sp;		/* an ensemble of sampled segment pairs (domain endpoints) */
  P7_TRACE       *tr;		/* reusable space for a trace of a domain                  */
  P7_TRACE       *gtr;		/* reusable space for a traceback of the entire target seq */

  /* Heuristic thresholds that control the region definition process */
  /* "rt" = "region threshold", for lack of better term  */
  float  rt1;   	/* controls when regions are called. mocc[i] post prob >= dt1 : triggers a region around i */
  float  rt2;		/* controls extent of regions. regions extended until mocc[i]-{b,e}occ[i] < dt2            */
  float  rt3;		/* controls when regions are flagged for split: if expected # of E preceding B is >= dt3   */

  /* Heuristic thresholds that control the stochastic traceback/clustering process */
  int    nsamples;	/* collect ensemble of this many stochastic traces */
  float  min_overlap;	/* 0.8 means >= 80% overlap of (smaller/larger) segment to link, both in seq and hmm            */
  int    of_smaller;	/* see above; TRUE means overlap denom is calc'ed wrt smaller segment; FALSE means larger       */
  int    max_diagdiff;	/* 4 means either start or endpoints of two segments must be within <=4 diagonals of each other */
  float  min_posterior;	/* 0.25 means a cluster must have >= 25% posterior prob in the sample to be reported            */
  float  min_endpointp;	/* 0.02 means choose widest endpoint with post prob of at least 2%                              */

  /* storage of the results; domain locations, scores, alignments          */
  P7_DOMAIN *dcl;
  int        ndom;	 /* number of domains defined, in the end.         */
  int        nalloc;     /* number of domain structures allocated in <dcl> */

  /* Additional results storage */
  float  nexpected;     /* posterior expected number of domains in the sequence (from posterior arrays) */
  int    nregions;	/* number of regions evaluated */
  int    nclustered;	/* number of regions evaluated by clustering ensemble of tracebacks */
  int    noverlaps;	/* number of envelopes defined in ensemble clustering that overlap w/ prev envelope */
  int    nenvelopes;	/* number of envelopes handed over for domain definition, null2, alignment, and scoring. */

} P7_DOMAINDEF;

/*****************************************************************
 * 11. P7_TOPHITS: ranking lists of top-scoring hits
 *****************************************************************/

#define p7_HITFLAGS_DEFAULT 0
#define p7_IS_INCLUDED      (1<<0)
#define p7_IS_REPORTED      (1<<1)
#define p7_IS_NEW           (1<<2)
#define p7_IS_DROPPED       (1<<3)
#define p7_IS_DUPLICATE     (1<<4)

/* Structure: P7_HIT
 *
 * Info about a high-scoring database hit, kept so we can output a
 * sorted list of high hits at the end.
 *
 * sqfrom and sqto are the coordinates that will be shown
 * in the results, not coords in arrays... therefore, reverse
 * complements have sqfrom > sqto
 */
typedef struct p7_hit_s {
  char   *name;			/* name of the target               (mandatory)           */
  char   *acc;			/* accession of the target          (optional; else NULL) */
  char   *desc;			/* description of the target        (optional; else NULL) */
  int    window_length;         /* for later use in e-value computation, when splitting long sequences */
  double sortkey;		/* number to sort by; big is better                       */

  float  score;			/* bit score of the sequence (all domains, w/ correction) */
  float  pre_score;		/* bit score of sequence before null2 correction          */
  float  sum_score;		/* bit score reconstructed from sum of domain envelopes   */

  double lnP;		        /* log(P-value) of the score               */
  double pre_lnP;		/* log(P-value) of the pre_score           */
  double sum_lnP;		/* log(P-value) of the sum_score           */

  float  nexpected;     /* posterior expected number of domains in the sequence (from posterior arrays) */
  int    nregions;	/* number of regions evaluated */
  int    nclustered;	/* number of regions evaluated by clustering ensemble of tracebacks */
  int    noverlaps;	/* number of envelopes defined in ensemble clustering that overlap w/ prev envelope */
  int    nenvelopes;	/* number of envelopes handed over for domain definition, null2, alignment, and scoring. */
  int    ndom;		/* total # of domains identified in this seq   */

  uint32_t flags;      	/* p7_IS_REPORTED | p7_IS_INCLUDED | p7_IS_NEW | p7_IS_DROPPED */
  int      nreported;	/* # of domains satisfying reporting thresholding  */
  int      nincluded;	/* # of domains satisfying inclusion thresholding */
  int      best_domain;	/* index of best-scoring domain in dcl */

  int64_t  seqidx;          /*unique identifier to track the database sequence from which this hit came*/
  int64_t  subseq_start; /*used to track which subsequence of a full_length target this hit came from, for purposes of removing duplicates */

  P7_DOMAIN *dcl;	/* domain coordinate list and alignment display */
  esl_pos_t  offset;	/* used in socket communications, in serialized communication: offset of P7_DOMAIN msg for this P7_HIT */
} P7_HIT;

/* Structure: P7_TOPHITS
 * merging when we prepare to output results. "hit" list is NULL and
 * unavailable until after we do a sort.
 */
typedef struct p7_tophits_s {
  P7_HIT **hit;         /* sorted pointer array                     */
  P7_HIT  *unsrt;	/* unsorted data storage                    */
  uint64_t Nalloc;	/* current allocation size                  */
  uint64_t N;		/* number of hits in list now               */
  uint64_t nreported;	/* number of hits that are reportable       */
  uint64_t nincluded;	/* number of hits that are includable       */
  int      is_sorted_by_sortkey; /* TRUE when hits sorted by sortkey and th->hit valid for all N hits */
  int      is_sorted_by_seqidx; /* TRUE when hits sorted by seq_idx, position, and th->hit valid for all N hits */
} P7_TOPHITS;

/*****************************************************************
 * 12. P7_SCOREDATA: data used in diagonal recovery and extension
 *****************************************************************/

enum p7_scoredatatype_e {
  p7_sd_std  = 0,
  p7_sd_fm   = 1,
};

/* This contains a compact representation of 8-bit bias-shifted scores for use in
 * diagonal recovery (standard SSV) and extension (standard and FM-SSV),
 * along with MAXL-associated prefix- and suffix-lengths, and optimal extensions
 * for FM-SSV.
 */
typedef struct p7_scoredata_s {
  int         type;
  int         M;
  union {//implicit (M+1)*K matrix, where M = # states, and K = # characters in alphabet
	uint8_t  *ssv_scores;    // this 2D array is used in the default nhmmer pipeline
	float    *ssv_scores_f;  // this 2D array is used in the FM-index based pipeline
  };
  float      *prefix_lengths;
  float      *suffix_lengths;
  float      *fwd_scores;
  float     **fwd_transitions;
  float     **opt_ext_fwd; // Used only for FM-index based pipeline
  float     **opt_ext_rev; // Used only for FM-index based pipeline
} P7_SCOREDATA;

/*****************************************************************
 * 13. P7_HMM_WINDOW: data used to track lists of sequence windows
 *****************************************************************/

typedef struct p7_hmm_window_s {
  float      score;
  float      null_sc;
  int32_t    id;          //sequence id of the database sequence hit
  int64_t    n;           //position in database sequence at which the diagonal/window starts
  int64_t    fm_n;        //position in the concatenated fm-index sequence at which the diagonal starts
  int32_t    length;      // length of the diagonal/window
  int16_t    k;           //position of the model at which the diagonal ends
  int64_t    target_len;  //length of the target sequence
  int8_t     complementarity;
  int8_t     used_to_extend;
} P7_HMM_WINDOW;

typedef struct p7_hmm_window_list_s {
  P7_HMM_WINDOW *windows;
  int       count;
  int       size;
} P7_HMM_WINDOWLIST;

/*****************************************************************
 * 14. The optimized implementation.
 *****************************************************************/
#if   defined (p7_IMPL_SSE)

/*** Start of inlined file: impl_sse.h ***/
#ifndef P7_IMPL_SSE_INCLUDED
#define P7_IMPL_SSE_INCLUDED

//#include "p7_config.h"

//#include "esl_alphabet.h"
//#include "esl_random.h"

#include <xmmintrin.h>    /* SSE  */
#include <emmintrin.h>    /* SSE2 */
#ifdef __SSE3__
#include <pmmintrin.h>   /* DENORMAL_MODE */
#endif
//#include "hmmer.h"

/* In calculating Q, the number of vectors we need in a row, we have
 * to make sure there's at least 2, or a striped implementation fails.
 */
#define p7O_NQB(M)   ( ESL_MAX(2, ((((M)-1) / 16) + 1)))   /* 16 uchars  */
#define p7O_NQW(M)   ( ESL_MAX(2, ((((M)-1) / 8)  + 1)))   /*  8 words   */
#define p7O_NQF(M)   ( ESL_MAX(2, ((((M)-1) / 4)  + 1)))   /*  4 floats  */

#define p7O_EXTRA_SB 17    /* see ssvfilter.c for explanation */

/*****************************************************************
 * 1. P7_OPROFILE: an optimized score profile
 *****************************************************************/
/* The OPROFILE is striped [Farrar07] and interleaved, as is the DP matrix.
 * For example, the layout of a profile for an M=14 model (xref J2/46):
 *
 * rsc[x] : striped blocks of M emissions, starting with q=0
 *                1     11     1      1
 *             1593   2604   371x   482x
 *
 * tsc:  grouped in order of accession in DP for 7 transition scores;
 *       starting at q=0 for all but the three transitions to M, which
 *       are rotated by -1 and rightshifted. DD's follow separately,
 *       starting at q=0.
 *
 *        {     1      1     1     1     1     1     1 }
 *        {  1593   x482  x482  x482  1593  1593  1593 }
 *        { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
 *
 *        {    11      1     1     1    11    11    11 }
 *        {  2604   1593  1593  1593  2604  2604  2604 }
 *        { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
 *
 *        {    1      11    11    11    1     1     1  }
 *        {  371x   2604  2604  2604  371x  371x  371x }
 *        { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
 *
 *        {    1      1     1     1     1     1     1  }
 *        {  482x   371x  371x  371x  482x  482x  482x }
 *        { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
 *
 *        {     1    11    1     1  }
 *        {  1593  2604  371x  482x }
 *        { [TDD] [TDD] [TDD] [TDD] }
 *
 */

#define p7O_NXSTATES  4    /* special states stored: ENJC                       */
#define p7O_NXTRANS   2         /* special states all have 2 transitions: move, loop */
#define p7O_NTRANS    8    /* 7 core transitions + BMk entry                    */
enum p7o_xstates_e      { p7O_E    = 0, p7O_N    = 1,  p7O_J  = 2,  p7O_C  = 3 };
enum p7o_xtransitions_e { p7O_MOVE = 0, p7O_LOOP = 1 };
enum p7o_tsc_e          { p7O_BM   = 0, p7O_MM   = 1,  p7O_IM = 2,  p7O_DM = 3, p7O_MD   = 4, p7O_MI   = 5,  p7O_II = 6,  p7O_DD = 7 };

typedef struct p7_oprofile_s {
  /* MSVFilter uses scaled, biased uchars: 16x unsigned byte vectors                 */
  __m128i **rbv;         /* match scores [x][q]: rm, rm[0] are allocated      */
  __m128i **sbv;         /* match scores for ssvfilter                        */
  uint8_t   tbm_b;    /* constant B->Mk cost:    scaled log 2/M(M+1)       */
  uint8_t   tec_b;    /* constant E->C  cost:    scaled log 0.5            */
  uint8_t   tjb_b;    /* constant NCJ move cost: scaled log 3/(L+3)        */
  float     scale_b;    /* typically 3 / log2: scores scale to 1/3 bits      */
  uint8_t   base_b;            /* typically +190: offset of uchar scores            */
  uint8_t   bias_b;    /* positive bias to emission scores, make them >=0   */

  /* ViterbiFilter uses scaled swords: 8x signed 16-bit integer vectors              */
  __m128i **rwv;    /* [x][q]: rw, rw[0] are allocated  [Kp][Q8]         */
  __m128i  *twv;    /* transition score blocks          [8*Q8]           */
  int16_t   xw[p7O_NXSTATES][p7O_NXTRANS]; /* NECJ state transition costs            */
  float     scale_w;            /* score units: typically 500 / log(2), 1/500 bits   */
  int16_t   base_w;             /* offset of sword scores: typically +12000          */
  int16_t   ddbound_w;    /* threshold precalculated for lazy DD evaluation    */
  float     ncj_roundoff;  /* missing precision on NN,CC,JJ after rounding      */

  /* Forward, Backward use IEEE754 single-precision floats: 4x vectors               */
  __m128 **rfv;         /* [x][q]:  rf, rf[0] are allocated [Kp][Q4]         */
  __m128  *tfv;          /* transition probability blocks    [8*Q4]           */
  float    xf[p7O_NXSTATES][p7O_NXTRANS]; /* NECJ transition costs                   */

  /* Our actual vector mallocs, before we align the memory                           */
  __m128i  *rbv_mem;
  __m128i  *sbv_mem;
  __m128i  *rwv_mem;
  __m128i  *twv_mem;
  __m128   *tfv_mem;
  __m128   *rfv_mem;

  /* Disk offset information for hmmpfam's fast model retrieval                      */
  off_t  offs[p7_NOFFSETS];     /* p7_{MFP}OFFSET, or -1                             */

  /* Disk offset bookkeeping for h3f:                                                */
  off_t  roff;                  /* record offset (start of record); -1 if none       */
  off_t  eoff;                  /* offset to last byte of record; -1 if unknown      */

  /* Information, annotation copied from parent profile:                             */
  char  *name;      /* unique name of model                              */
  char  *acc;      /* unique accession of model, or NULL                */
  char  *desc;                  /* brief (1-line) description of model, or NULL      */
  char  *rf;                    /* reference line           1..M; *ref=0: unused     */
  char  *mm;                    /* modelmask line           1..M; *ref=0: unused     */
  char  *cs;                    /* consensus structure line 1..M, *cs=0: unused      */
  char  *consensus;    /* consensus residues for ali display, 1..M          */
  float  evparam[p7_NEVPARAM];   /* parameters for determining E-values, or UNSET     */
  float  cutoff[p7_NCUTOFFS];   /* per-seq/per-dom bit cutoffs, or UNSET             */
  float  compo[p7_MAXABET];  /* per-model HMM filter composition, or UNSET        */
  const ESL_ALPHABET *abc;  /* copy of ptr to alphabet information               */

  /* Information about current configuration, size, allocation                       */
  int    L;      /* current configured target seq length              */
  int    M;      /* model length                                      */
  int    max_length;    /* upper bound on emitted sequence length            */
  int    allocM;    /* maximum model length currently allocated for      */
  int    allocQ4;    /* p7_NQF(allocM): alloc size for tf, rf             */
  int    allocQ8;    /* p7_NQW(allocM): alloc size for tw, rw             */
  int    allocQ16;    /* p7_NQB(allocM): alloc size for rb                 */
  int    mode;      /* currently must be p7_LOCAL                        */
  float  nj;      /* expected # of J's: 0 or 1, uni vs. multihit       */

  int    clone;                 /* this optimized profile structure is just a copy   */
								/* of another profile structre.  all pointers of     */
								/* this structure should not be freed.               */
} P7_OPROFILE;

typedef struct {
  int            count;       /* number of <P7_OPROFILE> objects in the block */
  int            listSize;    /* maximum number elements in the list          */
  P7_OPROFILE  **list;        /* array of <P7_OPROFILE> objects               */
} P7_OM_BLOCK;

/* retrieve match odds ratio [k][x]
 * this gets used in p7_alidisplay.c, when we're deciding if a residue is conserved or not */
static inline float
p7_oprofile_FGetEmission(const P7_OPROFILE *om, int k, int x)
{
  union { __m128 v; float p[4]; } u;
  int   Q = p7O_NQF(om->M);
  int   q = ((k-1) % Q);
  int   r = (k-1)/Q;
  u.v = om->rfv[x][q];
  return u.p[r];
}

/*****************************************************************
 * 2. P7_OMX: a one-row dynamic programming matrix
 *****************************************************************/

enum p7x_scells_e { p7X_M = 0, p7X_D = 1, p7X_I = 2 };
#define p7X_NSCELLS 3

/* Besides ENJBC states, we may also store a rescaling factor on each row  */
enum p7x_xcells_e { p7X_E = 0, p7X_N = 1, p7X_J = 2, p7X_B = 3, p7X_C = 4, p7X_SCALE = 5 };
#define p7X_NXCELLS 6

/*
 *
 * dpf[][]
 *    to access M(i,k) for i=0,1..L; k=1..M:  dpf[i][(k-1)/4 + p7X_M].element[(k-1)%4]
 *
 * xmx[] arrays for individual special states:
 *    xmx[ENJBC] = [0 1 2 3][4 5 6 7]..[L-2 L-1 L x]     XRQ >= (L/4)+1
 *    to access B[i] for example, for i=0..L:   xmx[B][i/4].x[i%4]  (quad i/4; element i%4).
 */
typedef struct p7_omx_s {
  int       M;      /* current actual model dimension                              */
  int       L;      /* current actual sequence dimension                           */

  /* The main dynamic programming matrix for M,D,I states                                      */
  __m128  **dpf;    /* striped DP matrix for [0,1..L][0..Q-1][MDI], float vectors  */
  __m128i **dpw;    /* striped DP matrix for [0,1..L][0..Q-1][MDI], sword vectors  */
  __m128i **dpb;    /* striped DP matrix for [0,1..L][0..Q-1] uchar vectors        */
  void     *dp_mem;    /* DP memory shared by <dpb>, <dpw>, <dpf>                     */
  int       allocR;    /* current allocated # rows in dp{uf}. allocR >= validR >= L+1 */
  int       validR;    /* current # of rows actually pointing at DP memory            */
  int       allocQ4;    /* current set row width in <dpf> quads:   allocQ4*4 >= M      */
  int       allocQ8;    /* current set row width in <dpw> octets:  allocQ8*8 >= M      */
  int       allocQ16;    /* current set row width in <dpb> 16-mers: allocQ16*16 >= M    */
  size_t    ncells;    /* current allocation size of <dp_mem>, in accessible cells    */

  /* The X states (for full,parser; or NULL, for scorer)                                       */
  float    *xmx;          /* logically [0.1..L][ENJBCS]; indexed [i*p7X_NXCELLS+s]       */
  void     *x_mem;    /* X memory before 16-byte alignment                           */
  int       allocXR;    /* # of rows allocated in each xmx[] array; allocXR >= L+1     */
  float     totscale;    /* log of the product of all scale factors (0.0 if unscaled)   */
  int       has_own_scales;  /* TRUE to use own scale factors; FALSE if scales provided     */

  /* Parsers,scorers only hold a row at a time, so to get them to dump full matrix, it
   * must be done during a DP calculation, after each row is calculated
   */
  int     debugging;    /* TRUE if we're in debugging mode                             */
  FILE   *dfp;      /* output stream for diagnostics                               */
} P7_OMX;

/* ?MXo(q) access macros work for either uchar or float, so long as you
 * init your "dp" to point to the appropriate array.
 */
#define MMXo(q)   (dp[(q) * p7X_NSCELLS + p7X_M])
#define DMXo(q)   (dp[(q) * p7X_NSCELLS + p7X_D])
#define IMXo(q)   (dp[(q) * p7X_NSCELLS + p7X_I])
#define XMXo(i,s) (xmx[(i) * p7X_NXCELLS + s])

/* and this version works with a ptr to the approp DP row. */
#define MMO(dp,q) ((dp)[(q) * p7X_NSCELLS + p7X_M])
#define DMO(dp,q) ((dp)[(q) * p7X_NSCELLS + p7X_D])
#define IMO(dp,q) ((dp)[(q) * p7X_NSCELLS + p7X_I])

static inline float
p7_omx_FGetMDI(const P7_OMX *ox, int s, int i, int k)
{
  union { __m128 v; float p[4]; } u;
  int   Q = p7O_NQF(ox->M);
  int   q = p7X_NSCELLS * ((k-1) % Q) + s;
  int   r = (k-1)/Q;
  u.v = ox->dpf[i][q];
  return u.p[r];
}

static inline void
p7_omx_FSetMDI(const P7_OMX *ox, int s, int i, int k, float val)
{
  union { __m128 v; float p[4]; } u;
  int   Q = p7O_NQF(ox->M);
  int   q = p7X_NSCELLS * ((k-1) % Q) + s;
  int   r = (k-1)/Q;

  u.v           = ox->dpf[i][q];
  u.p[r]        = val;
  ox->dpf[i][q] = u.v;
}

/*****************************************************************
 * 3. Declarations of the external API.
 *****************************************************************/

/* p7_omx.c */
P7_OMX      *p7_omx_Create(int allocM, int allocL, int allocXL);
int          p7_omx_GrowTo(P7_OMX *ox, int allocM, int allocL, int allocXL);
int          p7_omx_FDeconvert(P7_OMX *ox, P7_GMX *gx);
int          p7_omx_Reuse  (P7_OMX *ox);
void         p7_omx_Destroy(P7_OMX *ox);

int          p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox, int truefalse);
int          p7_omx_DumpMFRow(P7_OMX *ox, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
int          p7_omx_DumpVFRow(P7_OMX *ox, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
int          p7_omx_DumpFBRow(P7_OMX *ox, int logify, int rowi, int width, int precision, float xE, float xN, float xJ, float xB, float xC);

/* p7_oprofile.c */
P7_OPROFILE *p7_oprofile_Create(int M, const ESL_ALPHABET *abc);
int          p7_oprofile_IsLocal(const P7_OPROFILE *om);
void         p7_oprofile_Destroy(P7_OPROFILE *om);
size_t       p7_oprofile_Sizeof(P7_OPROFILE *om);
P7_OPROFILE *p7_oprofile_Copy(P7_OPROFILE *om);
P7_OPROFILE *p7_oprofile_Clone(const P7_OPROFILE *om);
int          p7_oprofile_UpdateFwdEmissionScores(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
int          p7_oprofile_UpdateVitEmissionScores(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
int          p7_oprofile_UpdateMSVEmissionScores(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);

int          p7_oprofile_Convert(const P7_PROFILE *gm, P7_OPROFILE *om);
int          p7_oprofile_ReconfigLength    (P7_OPROFILE *om, int L);
int          p7_oprofile_ReconfigMSVLength (P7_OPROFILE *om, int L);
int          p7_oprofile_ReconfigRestLength(P7_OPROFILE *om, int L);
int          p7_oprofile_ReconfigMultihit  (P7_OPROFILE *om, int L);
int          p7_oprofile_ReconfigUnihit    (P7_OPROFILE *om, int L);

int          p7_oprofile_Dump(FILE *fp, const P7_OPROFILE *om);
int          p7_oprofile_Sample(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, const P7_BG *bg, int M, int L,
			   P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_OPROFILE **ret_om);
int          p7_oprofile_Compare(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);
int          p7_profile_SameAsMF(const P7_OPROFILE *om, P7_PROFILE *gm);
int          p7_profile_SameAsVF(const P7_OPROFILE *om, P7_PROFILE *gm);

int          p7_oprofile_GetFwdTransitionArray(const P7_OPROFILE *om, int type, float *arr );
int          p7_oprofile_GetSSVEmissionScoreArray(const P7_OPROFILE *om, uint8_t *arr );
int          p7_oprofile_GetFwdEmissionScoreArray(const P7_OPROFILE *om, float *arr );
int          p7_oprofile_GetFwdEmissionArray(const P7_OPROFILE *om, P7_BG *bg, float *arr );

/* decoding.c */
int p7_Decoding      (const P7_OPROFILE *om, const P7_OMX *oxf,       P7_OMX *oxb, P7_OMX *pp);
int p7_DomainDecoding(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);

/* fwdback.c */
int p7_Forward       (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
int p7_ForwardParser (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
int p7_Backward      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
int p7_BackwardParser(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);

/* io.c */
int p7_oprofile_Write(FILE *ffp, FILE *pfp, P7_OPROFILE *om);
int p7_oprofile_ReadMSV (P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
int p7_oprofile_ReadInfoMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
int p7_oprofile_ReadBlockMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OM_BLOCK *hmmBlock);
int p7_oprofile_ReadRest(P7_HMMFILE *hfp, P7_OPROFILE *om);
int p7_oprofile_Position(P7_HMMFILE *hfp, off_t offset);

P7_OM_BLOCK *p7_oprofile_CreateBlock(int size);
void p7_oprofile_DestroyBlock(P7_OM_BLOCK *block);

/* ssvfilter.c */
int p7_SSVFilter    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);

/* msvfilter.c */
int p7_MSVFilter           (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
int p7_SSVFilter_longtarget(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *msvdata, P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist);

/* null2.c */
int p7_Null2_ByExpectation(const P7_OPROFILE *om, const P7_OMX *pp, float *null2);
int p7_Null2_ByTrace      (const P7_OPROFILE *om, const P7_TRACE *tr, int zstart, int zend, P7_OMX *wrk, float *null2);

/* optacc.c */
int p7_OptimalAccuracy(const P7_OPROFILE *om, const P7_OMX *pp,       P7_OMX *ox, float *ret_e);
int p7_OATrace        (const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr);

/* stotrace.c */
int p7_StochasticTrace(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr);

/* vitfilter.c */
int p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
int p7_ViterbiFilter_longtarget(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox,
										float filtersc, double P, P7_HMM_WINDOWLIST *windowlist);

/* vitscore.c */
int p7_ViterbiScore (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);

/*****************************************************************
 * 4. Implementation specific initialization
 *****************************************************************/
static inline void
impl_Init(void)
{
#ifdef HAVE_FLUSH_ZERO_MODE
  /* In order to avoid the performance penalty dealing with sub-normal
   * values in the floating point calculations, set the processor flag
   * so sub-normals are "flushed" immediately to zero.
   */
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

#ifdef _PMMINTRIN_H_INCLUDED
  /*
   * FLUSH_ZERO doesn't necessarily work in non-SIMD calculations
   * (yes on 64-bit, maybe not of 32-bit). This ensures that those
   * scalar calculations will agree across architectures.
   * (See TW notes  2012/0106_printf_underflow_bug/00NOTES for details)
   */
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
}
#endif /* P7_IMPL_SSE_INCLUDED */

/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Version 3.1b3; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 *
 * HMMER is distributed under the terms of the GNU General Public License
 * (GPLv3). See the LICENSE file for details.
 *****************************************************************/

/*
 * Currently (and this remains in flux as of 14 Dec 07) an optimized
 * implementation is required to provide an MSVFilter(),
 * ViterbiFilter() and a ForwardFilter() implementation. A call to
 * p7_oprofile_Convert() makes an optimized profile that works for
 * all filters.
 *
 * Any "Filter" returns a score may be an approximation (with
 * characterized or at least characterizable error), and which may
 * have limited upper range, such that high scores are returned as
 * eslINFINITY. Additionally, Filters might only work on local
 * alignment modes, because they are allowed to make assumptions about
 * the range of scores.
 *
 * Here, MSVFilter() and ViterbiFilter() are 8-bit lspace
 * implementations with limited precision and limited range (max 20
 * bits); ForwardFilter() is a pspace float implementation with
 * correct precision and limited range (max ~127 bits). Both require
 * local mode models.
 *
 * An optimized implementation may also provide other optimized
 * routines. It provides specialized Convert*() functions for these,
 * which may no-op (if the OPROFILE already suffices), or may
 * overwrite parts of the OPROFILE that Filters or other routines
 * might need. Therefore, after using a "bonus" function, a fresh
 * Convert() will be needed before a Filter() is called again. This
 * API is tentative.
 *
 * For example, here, ViterbiScore() is a 32-bit lspace float SSE
 * implementation of the Viterbi algorithm.
 *
 * A "Score" function might be an additional target for optimization,
 * for example. A "Score" function returns a correct score with full
 * floating-point precision and range, and works for any mode model.
 *
 * In the generic implementation, profile scores are 32-bit floating
 * point log-odds scores. In an optimized implementation, internally,
 * profile scores can be of any type, and may be in log space (lspace)
 * or probability space (pspace). (Calculations in probability space
 * are useful in the Forward algorithm, but always limit range.)  A
 * shorthand of "lspace uchar" means log-odds scores stored as
 * unsigned chars, for example; "pspace float" means odds ratios
 * stored as floats.
 *
 * A note on memory alignment: malloc() is required to return a
 * pointer "suitably aligned so that it may be aligned to a pointer of
 * any type of object" (C99 7.20.3). __m128 vectors are 128-bits wide,
 * so malloc() ought to return a pointer aligned on a 16-byte
 * boundary.  However, this is not the case for glibc, and apparently
 * other system libraries. Google turns up threads of arguments
 * between glibc and gcc developers over whose problem this is; this
 * argument has apparently not been resolved, and is of no help.
 * Here, we manually align the relevant pointers by overallocating in
 * *_mem with malloc, then arithmetically manipulating the address to
 * mask off (~0xf).
 */

/*** End of inlined file: impl_sse.h ***/


#elif defined (p7_IMPL_VMX)

/*** Start of inlined file: impl_vmx.h ***/
#ifndef P7_IMPL_VMX_INCLUDED
#define P7_IMPL_VMX_INCLUDED

#ifndef __APPLE_ALTIVEC__
#include <altivec.h>
#endif

//#include "p7_config.h"

//#include "esl_alphabet.h"
//#include "esl_random.h"

//#include "hmmer.h"

/* In calculating Q, the number of vectors we need in a row, we have
 * to make sure there's at least 2, or a striped implementation fails.
 */
#define p7O_NQB(M)   ( ESL_MAX(2, ((((M)-1) / 16) + 1)))   /* 16 uchars  */
#define p7O_NQW(M)   ( ESL_MAX(2, ((((M)-1) / 8)  + 1)))   /*  8 words   */
#define p7O_NQF(M)   ( ESL_MAX(2, ((((M)-1) / 4)  + 1)))   /*  4 floats  */

/*****************************************************************
 * 1. P7_OPROFILE: an optimized score profile
 *****************************************************************/
/* The OPROFILE is striped [Farrar07] and interleaved, as is the DP matrix.
 * For example, the layout of a profile for an M=14 model (xref J2/46):
 *
 * rsc[x] : striped blocks of M emissions, starting with q=0
 *                1     11     1      1
 *             1593   2604   371x   482x
 *
 * tsc:  grouped in order of accession in DP for 7 transition scores;
 *       starting at q=0 for all but the three transitions to M, which
 *       are rotated by -1 and rightshifted. DD's follow separately,
 *       starting at q=0.
 *
 *        {     1      1     1     1     1     1     1 }
 *        {  1593   x482  x482  x482  1593  1593  1593 }
 *        { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
 *
 *        {    11      1     1     1    11    11    11 }
 *        {  2604   1593  1593  1593  2604  2604  2604 }
 *        { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
 *
 *        {    1      11    11    11    1     1     1  }
 *        {  371x   2604  2604  2604  371x  371x  371x }
 *        { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
 *
 *        {    1      1     1     1     1     1     1  }
 *        {  482x   371x  371x  371x  482x  482x  482x }
 *        { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
 *
 *        {     1    11    1     1  }
 *        {  1593  2604  371x  482x }
 *        { [TDD] [TDD] [TDD] [TDD] }
 *
 */

#define p7O_NXSTATES  4		/* special states stored: ENJC                       */
#define p7O_NXTRANS   2         /* special states all have 2 transitions: move, loop */
#define p7O_NTRANS    8		/* 7 core transitions + BMk entry                    */
enum p7o_xstates_e      { p7O_E    = 0, p7O_N    = 1,  p7O_J  = 2,  p7O_C  = 3 };
enum p7o_xtransitions_e { p7O_MOVE = 0, p7O_LOOP = 1 };
enum p7o_tsc_e          { p7O_BM   = 0, p7O_MM   = 1,  p7O_IM = 2,  p7O_DM = 3, p7O_MD   = 4, p7O_MI   = 5,  p7O_II = 6,  p7O_DD = 7 };

typedef struct p7_oprofile_s {
  /* MSVFilter uses scaled, biased uchars: 16x unsigned byte vectors                 */
  vector unsigned char **rbv;   /* match scores [x][q]: rm, rm[0] are allocated      */
  uint8_t   tbm_b;		/* constant B->Mk cost:    scaled log 2/M(M+1)       */
  uint8_t   tec_b;		/* constant E->C  cost:    scaled log 0.5            */
  uint8_t   tjb_b;		/* constant NCJ move cost: scaled log 3/(L+3)        */
  float     scale_b;		/* typically 3 / log2: scores scale to 1/3 bits      */
  uint8_t   base_b;  	        /* typically +190: offset of uchar scores            */
  uint8_t   bias_b;		/* positive bias to emission scores, make them >=0   */

  /* ViterbiFilter uses scaled swords: 8x signed 16-bit integer vectors              */
  vector signed short **rwv;	/* [x][q]: rw, rw[0] are allocated  [Kp][Q8]         */
  vector signed short  *twv;	/* transition score blocks          [8*Q8]           */
  int16_t   xw[p7O_NXSTATES][p7O_NXTRANS]; /* NECJ state transition costs            */
  float     scale_w;            /* score units: typically 500 / log(2), 1/500 bits   */
  int16_t   base_w;             /* offset of sword scores: typically +12000          */
  int16_t   ddbound_w;		/* threshold precalculated for lazy DD evaluation    */
  float     ncj_roundoff;	/* missing precision on NN,CC,JJ after rounding      */

  /* Forward, Backward use IEEE754 single-precision floats: 4x vectors               */
  vector float **rfv;           /* [x][q]:  rf, rf[0] are allocated [Kp][Q4]         */
  vector float  *tfv;	    	/* transition probability blocks    [8*Q4]           */
  float    xf[p7O_NXSTATES][p7O_NXTRANS]; /* NECJ transition costs                   */

  /* Our actual vector mallocs, before we align the memory                           */
  vector unsigned char  *rbv_mem;
  vector signed short   *rwv_mem;
  vector signed short   *twv_mem;
  vector float          *tfv_mem;
  vector float          *rfv_mem;

  /* Disk offset information for hmmpfam's fast model retrieval                      */
  off_t  offs[p7_NOFFSETS];     /* p7_{MFP}OFFSET, or -1                             */

  /* Disk offset bookkeeping for h3f:                                                */
  off_t  roff;                  /* record offset (start of record); -1 if none       */
  off_t  eoff;                  /* offset to last byte of record; -1 if unknown      */

  /* Information, annotation copied from parent profile:                             */
  char  *name;			/* unique name of model                              */
  char  *acc;			/* unique accession of model, or NULL                */
  char  *desc;                  /* brief (1-line) description of model, or NULL      */
  char  *rf;                    /* reference line           1..M; *ref=0: unused     */
  char  *mm;                    /* modelmask line           1..M; *ref=0: unused     */
  char  *cs;                    /* consensus structure line 1..M, *cs=0: unused      */
  char  *consensus;		/* consensus residues for ali display, 1..M          */
  float  evparam[p7_NEVPARAM]; 	/* parameters for determining E-values, or UNSET     */
  float  cutoff[p7_NCUTOFFS]; 	/* per-seq/per-dom bit cutoffs, or UNSET             */
  float  compo[p7_MAXABET];	/* per-model HMM filter composition, or UNSET        */
  const ESL_ALPHABET *abc;	/* copy of ptr to alphabet information               */

  /* Information about current configuration, size, allocation                       */
  int    L;			/* current configured target seq length              */
  int    M;			/* model length                                      */
  int    max_length;		/* upper bound on emitted seq length                 */
  int    allocM;		/* maximum model length currently allocated for      */
  int    allocQ4;		/* p7_NQF(allocM): alloc size for tf, rf             */
  int    allocQ8;		/* p7_NQW(allocM): alloc size for tw, rw             */
  int    allocQ16;		/* p7_NQB(allocM): alloc size for rb                 */
  int    mode;			/* currently must be p7_LOCAL                        */
  float  nj;			/* expected # of J's: 0 or 1, uni vs. multihit       */

  int    clone;                 /* this optimized profile structure is just a copy   */
								/* of another profile structre.  all pointers of     */
								/* this structure should not be freed.               */
} P7_OPROFILE;

typedef struct {
  int            count;       /* number of <P7_OPROFILE> objects in the block */
  int            listSize;    /* maximum number elements in the list          */
  P7_OPROFILE  **list;        /* array of <P7_OPROFILE> objects               */
} P7_OM_BLOCK;

/* retrieve match odds ratio [k][x]
 * this gets used in p7_alidisplay.c, when we're deciding if a residue is conserved or not */
static inline float
p7_oprofile_FGetEmission(const P7_OPROFILE *om, int k, int x)
{
  union { vector float v; float p[4]; } u;
  int   Q = p7O_NQF(om->M);
  int   q = ((k-1) % Q);
  int   r = (k-1)/Q;
  u.v = om->rfv[x][q];
  return u.p[r];
}

/*****************************************************************
 * 2. P7_OMX: a one-row dynamic programming matrix
 *****************************************************************/

enum p7x_scells_e { p7X_M = 0, p7X_D = 1, p7X_I = 2 };
#define p7X_NSCELLS 3

/* Besides ENJBC states, we may also store a rescaling factor on each row  */
enum p7x_xcells_e { p7X_E = 0, p7X_N = 1, p7X_J = 2, p7X_B = 3, p7X_C = 4, p7X_SCALE = 5 };
#define p7X_NXCELLS 6

/*
 *
 * dpf[][]
 *    to access M(i,k) for i=0,1..L; k=1..M:  dpf[i][(k-1)/4 + p7X_M].element[(k-1)%4]
 *
 * xmx[] arrays for individual special states:
 *    xmx[ENJBC] = [0 1 2 3][4 5 6 7]..[L-2 L-1 L x]     XRQ >= (L/4)+1
 *    to access B[i] for example, for i=0..L:   xmx[B][i/4].x[i%4]  (quad i/4; element i%4).
 */
typedef struct p7_omx_s {
  int       M;			/* current actual model dimension                              */
  int       L;			/* current actual sequence dimension                           */

  /* The main dynamic programming matrix for M,D,I states                                      */
  vector float         **dpf;	/* striped DP matrix for [0,1..L][0..Q-1][MDI], float vectors  */
  vector signed short  **dpw;	/* striped DP matrix for [0,1..L][0..Q-1][MDI], sword vectors  */
  vector unsigned char **dpb;	/* striped DP matrix for [0,1..L][0..Q-1] uchar vectors        */
  void     *dp_mem;		/* DP memory shared by <dpb>, <dpw>, <dpf>                     */
  int       allocR;		/* current allocated # rows in dp{uf}. allocR >= validR >= L+1 */
  int       validR;		/* current # of rows actually pointing at DP memory            */
  int       allocQ4;		/* current set row width in <dpf> quads:   allocQ4*4 >= M      */
  int       allocQ8;		/* current set row width in <dpw> octets:  allocQ8*8 >= M      */
  int       allocQ16;		/* current set row width in <dpb> 16-mers: allocQ16*16 >= M    */
  size_t    ncells;		/* current allocation size of <dp_mem>, in accessible cells    */

  /* The X states (for full,parser; or NULL, for scorer)                                       */
  float    *xmx;        	/* logically [0.1..L][ENJBCS]; indexed [i*p7X_NXCELLS+s]       */
  void     *x_mem;		/* X memory before 16-byte alignment                           */
  int       allocXR;		/* # of rows allocated in each xmx[] array; allocXR >= L+1     */
  float     totscale;		/* log of the product of all scale factors (0.0 if unscaled)   */
  int       has_own_scales;	/* TRUE to use own scale factors; FALSE if scales provided     */

  /* Parsers,scorers only hold a row at a time, so to get them to dump full matrix, it
   * must be done during a DP calculation, after each row is calculated
   */
  int     debugging;		/* TRUE if we're in debugging mode                             */
  FILE   *dfp;			/* output stream for diagnostics                               */
} P7_OMX;

/* ?MXo(q) access macros work for either uchar or float, so long as you
 * init your "dp" to point to the appropriate array.
 */
#define MMXo(q)   (dp[(q) * p7X_NSCELLS + p7X_M])
#define DMXo(q)   (dp[(q) * p7X_NSCELLS + p7X_D])
#define IMXo(q)   (dp[(q) * p7X_NSCELLS + p7X_I])
#define XMXo(i,s) (xmx[(i) * p7X_NXCELLS + s])

/* and this version works with a ptr to the approp DP row. */
#define MMO(dp,q) ((dp)[(q) * p7X_NSCELLS + p7X_M])
#define DMO(dp,q) ((dp)[(q) * p7X_NSCELLS + p7X_D])
#define IMO(dp,q) ((dp)[(q) * p7X_NSCELLS + p7X_I])

static inline float
p7_omx_FGetMDI(const P7_OMX *ox, int s, int i, int k)
{
  union { vector float v; float p[4]; } u;
  int   Q = p7O_NQF(ox->M);
  int   q = p7X_NSCELLS * ((k-1) % Q) + s;
  int   r = (k-1)/Q;
  u.v = ox->dpf[i][q];
  return u.p[r];
}

static inline void
p7_omx_FSetMDI(const P7_OMX *ox, int s, int i, int k, float val)
{
  union { vector float v; float p[4]; } u;
  int   Q = p7O_NQF(ox->M);
  int   q = p7X_NSCELLS * ((k-1) % Q) + s;
  int   r = (k-1)/Q;

  u.v           = ox->dpf[i][q];
  u.p[r]        = val;
  ox->dpf[i][q] = u.v;
}

/*****************************************************************
 * 3. Declarations of the external API.
 *****************************************************************/

/* p7_omx.c */
P7_OMX      *p7_omx_Create(int allocM, int allocL, int allocXL);
int          p7_omx_GrowTo(P7_OMX *ox, int allocM, int allocL, int allocXL);
int          p7_omx_FDeconvert(P7_OMX *ox, P7_GMX *gx);
int          p7_omx_Reuse  (P7_OMX *ox);
void         p7_omx_Destroy(P7_OMX *ox);

int          p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox, int truefalse);
int          p7_omx_DumpMFRow(P7_OMX *ox, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
int          p7_omx_DumpVFRow(P7_OMX *ox, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
int          p7_omx_DumpFBRow(P7_OMX *ox, int logify, int rowi, int width, int precision, float xE, float xN, float xJ, float xB, float xC);

/* p7_oprofile.c */
P7_OPROFILE *p7_oprofile_Create(int M, const ESL_ALPHABET *abc);
int          p7_oprofile_IsLocal(const P7_OPROFILE *om);
void         p7_oprofile_Destroy(P7_OPROFILE *om);
size_t       p7_oprofile_Sizeof(P7_OPROFILE *om);
P7_OPROFILE *p7_oprofile_Copy(P7_OPROFILE *om);
P7_OPROFILE *p7_oprofile_Clone(const P7_OPROFILE *om);
int          p7_oprofile_UpdateFwdEmissionScores(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);

int          p7_oprofile_Convert(const P7_PROFILE *gm, P7_OPROFILE *om);
int          p7_oprofile_ReconfigLength    (P7_OPROFILE *om, int L);
int          p7_oprofile_ReconfigMSVLength (P7_OPROFILE *om, int L);
int          p7_oprofile_ReconfigRestLength(P7_OPROFILE *om, int L);
int          p7_oprofile_ReconfigMultihit  (P7_OPROFILE *om, int L);
int          p7_oprofile_ReconfigUnihit    (P7_OPROFILE *om, int L);

int          p7_oprofile_Dump(FILE *fp, const P7_OPROFILE *om);
int          p7_oprofile_Sample(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, const P7_BG *bg, int M, int L,
				       P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_OPROFILE **ret_om);
int          p7_oprofile_Compare(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);
int          p7_profile_SameAsMF(const P7_OPROFILE *om, P7_PROFILE *gm);
int          p7_profile_SameAsVF(const P7_OPROFILE *om, P7_PROFILE *gm);

int          p7_oprofile_GetFwdTransitionArray(const P7_OPROFILE *om, int type, float *arr );
int          p7_oprofile_GetSSVEmissionScoreArray(const P7_OPROFILE *om, uint8_t *arr );
int          p7_oprofile_GetFwdEmissionScoreArray(const P7_OPROFILE *om, float *arr );
int          p7_oprofile_GetFwdEmissionArray(const P7_OPROFILE *om, P7_BG *bg, float *arr );

/* decoding.c */
int p7_Decoding      (const P7_OPROFILE *om, const P7_OMX *oxf,       P7_OMX *oxb, P7_OMX *pp);
int p7_DomainDecoding(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);

/* fwdback.c */
int p7_Forward       (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
int p7_ForwardParser (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
int p7_Backward      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
int p7_BackwardParser(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);

/* io.c */
int p7_oprofile_Write(FILE *ffp, FILE *pfp, P7_OPROFILE *om);
int p7_oprofile_ReadMSV (P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
int p7_oprofile_ReadInfoMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
int p7_oprofile_ReadBlockMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OM_BLOCK *hmmBlock);
int p7_oprofile_ReadRest(P7_HMMFILE *hfp, P7_OPROFILE *om);
int p7_oprofile_Position(P7_HMMFILE *hfp, off_t offset);

P7_OM_BLOCK *p7_oprofile_CreateBlock(int size);
void p7_oprofile_DestroyBlock(P7_OM_BLOCK *block);

/* msvfilter.c */
int p7_MSVFilter    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
int p7_SSVFilter_longtarget(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *msvdata, P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist);

/* null2.c */
int p7_Null2_ByExpectation(const P7_OPROFILE *om, const P7_OMX *pp, float *null2);
int p7_Null2_ByTrace      (const P7_OPROFILE *om, const P7_TRACE *tr, int zstart, int zend, P7_OMX *wrk, float *null2);

/* optacc.c */
int p7_OptimalAccuracy(const P7_OPROFILE *om, const P7_OMX *pp,       P7_OMX *ox, float *ret_e);
int p7_OATrace        (const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr);

/* stotrace.c */
int p7_StochasticTrace(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox,
			      P7_TRACE *tr);

/* vitfilter.c */
int p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
int p7_ViterbiFilter_longtarget(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox,
							float filtersc, double P, P7_HMM_WINDOWLIST *windowlist);

/* vitscore.c */
int p7_ViterbiScore (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);

/*****************************************************************
 * 3. Implementation specific initialization
 *****************************************************************/
static inline void
impl_Init(void)
{
}

static inline void
impl_ThreadInit(void)
{
}

#endif /* P7_IMPL_VMX_INCLUDED */

/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Version 3.1b3; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 *
 * HMMER is distributed under the terms of the GNU General Public License
 * (GPLv3). See the LICENSE file for details.
 *****************************************************************/

/*
 * Currently (and this remains in flux as of 14 Dec 07) an optimized
 * implementation is required to provide an MSVFilter(),
 * ViterbiFilter() and a ForwardFilter() implementation. A call to
 * p7_oprofile_Convert() makes an optimized profile that works for
 * all filters.
 *
 * Any "Filter" returns a score may be an approximation (with
 * characterized or at least characterizable error), and which may
 * have limited upper range, such that high scores are returned as
 * eslINFINITY. Additionally, Filters might only work on local
 * alignment modes, because they are allowed to make assumptions about
 * the range of scores.
 *
 * Here, MSVFilter() and ViterbiFilter() are 8-bit lspace
 * implementations with limited precision and limited range (max 20
 * bits); ForwardFilter() is a pspace float implementation with
 * correct precision and limited range (max ~127 bits). Both require
 * local mode models.
 *
 * An optimized implementation may also provide other optimized
 * routines. It provides specialized Convert*() functions for these,
 * which may no-op (if the OPROFILE already suffices), or may
 * overwrite parts of the OPROFILE that Filters or other routines
 * might need. Therefore, after using a "bonus" function, a fresh
 * Convert() will be needed before a Filter() is called again. This
 * API is tentative.
 *
 * For example, here, ViterbiScore() is a 32-bit lspace float VMX
 * implementation of the Viterbi algorithm.
 *
 * A "Score" function might be an additional target for optimization,
 * for example. A "Score" function returns a correct score with full
 * floating-point precision and range, and works for any mode model.
 *
 * In the generic implementation, profile scores are 32-bit floating
 * point log-odds scores. In an optimized implementation, internally,
 * profile scores can be of any type, and may be in log space (lspace)
 * or probability space (pspace). (Calculations in probability space
 * are useful in the Forward algorithm, but always limit range.)  A
 * shorthand of "lspace uchar" means log-odds scores stored as
 * unsigned chars, for example; "pspace float" means odds ratios
 * stored as floats.
 *
 * A note on memory alignment: malloc() is required to return a
 * pointer "suitably aligned so that it may be aligned to a pointer of
 * any type of object" (C99 7.20.3). vectors are 128-bits wide,
 * so malloc() ought to return a pointer aligned on a 16-byte
 * boundary.  However, this is not the case for glibc, and apparently
 * other system libraries. Google turns up threads of arguments
 * between glibc and gcc developers over whose problem this is; this
 * argument has apparently not been resolved, and is of no help.
 * Here, we manually align the relevant pointers by overallocating in
 * *_mem with malloc, then arithmetically manipulating the address to
 * mask off (~0xf).
 */

/*** End of inlined file: impl_vmx.h ***/


#else

/*** Start of inlined file: impl_dummy.h ***/
#ifndef P7_IMPL_DUMMY_INCLUDED
#define P7_IMPL_DUMMY_INCLUDED

//#include "p7_config.h"

#include <math.h>		/* roundf() */

//#include "esl_alphabet.h"
//#include "esl_random.h"

//#include <pmmintrin.h>   /* DENORMAL_MODE */

//#include "hmmer.h"

/*****************************************************************
 * 1. P7_OPROFILE: an optimized score profile
 *****************************************************************/
typedef P7_PROFILE P7_OPROFILE;

#define p7O_NTRANS    8    /* 7 core transitions + BMk entry                    */
enum p7o_tsc_e          { p7O_BM   = 0, p7O_MM   = 1,  p7O_IM = 2,  p7O_DM = 3, p7O_MD   = 4, p7O_MI   = 5,  p7O_II = 6,  p7O_DD = 7 };

typedef struct {
  int            count;       /* number of <P7_OPROFILE> objects in the block */
  int            listSize;    /* maximum number elements in the list          */
  P7_OPROFILE  **list;        /* array of <P7_OPROFILE> objects               */
} P7_OM_BLOCK;

/* retrieve match odds ratio [k][x]
 * this gets used in p7_alidisplay.c, when we're deciding if a residue is conserved or not */
static inline float
p7_oprofile_FGetEmission(const P7_OPROFILE *om, int k, int x)
{
  return expf(p7P_MSC(om, k, x));
}

/*****************************************************************
 * 2. P7_OMX: a one-row dynamic programming matrix
 *****************************************************************/

typedef P7_GMX P7_OMX;

/*****************************************************************
 * 3. Declarations of the external API.
 *****************************************************************/

/* p7_omx.c */
P7_OMX      *p7_omx_Create(int allocM, int allocL, int allocXL);
int          p7_omx_GrowTo(P7_OMX *ox, int allocM, int allocL, int allocXL);
int          p7_omx_FDeconvert(P7_OMX *ox, P7_GMX *gx);
int          p7_omx_Reuse  (P7_OMX *ox);
void         p7_omx_Destroy(P7_OMX *ox);

int          p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox, int truefalse);
int          p7_omx_DumpMFRow(P7_OMX *ox, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
int          p7_omx_DumpVFRow(P7_OMX *ox, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
int          p7_omx_DumpFBRow(P7_OMX *ox, int logify, int rowi, int width, int precision, float xE, float xN, float xJ, float xB, float xC);

/* p7_oprofile.c */
P7_OPROFILE *p7_oprofile_Create(int M, const ESL_ALPHABET *abc);
int          p7_oprofile_IsLocal(const P7_OPROFILE *om);
void         p7_oprofile_Destroy(P7_OPROFILE *om);
size_t       p7_oprofile_Sizeof(P7_OPROFILE *om);
P7_OPROFILE *p7_oprofile_Copy(P7_OPROFILE *om);
P7_OPROFILE *p7_oprofile_Clone(P7_OPROFILE *om);
int          p7_oprofile_UpdateFwdEmissionScores(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);

int          p7_oprofile_Convert(const P7_PROFILE *gm, P7_OPROFILE *om);
int          p7_oprofile_ReconfigLength    (P7_OPROFILE *om, int L);
int          p7_oprofile_ReconfigMSVLength (P7_OPROFILE *om, int L);
int          p7_oprofile_ReconfigRestLength(P7_OPROFILE *om, int L);
int          p7_oprofile_ReconfigMultihit  (P7_OPROFILE *om, int L);
int          p7_oprofile_ReconfigUnihit    (P7_OPROFILE *om, int L);

int          p7_oprofile_Dump(FILE *fp, const P7_OPROFILE *om);
int          p7_oprofile_Sample(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, const P7_BG *bg, int M, int L,
				       P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_OPROFILE **ret_om);
int          p7_oprofile_Compare(P7_OPROFILE *om1, P7_OPROFILE *om2, float tol, char *errmsg);
int          p7_profile_SameAsMF(const P7_OPROFILE *om, P7_PROFILE *gm);
int          p7_profile_SameAsVF(const P7_OPROFILE *om, P7_PROFILE *gm);

int          p7_oprofile_GetFwdTransitionArray(const P7_OPROFILE *om, int type, float *arr );
int          p7_oprofile_GetSSVEmissionScoreArray(const P7_OPROFILE *om, uint8_t *arr );
int          p7_oprofile_GetFwdEmissionScoreArray(const P7_OPROFILE *om, float *arr );
int          p7_oprofile_GetFwdEmissionArray(const P7_OPROFILE *om, P7_BG *bg, float *arr );

/* decoding.c */
int p7_Decoding      (const P7_OPROFILE *om, const P7_OMX *oxf,       P7_OMX *oxb, P7_OMX *pp);
int p7_DomainDecoding(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);

/* fwdback.c */
int p7_Forward       (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
int p7_ForwardParser (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
int p7_Backward      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
int p7_BackwardParser(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);

/* io.c */
int p7_oprofile_Write(FILE *ffp, FILE *pfp, P7_OPROFILE *om);
int p7_oprofile_ReadMSV (P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om);
int p7_oprofile_ReadBlockMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OM_BLOCK *hmmBlock);
int p7_oprofile_ReadRest(P7_HMMFILE *hfp, P7_OPROFILE *om);

P7_OM_BLOCK *p7_oprofile_CreateBlock(int size);
void p7_oprofile_DestroyBlock(P7_OM_BLOCK *block);

/* msvfilter.c */
int p7_MSVFilter    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
int p7_SSVFilter_longtarget(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *msvdata, P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist);

/* null2.c */
int p7_Null2_ByExpectation(const P7_OPROFILE *om, P7_OMX *pp, float *null2);
int p7_Null2_ByTrace      (const P7_OPROFILE *om, const P7_TRACE *tr, int zstart, int zend, P7_OMX *wrk, float *null2);

/* optacc.c */
int p7_OptimalAccuracy(const P7_OPROFILE *om, const P7_OMX *pp,       P7_OMX *ox, float *ret_e);
int p7_OATrace        (const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr);

/* stotrace.c */
int p7_StochasticTrace(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox,
			      P7_TRACE *tr);

/* vitfilter.c */
int p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
int p7_ViterbiFilter_longtarget(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox,
							float filtersc, double P, P7_HMM_WINDOWLIST *windowlist);

/* vitscore.c */
int p7_ViterbiScore (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);

/*****************************************************************
 * 4. Implementation specific initialization
 *****************************************************************/
static inline void
impl_Init(void)
{
}

static inline void
impl_ThreadInit(void)
{
}

#endif /* P7_IMPL_DUMMY_INCLUDED */

/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Version 3.1b3; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 *
 * HMMER is distributed under the terms of the GNU General Public License
 * (GPLv3). See the LICENSE file for details.
 *****************************************************************/

/*
 * Currently (and this remains in flux as of 14 Dec 07) an optimized
 * implementation is required to provide an MSVFilter(),
 * ViterbiFilter() and a ForwardFilter() implementation. A call to
 * p7_oprofile_Convert() makes an optimized profile that works for
 * all filters.
 *
 * Any "Filter" returns a score may be an approximation (with
 * characterized or at least characterizable error), and which may
 * have limited upper range, such that high scores are returned as
 * eslINFINITY. Additionally, Filters might only work on local
 * alignment modes, because they are allowed to make assumptions about
 * the range of scores.
 *
 * Here, MSVFilter() and ViterbiFilter() are 8-bit lspace
 * implementations with limited precision and limited range (max 20
 * bits); ForwardFilter() is a pspace float implementation with
 * correct precision and limited range (max ~127 bits). Both require
 * local mode models.
 *
 * An optimized implementation may also provide other optimized
 * routines. It provides specialized Convert*() functions for these,
 * which may no-op (if the OPROFILE already suffices), or may
 * overwrite parts of the OPROFILE that Filters or other routines
 * might need. Therefore, after using a "bonus" function, a fresh
 * Convert() will be needed before a Filter() is called again. This
 * API is tentative.
 *
 * A "Score" function might be an additional target for optimization,
 * for example. A "Score" function returns a correct score with full
 * floating-point precision and range, and works for any mode model.
 *
 * In the generic implementation, profile scores are 32-bit floating
 * point log-odds scores. In an optimized implementation, internally,
 * profile scores can be of any type, and may be in log space (lspace)
 * or probability space (pspace). (Calculations in probability space
 * are useful in the Forward algorithm, but always limit range.)  A
 * shorthand of "lspace uchar" means log-odds scores stored as
 * unsigned chars, for example; "pspace float" means odds ratios
 * stored as floats.
 *
 * A note on memory alignment: malloc() is required to return a
 * pointer "suitably aligned so that it may be aligned to a pointer of
 * any type of object" (C99 7.20.3). __m128 vectors are 128-bits wide,
 * so malloc() ought to return a pointer aligned on a 16-byte
 * boundary.  However, this is not the case for glibc, and apparently
 * other system libraries. Google turns up threads of arguments
 * between glibc and gcc developers over whose problem this is; this
 * argument has apparently not been resolved, and is of no help.
 * Here, we manually align the relevant pointers by overallocating in
 * *_mem with malloc, then arithmetically manipulating the address to
 * mask off (~0xf).
 */

/*** End of inlined file: impl_dummy.h ***/


#endif

/*****************************************************************
 * 15. The FM-index acceleration to the SSV filter.  Only works for SSE
 *****************************************************************/
#define FM_MAX_LINE 256

/* Structure the 2D occ array into a single array.  "type" is either b or sb.
 * Note that one extra count value is required by RLE, one 4-byte int for
 * each superblock count vector, and one 2-byte short for each block count
 * vector. This is small overhead, even for a small alphabet like dna.
 */
#define FM_OCC_CNT( type, i, c)  ( occCnts_##type[(meta->alph_size)*(i) + (c)])

enum fm_alphabettypes_e {
  fm_DNA        = 0,  //acgt,  2 bit
  //fm_DNA_full   = 1,  //includes ambiguity codes, 4 bit.
  fm_AMINO      = 4,  // 5 bit
};
/*TODO: fm_DNA_full has currently been disabled because of problems with how the
 * FM index handles very long runs of the same character (in this case, Ns).
 * See wheelert/notebook/2013/12-11-FM-alphabet-speed notes on 12/12.
 */

enum fm_direction_e {
  fm_forward    = 0,
  fm_backward   = 1,
};

typedef struct fm_interval_s {
  int   lower;
  int   upper;
} FM_INTERVAL;

typedef struct fm_hit_s {
  uint64_t  start;
  uint32_t  block;
  int       direction;
  int       length;
  int       sortkey;
} FM_HIT;

typedef struct fm_ambiglist_s {
  FM_INTERVAL *ranges;
  uint32_t     count;
  uint32_t     size;
} FM_AMBIGLIST;

typedef struct fm_seqdata_s {

  uint32_t target_id;      // Which sequence in the target database did this segment come from (can be multiple segment per sequence, if a sequence has Ns)
  uint64_t target_start;   // The position in sequence {id} in the target database at which this sequence-block starts (usually 1, unless its a long sequence split out over multiple FMs)
  uint32_t fm_start;       // The position in the FM block at which this sequence begins
  uint32_t length;         // Length of this sequence segment  (usually the length of the target sequence, unless its a long sequence split out over multiple FMs)

  //meta data taken from the sequence this segment was taken from
  uint16_t name_length;
  uint16_t source_length;
  uint16_t acc_length;
  uint16_t desc_length;
  char     *name;
  char     *source;
  char     *acc;
  char     *desc;
} FM_SEQDATA;

typedef struct fm_metadata_s {
  uint8_t  fwd_only;
  uint8_t  alph_type;
  uint8_t  alph_size;
  uint8_t  charBits;
  uint32_t freq_SA; //frequency with which SA is sampled
  uint32_t freq_cnt_sb; //frequency with which full cumulative counts are captured
  uint32_t freq_cnt_b; //frequency with which intermittent counts are captured
  uint16_t block_count;
  uint32_t seq_count;
  uint64_t char_count; //total count of characters including those in and out of the alphabet
  char     *alph;
  char     *inv_alph;
  int      *compl_alph;
  FILE         *fp;
  FM_SEQDATA   *seq_data;
  FM_AMBIGLIST *ambig_list;
} FM_METADATA;

typedef struct fm_data_s {
  uint64_t N; //length of text
  uint32_t term_loc; // location in the BWT at which the '$' char is found (replaced in the sequence with 'a')
  uint32_t seq_offset;
  uint32_t ambig_offset;
  uint32_t seq_cnt;
  uint32_t ambig_cnt;
  uint32_t overlap; // number of bases at the beginning that overlap the FM-index for the preceding block
  uint8_t  *T;  //text corresponding to the BWT
  uint8_t  *BWT_mem;
  uint8_t  *BWT;
  uint32_t *SA; // sampled suffix array
  int64_t  *C; //the first position of each letter of the alphabet if all of T is sorted.  (signed, as I use that to keep tract of presence/absence)
  uint32_t *occCnts_sb;
  uint16_t *occCnts_b;
} FM_DATA;

typedef struct fm_dp_pair_s {
  uint16_t    pos;  // position of the diagonal in the model.
  float       score;
  float       max_score;
  uint8_t     score_peak_len; // how long was the diagonal when the most recent peak (within fm_drop_lim of the max score) was seen?
  uint8_t     consec_pos;
  uint8_t     max_consec_pos;
  uint8_t     consec_consensus;
  uint8_t     model_direction;
  uint8_t     complementarity;
} FM_DP_PAIR;

typedef struct fm_diag_s {
  uint64_t    n;  //position of the database sequence at which the diagonal starts
  union {
	double    sortkey;
	double    score;
  };
  uint16_t    k;  //position of the model at which the diagonal starts
  uint16_t    length;
  uint8_t     complementarity;
} FM_DIAG;

typedef struct fm_diaglist_s {
  FM_DIAG   *diags;
  int       count;
  int       size;
} FM_DIAGLIST;

/* Effectively global variables, to be initialized once in fm_initConfig(),
 * then passed around among threads to avoid recomputing them
 *
 * When allocated, must be 16-byte aligned, and all _m128i elements
 * must precede other types
 */
typedef struct {
#if   defined (p7_IMPL_SSE)
  /* mask arrays, and 16-byte-offsets into them */
  __m128i *fm_masks_mem;
  __m128i *fm_masks_v;
  __m128i *fm_reverse_masks_mem;
  __m128i *fm_reverse_masks_v;
  __m128i *fm_chars_mem;
  __m128i *fm_chars_v;

  /*various precomputed vectors*/
  __m128i fm_allones_v;
  __m128i fm_zeros_v;
  __m128i fm_neg128_v;
  __m128i fm_m0f;  //00 00 11 11
  __m128i fm_m01;  //01 01 01 01
  __m128i fm_m11;  //00 00 00 11

  /* no non-__m128i- elements above this line */
#endif //#if   defined (p7_IMPL_SSE)

  /*counter, to compute FM-index speed*/
  int occCallCnt;

  /*bounding cutoffs*/
  int max_depth;
  float drop_lim;  // 0.2 ; in seed, max drop in a run of length [fm_drop_max_len]
  int drop_max_len; // 4 ; maximum run length with score under (max - [fm_drop_lim])
  int consec_pos_req; //5
  int consensus_match_req; //10
  float score_density_req; //.5
  int ssv_length;
  float scthreshFM;
  float sc_thresh_ratio; //information content deficit,  actual_relent/target_relent

  /*pointer to FM-index metadata*/
  FM_METADATA *meta;

} FM_CFG;

#if   defined (p7_IMPL_SSE)
//used to convert from a byte array to an __m128i
typedef union {
		uint8_t bytes[16];
		__m128i m128;
		} byte_m128;

/* Gather the sum of all counts in a 16x8-bit element into a single 16-bit
 *  element of the register (the 0th element)
 *
 *  the _mm_sad_epu8  accumulates 8-bit counts into 16-bit counts:
 *      left 8 counts (64-bits) accumulate in counts_v[0],
 *      right 8 counts in counts_v[4]  (the other 6 16-bit ints are 0)
 *  the _mm_shuffle_epi32  flips the 4th int into the 0th slot
 */
#define FM_GATHER_8BIT_COUNTS( in_v, mid_v, out_v  ) do {\
	mid_v = _mm_sad_epu8 (in_v, cfg->fm_zeros_v);\
	tmp_v = _mm_shuffle_epi32(mid_v, _MM_SHUFFLE(1, 1, 1, 2));\
	out_v = _mm_add_epi16(mid_v, tmp_v);\
  } while (0)

/* Macro for SSE operations to turn 2-bit character values into 2-bit binary
 * (00 or 01) match/mismatch values representing occurrences of a character in a
 * 4-char-per-byte packed BWT.
 *
 * Typically followed by a call to FM_COUNT_SSE_4PACKED, possibly with a
 * mask in between to handle the case where we don't want to add over all
 * positions in the vector
 *
 * tmp_v and tmp2_v are used as temporary vectors throughout, and hold meaningless values
 * at the end
 *
 * xor(in_v, c_v)        : each 2-bit value will be 00 if a match, and non-0 if a mismatch
 * and(in_v, 01010101)   : look at the right bit of each 2-bit value,
 * srli(1)+and()         : look at the left bit of each 2-bit value,
 * or()                  : if either left bit or right bit is non-0, 01, else 00 (match is 00)
 *
 * subs()                : invert, so match is 01, mismatch is 00
 *
 */
#define FM_MATCH_2BIT(in_v, c_v, a_v, b_v, out_v) do {\
	a_v = _mm_xor_si128(in_v, c_v);\
	\
	b_v  = _mm_srli_epi16(a_v, 1);\
	a_v  = _mm_or_si128(a_v, b_v);\
	b_v  = _mm_and_si128(a_v, cfg->fm_m01);\
	\
	out_v  = _mm_subs_epi8(cfg->fm_m01,b_v);\
  } while (0)

/*Macro for SSE operations to count bits produced by FM_MATCH_SSE_4PACKED
 *
 * tmp_v and tmp2_v are used as temporary vectors throughout, and hold meaningless values
 * at the end
 *
 * then add up the 2-bit values:
 * srli(4)+add()         : left 4 bits shifted right, added to right 4 bits
 *
 * srli(2)+and(00000011) : left 2 bits (value 0..2) shifted right, masked, so no other bits active
 * and(00000011)         : right 2 bits (value 0..2) masked so no other bits active
 *
 * final 2 add()s        : tack current counts on to already-tabulated counts.
 */
#define FM_COUNT_2BIT(a_v, b_v, cnts_v) do {\
		b_v = _mm_srli_epi16(a_v, 4);\
		a_v  = _mm_add_epi16(a_v, b_v);\
		\
		b_v = _mm_srli_epi16(a_v, 2);\
		a_v  = _mm_and_si128(a_v,cfg->fm_m11);\
		b_v = _mm_and_si128(b_v,cfg->fm_m11);\
		\
		cnts_v = _mm_add_epi16(cnts_v, a_v);\
		cnts_v = _mm_add_epi16(cnts_v, b_v);\
  } while (0)

/* Macro for SSE operations that turns a vector of 4-bit character values into
 * 2 vectors representing matches. Each byte in the input vector consists of
 * a left half (4 bits) and a right half (4 bits). The 16 left-halves produce
 * one vector, which contains all-1s for bytes in which the left half matches
 * the c_v character (and 0s if it doesn't), while the 16 right-halves produce
 * the other vector, again with each byte either all-1s or all-0s.
 *
 * The expectation is that FM_COUNT_4BIT will be called after this, to
 * turn these binary values into sums over a series of vectors. The macros
 * are split up to allow one end or other to be trimmed in the case that
 * counting is not expected to include the full vector.
 *
 * srli(4)+and() : capture the left 4-bit value   (need the mask because 16-bit shift leaves garbage in left-4-bit chunks)
 * and()         : capture the right 4-bit value
 *
 * cmpeq()x2     : test if both left and right == c.  For each, if ==c , value = 11111111 (-1)
 */
#define FM_MATCH_4BIT(in_v, c_v, out1_v, out2_v) do {\
	out1_v    = _mm_srli_epi16(in_v, 4);\
	out2_v    = _mm_and_si128(in_v, cfg->fm_m0f);\
	out1_v    = _mm_and_si128(out1_v, cfg->fm_m0f);\
	\
	out1_v    = _mm_cmpeq_epi8(out1_v, c_v);\
	out2_v    = _mm_cmpeq_epi8(out2_v, c_v);\
  } while (0)

/* Macro for SSE operations that turns a vector of 4-bit character values into
 * 2 vectors representing matches. Each byte in the input vector consists of
 * a left half (4 bits) and a right half (4 bits). The 16 left-halves produce
 * one vector, which contains all-1s for bytes in which the left half is less than
 * the c_v character (and 0s if it doesn't), while the 16 right-halves produce
 * the other vector, again with each byte either all-1s or all-0s.
 *
 * The expectation is that FM_COUNT_4BIT will be called after this, to
 * turn these binary values into sums over a series of vectors. The macros
 * are split up to allow one end or other to be trimmed in the case that
 * counting is not expected to include the full vector.
 *
 * srli(4)+and() : capture the left 4-bit value   (need the mask because 16-bit shift leaves garbage in left-4-bit chunks)
 * and()         : capture the right 4-bit value
 *
 * cmplt()x2     : test if both left and right < c.  For each, if <c , value = 11111111 (-1)
 */
#define FM_LT_4BIT(in_v, c_v, out1_v, out2_v) do {\
	out1_v    = _mm_srli_epi16(in_v, 4);\
	out2_v    = _mm_and_si128(in_v, cfg->fm_m0f);\
	out1_v    = _mm_and_si128(out1_v, cfg->fm_m0f);\
	\
	out1_v    = _mm_cmplt_epi8(out1_v, c_v);\
	out2_v    = _mm_cmplt_epi8(out2_v, c_v);\
  } while (0)

/* Macro for SSE operations to add occurrence counts to the tally vector counts_v,
 * in the 4-bits-per-character case
 *
 * The expectation is that in[12]_v will contain bytes that are either
 *   00000000  =  0
 *  or
 *   11111111  = -1
 * so subtracting the value of the byte is the same as adding 0 or 1.
 */
#define FM_COUNT_4BIT(in1_v, in2_v, cnts_v) do {\
	cnts_v = _mm_subs_epi8(cnts_v, in1_v);\
	cnts_v = _mm_subs_epi8(cnts_v, in2_v);\
  } while (0)

#endif  // if  defined (p7_IMPL_SSE)

/*****************************************************************
 * 16. P7_PIPELINE: H3's accelerated seq/profile comparison pipeline
 *****************************************************************/

enum p7_pipemodes_e { p7_SEARCH_SEQS = 0, p7_SCAN_MODELS = 1 };
enum p7_zsetby_e    { p7_ZSETBY_NTARGETS = 0, p7_ZSETBY_OPTION = 1, p7_ZSETBY_FILEINFO = 2 };
enum p7_complementarity_e { p7_NOCOMPLEMENT    = 0, p7_COMPLEMENT   = 1 };

typedef struct p7_pipeline_s {
  /* Dynamic programming matrices                                           */
  P7_OMX     *oxf;		/* one-row Forward matrix, accel pipe       */
  P7_OMX     *oxb;		/* one-row Backward matrix, accel pipe      */
  P7_OMX     *fwd;		/* full Fwd matrix for domain envelopes     */
  P7_OMX     *bck;		/* full Bck matrix for domain envelopes     */

  /* Domain postprocessing                                                  */
  ESL_RANDOMNESS *r;		/* random number generator                  */
  int             do_reseeding; /* TRUE: reseed for reproducible results    */
  int             do_alignment_score_calc;
  P7_DOMAINDEF   *ddef;		/* domain definition workflow               */

  /* Reporting threshold settings                                           */
  int     by_E;		        /* TRUE to cut per-target report off by E   */
  double  E;	                /* per-target E-value threshold             */
  double  T;	                /* per-target bit score threshold           */
  int     dom_by_E;             /* TRUE to cut domain reporting off by E    */
  double  domE;	                /* domain E-value threshold                 */
  double  domT;	                /* domain bit score threshold               */
  int     use_bit_cutoffs;      /* (FALSE | p7H_GA | p7H_TC | p7H_NC)       */

  /* Inclusion threshold settings                                           */
  int     inc_by_E;		/* TRUE to threshold inclusion by E-values  */
  double  incE;			/* per-target inclusion E-value threshold   */
  double  incT;			/* per-target inclusion score threshold     */
  int     incdom_by_E;		/* TRUE to threshold domain inclusion by E  */
  double  incdomE;		/* per-domain inclusion E-value threshold   */
  double  incdomT;		/* per-domain inclusion E-value threshold   */

  /* Tracking search space sizes for E value calculations                   */
  double  Z;			/* eff # targs searched (per-target E-val)  */
  double  domZ;			/* eff # signific targs (per-domain E-val)  */
  enum p7_zsetby_e Z_setby;   	/* how Z was set                            */
  enum p7_zsetby_e domZ_setby;	/* how domZ was set                         */

  /* Threshold settings for pipeline                                        */
  int     do_max;	        /* TRUE to run in slow/max mode             */
  double  F1;		        /* MSV filter threshold                     */
  double  F2;		        /* Viterbi filter threshold                 */
  double  F3;		        /* uncorrected Forward filter threshold     */
  int     B1;               /* window length for biased-composition modifier - MSV*/
  int     B2;               /* window length for biased-composition modifier - Viterbi*/
  int     B3;               /* window length for biased-composition modifier - Forward*/
  int     do_biasfilter;	/* TRUE to use biased comp HMM filter       */
  int     do_null2;		/* TRUE to use null2 score corrections      */

  /* Accounting. (reduceable in threaded/MPI parallel version)              */
  uint64_t      nmodels;        /* # of HMMs searched                       */
  uint64_t      nseqs;	        /* # of sequences searched                  */
  uint64_t      nres;	        /* # of residues searched                   */
  uint64_t      nnodes;	        /* # of model nodes searched                */
  uint64_t      n_past_msv;	/* # comparisons that pass MSVFilter()      */
  uint64_t      n_past_bias;	/* # comparisons that pass bias filter      */
  uint64_t      n_past_vit;	/* # comparisons that pass ViterbiFilter()  */
  uint64_t      n_past_fwd;	/* # comparisons that pass ForwardFilter()  */
  uint64_t      n_output;	    /* # alignments that make it to the final output (used for nhmmer) */
  uint64_t      pos_past_msv;	/* # positions that pass MSVFilter()  (used for nhmmer) */
  uint64_t      pos_past_bias;	/* # positions that pass bias filter  (used for nhmmer) */
  uint64_t      pos_past_vit;	/* # positions that pass ViterbiFilter()  (used for nhmmer) */
  uint64_t      pos_past_fwd;	/* # positions that pass ForwardFilter()  (used for nhmmer) */
  uint64_t      pos_output;	    /* # positions that make it to the final output (used for nhmmer) */

  enum p7_pipemodes_e mode;    	/* p7_SCAN_MODELS | p7_SEARCH_SEQS          */
  int           long_targets;   /* TRUE if the target sequences are expected to be very long (e.g. dna chromosome search in nhmmer) */
  int           strands;         /*  p7_STRAND_TOPONLY  | p7_STRAND_BOTTOMONLY |  p7_STRAND_BOTH */
  int 		    	W;              /* window length for nhmmer scan - essentially maximum length of model that we expect to find*/
  int           block_length;   /* length of overlapping blocks read in the multi-threaded variant (default MAX_RESIDUE_COUNT) */

  int           show_accessions;/* TRUE to output accessions not names      */
  int           show_alignments;/* TRUE to output alignments (default)      */

  P7_HMMFILE   *hfp;		/* COPY of open HMM database (if scan mode) */
  char          errbuf[eslERRBUFSIZE];
} P7_PIPELINE;

/*****************************************************************
 * 17. P7_BUILDER: pipeline for new HMM construction
 *****************************************************************/

#define p7_DEFAULT_WINDOW_BETA  1e-7

enum p7_archchoice_e { p7_ARCH_FAST = 0, p7_ARCH_HAND = 1 };
enum p7_wgtchoice_e  { p7_WGT_NONE  = 0, p7_WGT_GIVEN = 1, p7_WGT_GSC    = 2, p7_WGT_PB       = 3, p7_WGT_BLOSUM = 4 };
enum p7_effnchoice_e { p7_EFFN_NONE = 0, p7_EFFN_SET  = 1, p7_EFFN_CLUST = 2, p7_EFFN_ENTROPY = 3, p7_EFFN_ENTROPY_EXP = 4 };

typedef struct p7_builder_s {
  /* Model architecture                                                                            */
  enum p7_archchoice_e arch_strategy;    /* choice of model architecture determination algorithm   */
  float                symfrac;	         /* residue occ thresh for fast architecture determination */
  float                fragthresh;	 /* if L <= fragthresh*alen, seq is called a fragment      */

  /* Relative sequence weights                                                                     */
  enum p7_wgtchoice_e  wgt_strategy;     /* choice of relative sequence weighting algorithm        */
  double               wid;		 /* %id threshold for BLOSUM relative weighting            */

  /* Effective sequence number                                                                     */
  enum p7_effnchoice_e effn_strategy;    /* choice of effective seq # determination algorithm      */
  double               re_target;	 /* rel entropy target for effn eweighting, if set; or -1.0*/
  double               esigma;		 /* min total rel ent parameter for effn entropy weights   */
  double               eid;		 /* %id threshold for effn clustering                      */
  double               eset;		 /* effective sequence number, if --eset; or -1.0          */

  /* Run-to-run variation due to random number generation                                          */
  ESL_RANDOMNESS      *r;	         /* RNG for E-value calibration simulations                */
  int                  do_reseeding;	 /* TRUE to reseed, making results reproducible            */

  /* E-value parameter calibration                                                                 */
  int                  EmL;            	 /* length of sequences generated for MSV fitting          */
  int                  EmN;	         /* # of sequences generated for MSV fitting               */
  int                  EvL;            	 /* length of sequences generated for Viterbi fitting      */
  int                  EvN;	         /* # of sequences generated for Viterbi fitting           */
  int                  EfL;	         /* length of sequences generated for Forward fitting      */
  int                  EfN;	         /* # of sequences generated for Forward fitting           */
  double               Eft;	         /* tail mass used for Forward fitting                     */

  /* Choice of prior                                                                               */
  P7_PRIOR            *prior;	         /* choice of prior when parameterizing from counts        */
  int                  max_insert_len;

  /* Optional: information used for parameterizing single sequence queries                         */
  ESL_SCOREMATRIX     *S;		 /* residue score matrix                                   */
  ESL_DMATRIX         *Q;	         /* Q->mx[a][b] = P(b|a) residue probabilities             */
  double               popen;         	 /* gap open probability                                   */
  double               pextend;          /* gap extend probability                                 */

  double               w_beta;    /*beta value used to compute W (window length)   */
  int                  w_len;     /*W (window length)  explicitly set */

  const ESL_ALPHABET  *abc;		 /* COPY of alphabet                                       */
  char errbuf[eslERRBUFSIZE];            /* informative message on model construction failure      */
} P7_BUILDER;

/*****************************************************************
 * 18. Routines in HMMER's exposed API.
 *****************************************************************/

/* build.c */
int p7_Handmodelmaker(ESL_MSA *msa,                P7_BUILDER *bld, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);
int p7_Fastmodelmaker(ESL_MSA *msa, float symfrac, P7_BUILDER *bld, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);

/* emit.c */
int p7_CoreEmit   (ESL_RANDOMNESS *r, const P7_HMM *hmm,                                        ESL_SQ *sq, P7_TRACE *tr);
int p7_ProfileEmit(ESL_RANDOMNESS *r, const P7_HMM *hmm, const P7_PROFILE *gm, const P7_BG *bg, ESL_SQ *sq, P7_TRACE *tr);
int p7_emit_SimpleConsensus(const P7_HMM *hmm, ESL_SQ *sq);
int p7_emit_FancyConsensus (const P7_HMM *hmm, float min_lower, float min_upper, ESL_SQ *sq);

/* errors.c */
void p7_Die (char *format, ...);
void p7_Fail(char *format, ...);

/* evalues.c */
int p7_Calibrate(P7_HMM *hmm, P7_BUILDER *cfg_b, ESL_RANDOMNESS **byp_rng, P7_BG **byp_bg, P7_PROFILE **byp_gm, P7_OPROFILE **byp_om);
int p7_Lambda(P7_HMM *hmm, P7_BG *bg, double *ret_lambda);
int p7_MSVMu     (ESL_RANDOMNESS *r, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda,               double *ret_mmu);
int p7_ViterbiMu (ESL_RANDOMNESS *r, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda,               double *ret_vmu);
int p7_Tau       (ESL_RANDOMNESS *r, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda, double tailp, double *ret_tau);

/* eweight.c */
int p7_EntropyWeight(const P7_HMM *hmm, const P7_BG *bg, const P7_PRIOR *pri, double infotarget, double *ret_Neff);

int p7_EntropyWeight_exp(const P7_HMM *hmm, const P7_BG *bg, const P7_PRIOR *pri, double etarget, double *ret_exp);
/* generic_decoding.c */
int p7_GDecoding      (const P7_PROFILE *gm, const P7_GMX *fwd,       P7_GMX *bck, P7_GMX *pp);
int p7_GDomainDecoding(const P7_PROFILE *gm, const P7_GMX *fwd, const P7_GMX *bck, P7_DOMAINDEF *ddef);

/* generic_fwdback.c */
int p7_GForward     (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *ret_sc);
int p7_GBackward    (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *ret_sc);
int p7_GHybrid      (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *opt_fwdscore, float *opt_hybscore);

/* generic_msv.c */
int p7_GMSV           (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *gx, float nu, float *ret_sc);
int p7_GMSV_longtarget(const ESL_DSQ *dsq, int L, P7_PROFILE *gm, P7_GMX *gx, float nu,  P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist);

/* generic_null2.c */
int p7_GNull2_ByExpectation(const P7_PROFILE *gm, P7_GMX *pp, float *null2);
int p7_GNull2_ByTrace      (const P7_PROFILE *gm, const P7_TRACE *tr, int zstart, int zend, P7_GMX *wrk, float *null2);

/* generic_optacc.c */
int p7_GOptimalAccuracy(const P7_PROFILE *gm, const P7_GMX *pp,       P7_GMX *gx, float *ret_e);
int p7_GOATrace        (const P7_PROFILE *gm, const P7_GMX *pp, const P7_GMX *gx, P7_TRACE *tr);

/* generic_stotrace.c */
int p7_GStochasticTrace(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *gx, P7_TRACE *tr);

/* generic_viterbi.c */
int p7_GViterbi     (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *ret_sc);

/* generic_vtrace.c */
int p7_GTrace       (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *gx, P7_TRACE *tr);
int p7_GViterbi_longtarget(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *gx,
					   float filtersc, double P, P7_HMM_WINDOWLIST *windowlist);

/* heatmap.c (evolving now, intend to move this to Easel in the future) */
double dmx_upper_max(ESL_DMATRIX *D);
double dmx_upper_min(ESL_DMATRIX *D);
double dmx_upper_element_sum(ESL_DMATRIX *D);
double dmx_upper_norm(ESL_DMATRIX *D);
int    dmx_Visualize(FILE *fp, ESL_DMATRIX *D, double min, double max);

/* hmmdutils.c */
void p7_openlog(const char *ident, int option, int facility);
void p7_syslog(int priority, const char *format, ...);
void p7_closelog(void);

/* hmmlogo.c */
float hmmlogo_maxHeight (P7_BG *bg);
int hmmlogo_RelativeEntropy_all      (P7_HMM *hmm, P7_BG *bg, float *rel_ents, float **probs, float **heights );
int hmmlogo_RelativeEntropy_above_bg (P7_HMM *hmm, P7_BG *bg, float *rel_ents, float **probs, float **heights );
int hmmlogo_ScoreHeights (P7_HMM *hmm, P7_BG *bg, float **heights );
int hmmlogo_IndelValues (P7_HMM *hmm, float *insert_P, float *insert_expL, float *delete_P );

/* hmmpgmd2msa.c */
int hmmpgmd2msa(void *data, P7_HMM *hmm, ESL_SQ *qsq,  int *incl, int incl_size, int *excl, int excl_size, ESL_MSA **ret_msa);

/* island.c */
int   p7_island_Viterbi(ESL_DSQ *dsq, int L, P7_PROFILE *gm, P7_GMX *mx, ESL_HISTOGRAM *h);

/* h2_io.c */
int   p7_h2io_WriteASCII(FILE *fp, P7_HMM *hmm);

/* hmmer.c */
void         p7_banner(FILE *fp, char *progname, char *banner);
ESL_GETOPTS *p7_CreateDefaultApp(ESL_OPTIONS *options, int nargs, int argc, char **argv, char *banner, char *usage);
int          p7_AminoFrequencies(float *f);

/* logsum.c */
int   p7_FLogsumInit(void);
float p7_FLogsum(float a, float b);
float p7_FLogsumError(float a, float b);
int   p7_ILogsumInit(void);
int   p7_ILogsum(int s1, int s2);

/* modelconfig.c */
int p7_ProfileConfig(const P7_HMM *hmm, const P7_BG *bg, P7_PROFILE *gm, int L, int mode);
int p7_ReconfigLength  (P7_PROFILE *gm, int L);
int p7_ReconfigMultihit(P7_PROFILE *gm, int L);
int p7_ReconfigUnihit  (P7_PROFILE *gm, int L);

/* modelstats.c */
double p7_MeanMatchInfo           (const P7_HMM *hmm, const P7_BG *bg);
double p7_MeanMatchEntropy        (const P7_HMM *hmm);
double p7_MeanMatchRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg);
double p7_MeanForwardScore        (const P7_HMM *hmm, const P7_BG *bg);
int    p7_MeanPositionRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg, double *ret_entropy);
int    p7_hmm_CompositionKLDist(P7_HMM *hmm, P7_BG *bg, float *ret_KL, float **opt_avp);

/* mpisupport.c */
#ifdef HAVE_MPI
int p7_hmm_MPISend(P7_HMM *hmm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
int p7_hmm_MPIPackSize(P7_HMM *hmm, MPI_Comm comm, int *ret_n);
int p7_hmm_MPIPack(P7_HMM *hmm, char *buf, int n, int *position, MPI_Comm comm);
int p7_hmm_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **abc, P7_HMM **ret_hmm);
int p7_hmm_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **abc, P7_HMM **ret_hmm);

int p7_profile_MPISend(P7_PROFILE *gm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
int p7_profile_MPIRecv(int source, int tag, MPI_Comm comm, const ESL_ALPHABET *abc, const P7_BG *bg,
			      char **buf, int *nalloc,  P7_PROFILE **ret_gm);

int p7_pipeline_MPISend(P7_PIPELINE *pli, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
int p7_pipeline_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_GETOPTS *go, P7_PIPELINE **ret_pli);

int p7_tophits_MPISend(P7_TOPHITS *th, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
int p7_tophits_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, P7_TOPHITS **ret_th);

int p7_oprofile_MPISend(P7_OPROFILE *om, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
int p7_oprofile_MPIPackSize(P7_OPROFILE *om, MPI_Comm comm, int *ret_n);
int p7_oprofile_MPIPack(P7_OPROFILE *om, char *buf, int n, int *pos, MPI_Comm comm);
int p7_oprofile_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **abc, P7_OPROFILE **ret_om);
int p7_oprofile_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **abc, P7_OPROFILE **ret_om);
#endif /*HAVE_MPI*/

/* tracealign.c */
int p7_tracealign_Seqs(ESL_SQ **sq,           P7_TRACE **tr, int nseq, int M,  int optflags, P7_HMM *hmm, ESL_MSA **ret_msa);
int p7_tracealign_MSA (const ESL_MSA *premsa, P7_TRACE **tr,           int M,  int optflags, ESL_MSA **ret_postmsa);
int p7_tracealign_computeTraces(P7_HMM *hmm, ESL_SQ  **sq, int offset, int N, P7_TRACE  **tr);
int p7_tracealign_getMSAandStats(P7_HMM *hmm, ESL_SQ  **sq, int N, ESL_MSA **ret_msa, float **ret_pp, float **ret_relent, float **ret_scores );

/* p7_alidisplay.c */
P7_ALIDISPLAY *p7_alidisplay_Create(const P7_TRACE *tr, int which, const P7_OPROFILE *om, const ESL_SQ *sq, const ESL_SQ *ntsq);
P7_ALIDISPLAY *p7_alidisplay_Clone(const P7_ALIDISPLAY *ad);
size_t         p7_alidisplay_Sizeof(const P7_ALIDISPLAY *ad);
int            p7_alidisplay_Serialize(P7_ALIDISPLAY *ad);
int            p7_alidisplay_Deserialize(P7_ALIDISPLAY *ad);
void           p7_alidisplay_Destroy(P7_ALIDISPLAY *ad);
char           p7_alidisplay_EncodePostProb(float p);
float          p7_alidisplay_DecodePostProb(char pc);
char           p7_alidisplay_EncodeAliPostProb(float p, float hi, float med, float lo);

int            p7_alidisplay_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth, P7_PIPELINE *pli);
int            p7_translated_alidisplay_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth, P7_PIPELINE *pli);
int            p7_nontranslated_alidisplay_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth, int show_accessions);

int            p7_alidisplay_Backconvert(const P7_ALIDISPLAY *ad, const ESL_ALPHABET *abc, ESL_SQ **ret_sq, P7_TRACE **ret_tr);
int            p7_alidisplay_Dump(FILE *fp, const P7_ALIDISPLAY *ad);
int            p7_alidisplay_Compare(const P7_ALIDISPLAY *ad1, const P7_ALIDISPLAY *ad2);

/* p7_bg.c */
P7_BG *p7_bg_Create(const ESL_ALPHABET *abc);
P7_BG *p7_bg_CreateUniform(const ESL_ALPHABET *abc);
P7_BG *p7_bg_Clone(const P7_BG *bg);
int    p7_bg_Dump(FILE *ofp, const P7_BG *bg);
void   p7_bg_Destroy(P7_BG *bg);
int    p7_bg_SetLength(P7_BG *bg, int L);
int    p7_bg_NullOne(const P7_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc);

int    p7_bg_Read(char *bgfile, P7_BG *bg, char *errbuf);
int    p7_bg_Write(FILE *fp, P7_BG *bg);

int    p7_bg_SetFilter  (P7_BG *bg, int M, const float *compo);
int    p7_bg_FilterScore(P7_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc);

/* p7_builder.c */
P7_BUILDER *p7_builder_Create(const ESL_GETOPTS *go, const ESL_ALPHABET *abc);
int         p7_builder_LoadScoreSystem(P7_BUILDER *bld, const char *matrix,                  double popen, double pextend, P7_BG *bg);
int         p7_builder_SetScoreSystem (P7_BUILDER *bld, const char *mxfile, const char *env, double popen, double pextend, P7_BG *bg);
void        p7_builder_Destroy(P7_BUILDER *bld);

int p7_Builder      (P7_BUILDER *bld, ESL_MSA *msa, P7_BG *bg, P7_HMM **opt_hmm, P7_TRACE ***opt_trarr, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om, ESL_MSA **opt_postmsa, FILE *seqweights_w_fp, FILE *seqweights_e_fp);
int p7_SingleBuilder(P7_BUILDER *bld, ESL_SQ *sq,   P7_BG *bg, P7_HMM **opt_hmm, P7_TRACE  **opt_tr,    P7_PROFILE **opt_gm, P7_OPROFILE **opt_om);
int p7_Builder_MaxLength      (P7_HMM *hmm, double emit_thresh);

/* p7_domaindef.c */
P7_DOMAINDEF *p7_domaindef_Create (ESL_RANDOMNESS *r);
int           p7_domaindef_Fetch  (P7_DOMAINDEF *ddef, int which, int *opt_i, int *opt_j, float *opt_sc, P7_ALIDISPLAY **opt_ad);
int           p7_domaindef_GrowTo (P7_DOMAINDEF *ddef, int L);
int           p7_domaindef_Reuse  (P7_DOMAINDEF *ddef);
int           p7_domaindef_DumpPosteriors(FILE *ofp, P7_DOMAINDEF *ddef);
void          p7_domaindef_Destroy(P7_DOMAINDEF *ddef);

int p7_domaindef_ByViterbi            (P7_PROFILE *gm, const ESL_SQ *sq, const ESL_SQ *ntsq, P7_GMX *gx1, P7_GMX *gx2, P7_DOMAINDEF *ddef);
int p7_domaindef_ByPosteriorHeuristics(const ESL_SQ *sq, const ESL_SQ *ntsq, P7_OPROFILE *om, P7_OMX *oxf, P7_OMX *oxb, P7_OMX *fwd, P7_OMX *bck,
				                                  P7_DOMAINDEF *ddef, P7_BG *bg, int long_target,
				                                  P7_BG *bg_tmp, float *scores_arr, float *fwd_emissions_arr);

/* p7_gmx.c */
P7_GMX *p7_gmx_Create (int allocM, int allocL);
int     p7_gmx_GrowTo (P7_GMX *gx, int allocM, int allocL);
size_t  p7_gmx_Sizeof (P7_GMX *gx);
int     p7_gmx_Reuse  (P7_GMX *gx);
void    p7_gmx_Destroy(P7_GMX *gx);
int     p7_gmx_Compare(P7_GMX *gx1, P7_GMX *gx2, float tolerance);
int     p7_gmx_Dump(FILE *fp, P7_GMX *gx, int flags);
int     p7_gmx_DumpWindow(FILE *fp, P7_GMX *gx, int istart, int iend, int kstart, int kend, int show_specials);

/* p7_hmm.c */
/*      1. The P7_HMM object: allocation, initialization, destruction. */
P7_HMM *p7_hmm_Create(int M, const ESL_ALPHABET *abc);
P7_HMM *p7_hmm_CreateShell(void);
int     p7_hmm_CreateBody(P7_HMM *hmm, int M, const ESL_ALPHABET *abc);
void    p7_hmm_Destroy(P7_HMM *hmm);
int     p7_hmm_CopyParameters(const P7_HMM *src, P7_HMM *dest);
P7_HMM *p7_hmm_Clone(const P7_HMM *hmm);
int     p7_hmm_Zero(P7_HMM *hmm);
char    p7_hmm_EncodeStatetype(char *typestring);
char   *p7_hmm_DecodeStatetype(char st);
/*      2. Convenience routines for setting fields in an HMM. */
int     p7_hmm_SetName       (P7_HMM *hmm, char *name);
int     p7_hmm_SetAccession  (P7_HMM *hmm, char *acc);
int     p7_hmm_SetDescription(P7_HMM *hmm, char *desc);
int     p7_hmm_AppendComlog  (P7_HMM *hmm, int argc, char **argv);
int     p7_hmm_SetCtime      (P7_HMM *hmm);
int     p7_hmm_SetComposition(P7_HMM *hmm);
int     p7_hmm_SetConsensus  (P7_HMM *hmm, ESL_SQ *sq);
/*      3. Renormalization and rescaling counts in core HMMs. */
int     p7_hmm_Scale      (P7_HMM *hmm, double scale);
int     p7_hmm_ScaleExponential(P7_HMM *hmm, double exp);
int     p7_hmm_Renormalize(P7_HMM *hmm);
/*      4. Debugging and development code. */
int     p7_hmm_Dump(FILE *fp, P7_HMM *hmm);
int     p7_hmm_Sample          (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm);
int     p7_hmm_SampleUngapped  (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm);
int     p7_hmm_SampleEnumerable(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm);
int     p7_hmm_SampleUniform   (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc,
				     float tmi, float tii, float tmd, float tdd,  P7_HMM **ret_hmm);
int     p7_hmm_Compare(P7_HMM *h1, P7_HMM *h2, float tol);
int     p7_hmm_Validate(P7_HMM *hmm, char *errbuf, float tol);
/*      5. Other routines in the API */
int     p7_hmm_CalculateOccupancy(const P7_HMM *hmm, float *mocc, float *iocc);

/* p7_hmmfile.c */
int  p7_hmmfile_OpenE    (char *filename, char *env, P7_HMMFILE **ret_hfp, char *errbuf);
int  p7_hmmfile_OpenENoDB(char *filename, char *env, P7_HMMFILE **ret_hfp, char *errbuf);
int  p7_hmmfile_Open     (char *filename, char *env, P7_HMMFILE **ret_hfp); /* deprecated */
int  p7_hmmfile_OpenNoDB (char *filename, char *env, P7_HMMFILE **ret_hfp); /* deprecated */
int  p7_hmmfile_OpenBuffer(char *buffer, int size, P7_HMMFILE **ret_hfp);
void p7_hmmfile_Close(P7_HMMFILE *hfp);
#ifdef HMMER_THREADS
int  p7_hmmfile_CreateLock(P7_HMMFILE *hfp);
#endif
int  p7_hmmfile_WriteBinary(FILE *fp, int format, P7_HMM *hmm);
int  p7_hmmfile_WriteASCII (FILE *fp, int format, P7_HMM *hmm);
int  p7_hmmfile_WriteToString (char **s, int format, P7_HMM *hmm);
int  p7_hmmfile_Read(P7_HMMFILE *hfp, ESL_ALPHABET **ret_abc,  P7_HMM **opt_hmm);
int  p7_hmmfile_PositionByKey(P7_HMMFILE *hfp, const char *key);
int  p7_hmmfile_Position(P7_HMMFILE *hfp, const off_t offset);

/* p7_hmmwindow.c */
int p7_hmmwindow_init (P7_HMM_WINDOWLIST *list);
P7_HMM_WINDOW *p7_hmmwindow_new (P7_HMM_WINDOWLIST *list, uint32_t id, uint32_t pos, uint32_t fm_pos, uint16_t k, uint32_t length, float score, uint8_t complementarity, uint32_t target_len);

/* p7_msvdata.c */
P7_SCOREDATA   *p7_hmm_ScoreDataCreate(P7_OPROFILE *om, P7_PROFILE *gm );
P7_SCOREDATA   *p7_hmm_ScoreDataClone(P7_SCOREDATA *src, int K);
int            p7_hmm_ScoreDataComputeRest(P7_OPROFILE *om, P7_SCOREDATA *data );
void           p7_hmm_ScoreDataDestroy( P7_SCOREDATA *data );
int            p7_hmm_initWindows (P7_HMM_WINDOWLIST *list);
P7_HMM_WINDOW *p7_hmm_newWindow (P7_HMM_WINDOWLIST *list, uint32_t id, uint32_t pos, uint32_t fm_pos, uint16_t k, uint32_t length, float score, uint8_t complementarity);

/* p7_null3.c */
void p7_null3_score(const ESL_ALPHABET *abc, const ESL_DSQ *dsq, P7_TRACE *tr, int start, int stop, P7_BG *bg, float *ret_sc);
void p7_null3_windowed_score(const ESL_ALPHABET *abc, const ESL_DSQ *dsq, int start, int stop, P7_BG *bg, float *ret_sc);

/* p7_pipeline.c */
P7_PIPELINE *p7_pipeline_Create(ESL_GETOPTS *go, int M_hint, int L_hint, int long_targets, enum p7_pipemodes_e mode);
int          p7_pipeline_Reuse  (P7_PIPELINE *pli);
void         p7_pipeline_Destroy(P7_PIPELINE *pli);
int          p7_pipeline_Merge  (P7_PIPELINE *p1, P7_PIPELINE *p2);

int p7_pli_ExtendAndMergeWindows (P7_OPROFILE *om, const P7_SCOREDATA *msvdata, P7_HMM_WINDOWLIST *windowlist, float pct_overlap);
int p7_pli_TargetReportable  (P7_PIPELINE *pli, float score,     double lnP);
int p7_pli_DomainReportable  (P7_PIPELINE *pli, float dom_score, double lnP);

int p7_pli_TargetIncludable  (P7_PIPELINE *pli, float score,     double lnP);
int p7_pli_DomainIncludable  (P7_PIPELINE *pli, float dom_score, double lnP);
int p7_pli_NewModel          (P7_PIPELINE *pli, const P7_OPROFILE *om, P7_BG *bg);
int p7_pli_NewModelThresholds(P7_PIPELINE *pli, const P7_OPROFILE *om);
int p7_pli_NewSeq            (P7_PIPELINE *pli, const ESL_SQ *sq);
int p7_Pipeline              (P7_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, const ESL_SQ *ntsq, P7_TOPHITS *th);
int p7_Pipeline_LongTarget   (P7_PIPELINE *pli, P7_OPROFILE *om, P7_SCOREDATA *data,
									 P7_BG *bg, P7_TOPHITS *hitlist, int64_t seqidx,
									 const ESL_SQ *sq, int complementarity,
									 const FM_DATA *fmf, const FM_DATA *fmb, FM_CFG *fm_cfg
/*                                     , ESL_STOPWATCH *ssv_watch_master
									 , ESL_STOPWATCH *postssv_watch_master
									 , ESL_STOPWATCH *watch_slave
*/
									 );

int p7_pli_Statistics(FILE *ofp, P7_PIPELINE *pli, ESL_STOPWATCH *w);

/* p7_prior.c */
P7_PRIOR  *p7_prior_CreateAmino(void);
P7_PRIOR  *p7_prior_CreateNucleic(void);
P7_PRIOR  *p7_prior_CreateLaplace(const ESL_ALPHABET *abc);
void       p7_prior_Destroy(P7_PRIOR *pri);

int        p7_ParameterEstimation(P7_HMM *hmm, const P7_PRIOR *pri);

/* p7_profile.c */
P7_PROFILE *p7_profile_Create(int M, const ESL_ALPHABET *abc);
P7_PROFILE *p7_profile_Clone(const P7_PROFILE *gm);
int         p7_profile_Copy(const P7_PROFILE *src, P7_PROFILE *dst);
int         p7_profile_SetNullEmissions(P7_PROFILE *gm);
int         p7_profile_Reuse(P7_PROFILE *gm);
size_t      p7_profile_Sizeof(P7_PROFILE *gm);
void        p7_profile_Destroy(P7_PROFILE *gm);
int         p7_profile_IsLocal(const P7_PROFILE *gm);
int         p7_profile_IsMultihit(const P7_PROFILE *gm);
int         p7_profile_GetT(const P7_PROFILE *gm, char st1, int k1,
				   char st2, int k2, float *ret_tsc);
int         p7_profile_Validate(const P7_PROFILE *gm, char *errbuf, float tol);
int         p7_profile_Compare(P7_PROFILE *gm1, P7_PROFILE *gm2, float tol);

/* p7_spensemble.c */
P7_SPENSEMBLE *p7_spensemble_Create(int init_n, int init_epc, int init_sigc);
int     p7_spensemble_Reuse(P7_SPENSEMBLE *sp);
int     p7_spensemble_Add(P7_SPENSEMBLE *sp, int sampleidx, int i, int j, int k, int m);
int     p7_spensemble_Cluster(P7_SPENSEMBLE *sp,
				     float min_overlap, int of_smaller, int max_diagdiff,
				     float min_posterior, float min_endpointp,
				     int *ret_nclusters);
int     p7_spensemble_GetClusterCoords(P7_SPENSEMBLE *sp, int which,
					      int *ret_i, int *ret_j, int *ret_k, int *ret_m, float *ret_p);
void    p7_spensemble_Destroy(P7_SPENSEMBLE *sp);

/* p7_tophits.c */
P7_TOPHITS *p7_tophits_Create(void);
int         p7_tophits_Grow(P7_TOPHITS *h);
int         p7_tophits_CreateNextHit(P7_TOPHITS *h, P7_HIT **ret_hit);
int         p7_tophits_Add(P7_TOPHITS *h,
				  char *name, char *acc, char *desc,
				  double sortkey,
				  float score,    double lnP,
				  float mothersc, double mother_lnP,
				  int sqfrom, int sqto, int sqlen,
				  int hmmfrom, int hmmto, int hmmlen,
				  int domidx, int ndom,
				  P7_ALIDISPLAY *ali);
int         p7_tophits_SortBySortkey(P7_TOPHITS *h);
int         p7_tophits_SortBySeqidxAndAlipos(P7_TOPHITS *h);
int         p7_tophits_SortByModelnameAndAlipos(P7_TOPHITS *h);

int         p7_tophits_Merge(P7_TOPHITS *h1, P7_TOPHITS *h2);
int         p7_tophits_GetMaxPositionLength(P7_TOPHITS *h);
int         p7_tophits_GetMaxNameLength(P7_TOPHITS *h);
int         p7_tophits_GetMaxAccessionLength(P7_TOPHITS *h);
int         p7_tophits_GetMaxShownLength(P7_TOPHITS *h);
void        p7_tophits_Destroy(P7_TOPHITS *h);

int p7_tophits_ComputeNhmmerEvalues(P7_TOPHITS *th, double N, int W);
int p7_tophits_RemoveDuplicates(P7_TOPHITS *th, int using_bit_cutoffs);
int p7_tophits_Threshold(P7_TOPHITS *th, P7_PIPELINE *pli);
int p7_tophits_CompareRanking(P7_TOPHITS *th, ESL_KEYHASH *kh, int *opt_nnew);
int p7_tophits_Targets(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw);
int p7_tophits_Domains(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw);

int p7_tophits_Alignment(const P7_TOPHITS *th, const ESL_ALPHABET *abc,
				ESL_SQ **inc_sqarr, P7_TRACE **inc_trarr, int inc_n, int optflags,
				ESL_MSA **ret_msa);
int p7_tophits_TabularTargets(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header);
int p7_tophits_TabularDomains(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header);
int p7_tophits_TabularXfam(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli);
int p7_tophits_TabularTail(FILE *ofp, const char *progname, enum p7_pipemodes_e pipemode,
				  const char *qfile, const char *tfile, const ESL_GETOPTS *go);
int p7_tophits_AliScores(FILE *ofp, char *qname, P7_TOPHITS *th );

/* p7_trace.c */
P7_TRACE *p7_trace_Create(void);
P7_TRACE *p7_trace_CreateWithPP(void);
int  p7_trace_Reuse(P7_TRACE *tr);
int  p7_trace_Grow(P7_TRACE *tr);
int  p7_trace_GrowIndex(P7_TRACE *tr);
int  p7_trace_GrowTo(P7_TRACE *tr, int N);
int  p7_trace_GrowIndexTo(P7_TRACE *tr, int ndom);
void p7_trace_Destroy(P7_TRACE *tr);
void p7_trace_DestroyArray(P7_TRACE **tr, int N);

int  p7_trace_GetDomainCount   (const P7_TRACE *tr, int *ret_ndom);
int  p7_trace_GetStateUseCounts(const P7_TRACE *tr, int *counts);
int  p7_trace_GetDomainCoords  (const P7_TRACE *tr, int which, int *ret_i1, int *ret_i2,
				       int *ret_k1, int *ret_k2);

int   p7_trace_Validate(const P7_TRACE *tr, const ESL_ALPHABET *abc, const ESL_DSQ *dsq, char *errbuf);
int   p7_trace_Dump(FILE *fp, const P7_TRACE *tr, const P7_PROFILE *gm, const ESL_DSQ *dsq);
int   p7_trace_Compare(P7_TRACE *tr1, P7_TRACE *tr2, float pptol);
int   p7_trace_Score(P7_TRACE *tr, ESL_DSQ *dsq, P7_PROFILE *gm, float *ret_sc);
int   p7_trace_SetPP(P7_TRACE *tr, const P7_GMX *pp);
float p7_trace_GetExpectedAccuracy(const P7_TRACE *tr);

int  p7_trace_Append(P7_TRACE *tr, char st, int k, int i);
int  p7_trace_AppendWithPP(P7_TRACE *tr, char st, int k, int i, float pp);
int  p7_trace_Reverse(P7_TRACE *tr);
int  p7_trace_Index(P7_TRACE *tr);

int  p7_trace_FauxFromMSA(ESL_MSA *msa, int *matassign, int optflags, P7_TRACE **tr);
int  p7_trace_Doctor(P7_TRACE *tr, int *opt_ndi, int *opt_nid);

int  p7_trace_Count(P7_HMM *hmm, ESL_DSQ *dsq, float wt, P7_TRACE *tr);

/* seqmodel.c */
int p7_Seqmodel(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int M, char *name,
		       ESL_DMATRIX *P, float *f, double popen, double pextend,
		       P7_HMM **ret_hmm);

/* fm_alphabet.c */
int fm_alphabetCreate (FM_METADATA *meta, uint8_t *alph_bits);
int fm_alphabetDestroy (FM_METADATA *meta);
int fm_reverseString (char *str, int N);
int fm_getComplement (char c, uint8_t alph_type);

/* fm_general.c */
uint64_t fm_computeSequenceOffset (const FM_DATA *fms, const FM_METADATA *meta, int block, uint64_t pos);
int fm_getOriginalPosition (const FM_DATA *fms, const FM_METADATA *meta, int fm_id, int length, int direction, uint64_t fm_pos,
									uint32_t *segment_id, uint64_t *seg_pos);
int fm_readFMmeta( FM_METADATA *meta);
int fm_FM_read( FM_DATA *fm, FM_METADATA *meta, int getAll );
void fm_FM_destroy ( FM_DATA *fm, int isMainFM);
uint8_t fm_getChar(uint8_t alph_type, int j, const uint8_t *B );
int fm_getSARangeReverse( const FM_DATA *fm, FM_CFG *cfg, char *query, char *inv_alph, FM_INTERVAL *interval);
int fm_getSARangeForward( const FM_DATA *fm, FM_CFG *cfg, char *query, char *inv_alph, FM_INTERVAL *interval);
int fm_configAlloc(FM_CFG **cfg);
int fm_configDestroy(FM_CFG *cfg);
int fm_metaDestroy(FM_METADATA *meta );
int fm_updateIntervalForward( const FM_DATA *fm, const FM_CFG *cfg, char c, FM_INTERVAL *interval_f, FM_INTERVAL *interval_bk);
int fm_updateIntervalReverse( const FM_DATA *fm, const FM_CFG *cfg, char c, FM_INTERVAL *interval);
int fm_initSeeds (FM_DIAGLIST *list) ;
FM_DIAG * fm_newSeed (FM_DIAGLIST *list);
int fm_initAmbiguityList (FM_AMBIGLIST *list);
int fm_addAmbiguityRange (FM_AMBIGLIST *list, uint32_t start, uint32_t stop);
int fm_convertRange2DSQ(const FM_DATA *fm, const FM_METADATA *meta, uint64_t first, int length, int complementarity, ESL_SQ *sq, int fix_ambiguities );
int fm_initConfigGeneric( FM_CFG *cfg, ESL_GETOPTS *go);

/* fm_ssv.c */
int p7_SSVFM_longlarget( P7_OPROFILE *om, float nu, P7_BG *bg, double F1,
					  const FM_DATA *fmf, const FM_DATA *fmb, FM_CFG *fm_cfg, const P7_SCOREDATA *ssvdata,
					  int strands, P7_HMM_WINDOWLIST *windowlist);

/* fm_sse.c */
int fm_configInit      (FM_CFG *cfg, ESL_GETOPTS *go);
int fm_getOccCount     (const FM_DATA *fm, const FM_CFG *cfg, int pos, uint8_t c);
int fm_getOccCountLT   (const FM_DATA *fm, const FM_CFG *cfg, int pos, uint8_t c, uint32_t *cnteq, uint32_t *cntlt);


/*** Start of inlined file: cachedb.h ***/
#ifndef P7_CACHEDB_INCLUDED
#define P7_CACHEDB_INCLUDED

typedef struct {
  char    *name;                   /* name; ("\0" if no name)               */
  ESL_DSQ *dsq;                    /* digitized sequence [1..n]             */
  int64_t  n;                      /* length of dsq                         */
  int64_t  idx;	                   /* ctr for this seq                      */
  uint64_t db_key;                 /* flag for included databases           */
  char    *desc;                   /* description                           */
} HMMER_SEQ;

typedef struct {
  uint32_t            count;       /* number of entries                     */
  uint32_t            K;           /* original number of entries            */
  HMMER_SEQ         **list;        /* list of sequences [0 .. count-1]      */
} SEQ_DB;

typedef struct {
  char               *name;        /* name of the seq database              */
  char               *id;          /* unique identifier string              */
  uint32_t            db_cnt;      /* number of sub databases               */
  SEQ_DB             *db;          /* list of databases [0 .. db_cnt-1]     */

  ESL_ALPHABET       *abc;         /* alphabet for database                 */

  uint32_t            count;       /* total number of sequences             */
  HMMER_SEQ          *list;        /* complete list of sequences (count)    */
  void               *residue_mem; /* memory holding the residues           */
  char               *header_mem;  /* memory holding the header strings     */

  uint64_t            res_size;    /* size of residue memory allocation     */
  uint64_t            hdr_size;    /* size of header memory allocation      */
} P7_SEQCACHE;

int    p7_seqcache_Open(char *seqfile, P7_SEQCACHE **ret_cache, char *errbuf);
void   p7_seqcache_Close(P7_SEQCACHE *cache);

#endif /*P7_CACHEDB_INCLUDED*/

/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Version 3.1b3; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 *
 * HMMER is distributed under the terms of the GNU General Public License
 * (GPLv3). See the LICENSE file for details.
 ************************************************************/

/*** End of inlined file: cachedb.h ***/


/*** Start of inlined file: p7_gbands.h ***/
#ifndef P7_GBANDS_INCLUDED
#define P7_GBANDS_INCLUDED

typedef struct {
  int     nseg;
  int     nrow;
  int     L;
  int     M;
  int64_t ncell;

  int *imem;
  int *kmem;

  int  segalloc;
  int  rowalloc;
} P7_GBANDS;

#define p7_GBANDS_NK 2		/* for IK banding. Or 6, for IKS banding */

P7_GBANDS *p7_gbands_Create  (void);
int        p7_gbands_Reuse   (P7_GBANDS *bnd);
int        p7_gbands_Append  (P7_GBANDS *bnd, int i, int ka, int kb);
int        p7_gbands_Prepend (P7_GBANDS *bnd, int i, int ka, int kb);
int        p7_gbands_Reverse (P7_GBANDS *bnd);
int        p7_gbands_GrowSegs(P7_GBANDS *bnd);
int        p7_gbands_GrowRows(P7_GBANDS *bnd);
void       p7_gbands_Destroy (P7_GBANDS *bnd);
int        p7_gbands_Dump(FILE *ofp, P7_GBANDS *bnd);

#endif /*P7_GBANDS_INCLUDED*/

/*** End of inlined file: p7_gbands.h ***/


/*** Start of inlined file: p7_gmxb.h ***/
#ifndef P7_GMXB_INCLUDED
#define P7_GMXB_INCLUDED

typedef struct {
  float     *dp;
  float     *xmx;
  P7_GBANDS *bnd;   /* a reference copy; caller remains responsible for free'ing banding */

  int64_t    dalloc;
  int        xalloc;
} P7_GMXB;

P7_GMXB *p7_gmxb_Create(P7_GBANDS *bnd);
int      p7_gmxb_Reinit(P7_GMXB *gxb, P7_GBANDS *bnd);
int      p7_gmxb_Reuse(P7_GMXB *gxb);
void     p7_gmxb_Destroy(P7_GMXB *gxb);
int      p7_gmxb_Dump(FILE *ofp, P7_GMXB *gxb, int flags);

int p7_GForwardBanded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *gxb, float *opt_sc);

#endif /*P7_GMXB_INCLUDED*/

/*** End of inlined file: p7_gmxb.h ***/


/*** Start of inlined file: p7_gmxchk.h ***/
#ifndef P7_GMXCHK_INCLUDED
#define P7_GMXCHK_INCLUDED

/*****************************************************************
 * 1. Exegesis: layout of rows in a checkpointed matrix.
 *****************************************************************/

/*
 * One P7_GMXCHK data structure is used for both Forward and Backward
 * computations on a target sequence. The Forward calculation is
 * checkpointed. The Backward calculation is linear memory in two
 * rows. The end result is a Forward score and a posterior-decoded set
 * of DP bands.
 *
 * In the diagram below, showing the row layout for the main matrix (MDI states):
 *   O = a checkpointed row;
 *   x = row that isn't checkpointed;
 *   * = boundary row 0, plus row(s) used for Backwards
 *
 *   i = index of residues in a target sequence of length L
 *   r = index of rows in the DP matrix, R0+R in total
 *
 *               |------------------------- L -------------------------------|
 *               |-----La----| |-Lb-| |-------------- Lc --------------------|
 * i =  .  .  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
 *      *  *  *  O  O  O  O  O  x  O  x  x  x  x  O  x  x  x  O  x  x  O  x  O
 * r =  0  1  2  3  4  5  6  7  .  8  .  .  .  .  9  .  .  . 10  .  . 11  . 12
 *      |--R0-|  |-----Ra----| |-Rb-| |-------------- Rc --------------------|
 *               |------------------------- R -------------------------------|
 *
 * There are four regions in the rows:
 *    region 0 (R0)                : boundary row 0, and Backwards' two rows
 *    region a ("all"; Ra)         : all rows are kept (no checkpointing)
 *    region b ("between"; Rb)     : partially checkpointed
 *    region c ("checkpointed; Rc) : fully checkpointed
 *
 * In region a, La = Rb
 * In region b, Rb = 0|1, Lb = 0..Rc+1
 *              more specificially: (Rb=0 && Lb=0) || (Rb=1 && 1 <= Lb <= Rc+1)
 * In region c, Lc = {{Rc+2} \choose {2}}-1 = (Rc+2)(Rc+1)/2 - 1
 *
 * In this example:
 *    R0 = 3
 *    Ra = 5  La = 5
 *    Rb = 1  La = 2
 *    Rc = 4  Lc = 14
 *
 * In checkpointed regions, we refer to "blocks", often indexed
 * <b>.  There are Rb+Rc blocks, and each block ends in a checkpointed
 * row. The "width" of each block, often called <w>, decrements from
 * Rc+1 down to 2 in the fully checkpointed region.
 *
 * The reason to mix checkpointing and non-checkpointing is that we
 * use as many rows as we can, given a set memory ceiling, to minimize
 * computation time.
 *
 * The special states (ENJBC) are kept in xmx for all rows 1..L, just
 * as in a normal (uncheckpointed) P7_GMX.
 */

/*****************************************************************
 * 2. Exegesis: layout of rows in a checkpointed matrix.
 *****************************************************************/

/* Layout of memory in a single DP row:
 *
 *  dpc:   [M  I  D] [M  I  D] [M  I  D]  ...  [M  I  D]  [E  N  JJ  J  B  CC  C]
 *    k:   |-- 0 --| |-- 1 --| |-- 2 --|  ...  |-- M --|
 *         |------------- (M+1)*p7G_NSCELLS -----------|  |---- p7GC_NXCELLS ---|
 *
 *  Row dp[r] = gxc->dp_mem+(r*allocW) = dpc
 *  Main state s={MID} at node k={0..M}: dpc[k*p7G_NSCELLS+s]
 *  Special state s={ENJBC,CC,JJ}:       dpc[(M+1)*p7G_NSCELLS+s]
 *
 *  We need to store "JJ" and "CC" states -- the partial path
 *  probabilities from C(i-1)->C and J(i-1)->J -- because the
 *  checkpointed implementation does not necessarily have access to
 *  values on row i-1 when it does posterior decoding. <p7G_NXCELLS>
 *  from <P7_GMX> is replaced by <p7GC_NXCELLS> in <P7_GMXCHK>, and
 *  <enum p7g_xcells_e> from 0..4 {ENJBC} with <p7gc_xcells_e> from
 *  0..6 {E,N,JJ,J,B,CC,C}.
 */

/*****************************************************************
 * 3. The P7_GMXCHK data structure.
 *****************************************************************/

/* p7GC_NXCELLS and p7gc_xcells_e
 *
 * For main states we share p7G_NSCELLS and p7G_{MID} with P7_GMX.
 * For special states, we replace <enum p7g_xcells_e> with an array
 * that inserts CC and JJ cells, which the checkpointed implementation
 * needs. (See note 2 above.) Note that the order of the p7GC_{X}
 * special states is not the same as p7G_{X}, so they should not be
 * mixed.
 */
enum p7gc_xcells_e {
  p7GC_E  = 0,
  p7GC_N  = 1,
  p7GC_JJ = 2,
  p7GC_J  = 3,
  p7GC_B  = 4,
  p7GC_CC = 5,
  p7GC_C  = 6
};
#define p7GC_NXCELLS 7

typedef struct p7_gmxchk_s {
  int      M;	        /* actual query model dimension of current comparison                 */
  int      L;	        /* actual target sequence dimension of current comparison             */
  int      R;	        /* actual # rows in current fwd matrix (<= Ra+Rb+Rc), excluding R0    */

  /* Checkpointed layout, mapping rows 1..R to residues 1..L:                                 */
  int      R0;	        /* # of extra rows: one for fwd[0] boundary, two for bck[prv,cur]     */
  int      Ra;	        /* # of rows used in "all" region (uncheckpointed)                    */
  int      Rb;	        /* # of rows in "between" region (one incomplete checkpoint segment)  */
  int      Rc;	        /* # of rows in "checkpointed" region                                 */
  int      La;	        /* residues 1..La are in "all" region                                 */
  int      Lb;      	/* residues La+1..La+Lb are in "between" region                       */
  int      Lc;	        /* residues La+Lb+1..La+Lb+Lc=L are in "checkpointed" region          */

  float   *dp_mem;	/* raw memory allocation, that dp[] rows point into                         */
  int      allocW;	/* allocated width/row, in cells ((M+1)*p7G_NSCELLS+p7G_NXCELLS) <= allocW) */
  int64_t  ncells;	/* total # of alloc'ed cells: ncells >= (validR)(allocW)                    */
  int64_t  ncell_limit;	/* recommended RAM limit on dp_mem; can temporarily exceed it               */

  float  **dp;		/* DP row pointers, dp[0..R0-1,R0..R0+R-1]. See note above for layout.      */
  int      allocR;	/* allocated size of dp[]. R+R0 <= R0+Ra+Rb+Rc <= validR <= allocR          */
  int      validR;	/* # of rows pointing at DP memory; may be < allocR after a GrowTo() call   */

  /* Info for debugging mode (conditionally compiled)                              */
#ifdef p7_DEBUGGING
  int      do_debugging;	/* TRUE if we're in debugging mode                 */
  FILE    *dfp;			/* output stream for debugging diagnostics         */
  int      dbg_width;		/* cell values in diagnostic output are fprintf'ed */
  int      dbg_precision;       /*     dfp, "%*.*f", dbg_width, dbg_precision, val */
  int      dbg_flags;		/* p7_DEFAULT | p7_HIDE_SPECIALS | p7_SHOW_LOG     */
#endif
} P7_GMXCHK;

#define MMR(p, k) ((p)[(k)* p7G_NSCELLS + p7G_M])
#define IMR(p, k) ((p)[(k)* p7G_NSCELLS + p7G_I])
#define DMR(p, k) ((p)[(k)* p7G_NSCELLS + p7G_D])
#define XMR(p, s) ((p)[(M+1)* p7G_NSCELLS + s])

/*****************************************************************
 * 4. Declarations of the p7_gmxchk API
 *****************************************************************/

P7_GMXCHK *p7_gmxchk_Create (int M, int L, int64_t ramlimit);
int        p7_gmxchk_GrowTo (P7_GMXCHK *gxc, int M, int L);
size_t     p7_gmxchk_Sizeof (const P7_GMXCHK *gxc);
int        p7_gmxchk_Reuse  (P7_GMXCHK *gxc);
void       p7_gmxchk_Destroy(P7_GMXCHK *gxc);

int        p7_gmxchk_Dump(FILE *ofp, P7_GMXCHK *gxc, int flags);
int        p7_gmxchk_SetDumpMode(P7_GMXCHK *gxc, FILE *ofp, int flags);
int        p7_gmxchk_DumpHeader(FILE *ofp, P7_GMXCHK *gxc,  int kstart, int kend, int flags);
int        p7_gmxchk_DumpRow(FILE *ofp, P7_GMXCHK *gxc, float *dpc, int i, int kstart, int kend, int flags);

/*
 * References:
 *    SRE:J8/109-112, Oct 2011: Implementation plan
 */

/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Version 3.1b3; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 *
 * HMMER is distributed under the terms of the GNU General Public License
 * (GPLv3). See the LICENSE file for details.
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

#endif /*P7_GMXCHK_INCLUDED*/

/*** End of inlined file: p7_gmxchk.h ***/


/*** Start of inlined file: p7_hmmcache.h ***/
#ifndef P7_HMMCACHE_INCLUDED
#define P7_HMMCACHE_INCLUDED


/*** Start of inlined file: esl_alphabet.h ***/
#ifndef eslALPHABET_INCLUDED
#define eslALPHABET_INCLUDED

#include <ctype.h>		/* isascii() */

/*** Start of inlined file: easel.h ***/
#ifndef eslEASEL_INCLUDED
#define eslEASEL_INCLUDED


/*** Start of inlined file: esl_config.h ***/
/* esl_config.h.in  [input to configure]
 *
 * System-dependent configuration of Easel, by autoconf.
 *
 * This file should be included in all Easel .c files before
 * anything else, because it may set #define's that control
 * behaviour of system includes and system libraries. An example
 * is large file support.
 *
 */
#ifndef eslCONFIG_INCLUDED
#define eslCONFIG_INCLUDED

/* Version info.
 */
#define EASEL_VERSION "0.43"
#define EASEL_DATE "July 2016"
#define EASEL_COPYRIGHT "Copyright (C) 2016 Howard Hughes Medical Institute"
#define EASEL_LICENSE "Freely distributed under a BSD open source license."

/* Large file support
 * Must precede any header file inclusion.
 */
/* #undef _FILE_OFFSET_BITS */
/* #undef _LARGE_FILES */
/* #undef _LARGEFILE_SOURCE */

/* Debugging verbosity (0=none;3=most verbose)
 */
#define eslDEBUGLEVEL 0

/* System headers
 */
#define HAVE_ENDIAN_H 1
#define HAVE_INTTYPES_H 1
#define HAVE_STDINT_H 1
#define HAVE_UNISTD_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_STRINGS_H 1

#define HAVE_SYS_PARAM_H 1
#define HAVE_SYS_SYSCTL_H 1

#define HAVE_EMMINTRIN_H 1
#define HAVE_PMMINTRIN_H 1
#define HAVE_XMMINTRIN_H 1

/* #undef HAVE_ALTIVEC_H */

/* Types
 */
/* #undef WORDS_BIGENDIAN */
/* #undef int8_t */
/* #undef int16_t */
/* #undef int32_t */
/* #undef int64_t */
/* #undef uint8_t */
/* #undef uint16_t */
/* #undef uint32_t */
/* #undef uint64_t */
/* #undef off_t */

/* Optional packages
 */
/* #undef HAVE_LIBGSL */

/* Optional parallel implementation support
 */
#define HAVE_SSE2 1
/* #undef HAVE_VMX */
/* #undef HAVE_MPI */
#define HAVE_PTHREAD 1

#define HAVE_SSE2_CAST 1

/* Programs */
#define HAVE_GZIP 1

/* Functions */
/* #undef HAVE_CHMOD */
#define HAVE_FSEEKO 1
#define HAVE_FSTAT 1
#define HAVE_GETCWD 1
#define HAVE_GETPID 1
#define HAVE_MKSTEMP 1
#define HAVE_POPEN 1
/* #undef HAVE_PUTENV */
#define HAVE_STAT 1
#define HAVE_STRCASECMP 1
#define HAVE_SYSCONF 1
#define HAVE_SYSCTL 1
#define HAVE_TIMES 1
/* #undef HAVE_ERFC */

/* #undef HAVE_FUNC_ATTRIBUTE_NORETURN */

/* Function behavior */
#define eslSTOPWATCH_HIGHRES

/*****************************************************************
 * Available augmentations.
 *
 * If you grab a single module from Easel to use it by itself,
 * leave all these #undef'd; you have no augmentations.
 *
 * If you grab additional Easel .c files, you can enable any
 * augmentations they provide to other modules by #defining the
 * modules you have below. Alternatively, you can -D them on
 * the compile line, as in cc -DeslAUGMENT_SSI -DeslAUGMENT_MSA.
 *
 * If you compile and install the complete Easel library, all of these
 * get #defined automatically by ./configure, plus the eslLIBRARY flag
 * which means the full library with all augmentations is
 * available. So, if you steal files from an installed library, just
 * set these all back to #undef (depending on which files you have).
 *****************************************************************/
#define eslLIBRARY 1

#ifndef eslLIBRARY
/* #undef eslAUGMENT_ALPHABET */
/* #undef eslAUGMENT_NCBI */
/* #undef eslAUGMENT_DMATRIX */
/* #undef eslAUGMENT_FILEPARSER */
/* #undef eslAUGMENT_GEV */
/* #undef eslAUGMENT_GUMBEL */
/* #undef eslAUGMENT_HISTOGRAM */
/* #undef eslAUGMENT_KEYHASH */
/* #undef eslAUGMENT_MINIMIZER */
/* #undef eslAUGMENT_MSA */
/* #undef eslAUGMENT_RANDOM */
/* #undef eslAUGMENT_RANDOMSEQ */
/* #undef eslAUGMENT_SSI */
/* #undef eslAUGMENT_STATS */
#endif

#ifdef eslLIBRARY
#define eslAUGMENT_ALPHABET
#define eslAUGMENT_NCBI
#define eslAUGMENT_DMATRIX
#define eslAUGMENT_FILEPARSER
#define eslAUGMENT_GEV
#define eslAUGMENT_GUMBEL
#define eslAUGMENT_HISTOGRAM
#define eslAUGMENT_KEYHASH
#define eslAUGMENT_MINIMIZER
#define eslAUGMENT_MSA
#define eslAUGMENT_RANDOM
#define eslAUGMENT_RANDOMSEQ
#define eslAUGMENT_SSI
#define eslAUGMENT_STATS
#endif

#endif /*eslCONFIG_INCLUDED*/

/*** End of inlined file: esl_config.h ***/

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
 *       extern void fatal(char *msg, ...) ESL_ANALYZER_NORETURN;
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

#endif /*eslEASEL_INCLUDED*/

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


/*** Start of inlined file: hmmer.h ***/
#ifndef P7_HMMERH_INCLUDED
#define P7_HMMERH_INCLUDED

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#ifdef HMMER_THREADS
#include <pthread.h>
#endif


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


/*** Start of inlined file: esl_histogram.h ***/
#ifndef eslHISTOGRAM_INCLUDED
#define eslHISTOGRAM_INCLUDED

#include <math.h>   /* floor() is in one of the macros */

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


/*** Start of inlined file: esl_msa.h ***/
#ifndef eslMSA_INCLUDED
#define eslMSA_INCLUDED

#include <stdio.h>


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
  int         max_ram;	        /* threshold in MB to trigger extern sort */

  char      **filenames;
  uint32_t   *fileformat;
  uint32_t   *bpl;
  uint32_t   *rpl;
  uint32_t    flen;		/* length of longest filename, inc '\0' */
  uint16_t    nfiles;		/* can store up to 2^15-1 (32767) files */

  ESL_PKEY   *pkeys;
  uint32_t    plen;	        /* length of longest pkey, including '\0'    */
  uint64_t    nprimary;		/* can store up to 2^63-1 = 9.2e18 keys      */
  char       *ptmpfile;		/* primary key tmpfile name, for extern sort */
  FILE       *ptmp;	        /* handle on open ptmpfile */

  ESL_SKEY   *skeys;
  uint32_t    slen;        	/* length of longest skey, including '\0' */
  uint64_t    nsecondary;
  char       *stmpfile;		/* secondary key tmpfile name, for extern sort */
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


/*** Start of inlined file: esl_sq.h ***/
#ifndef eslSQ_INCLUDED
#define eslSQ_INCLUDED

#ifdef eslAUGMENT_ALPHABET

#endif
#ifdef eslAUGMENT_MSA

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


/*** Start of inlined file: esl_stopwatch.h ***/
//#ifndef eslSTOPWATCH_INCLUDED
//#define eslSTOPWATCH_INCLUDED

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

//#endif /*eslSTOPWATCH_INCLUDED*/
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

/* Search modes. */
#define p7_NO_MODE   0
#define p7_LOCAL     1		/* multihit local:  "fs" mode   */
#define p7_GLOCAL    2		/* multihit glocal: "ls" mode   */
#define p7_UNILOCAL  3		/* unihit local: "sw" mode      */
#define p7_UNIGLOCAL 4		/* unihit glocal: "s" mode      */

#define p7_IsLocal(mode)  (mode == p7_LOCAL || mode == p7_UNILOCAL)
#define p7_IsMulti(mode)  (mode == p7_LOCAL || mode == p7_GLOCAL)

#define p7_NEVPARAM 6	/* number of statistical parameters stored in models                      */
#define p7_NCUTOFFS 6	/* number of Pfam score cutoffs stored in models                          */
#define p7_NOFFSETS 3	/* number of disk offsets stored in models for hmmscan's fast model input */
enum p7_evparams_e {    p7_MMU  = 0, p7_MLAMBDA = 1,     p7_VMU = 2,  p7_VLAMBDA = 3, p7_FTAU = 4, p7_FLAMBDA = 5 };
enum p7_cutoffs_e  {     p7_GA1 = 0,     p7_GA2 = 1,     p7_TC1 = 2,      p7_TC2 = 3,  p7_NC1 = 4,     p7_NC2 = 5 };
enum p7_offsets_e  { p7_MOFFSET = 0, p7_FOFFSET = 1, p7_POFFSET = 2 };

#define p7_EVPARAM_UNSET -99999.0f  /* if evparam[0] is unset, then all unset                         */
#define p7_CUTOFF_UNSET  -99999.0f  /* if cutoff[XX1] is unset, then cutoff[XX2] unset, XX={GA,TC,NC} */
#define p7_COMPO_UNSET   -1.0f      /* if compo[0] is unset, then all unset                           */

/* Option flags when creating multiple alignments with p7_tracealign_*() */
#define p7_DEFAULT             0
#define p7_DIGITIZE            (1<<0)
#define p7_ALL_CONSENSUS_COLS  (1<<1)
#define p7_TRIM                (1<<2)

/* Option flags when creating faux traces with p7_trace_FauxFromMSA() */
#define p7_MSA_COORDS	       (1<<0) /* default: i = unaligned seq residue coords     */

/* Which strand(s) should be searched */
enum p7_strands_e {    p7_STRAND_TOPONLY  = 0, p7_STRAND_BOTTOMONLY = 1,  p7_STRAND_BOTH = 2};

/*****************************************************************
 * 1. P7_HMM: a core model.
 *****************************************************************/

/* Bit flags used in <hmm->flags>: optional annotation in an HMM
 *
 * Flags marked with ! may not be changed nor used for other meanings,
 * because they're codes used by HMMER2 (and earlier) that must be
 * preserved for reverse compatibility with old HMMER files.
 *
 * Why use flags? (So I don't ask this question of myself again:)
 *   1. The way we allocate an HMM, we need to know if we're allocating
 *      M-width annotation fields (RF, CS, CA, MAP) before we read the
 *      annotation from the file.
 *   2. Historically, H2 used flags, so we still need to read H2 flags
 *      for backwards compatibility; so we may as well keep using them.
 */
#define p7H_HASBITS (1<<0)    /* obsolete (was: model has log-odds scores)       !*/
#define p7H_DESC    (1<<1)    /* description exists (legacy; xref SRE:J5/114)    !*/
#define p7H_RF      (1<<2)    /* #RF annotation available                        !*/
#define p7H_CS      (1<<3)    /* #CS annotation available                        !*/
#define p7H_XRAY    (1<<4)    /* obsolete (was: structural data available)       !*/
#define p7H_HASPROB (1<<5)    /* obsolete (was: model in probability form)       !*/
#define p7H_HASDNA  (1<<6)    /* obsolete (was: protein HMM->DNA seq params set) !*/
#define p7H_STATS   (1<<7)    /* model has E-value statistics calibrated         !*/
#define p7H_MAP     (1<<8)    /* alignment map is available                      !*/
#define p7H_ACC     (1<<9)    /* accession is available (legacy; xref SRE:J5/114)!*/
#define p7H_GA      (1<<10)   /* gathering thresholds available                  !*/
#define p7H_TC      (1<<11)   /* trusted cutoffs available                       !*/
#define p7H_NC      (1<<12)   /* noise cutoffs available                         !*/
#define p7H_CA      (1<<13)   /* surface accessibilities available               !*/
#define p7H_COMPO   (1<<14)   /* model-specific residue composition available     */
#define p7H_CHKSUM  (1<<15)   /* model has an alignment checksum                  */
#define p7H_CONS    (1<<16)   /* consensus residue line available                 */
#define p7H_MMASK   (1<<17)   /* #MM annotation available                        !*/

/* Indices of Plan7 main model state transitions, hmm->t[k][] */
enum p7h_transitions_e {
  p7H_MM = 0,
  p7H_MI = 1,
  p7H_MD = 2,
  p7H_IM = 3,
  p7H_II = 4,
  p7H_DM = 5,
  p7H_DD = 6
};
#define p7H_NTRANSITIONS 7

/* How the hmm->t[k] vector is interpreted as separate probability vectors. */
#define P7H_TMAT(hmm, k) ((hmm)->t[k])
#define P7H_TINS(hmm, k) ((hmm)->t[k]+3)
#define P7H_TDEL(hmm, k) ((hmm)->t[k]+5)
#define p7H_NTMAT 3
#define p7H_NTDEL 2
#define p7H_NTINS 2

/* Some notes:
 *   0. The model might be either in counts or probability form.
 *   1. t[0] is special: t[0][TMM,TMI,TMD] are the begin->M_1,I_0,D_1 entry probabilities,
 *      t[0][TIM,TII] are the I_0 transitions, and delete state 0 doesn't
 *      exist. Therefore D[0] transitions and mat[0] emissions are unused.
 *      To simplify some normalization code, we adopt a convention that these are set
 *      to valid probability distributions: 1.0 for t[0][TDM] and mat[0][0],
 *      and 0 for the rest.
 *   2. t[M] is also special: TMD and TDD are 0 because there is no next delete state;
 *      TDM is therefore 1.0 by definition. TMM and TDM are interpreted as the
 *      M->E and D->E end transitions. t[M][TDM] must be 1.0, therefore.
 */
typedef struct p7_hmm_s {
  /*::cexcerpt::plan7_core::begin::*/
  int     M;                    /* length of the model (# nodes)                           */
  float **t;                    /* transition prob's. t[(0),1..M][0..p7H_NTRANSITIONS-1]   */
  float **mat;                  /* match emissions.  mat[1..M][0..K-1]                     */
  float **ins;                  /* insert emissions. ins[1..M][0..K-1]                     */
  /*::cexcerpt::plan7_core::end::*/

  /* Annotation. Everything but <name> is optional. Flags are set when
   * optional values are set. All the char *'s are proper nul-terminated
   * strings, not just arrays. (hmm->map is an int array).
   */
  char    *name;                 /* name of the model                     (mandatory)      */ /* String, \0-terminated   */
  char    *acc;	                 /* accession number of model (Pfam)      (optional: NULL) */ /* String, \0-terminated   */
  char    *desc;                 /* brief (1-line) description of model   (optional: NULL) */ /* String, \0-terminated   */
  char    *rf;                   /* reference line from alignment 1..M    (p7H_RF)         */ /* String; 0=' ', M+1='\0' */
  char    *mm;                   /* model mask line from alignment 1..M   (p7H_MM)         */ /* String; 0=' ', M+1='\0' */
  char    *consensus;		         /* consensus residue line        1..M    (p7H_CONS)       */ /* String; 0=' ', M+1='\0' */
  char    *cs;                   /* consensus structure line      1..M    (p7H_CS)         */ /* String; 0=' ', M+1='\0' */
  char    *ca;	                 /* consensus accessibility line  1..M    (p7H_CA)         */ /* String; 0=' ', M+1='\0' */

  char    *comlog;               /* command line(s) that built model      (optional: NULL) */ /* String, \0-terminated   */
  int      nseq;	         /* number of training sequences          (optional: -1)   */
  float    eff_nseq;             /* effective number of seqs (<= nseq)    (optional: -1)   */
  int	   max_length;           /* upper bound length, all but 1e-7 prob (optional: -1)   */
  char    *ctime;	         /* creation date                         (optional: NULL) */
  int     *map;	                 /* map of alignment cols onto model 1..M (p7H_MAP)        */ /* Array; map[0]=0 */
  uint32_t checksum;             /* checksum of training sequences        (p7H_CHKSUM)     */
  float    evparam[p7_NEVPARAM]; /* E-value params                        (p7H_STATS)      */
  float    cutoff[p7_NCUTOFFS];  /* Pfam score cutoffs                    (p7H_{GA,TC,NC}) */
  float    compo[p7_MAXABET];    /* model bg residue comp                 (p7H_COMPO)      */

  off_t    offset;               /* HMM record offset on disk                              */
  const ESL_ALPHABET *abc;       /* ptr to alphabet info (hmm->abc->K is alphabet size)    */
  int      flags;                /* status flags                                           */
} P7_HMM;

/*****************************************************************
 * 2. P7_PROFILE: a scoring profile, and its implicit model.
 *****************************************************************/

/* Indices for special state types in the length model, gm->xsc[x][]
 */
enum p7p_xstates_e {
  p7P_E = 0,
  p7P_N = 1,
  p7P_J = 2,
  p7P_C = 3
};
#define p7P_NXSTATES 4

/* Indices for transitions from the length modeling scores gm->xsc[][x]
 */
enum p7p_xtransitions_e {
  p7P_LOOP = 0,
  p7P_MOVE = 1
};
#define p7P_NXTRANS 2

/* Indices for transition scores gm->tsc[k][] */
/* order is optimized for dynamic programming */
enum p7p_tsc_e {
  p7P_MM = 0,
  p7P_IM = 1,
  p7P_DM = 2,
  p7P_BM = 3,
  p7P_MD = 4,
  p7P_DD = 5,
  p7P_MI = 6,
  p7P_II = 7,
};
#define p7P_NTRANS 8

/* Indices for residue emission score vectors
 */
enum p7p_rsc_e {
  p7P_MSC = 0,
  p7P_ISC = 1
};
#define p7P_NR 2

/* Accessing transition, emission scores */
/* _BM is specially stored off-by-one: [k-1][p7P_BM] is score for entering at Mk */
#define p7P_TSC(gm, k, s) ((gm)->tsc[(k) * p7P_NTRANS + (s)])
#define p7P_MSC(gm, k, x) ((gm)->rsc[x][(k) * p7P_NR + p7P_MSC])
#define p7P_ISC(gm, k, x) ((gm)->rsc[x][(k) * p7P_NR + p7P_ISC])

typedef struct p7_profile_s {
  float  *tsc;          /* transitions  [0.1..M-1][0..p7P_NTRANS-1], hand-indexed  */
  float **rsc;          /* emissions [0..Kp-1][0.1..M][p7P_NR], hand-indexed       */
  float   xsc[p7P_NXSTATES][p7P_NXTRANS]; /* special transitions [NECJ][LOOP,MOVE] */

  int     mode;        	/* configured algorithm mode (e.g. p7_LOCAL)               */
  int     L;		/* current configured target seq length                    */
  int     allocM;	/* max # of nodes allocated in this structure              */
  int     M;		/* number of nodes in the model                            */
  int     max_length;	/* calculated upper bound on emitted seq length            */
  float   nj;		/* expected # of uses of J; precalculated from loop config */

  /* Info, most of which is a copy from parent HMM:                                       */
  char  *name;			/* unique name of model                                   */
  char  *acc;			/* unique accession of model, or NULL                     */
  char  *desc;                  /* brief (1-line) description of model, or NULL           */
  char  *rf;                    /* reference line from alignment 1..M; *rf=0 means unused */
  char  *mm;                    /* modelmask line           1..M; *ref=0: unused     */
  char  *cs;                    /* consensus structure line      1..M, *cs=0 means unused */
  char  *consensus;		/* consensus residues to display in alignments, 1..M      */
  float  evparam[p7_NEVPARAM]; 	/* parameters for determining E-values, or UNSET          */
  float  cutoff[p7_NCUTOFFS]; 	/* per-seq/per-domain bit score cutoffs, or UNSET         */
  float  compo[p7_MAXABET];	/* per-model HMM filter composition, or UNSET             */

  /* Disk offset information for hmmpfam's fast model retrieval                           */
  off_t  offs[p7_NOFFSETS];     /* p7_{MFP}OFFSET, or -1                                  */

  off_t  roff;                  /* record offset (start of record); -1 if none            */
  off_t  eoff;                  /* offset to last byte of record; -1 if unknown           */

  const ESL_ALPHABET *abc;	/* copy of pointer to appropriate alphabet                */
} P7_PROFILE;

/*****************************************************************
 * 3. P7_BG: a null (background) model.
 *****************************************************************/

/* This really contains three different things:
 *
 *   - the "null1" model, a one-state HMM consisting of background
 *     frequencies <f> and a parameter <p1> for a target-length
 *     dependent geometric;
 *
 *   - the "bias filter" <fhmm> a two-state HMM composed from null1's
 *     background <f> and the model's mean composition <compo>. This
 *     model is constructed dynamically, every time a new profile is
 *     considered;
 *
 *   - a single term <omega> that's needed by the "null2" model to set
 *     a balance between the null1 and null2 scoring terms.  The null2
 *     model is otherwise defined by construction, in p7_domaindef.c.
 *
 * Someday we might pull this apart into two or three separate
 * objects.
 */
typedef struct p7_bg_s {
  float   *f;		/* null1 background residue frequencies [0..K-1]: set at initialization    */
  float    p1;		/* null1's transition prob: p7_bg_SetLength() sets this from target seq L  */

  ESL_HMM *fhmm;	/* bias filter: p7_bg_SetFilter() sets this, from model's mean composition */

  float    omega;	/* the "prior" on null2/null3: set at initialization (one omega for both null types)  */

  const ESL_ALPHABET *abc;	/* reference to alphabet in use: set at initialization             */
} P7_BG;

/*****************************************************************
 * 4. P7_TRACE:  a traceback (alignment of seq to profile).
 *****************************************************************/

/* Traceback structure for alignment of a model to a sequence.
 *
 * A traceback only makes sense in a triplet (tr, gm, dsq), for a
 * given profile or HMM (with nodes 1..M) and a given digital sequence
 * (with positions 1..L).
 *
 * A traceback may be relative to a profile (usually) or to a core
 * model (as a special case in model construction; see build.c). You
 * can tell the difference by looking at the first statetype,
 * tr->st[0]; if it's a p7T_S, it's for a profile, and if it's p7T_B,
 * it's for a core model.
 *
 * A "profile" trace uniquely has S,N,C,T,J states and their
 * transitions; it also can have B->Mk and Mk->E internal entry/exit
 * transitions for local alignments. It may not contain X states.
 *
 * A "core" trace may contain I0, IM, and D1 states and their
 * transitions. A "core" trace can also have B->X->{MDI}k and
 * {MDI}k->X->E transitions as a special hack in a build procedure, to
 * deal with the case of a local alignment fragment implied by an
 * input alignment, which is "impossible" for a core model.
 * X "states" only appear in core traces, and only at these
 * entry/exit places; some code depends on this.
 *
 * A profile's N,C,J states emit on transition, not on state, so a
 * path of N emits 0 residues, NN emits 1 residue, NNN emits 2
 * residues, and so on. By convention, the trace always associates an
 * emission-on-transition with the trailing (destination) state, so
 * the first N, C, or J is stored in a trace as a nonemitter (i=0).
 *
 * A i coords in a traceback are usually 1..L with respect to an
 * unaligned digital target sequence, but in the special case of
 * traces faked from existing MSAs (as in hmmbuild), the coords may
 * be 1..alen relative to an MSA's columns.
 */

/* State types */
enum p7t_statetype_e {
  p7T_BOGUS =  0,
  p7T_M     =  1,
  p7T_D     =  2,
  p7T_I     =  3,
  p7T_S     =  4,
  p7T_N     =  5,
  p7T_B     =  6,
  p7T_E     =  7,
  p7T_C     =  8,
  p7T_T     =  9,
  p7T_J     = 10,
  p7T_X     = 11, 	/* missing data: used esp. for local entry/exits */
};
#define p7T_NSTATETYPES 12

typedef struct p7_trace_s {
  int    N;		/* length of traceback                       */
  int    nalloc;        /* allocated length of traceback             */
  char  *st;		/* state type code                   [0..N-1]*/
  int   *k;		/* node index; 1..M if M,D,I; else 0 [0..N-1]*/
  int   *i;		/* pos emitted in dsq, 1..L; else 0  [0..N-1]*/
  float *pp;		/* posterior prob of x_i; else 0     [0..N-1]*/
  int    M;		/* model length M (maximum k)                */
  int    L;		/* sequence length L (maximum i)             */

  /* The following section is data generated by "indexing" a trace's domains */
  int   ndom;		/* number of domains in trace (= # of B or E states) */
  int  *tfrom,   *tto;	/* locations of B/E states in trace (0..tr->N-1)     */
  int  *sqfrom,  *sqto;	/* first/last M-emitted residue on sequence (1..L)   */
  int  *hmmfrom, *hmmto;/* first/last M state on model (1..M)                */
  int   ndomalloc;	/* current allocated size of these stacks            */

} P7_TRACE;

/*****************************************************************
 * 5. P7_HMMFILE:  an HMM save file or database, open for reading.
 *****************************************************************/

/* These tags need to be in temporal order, so we can do tests
 * like "if (format >= p7_HMMFILE_3b) ..."
 */
enum p7_hmmfile_formats_e {
  p7_HMMFILE_20 = 0,
  p7_HMMFILE_3a = 1,
  p7_HMMFILE_3b = 2,
  p7_HMMFILE_3c = 3,
  p7_HMMFILE_3d = 4,
  p7_HMMFILE_3e = 5,
  p7_HMMFILE_3f = 6,
};

typedef struct p7_hmmfile_s {
  FILE         *f;		 /* pointer to stream for reading models                 */
  char         *fname;	         /* (fully qualified) name of the HMM file; [STDIN] if - */
  ESL_SSI      *ssi;		 /* open SSI index for model file <f>; NULL if none.     */

  int           do_gzip;	/* TRUE if f is "gzip -dc |" (will pclose(f))           */
  int           do_stdin;       /* TRUE if f is stdin (won't close f)                   */
  int           newly_opened;	/* TRUE if we just opened the stream (and parsed magic) */
  int           is_pressed;	/* TRUE if a pressed HMM database file (Pfam or equiv)  */

  int            format;	/* HMM file format code */
  int           (*parser)(struct p7_hmmfile_s *, ESL_ALPHABET **, P7_HMM **);
  ESL_FILEPARSER *efp;

  /* If <is_pressed>, we can read optimized profiles directly, via:  */
  FILE         *ffp;		/* MSV part of the optimized profile */
  FILE         *pfp;		/* rest of the optimized profile     */

#ifdef HMMER_THREADS
  int              syncRead;
  pthread_mutex_t  readMutex;
#endif

  char          errbuf[eslERRBUFSIZE];
} P7_HMMFILE;

/* note on <fname>, above:
 * this is the actual name of the HMM file being read.
 *
 * The way p7_hmmfile_Open() works, it will preferentially look for
 * hmmpress'ed binary files. If you open "foo", it will first try to
 * open "foo.h3m" and <fname> will be "foo.h3m". "foo" does not even
 * have to exist. If a parsing error occurs, you want <fname> to
 * be "foo.h3m", so error messages report blame correctly.
 * In the special case of reading from stdin, <fname> is "[STDIN]".
 */

/*****************************************************************
 * 6. P7_GMX: a "generic" dynamic programming matrix
 *****************************************************************/

enum p7g_scells_e {
  p7G_M = 0,
  p7G_I = 1,
  p7G_D = 2,
};
#define p7G_NSCELLS 3

enum p7g_xcells_e {
  p7G_E  = 0,
  p7G_N  = 1,
  p7G_J  = 2,
  p7G_B  = 3,
  p7G_C  = 4
};
#define p7G_NXCELLS 5

typedef struct p7_gmx_s {
  int  M;		/* actual model dimension (model 1..M)    */
  int  L;		/* actual sequence dimension (seq 1..L)   */

  int      allocR;      /* current allocated # of rows : L+1 <= validR <= allocR                */
  int      validR;	/* # of rows actually pointing at DP memory                             */
  int      allocW;	/* current set row width :  M+1 <= allocW                               */
  uint64_t ncells;	/* total # of allocated cells in 2D matrix : ncells >= (validR)(allocW) */

  float **dp;           /* logically [0.1..L][0.1..M][0..p7G_NSCELLS-1]; indexed [i][k*p7G_NSCELLS+s] */
  float  *xmx;          /* logically [0.1..L][0..p7G_NXCELLS-1]; indexed [i*p7G_NXCELLS+s]            */

  float  *dp_mem;
} P7_GMX;

/* Macros below implement indexing idioms for generic DP routines.
 * They require the following setup, for profile <gm> and matrix <gx>:
 *   float const *tsc = gm->tsc;
 *   float      **dp  = gx->dp;
 *   float       *xmx = gx->xmx;
 * and for each row i (target residue x_i in digital seq <dsq>):
 *   float const *rsc = gm->rsc[dsq[i]];
 */
#define MMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_M])
#define IMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_I])
#define DMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_D])
#define XMX(i,s) (xmx[(i) * p7G_NXCELLS + (s)])

#define TSC(s,k) (tsc[(k) * p7P_NTRANS + (s)])
#define MSC(k)   (rsc[(k) * p7P_NR     + p7P_MSC])
#define ISC(k)   (rsc[(k) * p7P_NR     + p7P_ISC])

/* Flags that control P7_GMX debugging dumps */
#define p7_HIDE_SPECIALS (1<<0)
#define p7_SHOW_LOG      (1<<1)

/*****************************************************************
 * 7. P7_PRIOR: mixture Dirichlet prior for profile HMMs
 *****************************************************************/

typedef struct p7_prior_s {
  ESL_MIXDCHLET *tm;		/*  match transitions */
  ESL_MIXDCHLET *ti;		/* insert transitions */
  ESL_MIXDCHLET *td;		/* delete transitions */
  ESL_MIXDCHLET *em;		/*  match emissions   */
  ESL_MIXDCHLET *ei;		/* insert emissions   */
} P7_PRIOR;

/*****************************************************************
 * 8. P7_SPENSEMBLE: segment pair ensembles for domain locations
 *****************************************************************/

/* struct p7_spcoord_s:
 *    a coord quad defining a segment pair.
 */
struct p7_spcoord_s {
  int idx; 	/* backreference index: which trace a seg came from, or which cluster a domain came from */
  int i, j;	/* start,end in a target sequence (1..L)  */
  int k, m;     /* start,end in a query model (1..M)      */
  float prob;	/* posterior probability of segment       */
};

/* Structure: P7_SPENSEMBLE
 *
 * Collection and clustering of an ensemble of sampled segment pairs,
 * in order to define domain locations using their posterior
 * probability distribution (as opposed to Viterbi MAP tracebacks).
 */
typedef struct p7_spensemble_s {
  /* Section 1: a collected ensemble of segment pairs                                       */
  int                  nsamples;    /* number of sampled traces                             */
  struct p7_spcoord_s *sp;	    /* array of sampled seg pairs; [0..n-1]                 */
  int                  nalloc;	    /* allocated size of <sp>                               */
  int                  n;	    /* number of seg pairs in <sp>                          */

  /* Section 2: then the ensemble is clustered by single-linkage clustering                 */
  int *workspace;                   /* temp space for Easel SLC algorithm: 2*n              */
  int *assignment;                  /* each seg pair's cluster index: [0..n-1] = (0..nc-1)  */
  int  nc;	                    /* number of different clusters                         */

  /* Section 3: then endpoint distribution is examined within each large cluster            */
  int *epc;	                    /* array counting frequency of each endpoint            */
  int  epc_alloc;	            /* allocated width of <epc>                             */

  /* Section 4: finally each large cluster is resolved into domain coords                   */
  struct p7_spcoord_s *sigc;	    /* array of coords for each domain, [0..nsigc-1]        */
  int                  nsigc;	    /* number of "significant" clusters, domains            */
  int                  nsigc_alloc; /* current allocated max for nsigc                      */
} P7_SPENSEMBLE;

/*****************************************************************
 * 9. P7_ALIDISPLAY: an alignment formatted for printing
 *****************************************************************/

/* Structure: P7_ALIDISPLAY
 *
 * Alignment of a sequence domain to an HMM, formatted for printing.
 *
 * For an alignment of L residues and names C chars long, requires
 * 6L + 2C + 30 bytes; for typical case of L=100,C=10, that's
 * <0.7 Kb.
 */
typedef struct p7_alidisplay_s {
  char *rfline;                 /* reference coord info; or NULL        */
  char *mmline;                 /* modelmask coord info; or NULL        */
  char *csline;                 /* consensus structure info; or NULL    */
  char *model;                  /* aligned query consensus sequence     */
  char *mline;                  /* "identities", conservation +'s, etc. */
  char *aseq;                   /* aligned target sequence              */
  char *ntseq;                  /* nucleotide target sequence if nhmmscant */
  char *ppline;			        /* posterior prob annotation; or NULL   */
  int   N;			            /* length of strings                    */

  char *hmmname;		/* name of HMM                          */
  char *hmmacc;			/* accession of HMM; or [0]='\0'        */
  char *hmmdesc;		/* description of HMM; or [0]='\0'      */
  int   hmmfrom;		/* start position on HMM (1..M, or -1)  */
  int   hmmto;			/* end position on HMM (1..M, or -1)    */
  int   M;			/* length of model                      */

  char *sqname;			/* name of target sequence              */
  char *sqacc;			/* accession of target seq; or [0]='\0' */
  char *sqdesc;			/* description of targ seq; or [0]='\0' */
  long  sqfrom;			/* start position on sequence (1..L)    */
  long  sqto;		    /* end position on sequence   (1..L)    */
  long  L;			/* length of sequence                   */

  int   memsize;                /* size of allocated block of memory    */
  char *mem;			/* memory used for the char data above  */
} P7_ALIDISPLAY;

/*****************************************************************
 * 10. P7_DOMAINDEF: reusably managing workflow in defining domains
 *****************************************************************/

typedef struct p7_dom_s {
  int            ienv, jenv;
  int            iali, jali;
  int            iorf, jorf; /*Used in translated search to capture the range in the DNA sequence of the ORF containing the match to a protein query */
  float          envsc;  	/* Forward score in envelope ienv..jenv; NATS; without null2 correction       */
  float          domcorrection;	/* null2 score when calculating a per-domain score; NATS                      */
  float          dombias;	/* FLogsum(0, log(bg->omega) + domcorrection): null2 score contribution; NATS */
  float          oasc;		/* optimal accuracy score (units: expected # residues correctly aligned)      */
  float          bitscore;	/* overall score in BITS, null corrected, if this were the only domain in seq */
  double         lnP;	        /* log(P-value) of the bitscore                                               */
  int            is_reported;	/* TRUE if domain meets reporting thresholds                                  */
  int            is_included;	/* TRUE if domain meets inclusion thresholds                                  */
  float         *scores_per_pos; /* score in BITS that each position in the alignment contributes to an overall viterbi score */
  P7_ALIDISPLAY *ad;
} P7_DOMAIN;

/* Structure: P7_DOMAINDEF
 *
 * This is a container for all the necessary information for domain
 * definition procedures in <p7_domaindef.c>, including a bunch of
 * heuristic thresholds. The structure is reusable to minimize the
 * number of allocation/free cycles that need to be done when
 * processing a large number of sequences. You create the structure
 * with <p7_domaindef_Create()>; after you're done with defining
 * domains on a sequence, you call <p7_domaindef_Reuse()> before using
 * it on the next sequence; and when you're completely done, you free
 * it with <p7_domaindef_Destroy()>. All memory management is handled
 * internally; you don't need to reallocate anything yourself.
 */
typedef struct p7_domaindef_s {
  /* for posteriors of being in a domain, B, E */
  float *mocc;			/* mocc[i=1..L] = prob that i is emitted by core model (is in a domain)       */
  float *btot; 			/* btot[i=1..L] = cumulative expected times that domain starts at or before i */
  float *etot;			/* etot[i=1..L] = cumulative expected times that domain ends at or before i   */
  int    L;
  int    Lalloc;

  /* the ad hoc null2 model: 1..L nat scores for each residue, log f'(x_i) / f(x_i) */
  float *n2sc;

  /* rng and reusable memory for stochastic tracebacks */
  ESL_RANDOMNESS *r;		/* random number generator                                 */
  int             do_reseeding;	/* TRUE to reset the RNG, make results reproducible        */
  P7_SPENSEMBLE  *sp;		/* an ensemble of sampled segment pairs (domain endpoints) */
  P7_TRACE       *tr;		/* reusable space for a trace of a domain                  */
  P7_TRACE       *gtr;		/* reusable space for a traceback of the entire target seq */

  /* Heuristic thresholds that control the region definition process */
  /* "rt" = "region threshold", for lack of better term  */
  float  rt1;   	/* controls when regions are called. mocc[i] post prob >= dt1 : triggers a region around i */
  float  rt2;		/* controls extent of regions. regions extended until mocc[i]-{b,e}occ[i] < dt2            */
  float  rt3;		/* controls when regions are flagged for split: if expected # of E preceding B is >= dt3   */

  /* Heuristic thresholds that control the stochastic traceback/clustering process */
  int    nsamples;	/* collect ensemble of this many stochastic traces */
  float  min_overlap;	/* 0.8 means >= 80% overlap of (smaller/larger) segment to link, both in seq and hmm            */
  int    of_smaller;	/* see above; TRUE means overlap denom is calc'ed wrt smaller segment; FALSE means larger       */
  int    max_diagdiff;	/* 4 means either start or endpoints of two segments must be within <=4 diagonals of each other */
  float  min_posterior;	/* 0.25 means a cluster must have >= 25% posterior prob in the sample to be reported            */
  float  min_endpointp;	/* 0.02 means choose widest endpoint with post prob of at least 2%                              */

  /* storage of the results; domain locations, scores, alignments          */
  P7_DOMAIN *dcl;
  int        ndom;	 /* number of domains defined, in the end.         */
  int        nalloc;     /* number of domain structures allocated in <dcl> */

  /* Additional results storage */
  float  nexpected;     /* posterior expected number of domains in the sequence (from posterior arrays) */
  int    nregions;	/* number of regions evaluated */
  int    nclustered;	/* number of regions evaluated by clustering ensemble of tracebacks */
  int    noverlaps;	/* number of envelopes defined in ensemble clustering that overlap w/ prev envelope */
  int    nenvelopes;	/* number of envelopes handed over for domain definition, null2, alignment, and scoring. */

} P7_DOMAINDEF;

/*****************************************************************
 * 11. P7_TOPHITS: ranking lists of top-scoring hits
 *****************************************************************/

#define p7_HITFLAGS_DEFAULT 0
#define p7_IS_INCLUDED      (1<<0)
#define p7_IS_REPORTED      (1<<1)
#define p7_IS_NEW           (1<<2)
#define p7_IS_DROPPED       (1<<3)
#define p7_IS_DUPLICATE     (1<<4)

/* Structure: P7_HIT
 *
 * Info about a high-scoring database hit, kept so we can output a
 * sorted list of high hits at the end.
 *
 * sqfrom and sqto are the coordinates that will be shown
 * in the results, not coords in arrays... therefore, reverse
 * complements have sqfrom > sqto
 */
typedef struct p7_hit_s {
  char   *name;			/* name of the target               (mandatory)           */
  char   *acc;			/* accession of the target          (optional; else NULL) */
  char   *desc;			/* description of the target        (optional; else NULL) */
  int    window_length;         /* for later use in e-value computation, when splitting long sequences */
  double sortkey;		/* number to sort by; big is better                       */

  float  score;			/* bit score of the sequence (all domains, w/ correction) */
  float  pre_score;		/* bit score of sequence before null2 correction          */
  float  sum_score;		/* bit score reconstructed from sum of domain envelopes   */

  double lnP;		        /* log(P-value) of the score               */
  double pre_lnP;		/* log(P-value) of the pre_score           */
  double sum_lnP;		/* log(P-value) of the sum_score           */

  float  nexpected;     /* posterior expected number of domains in the sequence (from posterior arrays) */
  int    nregions;	/* number of regions evaluated */
  int    nclustered;	/* number of regions evaluated by clustering ensemble of tracebacks */
  int    noverlaps;	/* number of envelopes defined in ensemble clustering that overlap w/ prev envelope */
  int    nenvelopes;	/* number of envelopes handed over for domain definition, null2, alignment, and scoring. */
  int    ndom;		/* total # of domains identified in this seq   */

  uint32_t flags;      	/* p7_IS_REPORTED | p7_IS_INCLUDED | p7_IS_NEW | p7_IS_DROPPED */
  int      nreported;	/* # of domains satisfying reporting thresholding  */
  int      nincluded;	/* # of domains satisfying inclusion thresholding */
  int      best_domain;	/* index of best-scoring domain in dcl */

  int64_t  seqidx;          /*unique identifier to track the database sequence from which this hit came*/
  int64_t  subseq_start; /*used to track which subsequence of a full_length target this hit came from, for purposes of removing duplicates */

  P7_DOMAIN *dcl;	/* domain coordinate list and alignment display */
  esl_pos_t  offset;	/* used in socket communications, in serialized communication: offset of P7_DOMAIN msg for this P7_HIT */
} P7_HIT;

/* Structure: P7_TOPHITS
 * merging when we prepare to output results. "hit" list is NULL and
 * unavailable until after we do a sort.
 */
typedef struct p7_tophits_s {
  P7_HIT **hit;         /* sorted pointer array                     */
  P7_HIT  *unsrt;	/* unsorted data storage                    */
  uint64_t Nalloc;	/* current allocation size                  */
  uint64_t N;		/* number of hits in list now               */
  uint64_t nreported;	/* number of hits that are reportable       */
  uint64_t nincluded;	/* number of hits that are includable       */
  int      is_sorted_by_sortkey; /* TRUE when hits sorted by sortkey and th->hit valid for all N hits */
  int      is_sorted_by_seqidx; /* TRUE when hits sorted by seq_idx, position, and th->hit valid for all N hits */
} P7_TOPHITS;

/*****************************************************************
 * 12. P7_SCOREDATA: data used in diagonal recovery and extension
 *****************************************************************/

enum p7_scoredatatype_e {
  p7_sd_std  = 0,
  p7_sd_fm   = 1,
};

/* This contains a compact representation of 8-bit bias-shifted scores for use in
 * diagonal recovery (standard SSV) and extension (standard and FM-SSV),
 * along with MAXL-associated prefix- and suffix-lengths, and optimal extensions
 * for FM-SSV.
 */
typedef struct p7_scoredata_s {
  int         type;
  int         M;
  union {//implicit (M+1)*K matrix, where M = # states, and K = # characters in alphabet
	uint8_t  *ssv_scores;    // this 2D array is used in the default nhmmer pipeline
	float    *ssv_scores_f;  // this 2D array is used in the FM-index based pipeline
  };
  float      *prefix_lengths;
  float      *suffix_lengths;
  float      *fwd_scores;
  float     **fwd_transitions;
  float     **opt_ext_fwd; // Used only for FM-index based pipeline
  float     **opt_ext_rev; // Used only for FM-index based pipeline
} P7_SCOREDATA;

/*****************************************************************
 * 13. P7_HMM_WINDOW: data used to track lists of sequence windows
 *****************************************************************/

typedef struct p7_hmm_window_s {
  float      score;
  float      null_sc;
  int32_t    id;          //sequence id of the database sequence hit
  int64_t    n;           //position in database sequence at which the diagonal/window starts
  int64_t    fm_n;        //position in the concatenated fm-index sequence at which the diagonal starts
  int32_t    length;      // length of the diagonal/window
  int16_t    k;           //position of the model at which the diagonal ends
  int64_t    target_len;  //length of the target sequence
  int8_t     complementarity;
  int8_t     used_to_extend;
} P7_HMM_WINDOW;

typedef struct p7_hmm_window_list_s {
  P7_HMM_WINDOW *windows;
  int       count;
  int       size;
} P7_HMM_WINDOWLIST;

/*****************************************************************
 * 14. The optimized implementation.
 *****************************************************************/
#if   defined (p7_IMPL_SSE)

#elif defined (p7_IMPL_VMX)

#else

#endif

/*****************************************************************
 * 15. The FM-index acceleration to the SSV filter.  Only works for SSE
 *****************************************************************/
#define FM_MAX_LINE 256

/* Structure the 2D occ array into a single array.  "type" is either b or sb.
 * Note that one extra count value is required by RLE, one 4-byte int for
 * each superblock count vector, and one 2-byte short for each block count
 * vector. This is small overhead, even for a small alphabet like dna.
 */
#define FM_OCC_CNT( type, i, c)  ( occCnts_##type[(meta->alph_size)*(i) + (c)])

enum fm_alphabettypes_e {
  fm_DNA        = 0,  //acgt,  2 bit
  //fm_DNA_full   = 1,  //includes ambiguity codes, 4 bit.
  fm_AMINO      = 4,  // 5 bit
};
/*TODO: fm_DNA_full has currently been disabled because of problems with how the
 * FM index handles very long runs of the same character (in this case, Ns).
 * See wheelert/notebook/2013/12-11-FM-alphabet-speed notes on 12/12.
 */

enum fm_direction_e {
  fm_forward    = 0,
  fm_backward   = 1,
};

typedef struct fm_interval_s {
  int   lower;
  int   upper;
} FM_INTERVAL;

typedef struct fm_hit_s {
  uint64_t  start;
  uint32_t  block;
  int       direction;
  int       length;
  int       sortkey;
} FM_HIT;

typedef struct fm_ambiglist_s {
  FM_INTERVAL *ranges;
  uint32_t     count;
  uint32_t     size;
} FM_AMBIGLIST;

typedef struct fm_seqdata_s {

  uint32_t target_id;      // Which sequence in the target database did this segment come from (can be multiple segment per sequence, if a sequence has Ns)
  uint64_t target_start;   // The position in sequence {id} in the target database at which this sequence-block starts (usually 1, unless its a long sequence split out over multiple FMs)
  uint32_t fm_start;       // The position in the FM block at which this sequence begins
  uint32_t length;         // Length of this sequence segment  (usually the length of the target sequence, unless its a long sequence split out over multiple FMs)

  //meta data taken from the sequence this segment was taken from
  uint16_t name_length;
  uint16_t source_length;
  uint16_t acc_length;
  uint16_t desc_length;
  char     *name;
  char     *source;
  char     *acc;
  char     *desc;
} FM_SEQDATA;

typedef struct fm_metadata_s {
  uint8_t  fwd_only;
  uint8_t  alph_type;
  uint8_t  alph_size;
  uint8_t  charBits;
  uint32_t freq_SA; //frequency with which SA is sampled
  uint32_t freq_cnt_sb; //frequency with which full cumulative counts are captured
  uint32_t freq_cnt_b; //frequency with which intermittent counts are captured
  uint16_t block_count;
  uint32_t seq_count;
  uint64_t char_count; //total count of characters including those in and out of the alphabet
  char     *alph;
  char     *inv_alph;
  int      *compl_alph;
  FILE         *fp;
  FM_SEQDATA   *seq_data;
  FM_AMBIGLIST *ambig_list;
} FM_METADATA;

typedef struct fm_data_s {
  uint64_t N; //length of text
  uint32_t term_loc; // location in the BWT at which the '$' char is found (replaced in the sequence with 'a')
  uint32_t seq_offset;
  uint32_t ambig_offset;
  uint32_t seq_cnt;
  uint32_t ambig_cnt;
  uint32_t overlap; // number of bases at the beginning that overlap the FM-index for the preceding block
  uint8_t  *T;  //text corresponding to the BWT
  uint8_t  *BWT_mem;
  uint8_t  *BWT;
  uint32_t *SA; // sampled suffix array
  int64_t  *C; //the first position of each letter of the alphabet if all of T is sorted.  (signed, as I use that to keep tract of presence/absence)
  uint32_t *occCnts_sb;
  uint16_t *occCnts_b;
} FM_DATA;

typedef struct fm_dp_pair_s {
  uint16_t    pos;  // position of the diagonal in the model.
  float       score;
  float       max_score;
  uint8_t     score_peak_len; // how long was the diagonal when the most recent peak (within fm_drop_lim of the max score) was seen?
  uint8_t     consec_pos;
  uint8_t     max_consec_pos;
  uint8_t     consec_consensus;
  uint8_t     model_direction;
  uint8_t     complementarity;
} FM_DP_PAIR;

typedef struct fm_diag_s {
  uint64_t    n;  //position of the database sequence at which the diagonal starts
  union {
	double    sortkey;
	double    score;
  };
  uint16_t    k;  //position of the model at which the diagonal starts
  uint16_t    length;
  uint8_t     complementarity;
} FM_DIAG;

typedef struct fm_diaglist_s {
  FM_DIAG   *diags;
  int       count;
  int       size;
} FM_DIAGLIST;

/* Effectively global variables, to be initialized once in fm_initConfig(),
 * then passed around among threads to avoid recomputing them
 *
 * When allocated, must be 16-byte aligned, and all _m128i elements
 * must precede other types
 */
typedef struct {
#if   defined (p7_IMPL_SSE)
  /* mask arrays, and 16-byte-offsets into them */
  __m128i *fm_masks_mem;
  __m128i *fm_masks_v;
  __m128i *fm_reverse_masks_mem;
  __m128i *fm_reverse_masks_v;
  __m128i *fm_chars_mem;
  __m128i *fm_chars_v;

  /*various precomputed vectors*/
  __m128i fm_allones_v;
  __m128i fm_zeros_v;
  __m128i fm_neg128_v;
  __m128i fm_m0f;  //00 00 11 11
  __m128i fm_m01;  //01 01 01 01
  __m128i fm_m11;  //00 00 00 11

  /* no non-__m128i- elements above this line */
#endif //#if   defined (p7_IMPL_SSE)

  /*counter, to compute FM-index speed*/
  int occCallCnt;

  /*bounding cutoffs*/
  int max_depth;
  float drop_lim;  // 0.2 ; in seed, max drop in a run of length [fm_drop_max_len]
  int drop_max_len; // 4 ; maximum run length with score under (max - [fm_drop_lim])
  int consec_pos_req; //5
  int consensus_match_req; //10
  float score_density_req; //.5
  int ssv_length;
  float scthreshFM;
  float sc_thresh_ratio; //information content deficit,  actual_relent/target_relent

  /*pointer to FM-index metadata*/
  FM_METADATA *meta;

} FM_CFG;

#if   defined (p7_IMPL_SSE)
//used to convert from a byte array to an __m128i
typedef union {
		uint8_t bytes[16];
		__m128i m128;
		} byte_m128;

/* Gather the sum of all counts in a 16x8-bit element into a single 16-bit
 *  element of the register (the 0th element)
 *
 *  the _mm_sad_epu8  accumulates 8-bit counts into 16-bit counts:
 *      left 8 counts (64-bits) accumulate in counts_v[0],
 *      right 8 counts in counts_v[4]  (the other 6 16-bit ints are 0)
 *  the _mm_shuffle_epi32  flips the 4th int into the 0th slot
 */
#define FM_GATHER_8BIT_COUNTS( in_v, mid_v, out_v  ) do {\
	mid_v = _mm_sad_epu8 (in_v, cfg->fm_zeros_v);\
	tmp_v = _mm_shuffle_epi32(mid_v, _MM_SHUFFLE(1, 1, 1, 2));\
	out_v = _mm_add_epi16(mid_v, tmp_v);\
  } while (0)

/* Macro for SSE operations to turn 2-bit character values into 2-bit binary
 * (00 or 01) match/mismatch values representing occurrences of a character in a
 * 4-char-per-byte packed BWT.
 *
 * Typically followed by a call to FM_COUNT_SSE_4PACKED, possibly with a
 * mask in between to handle the case where we don't want to add over all
 * positions in the vector
 *
 * tmp_v and tmp2_v are used as temporary vectors throughout, and hold meaningless values
 * at the end
 *
 * xor(in_v, c_v)        : each 2-bit value will be 00 if a match, and non-0 if a mismatch
 * and(in_v, 01010101)   : look at the right bit of each 2-bit value,
 * srli(1)+and()         : look at the left bit of each 2-bit value,
 * or()                  : if either left bit or right bit is non-0, 01, else 00 (match is 00)
 *
 * subs()                : invert, so match is 01, mismatch is 00
 *
 */
#define FM_MATCH_2BIT(in_v, c_v, a_v, b_v, out_v) do {\
	a_v = _mm_xor_si128(in_v, c_v);\
	\
	b_v  = _mm_srli_epi16(a_v, 1);\
	a_v  = _mm_or_si128(a_v, b_v);\
	b_v  = _mm_and_si128(a_v, cfg->fm_m01);\
	\
	out_v  = _mm_subs_epi8(cfg->fm_m01,b_v);\
  } while (0)

/*Macro for SSE operations to count bits produced by FM_MATCH_SSE_4PACKED
 *
 * tmp_v and tmp2_v are used as temporary vectors throughout, and hold meaningless values
 * at the end
 *
 * then add up the 2-bit values:
 * srli(4)+add()         : left 4 bits shifted right, added to right 4 bits
 *
 * srli(2)+and(00000011) : left 2 bits (value 0..2) shifted right, masked, so no other bits active
 * and(00000011)         : right 2 bits (value 0..2) masked so no other bits active
 *
 * final 2 add()s        : tack current counts on to already-tabulated counts.
 */
#define FM_COUNT_2BIT(a_v, b_v, cnts_v) do {\
		b_v = _mm_srli_epi16(a_v, 4);\
		a_v  = _mm_add_epi16(a_v, b_v);\
		\
		b_v = _mm_srli_epi16(a_v, 2);\
		a_v  = _mm_and_si128(a_v,cfg->fm_m11);\
		b_v = _mm_and_si128(b_v,cfg->fm_m11);\
		\
		cnts_v = _mm_add_epi16(cnts_v, a_v);\
		cnts_v = _mm_add_epi16(cnts_v, b_v);\
  } while (0)

/* Macro for SSE operations that turns a vector of 4-bit character values into
 * 2 vectors representing matches. Each byte in the input vector consists of
 * a left half (4 bits) and a right half (4 bits). The 16 left-halves produce
 * one vector, which contains all-1s for bytes in which the left half matches
 * the c_v character (and 0s if it doesn't), while the 16 right-halves produce
 * the other vector, again with each byte either all-1s or all-0s.
 *
 * The expectation is that FM_COUNT_4BIT will be called after this, to
 * turn these binary values into sums over a series of vectors. The macros
 * are split up to allow one end or other to be trimmed in the case that
 * counting is not expected to include the full vector.
 *
 * srli(4)+and() : capture the left 4-bit value   (need the mask because 16-bit shift leaves garbage in left-4-bit chunks)
 * and()         : capture the right 4-bit value
 *
 * cmpeq()x2     : test if both left and right == c.  For each, if ==c , value = 11111111 (-1)
 */
#define FM_MATCH_4BIT(in_v, c_v, out1_v, out2_v) do {\
	out1_v    = _mm_srli_epi16(in_v, 4);\
	out2_v    = _mm_and_si128(in_v, cfg->fm_m0f);\
	out1_v    = _mm_and_si128(out1_v, cfg->fm_m0f);\
	\
	out1_v    = _mm_cmpeq_epi8(out1_v, c_v);\
	out2_v    = _mm_cmpeq_epi8(out2_v, c_v);\
  } while (0)

/* Macro for SSE operations that turns a vector of 4-bit character values into
 * 2 vectors representing matches. Each byte in the input vector consists of
 * a left half (4 bits) and a right half (4 bits). The 16 left-halves produce
 * one vector, which contains all-1s for bytes in which the left half is less than
 * the c_v character (and 0s if it doesn't), while the 16 right-halves produce
 * the other vector, again with each byte either all-1s or all-0s.
 *
 * The expectation is that FM_COUNT_4BIT will be called after this, to
 * turn these binary values into sums over a series of vectors. The macros
 * are split up to allow one end or other to be trimmed in the case that
 * counting is not expected to include the full vector.
 *
 * srli(4)+and() : capture the left 4-bit value   (need the mask because 16-bit shift leaves garbage in left-4-bit chunks)
 * and()         : capture the right 4-bit value
 *
 * cmplt()x2     : test if both left and right < c.  For each, if <c , value = 11111111 (-1)
 */
#define FM_LT_4BIT(in_v, c_v, out1_v, out2_v) do {\
	out1_v    = _mm_srli_epi16(in_v, 4);\
	out2_v    = _mm_and_si128(in_v, cfg->fm_m0f);\
	out1_v    = _mm_and_si128(out1_v, cfg->fm_m0f);\
	\
	out1_v    = _mm_cmplt_epi8(out1_v, c_v);\
	out2_v    = _mm_cmplt_epi8(out2_v, c_v);\
  } while (0)

/* Macro for SSE operations to add occurrence counts to the tally vector counts_v,
 * in the 4-bits-per-character case
 *
 * The expectation is that in[12]_v will contain bytes that are either
 *   00000000  =  0
 *  or
 *   11111111  = -1
 * so subtracting the value of the byte is the same as adding 0 or 1.
 */
#define FM_COUNT_4BIT(in1_v, in2_v, cnts_v) do {\
	cnts_v = _mm_subs_epi8(cnts_v, in1_v);\
	cnts_v = _mm_subs_epi8(cnts_v, in2_v);\
  } while (0)

#endif  // if  defined (p7_IMPL_SSE)

/*****************************************************************
 * 16. P7_PIPELINE: H3's accelerated seq/profile comparison pipeline
 *****************************************************************/

enum p7_pipemodes_e { p7_SEARCH_SEQS = 0, p7_SCAN_MODELS = 1 };
enum p7_zsetby_e    { p7_ZSETBY_NTARGETS = 0, p7_ZSETBY_OPTION = 1, p7_ZSETBY_FILEINFO = 2 };
enum p7_complementarity_e { p7_NOCOMPLEMENT    = 0, p7_COMPLEMENT   = 1 };

typedef struct p7_pipeline_s {
  /* Dynamic programming matrices                                           */
  P7_OMX     *oxf;		/* one-row Forward matrix, accel pipe       */
  P7_OMX     *oxb;		/* one-row Backward matrix, accel pipe      */
  P7_OMX     *fwd;		/* full Fwd matrix for domain envelopes     */
  P7_OMX     *bck;		/* full Bck matrix for domain envelopes     */

  /* Domain postprocessing                                                  */
  ESL_RANDOMNESS *r;		/* random number generator                  */
  int             do_reseeding; /* TRUE: reseed for reproducible results    */
  int             do_alignment_score_calc;
  P7_DOMAINDEF   *ddef;		/* domain definition workflow               */

  /* Reporting threshold settings                                           */
  int     by_E;		        /* TRUE to cut per-target report off by E   */
  double  E;	                /* per-target E-value threshold             */
  double  T;	                /* per-target bit score threshold           */
  int     dom_by_E;             /* TRUE to cut domain reporting off by E    */
  double  domE;	                /* domain E-value threshold                 */
  double  domT;	                /* domain bit score threshold               */
  int     use_bit_cutoffs;      /* (FALSE | p7H_GA | p7H_TC | p7H_NC)       */

  /* Inclusion threshold settings                                           */
  int     inc_by_E;		/* TRUE to threshold inclusion by E-values  */
  double  incE;			/* per-target inclusion E-value threshold   */
  double  incT;			/* per-target inclusion score threshold     */
  int     incdom_by_E;		/* TRUE to threshold domain inclusion by E  */
  double  incdomE;		/* per-domain inclusion E-value threshold   */
  double  incdomT;		/* per-domain inclusion E-value threshold   */

  /* Tracking search space sizes for E value calculations                   */
  double  Z;			/* eff # targs searched (per-target E-val)  */
  double  domZ;			/* eff # signific targs (per-domain E-val)  */
  enum p7_zsetby_e Z_setby;   	/* how Z was set                            */
  enum p7_zsetby_e domZ_setby;	/* how domZ was set                         */

  /* Threshold settings for pipeline                                        */
  int     do_max;	        /* TRUE to run in slow/max mode             */
  double  F1;		        /* MSV filter threshold                     */
  double  F2;		        /* Viterbi filter threshold                 */
  double  F3;		        /* uncorrected Forward filter threshold     */
  int     B1;               /* window length for biased-composition modifier - MSV*/
  int     B2;               /* window length for biased-composition modifier - Viterbi*/
  int     B3;               /* window length for biased-composition modifier - Forward*/
  int     do_biasfilter;	/* TRUE to use biased comp HMM filter       */
  int     do_null2;		/* TRUE to use null2 score corrections      */

  /* Accounting. (reduceable in threaded/MPI parallel version)              */
  uint64_t      nmodels;        /* # of HMMs searched                       */
  uint64_t      nseqs;	        /* # of sequences searched                  */
  uint64_t      nres;	        /* # of residues searched                   */
  uint64_t      nnodes;	        /* # of model nodes searched                */
  uint64_t      n_past_msv;	/* # comparisons that pass MSVFilter()      */
  uint64_t      n_past_bias;	/* # comparisons that pass bias filter      */
  uint64_t      n_past_vit;	/* # comparisons that pass ViterbiFilter()  */
  uint64_t      n_past_fwd;	/* # comparisons that pass ForwardFilter()  */
  uint64_t      n_output;	    /* # alignments that make it to the final output (used for nhmmer) */
  uint64_t      pos_past_msv;	/* # positions that pass MSVFilter()  (used for nhmmer) */
  uint64_t      pos_past_bias;	/* # positions that pass bias filter  (used for nhmmer) */
  uint64_t      pos_past_vit;	/* # positions that pass ViterbiFilter()  (used for nhmmer) */
  uint64_t      pos_past_fwd;	/* # positions that pass ForwardFilter()  (used for nhmmer) */
  uint64_t      pos_output;	    /* # positions that make it to the final output (used for nhmmer) */

  enum p7_pipemodes_e mode;    	/* p7_SCAN_MODELS | p7_SEARCH_SEQS          */
  int           long_targets;   /* TRUE if the target sequences are expected to be very long (e.g. dna chromosome search in nhmmer) */
  int           strands;         /*  p7_STRAND_TOPONLY  | p7_STRAND_BOTTOMONLY |  p7_STRAND_BOTH */
  int 		    	W;              /* window length for nhmmer scan - essentially maximum length of model that we expect to find*/
  int           block_length;   /* length of overlapping blocks read in the multi-threaded variant (default MAX_RESIDUE_COUNT) */

  int           show_accessions;/* TRUE to output accessions not names      */
  int           show_alignments;/* TRUE to output alignments (default)      */

  P7_HMMFILE   *hfp;		/* COPY of open HMM database (if scan mode) */
  char          errbuf[eslERRBUFSIZE];
} P7_PIPELINE;

/*****************************************************************
 * 17. P7_BUILDER: pipeline for new HMM construction
 *****************************************************************/

#define p7_DEFAULT_WINDOW_BETA  1e-7

enum p7_archchoice_e { p7_ARCH_FAST = 0, p7_ARCH_HAND = 1 };
enum p7_wgtchoice_e  { p7_WGT_NONE  = 0, p7_WGT_GIVEN = 1, p7_WGT_GSC    = 2, p7_WGT_PB       = 3, p7_WGT_BLOSUM = 4 };
enum p7_effnchoice_e { p7_EFFN_NONE = 0, p7_EFFN_SET  = 1, p7_EFFN_CLUST = 2, p7_EFFN_ENTROPY = 3, p7_EFFN_ENTROPY_EXP = 4 };

typedef struct p7_builder_s {
  /* Model architecture                                                                            */
  enum p7_archchoice_e arch_strategy;    /* choice of model architecture determination algorithm   */
  float                symfrac;	         /* residue occ thresh for fast architecture determination */
  float                fragthresh;	 /* if L <= fragthresh*alen, seq is called a fragment      */

  /* Relative sequence weights                                                                     */
  enum p7_wgtchoice_e  wgt_strategy;     /* choice of relative sequence weighting algorithm        */
  double               wid;		 /* %id threshold for BLOSUM relative weighting            */

  /* Effective sequence number                                                                     */
  enum p7_effnchoice_e effn_strategy;    /* choice of effective seq # determination algorithm      */
  double               re_target;	 /* rel entropy target for effn eweighting, if set; or -1.0*/
  double               esigma;		 /* min total rel ent parameter for effn entropy weights   */
  double               eid;		 /* %id threshold for effn clustering                      */
  double               eset;		 /* effective sequence number, if --eset; or -1.0          */

  /* Run-to-run variation due to random number generation                                          */
  ESL_RANDOMNESS      *r;	         /* RNG for E-value calibration simulations                */
  int                  do_reseeding;	 /* TRUE to reseed, making results reproducible            */

  /* E-value parameter calibration                                                                 */
  int                  EmL;            	 /* length of sequences generated for MSV fitting          */
  int                  EmN;	         /* # of sequences generated for MSV fitting               */
  int                  EvL;            	 /* length of sequences generated for Viterbi fitting      */
  int                  EvN;	         /* # of sequences generated for Viterbi fitting           */
  int                  EfL;	         /* length of sequences generated for Forward fitting      */
  int                  EfN;	         /* # of sequences generated for Forward fitting           */
  double               Eft;	         /* tail mass used for Forward fitting                     */

  /* Choice of prior                                                                               */
  P7_PRIOR            *prior;	         /* choice of prior when parameterizing from counts        */
  int                  max_insert_len;

  /* Optional: information used for parameterizing single sequence queries                         */
  ESL_SCOREMATRIX     *S;		 /* residue score matrix                                   */
  ESL_DMATRIX         *Q;	         /* Q->mx[a][b] = P(b|a) residue probabilities             */
  double               popen;         	 /* gap open probability                                   */
  double               pextend;          /* gap extend probability                                 */

  double               w_beta;    /*beta value used to compute W (window length)   */
  int                  w_len;     /*W (window length)  explicitly set */

  const ESL_ALPHABET  *abc;		 /* COPY of alphabet                                       */
  char errbuf[eslERRBUFSIZE];            /* informative message on model construction failure      */
} P7_BUILDER;

/*****************************************************************
 * 18. Routines in HMMER's exposed API.
 *****************************************************************/

/* build.c */
int p7_Handmodelmaker(ESL_MSA *msa,                P7_BUILDER *bld, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);
int p7_Fastmodelmaker(ESL_MSA *msa, float symfrac, P7_BUILDER *bld, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);

/* emit.c */
int p7_CoreEmit   (ESL_RANDOMNESS *r, const P7_HMM *hmm,                                        ESL_SQ *sq, P7_TRACE *tr);
int p7_ProfileEmit(ESL_RANDOMNESS *r, const P7_HMM *hmm, const P7_PROFILE *gm, const P7_BG *bg, ESL_SQ *sq, P7_TRACE *tr);
int p7_emit_SimpleConsensus(const P7_HMM *hmm, ESL_SQ *sq);
int p7_emit_FancyConsensus (const P7_HMM *hmm, float min_lower, float min_upper, ESL_SQ *sq);

/* errors.c */
void p7_Die (char *format, ...);
void p7_Fail(char *format, ...);

/* evalues.c */
int p7_Calibrate(P7_HMM *hmm, P7_BUILDER *cfg_b, ESL_RANDOMNESS **byp_rng, P7_BG **byp_bg, P7_PROFILE **byp_gm, P7_OPROFILE **byp_om);
int p7_Lambda(P7_HMM *hmm, P7_BG *bg, double *ret_lambda);
int p7_MSVMu     (ESL_RANDOMNESS *r, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda,               double *ret_mmu);
int p7_ViterbiMu (ESL_RANDOMNESS *r, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda,               double *ret_vmu);
int p7_Tau       (ESL_RANDOMNESS *r, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda, double tailp, double *ret_tau);

/* eweight.c */
int p7_EntropyWeight(const P7_HMM *hmm, const P7_BG *bg, const P7_PRIOR *pri, double infotarget, double *ret_Neff);

int p7_EntropyWeight_exp(const P7_HMM *hmm, const P7_BG *bg, const P7_PRIOR *pri, double etarget, double *ret_exp);
/* generic_decoding.c */
int p7_GDecoding      (const P7_PROFILE *gm, const P7_GMX *fwd,       P7_GMX *bck, P7_GMX *pp);
int p7_GDomainDecoding(const P7_PROFILE *gm, const P7_GMX *fwd, const P7_GMX *bck, P7_DOMAINDEF *ddef);

/* generic_fwdback.c */
int p7_GForward     (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *ret_sc);
int p7_GBackward    (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *ret_sc);
int p7_GHybrid      (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *opt_fwdscore, float *opt_hybscore);

/* generic_msv.c */
int p7_GMSV           (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *gx, float nu, float *ret_sc);
int p7_GMSV_longtarget(const ESL_DSQ *dsq, int L, P7_PROFILE *gm, P7_GMX *gx, float nu,  P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist);

/* generic_null2.c */
int p7_GNull2_ByExpectation(const P7_PROFILE *gm, P7_GMX *pp, float *null2);
int p7_GNull2_ByTrace      (const P7_PROFILE *gm, const P7_TRACE *tr, int zstart, int zend, P7_GMX *wrk, float *null2);

/* generic_optacc.c */
int p7_GOptimalAccuracy(const P7_PROFILE *gm, const P7_GMX *pp,       P7_GMX *gx, float *ret_e);
int p7_GOATrace        (const P7_PROFILE *gm, const P7_GMX *pp, const P7_GMX *gx, P7_TRACE *tr);

/* generic_stotrace.c */
int p7_GStochasticTrace(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *gx, P7_TRACE *tr);

/* generic_viterbi.c */
int p7_GViterbi     (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *ret_sc);

/* generic_vtrace.c */
int p7_GTrace       (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *gx, P7_TRACE *tr);
int p7_GViterbi_longtarget(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *gx,
					   float filtersc, double P, P7_HMM_WINDOWLIST *windowlist);

/* heatmap.c (evolving now, intend to move this to Easel in the future) */
double dmx_upper_max(ESL_DMATRIX *D);
double dmx_upper_min(ESL_DMATRIX *D);
double dmx_upper_element_sum(ESL_DMATRIX *D);
double dmx_upper_norm(ESL_DMATRIX *D);
int    dmx_Visualize(FILE *fp, ESL_DMATRIX *D, double min, double max);

/* hmmdutils.c */
void p7_openlog(const char *ident, int option, int facility);
void p7_syslog(int priority, const char *format, ...);
void p7_closelog(void);

/* hmmlogo.c */
float hmmlogo_maxHeight (P7_BG *bg);
int hmmlogo_RelativeEntropy_all      (P7_HMM *hmm, P7_BG *bg, float *rel_ents, float **probs, float **heights );
int hmmlogo_RelativeEntropy_above_bg (P7_HMM *hmm, P7_BG *bg, float *rel_ents, float **probs, float **heights );
int hmmlogo_ScoreHeights (P7_HMM *hmm, P7_BG *bg, float **heights );
int hmmlogo_IndelValues (P7_HMM *hmm, float *insert_P, float *insert_expL, float *delete_P );

/* hmmpgmd2msa.c */
int hmmpgmd2msa(void *data, P7_HMM *hmm, ESL_SQ *qsq,  int *incl, int incl_size, int *excl, int excl_size, ESL_MSA **ret_msa);

/* island.c */
int   p7_island_Viterbi(ESL_DSQ *dsq, int L, P7_PROFILE *gm, P7_GMX *mx, ESL_HISTOGRAM *h);

/* h2_io.c */
int   p7_h2io_WriteASCII(FILE *fp, P7_HMM *hmm);

/* hmmer.c */
void         p7_banner(FILE *fp, char *progname, char *banner);
ESL_GETOPTS *p7_CreateDefaultApp(ESL_OPTIONS *options, int nargs, int argc, char **argv, char *banner, char *usage);
int          p7_AminoFrequencies(float *f);

/* logsum.c */
int   p7_FLogsumInit(void);
float p7_FLogsum(float a, float b);
float p7_FLogsumError(float a, float b);
int   p7_ILogsumInit(void);
int   p7_ILogsum(int s1, int s2);

/* modelconfig.c */
int p7_ProfileConfig(const P7_HMM *hmm, const P7_BG *bg, P7_PROFILE *gm, int L, int mode);
int p7_ReconfigLength  (P7_PROFILE *gm, int L);
int p7_ReconfigMultihit(P7_PROFILE *gm, int L);
int p7_ReconfigUnihit  (P7_PROFILE *gm, int L);

/* modelstats.c */
double p7_MeanMatchInfo           (const P7_HMM *hmm, const P7_BG *bg);
double p7_MeanMatchEntropy        (const P7_HMM *hmm);
double p7_MeanMatchRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg);
double p7_MeanForwardScore        (const P7_HMM *hmm, const P7_BG *bg);
int    p7_MeanPositionRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg, double *ret_entropy);
int    p7_hmm_CompositionKLDist(P7_HMM *hmm, P7_BG *bg, float *ret_KL, float **opt_avp);

/* mpisupport.c */
#ifdef HAVE_MPI
int p7_hmm_MPISend(P7_HMM *hmm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
int p7_hmm_MPIPackSize(P7_HMM *hmm, MPI_Comm comm, int *ret_n);
int p7_hmm_MPIPack(P7_HMM *hmm, char *buf, int n, int *position, MPI_Comm comm);
int p7_hmm_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **abc, P7_HMM **ret_hmm);
int p7_hmm_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **abc, P7_HMM **ret_hmm);

int p7_profile_MPISend(P7_PROFILE *gm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
int p7_profile_MPIRecv(int source, int tag, MPI_Comm comm, const ESL_ALPHABET *abc, const P7_BG *bg,
			      char **buf, int *nalloc,  P7_PROFILE **ret_gm);

int p7_pipeline_MPISend(P7_PIPELINE *pli, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
int p7_pipeline_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_GETOPTS *go, P7_PIPELINE **ret_pli);

int p7_tophits_MPISend(P7_TOPHITS *th, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
int p7_tophits_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, P7_TOPHITS **ret_th);

int p7_oprofile_MPISend(P7_OPROFILE *om, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
int p7_oprofile_MPIPackSize(P7_OPROFILE *om, MPI_Comm comm, int *ret_n);
int p7_oprofile_MPIPack(P7_OPROFILE *om, char *buf, int n, int *pos, MPI_Comm comm);
int p7_oprofile_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **abc, P7_OPROFILE **ret_om);
int p7_oprofile_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **abc, P7_OPROFILE **ret_om);
#endif /*HAVE_MPI*/

/* tracealign.c */
int p7_tracealign_Seqs(ESL_SQ **sq,           P7_TRACE **tr, int nseq, int M,  int optflags, P7_HMM *hmm, ESL_MSA **ret_msa);
int p7_tracealign_MSA (const ESL_MSA *premsa, P7_TRACE **tr,           int M,  int optflags, ESL_MSA **ret_postmsa);
int p7_tracealign_computeTraces(P7_HMM *hmm, ESL_SQ  **sq, int offset, int N, P7_TRACE  **tr);
int p7_tracealign_getMSAandStats(P7_HMM *hmm, ESL_SQ  **sq, int N, ESL_MSA **ret_msa, float **ret_pp, float **ret_relent, float **ret_scores );

/* p7_alidisplay.c */
P7_ALIDISPLAY *p7_alidisplay_Create(const P7_TRACE *tr, int which, const P7_OPROFILE *om, const ESL_SQ *sq, const ESL_SQ *ntsq);
P7_ALIDISPLAY *p7_alidisplay_Clone(const P7_ALIDISPLAY *ad);
size_t         p7_alidisplay_Sizeof(const P7_ALIDISPLAY *ad);
int            p7_alidisplay_Serialize(P7_ALIDISPLAY *ad);
int            p7_alidisplay_Deserialize(P7_ALIDISPLAY *ad);
void           p7_alidisplay_Destroy(P7_ALIDISPLAY *ad);
char           p7_alidisplay_EncodePostProb(float p);
float          p7_alidisplay_DecodePostProb(char pc);
char           p7_alidisplay_EncodeAliPostProb(float p, float hi, float med, float lo);

int            p7_alidisplay_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth, P7_PIPELINE *pli);
int            p7_translated_alidisplay_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth, P7_PIPELINE *pli);
int            p7_nontranslated_alidisplay_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth, int show_accessions);

int            p7_alidisplay_Backconvert(const P7_ALIDISPLAY *ad, const ESL_ALPHABET *abc, ESL_SQ **ret_sq, P7_TRACE **ret_tr);
int            p7_alidisplay_Dump(FILE *fp, const P7_ALIDISPLAY *ad);
int            p7_alidisplay_Compare(const P7_ALIDISPLAY *ad1, const P7_ALIDISPLAY *ad2);

/* p7_bg.c */
P7_BG *p7_bg_Create(const ESL_ALPHABET *abc);
P7_BG *p7_bg_CreateUniform(const ESL_ALPHABET *abc);
P7_BG *p7_bg_Clone(const P7_BG *bg);
int    p7_bg_Dump(FILE *ofp, const P7_BG *bg);
void   p7_bg_Destroy(P7_BG *bg);
int    p7_bg_SetLength(P7_BG *bg, int L);
int    p7_bg_NullOne(const P7_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc);

int    p7_bg_Read(char *bgfile, P7_BG *bg, char *errbuf);
int    p7_bg_Write(FILE *fp, P7_BG *bg);

int    p7_bg_SetFilter  (P7_BG *bg, int M, const float *compo);
int    p7_bg_FilterScore(P7_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc);

/* p7_builder.c */
P7_BUILDER *p7_builder_Create(const ESL_GETOPTS *go, const ESL_ALPHABET *abc);
int         p7_builder_LoadScoreSystem(P7_BUILDER *bld, const char *matrix,                  double popen, double pextend, P7_BG *bg);
int         p7_builder_SetScoreSystem (P7_BUILDER *bld, const char *mxfile, const char *env, double popen, double pextend, P7_BG *bg);
void        p7_builder_Destroy(P7_BUILDER *bld);

int p7_Builder      (P7_BUILDER *bld, ESL_MSA *msa, P7_BG *bg, P7_HMM **opt_hmm, P7_TRACE ***opt_trarr, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om, ESL_MSA **opt_postmsa, FILE *seqweights_w_fp, FILE *seqweights_e_fp);
int p7_SingleBuilder(P7_BUILDER *bld, ESL_SQ *sq,   P7_BG *bg, P7_HMM **opt_hmm, P7_TRACE  **opt_tr,    P7_PROFILE **opt_gm, P7_OPROFILE **opt_om);
int p7_Builder_MaxLength      (P7_HMM *hmm, double emit_thresh);

/* p7_domaindef.c */
P7_DOMAINDEF *p7_domaindef_Create (ESL_RANDOMNESS *r);
int           p7_domaindef_Fetch  (P7_DOMAINDEF *ddef, int which, int *opt_i, int *opt_j, float *opt_sc, P7_ALIDISPLAY **opt_ad);
int           p7_domaindef_GrowTo (P7_DOMAINDEF *ddef, int L);
int           p7_domaindef_Reuse  (P7_DOMAINDEF *ddef);
int           p7_domaindef_DumpPosteriors(FILE *ofp, P7_DOMAINDEF *ddef);
void          p7_domaindef_Destroy(P7_DOMAINDEF *ddef);

int p7_domaindef_ByViterbi            (P7_PROFILE *gm, const ESL_SQ *sq, const ESL_SQ *ntsq, P7_GMX *gx1, P7_GMX *gx2, P7_DOMAINDEF *ddef);
int p7_domaindef_ByPosteriorHeuristics(const ESL_SQ *sq, const ESL_SQ *ntsq, P7_OPROFILE *om, P7_OMX *oxf, P7_OMX *oxb, P7_OMX *fwd, P7_OMX *bck,
				                                  P7_DOMAINDEF *ddef, P7_BG *bg, int long_target,
				                                  P7_BG *bg_tmp, float *scores_arr, float *fwd_emissions_arr);

/* p7_gmx.c */
P7_GMX *p7_gmx_Create (int allocM, int allocL);
int     p7_gmx_GrowTo (P7_GMX *gx, int allocM, int allocL);
size_t  p7_gmx_Sizeof (P7_GMX *gx);
int     p7_gmx_Reuse  (P7_GMX *gx);
void    p7_gmx_Destroy(P7_GMX *gx);
int     p7_gmx_Compare(P7_GMX *gx1, P7_GMX *gx2, float tolerance);
int     p7_gmx_Dump(FILE *fp, P7_GMX *gx, int flags);
int     p7_gmx_DumpWindow(FILE *fp, P7_GMX *gx, int istart, int iend, int kstart, int kend, int show_specials);

/* p7_hmm.c */
/*      1. The P7_HMM object: allocation, initialization, destruction. */
P7_HMM *p7_hmm_Create(int M, const ESL_ALPHABET *abc);
P7_HMM *p7_hmm_CreateShell(void);
int     p7_hmm_CreateBody(P7_HMM *hmm, int M, const ESL_ALPHABET *abc);
void    p7_hmm_Destroy(P7_HMM *hmm);
int     p7_hmm_CopyParameters(const P7_HMM *src, P7_HMM *dest);
P7_HMM *p7_hmm_Clone(const P7_HMM *hmm);
int     p7_hmm_Zero(P7_HMM *hmm);
char    p7_hmm_EncodeStatetype(char *typestring);
char   *p7_hmm_DecodeStatetype(char st);
/*      2. Convenience routines for setting fields in an HMM. */
int     p7_hmm_SetName       (P7_HMM *hmm, char *name);
int     p7_hmm_SetAccession  (P7_HMM *hmm, char *acc);
int     p7_hmm_SetDescription(P7_HMM *hmm, char *desc);
int     p7_hmm_AppendComlog  (P7_HMM *hmm, int argc, char **argv);
int     p7_hmm_SetCtime      (P7_HMM *hmm);
int     p7_hmm_SetComposition(P7_HMM *hmm);
int     p7_hmm_SetConsensus  (P7_HMM *hmm, ESL_SQ *sq);
/*      3. Renormalization and rescaling counts in core HMMs. */
int     p7_hmm_Scale      (P7_HMM *hmm, double scale);
int     p7_hmm_ScaleExponential(P7_HMM *hmm, double exp);
int     p7_hmm_Renormalize(P7_HMM *hmm);
/*      4. Debugging and development code. */
int     p7_hmm_Dump(FILE *fp, P7_HMM *hmm);
int     p7_hmm_Sample          (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm);
int     p7_hmm_SampleUngapped  (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm);
int     p7_hmm_SampleEnumerable(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm);
int     p7_hmm_SampleUniform   (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc,
				     float tmi, float tii, float tmd, float tdd,  P7_HMM **ret_hmm);
int     p7_hmm_Compare(P7_HMM *h1, P7_HMM *h2, float tol);
int     p7_hmm_Validate(P7_HMM *hmm, char *errbuf, float tol);
/*      5. Other routines in the API */
int     p7_hmm_CalculateOccupancy(const P7_HMM *hmm, float *mocc, float *iocc);

/* p7_hmmfile.c */
int  p7_hmmfile_OpenE    (char *filename, char *env, P7_HMMFILE **ret_hfp, char *errbuf);
int  p7_hmmfile_OpenENoDB(char *filename, char *env, P7_HMMFILE **ret_hfp, char *errbuf);
int  p7_hmmfile_Open     (char *filename, char *env, P7_HMMFILE **ret_hfp); /* deprecated */
int  p7_hmmfile_OpenNoDB (char *filename, char *env, P7_HMMFILE **ret_hfp); /* deprecated */
int  p7_hmmfile_OpenBuffer(char *buffer, int size, P7_HMMFILE **ret_hfp);
void p7_hmmfile_Close(P7_HMMFILE *hfp);
#ifdef HMMER_THREADS
int  p7_hmmfile_CreateLock(P7_HMMFILE *hfp);
#endif
int  p7_hmmfile_WriteBinary(FILE *fp, int format, P7_HMM *hmm);
int  p7_hmmfile_WriteASCII (FILE *fp, int format, P7_HMM *hmm);
int  p7_hmmfile_WriteToString (char **s, int format, P7_HMM *hmm);
int  p7_hmmfile_Read(P7_HMMFILE *hfp, ESL_ALPHABET **ret_abc,  P7_HMM **opt_hmm);
int  p7_hmmfile_PositionByKey(P7_HMMFILE *hfp, const char *key);
int  p7_hmmfile_Position(P7_HMMFILE *hfp, const off_t offset);

/* p7_hmmwindow.c */
int p7_hmmwindow_init (P7_HMM_WINDOWLIST *list);
P7_HMM_WINDOW *p7_hmmwindow_new (P7_HMM_WINDOWLIST *list, uint32_t id, uint32_t pos, uint32_t fm_pos, uint16_t k, uint32_t length, float score, uint8_t complementarity, uint32_t target_len);

/* p7_msvdata.c */
P7_SCOREDATA   *p7_hmm_ScoreDataCreate(P7_OPROFILE *om, P7_PROFILE *gm );
P7_SCOREDATA   *p7_hmm_ScoreDataClone(P7_SCOREDATA *src, int K);
int            p7_hmm_ScoreDataComputeRest(P7_OPROFILE *om, P7_SCOREDATA *data );
void           p7_hmm_ScoreDataDestroy( P7_SCOREDATA *data );
int            p7_hmm_initWindows (P7_HMM_WINDOWLIST *list);
P7_HMM_WINDOW *p7_hmm_newWindow (P7_HMM_WINDOWLIST *list, uint32_t id, uint32_t pos, uint32_t fm_pos, uint16_t k, uint32_t length, float score, uint8_t complementarity);

/* p7_null3.c */
void p7_null3_score(const ESL_ALPHABET *abc, const ESL_DSQ *dsq, P7_TRACE *tr, int start, int stop, P7_BG *bg, float *ret_sc);
void p7_null3_windowed_score(const ESL_ALPHABET *abc, const ESL_DSQ *dsq, int start, int stop, P7_BG *bg, float *ret_sc);

/* p7_pipeline.c */
P7_PIPELINE *p7_pipeline_Create(ESL_GETOPTS *go, int M_hint, int L_hint, int long_targets, enum p7_pipemodes_e mode);
int          p7_pipeline_Reuse  (P7_PIPELINE *pli);
void         p7_pipeline_Destroy(P7_PIPELINE *pli);
int          p7_pipeline_Merge  (P7_PIPELINE *p1, P7_PIPELINE *p2);

int p7_pli_ExtendAndMergeWindows (P7_OPROFILE *om, const P7_SCOREDATA *msvdata, P7_HMM_WINDOWLIST *windowlist, float pct_overlap);
int p7_pli_TargetReportable  (P7_PIPELINE *pli, float score,     double lnP);
int p7_pli_DomainReportable  (P7_PIPELINE *pli, float dom_score, double lnP);

int p7_pli_TargetIncludable  (P7_PIPELINE *pli, float score,     double lnP);
int p7_pli_DomainIncludable  (P7_PIPELINE *pli, float dom_score, double lnP);
int p7_pli_NewModel          (P7_PIPELINE *pli, const P7_OPROFILE *om, P7_BG *bg);
int p7_pli_NewModelThresholds(P7_PIPELINE *pli, const P7_OPROFILE *om);
int p7_pli_NewSeq            (P7_PIPELINE *pli, const ESL_SQ *sq);
int p7_Pipeline              (P7_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, const ESL_SQ *ntsq, P7_TOPHITS *th);
int p7_Pipeline_LongTarget   (P7_PIPELINE *pli, P7_OPROFILE *om, P7_SCOREDATA *data,
									 P7_BG *bg, P7_TOPHITS *hitlist, int64_t seqidx,
									 const ESL_SQ *sq, int complementarity,
									 const FM_DATA *fmf, const FM_DATA *fmb, FM_CFG *fm_cfg
/*                                     , ESL_STOPWATCH *ssv_watch_master
									 , ESL_STOPWATCH *postssv_watch_master
									 , ESL_STOPWATCH *watch_slave
*/
									 );

int p7_pli_Statistics(FILE *ofp, P7_PIPELINE *pli, ESL_STOPWATCH *w);

/* p7_prior.c */
P7_PRIOR  *p7_prior_CreateAmino(void);
P7_PRIOR  *p7_prior_CreateNucleic(void);
P7_PRIOR  *p7_prior_CreateLaplace(const ESL_ALPHABET *abc);
void       p7_prior_Destroy(P7_PRIOR *pri);

int        p7_ParameterEstimation(P7_HMM *hmm, const P7_PRIOR *pri);

/* p7_profile.c */
P7_PROFILE *p7_profile_Create(int M, const ESL_ALPHABET *abc);
P7_PROFILE *p7_profile_Clone(const P7_PROFILE *gm);
int         p7_profile_Copy(const P7_PROFILE *src, P7_PROFILE *dst);
int         p7_profile_SetNullEmissions(P7_PROFILE *gm);
int         p7_profile_Reuse(P7_PROFILE *gm);
size_t      p7_profile_Sizeof(P7_PROFILE *gm);
void        p7_profile_Destroy(P7_PROFILE *gm);
int         p7_profile_IsLocal(const P7_PROFILE *gm);
int         p7_profile_IsMultihit(const P7_PROFILE *gm);
int         p7_profile_GetT(const P7_PROFILE *gm, char st1, int k1,
				   char st2, int k2, float *ret_tsc);
int         p7_profile_Validate(const P7_PROFILE *gm, char *errbuf, float tol);
int         p7_profile_Compare(P7_PROFILE *gm1, P7_PROFILE *gm2, float tol);

/* p7_spensemble.c */
P7_SPENSEMBLE *p7_spensemble_Create(int init_n, int init_epc, int init_sigc);
int     p7_spensemble_Reuse(P7_SPENSEMBLE *sp);
int     p7_spensemble_Add(P7_SPENSEMBLE *sp, int sampleidx, int i, int j, int k, int m);
int     p7_spensemble_Cluster(P7_SPENSEMBLE *sp,
				     float min_overlap, int of_smaller, int max_diagdiff,
				     float min_posterior, float min_endpointp,
				     int *ret_nclusters);
int     p7_spensemble_GetClusterCoords(P7_SPENSEMBLE *sp, int which,
					      int *ret_i, int *ret_j, int *ret_k, int *ret_m, float *ret_p);
void    p7_spensemble_Destroy(P7_SPENSEMBLE *sp);

/* p7_tophits.c */
P7_TOPHITS *p7_tophits_Create(void);
int         p7_tophits_Grow(P7_TOPHITS *h);
int         p7_tophits_CreateNextHit(P7_TOPHITS *h, P7_HIT **ret_hit);
int         p7_tophits_Add(P7_TOPHITS *h,
				  char *name, char *acc, char *desc,
				  double sortkey,
				  float score,    double lnP,
				  float mothersc, double mother_lnP,
				  int sqfrom, int sqto, int sqlen,
				  int hmmfrom, int hmmto, int hmmlen,
				  int domidx, int ndom,
				  P7_ALIDISPLAY *ali);
int         p7_tophits_SortBySortkey(P7_TOPHITS *h);
int         p7_tophits_SortBySeqidxAndAlipos(P7_TOPHITS *h);
int         p7_tophits_SortByModelnameAndAlipos(P7_TOPHITS *h);

int         p7_tophits_Merge(P7_TOPHITS *h1, P7_TOPHITS *h2);
int         p7_tophits_GetMaxPositionLength(P7_TOPHITS *h);
int         p7_tophits_GetMaxNameLength(P7_TOPHITS *h);
int         p7_tophits_GetMaxAccessionLength(P7_TOPHITS *h);
int         p7_tophits_GetMaxShownLength(P7_TOPHITS *h);
void        p7_tophits_Destroy(P7_TOPHITS *h);

int p7_tophits_ComputeNhmmerEvalues(P7_TOPHITS *th, double N, int W);
int p7_tophits_RemoveDuplicates(P7_TOPHITS *th, int using_bit_cutoffs);
int p7_tophits_Threshold(P7_TOPHITS *th, P7_PIPELINE *pli);
int p7_tophits_CompareRanking(P7_TOPHITS *th, ESL_KEYHASH *kh, int *opt_nnew);
int p7_tophits_Targets(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw);
int p7_tophits_Domains(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw);

int p7_tophits_Alignment(const P7_TOPHITS *th, const ESL_ALPHABET *abc,
				ESL_SQ **inc_sqarr, P7_TRACE **inc_trarr, int inc_n, int optflags,
				ESL_MSA **ret_msa);
int p7_tophits_TabularTargets(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header);
int p7_tophits_TabularDomains(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header);
int p7_tophits_TabularXfam(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli);
int p7_tophits_TabularTail(FILE *ofp, const char *progname, enum p7_pipemodes_e pipemode,
				  const char *qfile, const char *tfile, const ESL_GETOPTS *go);
int p7_tophits_AliScores(FILE *ofp, char *qname, P7_TOPHITS *th );

/* p7_trace.c */
P7_TRACE *p7_trace_Create(void);
P7_TRACE *p7_trace_CreateWithPP(void);
int  p7_trace_Reuse(P7_TRACE *tr);
int  p7_trace_Grow(P7_TRACE *tr);
int  p7_trace_GrowIndex(P7_TRACE *tr);
int  p7_trace_GrowTo(P7_TRACE *tr, int N);
int  p7_trace_GrowIndexTo(P7_TRACE *tr, int ndom);
void p7_trace_Destroy(P7_TRACE *tr);
void p7_trace_DestroyArray(P7_TRACE **tr, int N);

int  p7_trace_GetDomainCount   (const P7_TRACE *tr, int *ret_ndom);
int  p7_trace_GetStateUseCounts(const P7_TRACE *tr, int *counts);
int  p7_trace_GetDomainCoords  (const P7_TRACE *tr, int which, int *ret_i1, int *ret_i2,
				       int *ret_k1, int *ret_k2);

int   p7_trace_Validate(const P7_TRACE *tr, const ESL_ALPHABET *abc, const ESL_DSQ *dsq, char *errbuf);
int   p7_trace_Dump(FILE *fp, const P7_TRACE *tr, const P7_PROFILE *gm, const ESL_DSQ *dsq);
int   p7_trace_Compare(P7_TRACE *tr1, P7_TRACE *tr2, float pptol);
int   p7_trace_Score(P7_TRACE *tr, ESL_DSQ *dsq, P7_PROFILE *gm, float *ret_sc);
int   p7_trace_SetPP(P7_TRACE *tr, const P7_GMX *pp);
float p7_trace_GetExpectedAccuracy(const P7_TRACE *tr);

int  p7_trace_Append(P7_TRACE *tr, char st, int k, int i);
int  p7_trace_AppendWithPP(P7_TRACE *tr, char st, int k, int i, float pp);
int  p7_trace_Reverse(P7_TRACE *tr);
int  p7_trace_Index(P7_TRACE *tr);

int  p7_trace_FauxFromMSA(ESL_MSA *msa, int *matassign, int optflags, P7_TRACE **tr);
int  p7_trace_Doctor(P7_TRACE *tr, int *opt_ndi, int *opt_nid);

int  p7_trace_Count(P7_HMM *hmm, ESL_DSQ *dsq, float wt, P7_TRACE *tr);

/* seqmodel.c */
int p7_Seqmodel(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int M, char *name,
		       ESL_DMATRIX *P, float *f, double popen, double pextend,
		       P7_HMM **ret_hmm);

/* fm_alphabet.c */
int fm_alphabetCreate (FM_METADATA *meta, uint8_t *alph_bits);
int fm_alphabetDestroy (FM_METADATA *meta);
int fm_reverseString (char *str, int N);
int fm_getComplement (char c, uint8_t alph_type);

/* fm_general.c */
uint64_t fm_computeSequenceOffset (const FM_DATA *fms, const FM_METADATA *meta, int block, uint64_t pos);
int fm_getOriginalPosition (const FM_DATA *fms, const FM_METADATA *meta, int fm_id, int length, int direction, uint64_t fm_pos,
									uint32_t *segment_id, uint64_t *seg_pos);
int fm_readFMmeta( FM_METADATA *meta);
int fm_FM_read( FM_DATA *fm, FM_METADATA *meta, int getAll );
void fm_FM_destroy ( FM_DATA *fm, int isMainFM);
uint8_t fm_getChar(uint8_t alph_type, int j, const uint8_t *B );
int fm_getSARangeReverse( const FM_DATA *fm, FM_CFG *cfg, char *query, char *inv_alph, FM_INTERVAL *interval);
int fm_getSARangeForward( const FM_DATA *fm, FM_CFG *cfg, char *query, char *inv_alph, FM_INTERVAL *interval);
int fm_configAlloc(FM_CFG **cfg);
int fm_configDestroy(FM_CFG *cfg);
int fm_metaDestroy(FM_METADATA *meta );
int fm_updateIntervalForward( const FM_DATA *fm, const FM_CFG *cfg, char c, FM_INTERVAL *interval_f, FM_INTERVAL *interval_bk);
int fm_updateIntervalReverse( const FM_DATA *fm, const FM_CFG *cfg, char c, FM_INTERVAL *interval);
int fm_initSeeds (FM_DIAGLIST *list) ;
FM_DIAG * fm_newSeed (FM_DIAGLIST *list);
int fm_initAmbiguityList (FM_AMBIGLIST *list);
int fm_addAmbiguityRange (FM_AMBIGLIST *list, uint32_t start, uint32_t stop);
int fm_convertRange2DSQ(const FM_DATA *fm, const FM_METADATA *meta, uint64_t first, int length, int complementarity, ESL_SQ *sq, int fix_ambiguities );
int fm_initConfigGeneric( FM_CFG *cfg, ESL_GETOPTS *go);

/* fm_ssv.c */
int p7_SSVFM_longlarget( P7_OPROFILE *om, float nu, P7_BG *bg, double F1,
					  const FM_DATA *fmf, const FM_DATA *fmb, FM_CFG *fm_cfg, const P7_SCOREDATA *ssvdata,
					  int strands, P7_HMM_WINDOWLIST *windowlist);

/* fm_sse.c */
int fm_configInit      (FM_CFG *cfg, ESL_GETOPTS *go);
int fm_getOccCount     (const FM_DATA *fm, const FM_CFG *cfg, int pos, uint8_t c);
int fm_getOccCountLT   (const FM_DATA *fm, const FM_CFG *cfg, int pos, uint8_t c, uint32_t *cnteq, uint32_t *cntlt);

#endif /*P7_HMMERH_INCLUDED*/

/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Version 3.1b3; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 *
 * HMMER is distributed under the terms of the GNU General Public License
 * (GPLv3). See the LICENSE file for details.
 *
 * SVN $Id$
 * SVN $URL$
 ************************************************************/

/*** End of inlined file: hmmer.h ***/

typedef struct {
  char               *name;        /* name of the hmm database              */
  ESL_ALPHABET       *abc;         /* alphabet for database                 */

  P7_OPROFILE       **list;        /* list of profiles [0 .. n-1]           */
  uint32_t            lalloc;	   /* allocated length of <list>            */
  uint32_t            n;           /* number of entries in <list>           */
} P7_HMMCACHE;

int    p7_hmmcache_Open (char *hmmfile, P7_HMMCACHE **ret_cache, char *errbuf);
size_t p7_hmmcache_Sizeof         (P7_HMMCACHE *cache);
int    p7_hmmcache_SetNumericNames(P7_HMMCACHE *cache);
void   p7_hmmcache_Close          (P7_HMMCACHE *cache);

#endif /*P7_HMMCACHE_INCLUDED*/

/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Version 3.1b3; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 *
 * HMMER is distributed under the terms of the GNU General Public License
 * (GPLv3). See the LICENSE file for details.
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/*** End of inlined file: p7_hmmcache.h ***/

/*
 * SVN $URL$
 * SVN $Id$
 */
#ifndef P7_HMMPGMD_INCLUDED
#define P7_HMMPGMD_INCLUDED


typedef struct {
  uint32_t   status;            /* error status                             */
  uint64_t   msg_size;          /* size of the next packet.  if status not  */
                                /* zero, the length is for the error string */
                                /* otherwise it is the length of the data   */
                                /* to follow.                               */
} HMMD_SEARCH_STATUS;

typedef struct {
  double     elapsed;         	/* elapsed time, seconds                    */
  double     user;            	/* CPU time, seconds                        */
  double     sys;             	/* system time, seconds                     */

  double  Z;			/* eff # targs searched (per-target E-val)  */
  double  domZ;			/* eff # signific targs (per-domain E-val)  */
  enum p7_zsetby_e Z_setby;   	/* how Z was set                            */
  enum p7_zsetby_e domZ_setby;	/* how domZ was set                         */

  uint64_t   nmodels;         	/* # of HMMs searched                       */
  uint64_t   nseqs;           	/* # of sequences searched                  */
  uint64_t   n_past_msv;      	/* # comparisons that pass MSVFilter()      */
  uint64_t   n_past_bias;     	/* # comparisons that pass bias filter      */
  uint64_t   n_past_vit;      	/* # comparisons that pass ViterbiFilter()  */
  uint64_t   n_past_fwd;      	/* # comparisons that pass ForwardFilter()  */

  uint64_t   nhits;           	/* number of hits in list now               */
  uint64_t   nreported;       	/* number of hits that are reportable       */
  uint64_t   nincluded;       	/* number of hits that are includable       */
} HMMD_SEARCH_STATS;

#define HMMD_SEQUENCE   101
#define HMMD_HMM        102

/* commands between master and worker */
#define HMMD_CMD_SEARCH     10001
#define HMMD_CMD_SCAN       10002
#define HMMD_CMD_INIT       10003
#define HMMD_CMD_SHUTDOWN   10004
#define HMMD_CMD_RESET      10005

#define MAX_INIT_DESC 32

/* HMMD_CMD_SEARCH or HMMD_CMD_SCAN */
typedef struct {
  uint32_t    db_inx;               /* database index to search                 */
  uint32_t    db_type;              /* database type to search                  */
  uint32_t    inx;                  /* index to begin search                    */
  uint32_t    cnt;                  /* number of sequences to search            */
  uint32_t    query_type;           /* sequence / hmm                           */
  uint32_t    query_length;         /* length of the query data                 */
  uint32_t    opts_length;          /* length of the options string             */
  char        data[1];              /* search data                              */
} HMMD_SEARCH_CMD;

/* HMMD_CMD_INIT */
typedef struct {
  char        sid[MAX_INIT_DESC];   /* unique id for sequence database          */
  char        hid[MAX_INIT_DESC];   /* unique id for hmm database               */
  uint32_t    seqdb_off;            /* offset to seq database name, 0 if none   */
  uint32_t    hmmdb_off;            /* offset to hmm database name, 0 if none   */
  uint32_t    db_cnt;               /* total number of sequence databases       */
  uint32_t    seq_cnt;              /* sequences in database                    */
  uint32_t    hmm_cnt;              /* total number hmm databases               */
  uint32_t    model_cnt;            /* models in hmm database                   */
  char        data[1];              /* string data                              */
} HMMD_INIT_CMD;

/* HMMD_CMD_RESET */
typedef struct {
  char        ip_addr[1];           /* ip address                               */
} HMMD_INIT_RESET;

/* HMMD_HEADER */
typedef struct {
  uint32_t   length;                /* message length                           */
  uint32_t   command;               /* message type                             */
  uint32_t    status;               /* 0 - success                              */
} HMMD_HEADER;

typedef struct {
  HMMD_HEADER hdr;                  /* length and type of message               */
  union {
    HMMD_INIT_CMD   init;
    HMMD_SEARCH_CMD srch;
    HMMD_INIT_RESET reset;
  };
} HMMD_COMMAND;

#define MSG_SIZE(x) (sizeof(HMMD_HEADER) + ((HMMD_HEADER *)(x))->length)

size_t writen(int fd, const void *vptr, size_t n);
size_t readn(int fd, void *vptr, size_t n);

typedef struct queue_data_s {
  uint32_t       cmd_type;    /* type of command to perform     */
  uint32_t       query_type;  /* type of the query              */
  P7_HMM        *hmm;         /* query HMM                      */
  ESL_SQ        *seq;         /* query sequence                 */
  ESL_ALPHABET  *abc;         /* digital alphabet               */
  ESL_GETOPTS   *opts;        /* search specific options        */
  HMMD_COMMAND  *cmd;         /* workers search command         */

  int            sock;        /* socket descriptor of client    */
  char           ip_addr[64];

  int            dbx;         /* database index to search       */
  int            inx;         /* sequence index to start search */
  int            cnt;         /* number of sequences to search  */

} QUEUE_DATA;


typedef struct {
  int       N;       /* number of ranges */
  uint32_t *starts;  /* 0..N-1  start positions */
  uint32_t *ends;    /* 0..N-1  start positions */
} RANGE_LIST;

extern void free_QueueData(QUEUE_DATA *data);
extern int  hmmpgmd_IsWithinRanges (int64_t sq_idx, RANGE_LIST *list );
extern int  hmmpgmd_GetRanges (RANGE_LIST *list, char *rangestr);

extern int  process_searchopts(int fd, char *cmdstr, ESL_GETOPTS **ret_opts);

extern void worker_process(ESL_GETOPTS *go);
extern void master_process(ESL_GETOPTS *go);

#define LOG_FATAL_MSG(str, err) {                                               \
    p7_syslog(LOG_CRIT,"[%s:%d] - %s error %d - %s\n", __FILE__, __LINE__, str, err, strerror(err)); \
    exit(0); \
  }

#endif /*P7_HMMPGMD_INCLUDED*/

/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Version 3.1b3; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * HMMER is distributed under the terms of the GNU General Public License
 * (GPLv3). See the LICENSE file for details.
 ************************************************************/

#endif /*P7_HMMERH_INCLUDED*/

/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Version 3.1b3; July 2016
 * Copyright (C) 2016 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 *
 * HMMER is distributed under the terms of the GNU General Public License
 * (GPLv3). See the LICENSE file for details.
 *
 * SVN $Id$
 * SVN $URL$
 ************************************************************/


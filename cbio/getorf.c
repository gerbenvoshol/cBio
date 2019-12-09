/* @source getorf application
**
** Finds and extracts open reading frames (ORFs)
**
** @author Copyright (C) Gary Williams (gwilliam@hgmp.mrc.ac.uk)
** @author Copyright (C) Gerben Voshol (gpvoshol@gmail.com)
** @@
**
** This program is free software; you can redistribute it and/or
** modify it under the terms of the GNU General Public License
** as published by the Free Software Foundation; either version 2
** of the License, or (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
******************************************************************************/

#include "cbio.h"
#include <math.h> // For exp()



static void getorf_WriteORF(const AjPSeq seq, ajint len, ajint seqlen,
                            AjBool sense, ajint find, ajint *orf_no,
                            ajint start, ajint pos, const AjPStr str,
                            AjPSeqout seqout, ajint around, AjBool stop);

static void getorf_AppORF(ajint find, AjPStr *str,
                          const char *chrseq, ajint pos,
                          char aa);

static void getorf_FindORFs(const AjPSeq seq, ajint len, const AjPTrn trnTable,
                            ajuint minsize, ajuint maxsize, AjPSeqout seqout,
                            AjBool sense, AjBool circular, ajint find,
                            ajint *orf_no, AjBool methionine, AjBool stop, ajint around);




/* types of control codon */
#define START 1
#define STOP -1

/* types of ORF to find */
#define P_STOP2STOP      0
#define P_START2STOP     1
#define N_STOP2STOP      2
#define N_START2STOP     3
#define AROUND_START     4
#define AROUND_INIT_STOP 5
#define AROUND_END_STOP  6




/* @prog getorf ***************************************************************
**
** Finds and extracts open reading frames (ORFs)
**
******************************************************************************/

int main(int argc, char **argv)
{

	AjPSeqall seqall;
	AjPSeqout seqout;
	AjPStr tablestr;
	ajint table;
	ajuint minsize;
	ajuint maxsize;
	AjPStr findstr;
	ajint find;
	AjBool methionine;
	AjBool stop;
	AjBool circular;
	AjBool reverse;
	ajint around;

	AjPSeq seq = NULL;
	AjPTrn trnTable;
	AjPStr sseq = NULL;	/* sequence string */

	/* ORF number to append to name of sequence to create unique name */
	ajint orf_no;


	AjBool sense;	/* ajTrue = forward sense */
	ajint len;

	embInit("getorf", argc, argv);

	seqout     = ajAcdGetSeqoutall("outseq");
	seqall     = ajAcdGetSeqall("sequence");
	tablestr   = ajAcdGetListSingle("table");
	minsize    = ajAcdGetInt("minsize");
	maxsize    = ajAcdGetInt("maxsize");
	findstr    = ajAcdGetListSingle("find");
	methionine = ajAcdGetBoolean("methionine");
	stop       = ajAcdGetBoolean("stop");
	circular   = ajAcdGetBoolean("circular");
	reverse    = ajAcdGetBoolean("reverse");
	around     = ajAcdGetInt("flanking");


	/* initialise the translation table */
	ajStrToInt(tablestr, &table);
	trnTable = ajTrnNewI(table);

	/* what sort of ORF are we looking for */
	ajStrToInt(findstr, &find);

	/*
	** get the minimum size converted to protein length if storing
	** protein sequences
	*/
	if (find == P_STOP2STOP || find == P_START2STOP || find == AROUND_START) {
		minsize /= 3;
		maxsize /= 3;
	}

	while (ajSeqallNext(seqall, &seq)) {
		ajSeqTrim(seq);

		orf_no = 1;		   /* number of the next ORF */
		sense = ajTrue;		   /* forward sense initially */

		/* get the length of the sequence */
		len = ajSeqGetLen(seq);

		/*
		** if the sequence is circular, append it to itself to triple its
		**  length so can deal easily with wrapped ORFs, but don't update
		** len
		*/
		if (circular) {
			ajStrAssignS(&sseq, ajSeqGetSeqS(seq));
			ajStrAppendS(&sseq, ajSeqGetSeqS(seq));
			ajStrAppendS(&sseq, ajSeqGetSeqS(seq));
			ajSeqAssignSeqS(seq, sseq);
		}

		/* find the ORFs */
		getorf_FindORFs(seq, len, trnTable, minsize, maxsize, seqout, sense,
		                circular, find, &orf_no, methionine, stop, around);

		/* now reverse complement the sequence and do it again */
		if (reverse) {
			sense = ajFalse;
			ajSeqReverseForce(seq);
			getorf_FindORFs(seq, len, trnTable, minsize, maxsize, seqout, sense,
			                circular, find, &orf_no, methionine, stop, 
			                around);
		}
	}

	ajSeqoutClose(seqout);
	ajTrnDel(&trnTable);

	ajSeqallDel(&seqall);
	ajSeqDel(&seq);
	ajStrDel(&sseq);
	ajSeqoutDel(&seqout);
	ajStrDel(&tablestr);
	ajStrDel(&findstr);

	embExit();

	return 0;
}




/* @funcstatic getorf_FindORFs ************************************************
**
** finds all orfs in the current sense and writes them out
**
** @param [r] seq [const AjPSeq] Nucleotide sequence
** @param [r] len [ajint] Sequence length
** @param [r] trnTable [const AjPTrn] Translation table
** @param [r] minsize [ajuint] Minimum size ORF to find
** @param [r] maxsize [ajuint] Maximum size ORF to find
** @param [u] seqout [AjPSeqout] Sequence output object
** @param [r] sense [AjBool] Forward (sense) strand if true, else reverse strand
** @param [r] circular [AjBool] True if sequence is circular
** @param [r] find [ajint] Find code (see main program)
** @param [w] orf_no [ajint*] ORF counter updated
** @param [r] methionine [AjBool] If true report all start codons as 'M'
** @param [r] stop [AjBool] If true report stop codons
** @param [r] around [ajint] Number of bases around start/stop codons
** @@
******************************************************************************/

static void getorf_FindORFs(const AjPSeq seq, ajint len, const AjPTrn trnTable,
                            ajuint minsize, ajuint maxsize, AjPSeqout seqout,
                            AjBool sense, AjBool circular, ajint find,
                            ajint *orf_no, AjBool methionine, AjBool stop, ajint around)
{
	AjBool ORF[3];			/* true if found an ORF */
	AjBool LASTORF[3];		 /* true if hit the end of an ORF past
				    the end on the genome in this
				    frame */
	AjBool GOTSTOP[3];		 /* true if found a STOP in a circular
				    genome's frame when
				    find = P_STOP2STOP or
				    N_STOP2STOP */
	ajint start[3];		  /* possible starting position of the
				     three frames */
	ajint pos;
	ajint codon;
	char aa;
	ajint frame;
	AjPStr newstr[3];		 /* strings of the three frames of ORF
				    sequences that we are growing */
	AjPSeq pep = NULL;
	ajint i;

	ajint seqlen;
	const char *chrseq;

	seqlen = ajSeqGetLen(seq);
	chrseq = ajSeqGetSeqC(seq);

	/* initialise the ORF sequences */
	newstr[0] = NULL;
	newstr[1] = NULL;
	newstr[2] = NULL;

	/*
	** initialise flags for found the last ORF past the end of a circular
	** genome
	*/
	LASTORF[0] = ajFalse;
	LASTORF[1] = ajFalse;
	LASTORF[2] = ajFalse;

	/* initialise flags for found at least one STOP codon in a frame */
	GOTSTOP[0] = ajFalse;
	GOTSTOP[1] = ajFalse;
	GOTSTOP[2] = ajFalse;

	if (circular || find == P_START2STOP || find == N_START2STOP ||
	        find == AROUND_START) {
		ORF[0] = ajFalse;
		ORF[1] = ajFalse;
		ORF[2] = ajFalse;
	} else {
		/*
		** assume already in a ORF so we get ORFs at the start of the
		** sequence
		*/
		ORF[0] = ajTrue;
		ORF[1] = ajTrue;
		ORF[2] = ajTrue;
		start[0] = 0;
		start[1] = 1;
		start[2] = 2;
	}

	for (pos = 0; pos < seqlen - 2; pos++) {
		codon = ajTrnCodonstrTypeC(trnTable, &chrseq[pos], &aa);
		frame = pos % 3;
		ajDebug("len=%d, Pos=%d, Frame=%d start/stop=%d, aa=%c '%c%c%c'\n",
		        len, pos, frame, codon, aa,
		        chrseq[pos], chrseq[pos + 1], chrseq[pos + 2]);

		/* don't want to find extra ORFs when already been round circ */
		if (LASTORF[frame]) {
			continue;
		}

		if (find == P_STOP2STOP || find == N_STOP2STOP ||
		        find == AROUND_INIT_STOP || find == AROUND_END_STOP) {
			/* look for stop codon to begin reporting ORF */
			/* note that there was at least one STOP in a circular genome */
			if (codon == STOP) {
				GOTSTOP[frame] = ajTrue;
			}

			/* write details if a STOP is hit or the end of the sequence */
			if (codon == STOP || pos >= seqlen - 5) {

				/*
				** End of the sequence? If so, append any
				** last codon to the sequence - otherwise, ignore the STOP
				** codon
				*/
				if (codon != STOP)
					getorf_AppORF(find, &newstr[frame], chrseq, pos,
					              aa);

				/* Already have a sequence to write out? */
				if (ORF[frame]) {
					if (ajStrGetLen(newstr[frame]) >= minsize &&
					        ajStrGetLen(newstr[frame]) <= maxsize) {
						/* create a new sequence */
						if (codon == STOP)
							getorf_WriteORF(seq, len, seqlen, sense,
							                find, orf_no, start[frame],
							                pos - 1, newstr[frame],
							                seqout, around, stop);
						else
							getorf_WriteORF(seq, len, seqlen, sense,
							                find, orf_no, start[frame],
							                pos + 2, newstr[frame],
							                seqout, around, stop);
					}

					ajStrSetClear(&newstr[frame]);
				}

				/*
				** if its a circular genome and the STOP codon hits past
				** the end of the genome in all frames, then break
				*/
				if (circular && pos >= len) {
					ORF[frame] = ajFalse; /* past the end of the genome */
					LASTORF[frame] = ajTrue; /* finished getting ORFs */
					if (LASTORF[0] && LASTORF[1] && LASTORF[2]) {
						break;
					}
				} else {
					/*
					** hit a STOP, therefore a potential ORF to write
					** out next time, even if the genome is circular
					*/
					ORF[frame]   = ajTrue;
					start[frame] = pos + 3; /* next start of the ORF */
				}

			} else if (ORF[frame])
				/* append sequence to newstr if in an ORF */
			{
				getorf_AppORF(find, &newstr[frame], chrseq, pos, aa);
			}
		} else { /* Look for start: P_START2STOP N_START2STOP AROUND_START */

			if (codon == START && !ORF[frame]) {
				/* not in a ORF already and found a START */
				if (pos < len) {
					/*
					**  reset the newstr to zero length to enable
					**  storing the ORF for this
					*/
					ajStrSetClear(&newstr[frame]);
					ORF[frame] = ajTrue; /* now in an ORF */
					start[frame] = pos;	/* start of the ORF for this frame */
					if (methionine)
						getorf_AppORF(find, &newstr[frame], chrseq,
						              pos, 'M');
					else
						getorf_AppORF(find, &newstr[frame], chrseq,
						              pos, aa);
				}
			} else if (codon == STOP) {
				/* hit a STOP */

				/* Already have a sequence to write out? */
				if (ORF[frame]) {
					ORF[frame] = ajFalse; /* not in an ORF */

					if (ajStrGetLen(newstr[frame]) >= minsize &&
					        ajStrGetLen(newstr[frame]) <= maxsize) {
						/* create a new sequence */
						getorf_WriteORF(seq, len, seqlen, sense,
						                find, orf_no, start[frame],
						                pos - 1, newstr[frame],
						                seqout, around, stop);
					}
				}

				/*
				** if a circular genome and hit the STOP past
				** the end of the genome in all frames, then break
				*/
				if (circular && pos >= len) {
					LASTORF[frame] = ajTrue; /* finished getting ORFs */
					if (LASTORF[0] && LASTORF[1] && LASTORF[2]) {
						break;
					}
				}

				ajStrSetClear(&newstr[frame]);
			} else if (pos >= seqlen - 5) {
				/* hit the end of the sequence  without a stop */

				/* Already have a sequence to write out? */
				if (ORF[frame]) {
					ORF[frame] = ajFalse; /* not in an ORF */

					/*
					** End of the sequence? If so, append any
					** last codon to the sequence - otherwise, ignore the
					** STOP codon
					*/
					if (pos >= seqlen - 5 && pos < seqlen - 2)
						getorf_AppORF(find, &newstr[frame], chrseq,
						              pos, aa);

					if (ajStrGetLen(newstr[frame]) >= minsize &&
					        ajStrGetLen(newstr[frame]) <= maxsize) {
						/* create a new sequence */
						getorf_WriteORF(seq, len, seqlen, sense,
						                find, orf_no, start[frame],
						                pos + 2, newstr[frame],
						                seqout, around, stop);
					}
				}

				/*
				** if a circular genome and hit the STOP past
				** the end of the genome in all frames, then break
				*/
				if (circular && pos >= len) {
					LASTORF[frame] = ajTrue; /* finished getting ORFs */
					if (LASTORF[0] && LASTORF[1] && LASTORF[2]) {
						break;
					}
				}

				ajStrSetClear(&newstr[frame]);
			} else if (ORF[frame])
				getorf_AppORF(find, &newstr[frame], chrseq, pos,
				              aa);

		}
	}

	/*
	** Currently miss reporting a STOP-to-STOP ORF that is
	** the full length of a circular genome when there are no STOP codons in
	** that frame
	*/
	if ((find == P_STOP2STOP || find == N_STOP2STOP) && circular) {
		if (!GOTSTOP[0]) {
			/* translate frame 1 into pep */
			pep = ajTrnSeqOrig(trnTable, seq, 1);
			if (ajSeqGetLen(pep) >= minsize &&
			        ajSeqGetLen(pep) <= maxsize)
				getorf_WriteORF(seq, len, seqlen, sense, find, orf_no,
				                0, seqlen - 1, ajSeqGetSeqS(pep), seqout,
				                around, stop);
			ajSeqDel(&pep);
		}

		if (!GOTSTOP[1]) {
			/* translate frame 2 into pep */
			pep = ajTrnSeqOrig(trnTable, seq, 2);
			if (ajSeqGetLen(pep) >= minsize &&
			        ajSeqGetLen(pep) <= maxsize)
				getorf_WriteORF(seq, len, seqlen, sense, find, orf_no,
				                1, seqlen - 1, ajSeqGetSeqS(pep), seqout,
				                around, stop);
			ajSeqDel(&pep);
		}

		if (!GOTSTOP[2]) {
			/* translate frame 3 into pep */
			pep = ajTrnSeqOrig(trnTable, seq, 3);
			if (ajSeqGetLen(pep) >= minsize &&
			        ajSeqGetLen(pep) >= maxsize)
				getorf_WriteORF(seq, len, seqlen, sense, find, orf_no,
				                2, seqlen - 1, ajSeqGetSeqS(pep), seqout,
				                around, stop);
			ajSeqDel(&pep);
		}
	}

	for (i = 0; i < 3; ++i) {
		ajStrDel(&newstr[i]);
	}

	return;
}

//ajFmtStr

double sORFscore(AjPStr sequence)
{
    int hash;        // Index of feature (hash value)
    int alph_sz = 4; // Size of the alphabet (=A, T, C, G (4))
	int kmer_sz = 4; // number of nucleotides in a k-mer
	int step = 3;    // Step size

	// Allocate space to hold the feature vector plus 1 for unknown base pairs (N, Y, S, etc.)
	int kmer_cnt_sz = pow(alph_sz, kmer_sz);
	double *kmer_cnt = calloc(kmer_cnt_sz + 1, sizeof(double));

	// self.letters = ['A', 'T', 'C', 'G']
    int i, n;

   	int total_cnt = 0;
   	int limit = sequence->Len - kmer_sz;
    for (i = 0; i < limit; i += step) {
        /* calculate the lookup value for this k-mer */
        hash = 0;
        for (n = 0; n < kmer_sz; n++) {
            hash *= alph_sz;
            switch (sequence->Ptr[i + n]) {
                case 'A': 
                case 'a':
                	//hash += 0; NoOp
                    break;
                case 'T': 
                case 't': 
                	hash += 1;
                	break;
                case 'C': 
                case 'c': 
                    hash += 2; 
                    break;
                case 'G': 
                case 'g': 
                    hash += 3; 
                    break;
                default: 
                    hash = kmer_cnt_sz; 
                    break;
            }
            // encountered a unknown base pair (N, Y, S, etc.)
            if (hash == (kmer_cnt_sz))
                break;
        }
        kmer_cnt[hash]++;
        total_cnt++;
    }

    // To normalize for length divide the kmer count by the total kmers found
    for (i = 0; i < kmer_cnt_sz; i++) {
    	kmer_cnt[i] /= total_cnt;
    }

    
    double product = 0.0;

    // Logistic regression values from MiPepid
    // MiPepid: MicroPeptide identification tool using machine learning
    // Zhu and Gribskov BMC Bioinformatics (2019) 20:559
    // https://doi.org/10.1186/s12859-019-3033-9
	//double bias = 3.89456584;
    //double LGR[256] = {-5.57759579e+00, -1.05859130e+02, -4.04099724e+01, -4.73003909e+01, -1.15944023e+01, -1.35499324e+02, -6.96921599e+01, -5.05889625e+01, 1.06909132e+01, -1.13998582e+02, -4.36024991e+01, -6.08945003e+01, 3.32641905e+01, -9.52333770e+01, -4.25772483e+01, -4.01171853e+01, -8.74991254e+00, -1.28555437e+02, -6.63367327e+01, -4.85316525e+01, -7.27756568e+00, -1.31792063e+02, -5.95844301e+01, -4.74150789e+01, 5.83648035e+00, -1.19405397e+02, -4.75617470e+01, -2.09248189e+01, -6.45712491e+00, -1.17808779e+02, -6.12794174e+01, -4.39675011e+01, 1.01396323e+01, -1.14291203e+02, -5.12373807e+01, -5.56194754e+01, 1.00624046e+01, -1.08579545e+02, -5.75149113e+01, -4.85631936e+01, 2.00764723e+01, -1.23102599e+02, -3.97539501e+01, -3.15539007e+01, 1.18977732e+01, -9.15854740e+01, -4.00769470e+01, -3.79413458e+01, 3.75096913e+00, -1.06887052e+02, -4.98366159e+01, -4.20672905e+01, -2.23796586e+01, -1.07964664e+02, -6.21945424e+01, -6.48052671e+01, -1.12948130e+01, -9.65093727e+01, -5.58185035e+01, -3.37525318e+01, 5.98103136e+00, -9.96903732e+01, -6.33551089e+01, -7.55357475e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.15003833e+02, 2.19461057e-01, 8.03316537e+01, 6.25697377e+01, 1.43243972e+02, 1.14703945e+01, 8.63523395e+01, 8.25277112e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.18865768e+02, -7.75475006e+00, 5.41687205e+01, 4.44312565e+01, 1.12756381e+02, -5.02439599e+00, 7.28463421e+01, 6.38045711e+01, 1.29379677e+02, -4.58037434e+00, 5.53781537e+01, 9.46703031e+01, 1.11493114e+02, -3.51149548e+00, 5.87383416e+01, 6.28688072e+01, 1.22052654e+02, -1.00137602e+01, 5.84068672e+01, 5.51824638e+01, 1.11169944e+02, 1.65807172e+01, 3.87404617e+01, 5.94441954e+01, 1.31421290e+02, -1.27773744e+01, 4.91866494e+01, 6.43869315e+01, 1.21636222e+02, 3.93777797e+01, 6.92523571e+01, 4.05135279e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.06669167e+02, -1.84793880e+01, 3.47484935e+01, 4.72586934e+01, 1.12854641e+02, 8.59001107e+00, 5.48507321e+01, 3.63704528e+01, 1.14273590e+02, 1.14382809e+00, 5.96244323e+01, 4.20256539e+01, 5.36960710e+01, -6.10556143e+01, -3.08631437e+00, 3.79123320e+00, 3.50447529e+01, -6.50109528e+01, 9.75187382e+00, -1.45064740e+00, 5.28125614e+01, -6.21375842e+01, -5.16742564e+00, -6.09211764e-01, 7.05775912e+01, -5.24637963e+01, 6.37673108e+00, 8.11037559e+00, 4.47679038e+01, -5.81192484e+01, -1.33135841e+01, 6.95529030e+00, 4.95470920e+01, -6.48223723e+01, -8.70891870e+00, -1.90696891e+01, 6.08267031e+01, -5.98868339e+01, -1.07505246e+01, 5.42047226e-01, 5.28722573e+01, -3.37377870e+01, -2.15058102e+00, 5.62076932e+00, 7.42038831e+01, -7.30633343e+01, -1.03414191e+00, 4.03259620e-01, 6.25798033e+01, -6.20342963e+01, -2.48348527e+00, -1.28989162e+00, 4.67894663e+01, -5.85141036e+01, -7.96355984e+00, 8.79433715e+00, 5.69606600e+01, -6.63854369e+01, 3.27504433e+01, 9.62616425e-01, 5.62602889e+01, -6.52692461e+01, 8.62970241e+00, -8.18086603e-01, 7.54924311e+01, -5.47740270e+01, 2.00177003e+01, 9.08236915e+00, 6.93157465e+01, -5.19226615e+01, 4.08674764e+00, -1.44607527e+01, 8.05478099e+01, -6.51419322e+00, 4.40737781e+01, 3.15309602e+00, 6.00953853e+01, -5.66429293e+01, 1.35737043e+01, 9.15312420e+00, 6.59706127e+01, -6.62643170e+01, 1.84054952e+01, 4.27073472e-01, 5.33166663e+01, -4.97146256e+01, 2.56342582e+01, 1.90745293e+01, 6.49879744e+01, -5.90936211e+01, 1.92314208e+00, 1.59903275e+01, 6.15077598e+01, -5.70131969e+01, -9.08418837e+00, -2.28102018e+00, 5.35815571e+01, -7.11274886e+01, -3.37726880e+00, -4.75182346e+00, 6.47758366e+01, -6.90522469e+01, -5.16711990e-02, 3.48228733e+00, 4.91129306e+01, -6.29615636e+01, -1.69272582e+01, -8.06704690e+00, 4.44769390e+01, -7.25176651e+01, -1.38187706e+01, 7.99919964e+00, 5.80777355e+01, -6.39614145e+01, 2.70245140e+01, -5.80524646e+00, 5.12273533e+01, -5.99319562e+01, 2.08131354e+01, 3.44488487e+01, 7.76892172e+01, -7.27316393e+01, 1.77333630e+01, 3.55487592e+01, 5.80743102e+01, -5.34346189e+01, 1.43489816e+01, -2.44180727e-01, 4.13457278e+01, -4.81705250e+01, -5.30251955e+00, 8.11774427e+00, 4.92675520e+01, -5.87421900e+01, 2.33518040e+01, 4.45204265e+00, 4.75086432e+01, -5.84129208e+01, -8.79602240e+00, -3.60929319e+00};
    
    double bias = 0.0;
    // Used liblinear to retrain based on the original dataset the following:
    // L2-regularized logistic regression with a cost of 512
    // command used: ./train -s 0 -c 512 -v 10 mipepid.train model
    // Cross Validation Accuracy = 95.6278% 
    // F1 score = 0.978 versus 0.964 of the MiPepid paper based on a threshold of 0.6
    double LGR[256] = {0.44369950102579697, -38.640531070274335, -6.9905512052639462, -8.1980998153248201, 1.2067632474151635, -52.511824021803591, -17.158555386449834, -11.233992371339454, 9.7430886557140912, -42.959087268717866, -7.0196589793172262, -18.469656457152048, 25.283208467763796, -34.683833344665942, -6.166713736589692, -8.7797145892706574, -1.2428184105999529, -46.745790449655587, -33.723307994700761, -19.275995222650575, -2.4535021293875126, -49.233455459331388, -20.904067185763186, -11.283694812438748, 13.206700626509461, -45.559231918455211, -13.853626827592064, 0.59741923182191936, -16.097516618615021, -64.419415641175135, -42.061409083122875, -33.379880180665353, 7.3647757542136532, -40.166325434255242, -16.982014921700227, -14.600719923419442, 6.6063152682808717, -42.203622865896456, -19.127272346347642, -14.651789102107445, 13.307850344651522, -44.586705640412951, -9.6260808654394552, 2.3400864004443029, 14.111897671047149, -34.389550672434432, -10.02656824295164, -3.7618333553883159, 3.7125352961095661, -36.554481035912566, -4.7199136148494842, -13.722127637016218, 0.24667449612925618, -43.639419193129314, -17.811765333432884, -23.310937461158712, -6.4440570042804799, -35.779800841042317, -17.342161122297771, -0.44241175044745534, 7.7236511416724669, -37.348432943233995, -21.84988316175837, -27.591457652391149, 0, 0, 0, 0, 46.342182263952523, 3.4709161421067609, 29.754220546797697, 28.71186423544847, 60.179125631496426, 15.059449115226972, 42.558055513207151, 34.538354144173361, 0, 0, 0, 0, 52.406197562251712, 2.7490634019177089, 19.794112383877831, 20.260354058887003, 46.841984090052897, -6.5164768771746022, 31.519615173371953, 31.286075446700561, 54.878014514798153, 3.4496380568785385, 24.545863401770355, 45.334886696657371, 47.383173524854698, 0.53222203823961234, 28.200723384284181, 26.724165260193118, 58.005971130912442, -4.2595606753837068, 24.362905596169991, 24.820630529247616, 46.600411613133289, 7.3809723666323368, 12.470428787411313, 22.580806531402601, 56.166204472548358, -3.9543067455853582, 19.990098113146686, 28.279553629040166, 58.549571966323562, 16.247336508278597, 21.745664602161956, 17.567484033390681, 0, 0, 0, 0, 37.201815400204417, -9.6072565452460843, 10.699719097136253, 21.980042899494567, 46.803190193787657, 0.67893699295842524, 22.609221967460169, 9.2675211740251147, 50.720636832509221, 3.7433625487037414, 25.395922032162275, 16.651456557095582, 22.220123356219577, -23.130781624815938, 4.4060385172796499, 3.0606345058295394, 8.5270796357489171, -29.934346798901483, 0.52439258274557932, -6.3537825017816241, 15.841290589217525, -24.599892729528253, 1.0405342753893116, 0.9372998386992013, 35.135996514435817, -18.523546125929006, 10.096881670723178, 12.058382517102165, 13.383337922355997, -21.562421862012528, -4.3616295731256205, -0.9903205553275114, 20.219338468197019, -29.01510082308894, -0.20835069860145655, -9.5126412927381487, 24.222444637628602, -20.959369796907094, -0.49353570190420482, 7.0232155206592441, 22.969458076412248, -16.927028116070456, 7.7596815258729919, 10.874073106456041, 31.189718972381872, -27.570731584881909, 10.267400024451032, 0.87733971221255169, 21.021874328012107, -19.045454137831225, -1.837755060785309, -1.2038922776540957, 24.160537033373192, -22.981346382674189, -2.6671667614260195, 10.335428723256676, 28.354909741754984, -21.269375901922484, 29.85443063730003, 10.09394312637809, 28.050861384921703, -14.121513688062187, 16.83162148435968, 7.3849637596243323, 27.87254620101475, -18.93037210719476, 16.93324698778688, 7.496143636437643, 26.643196676017165, -10.905677942262299, -0.48414387489861288, -7.0339847626203467, 41.4285513872343, 6.1501211777189377, 22.111898686991616, 14.177373011728013, 28.376292856481808, -19.586781273822119, 10.599342699308995, 10.891308973310288, 27.644605757413846, -21.738798665235961, 14.015838542877752, 5.8312666873623691, 26.23079071046595, -18.194023386995287, 21.485604539572915, 15.313300978420806, 30.596706162375007, -19.645226790520443, 6.1606043257207403, 16.34473398620074, 27.402321740386977, -16.751200388780362, -8.0150254975549995, -3.160252475394437, 20.136690140335492, -27.442685504317243, 3.8236985471025751, -5.8403274808791261, 27.808219969210537, -29.8868403498894, 2.5521382594900492, 3.5672468109250515, 19.495310963038918, -23.070817802717013, -5.9382160305065064, -1.5342390429249779, 18.326951099085857, -30.867855485694868, -5.3679264119862973, 10.991749396393052, 18.600217008634296, -23.596716649356473, 14.420950685023399, 0.94670631445065956, 25.804209093279155, -23.375779441165545, 14.690922663033509, 22.707517461195142, 32.229845316103841, -17.502166971475525, 6.4977335921150736, 21.481499059656986, 29.128412915874605, -20.471075491956505, 15.341772500193128, 0.66507076672007115, 17.52944064276133, -19.52417738376295, 2.8108344094416258, 5.4622599970819641, 22.07538323716204, -15.781613670543432, 10.348241257277634, 11.898889015584077, 25.951335140132802, -25.618472728389754, -0.73698164819825318, 3.070321878350784};
    for (i = 0; i < kmer_cnt_sz; i++) {
    	product += LGR[i] * kmer_cnt[i];
    }
    double score = 1/(1+ exp(-product - bias));

    free(kmer_cnt);

    return score;
}


/* @funcstatic getorf_WriteORF ************************************************
**
** Undocumented.
**
** @param [r] seq [const AjPSeq] Undocumented
** @param [r] len [ajint] Undocumented
** @param [r] seqlen [ajint] Undocumented
** @param [r] sense [AjBool] Undocumented
** @param [r] find [ajint] Undocumented
** @param [w] orf_no [ajint*] Undocumented
** @param [r] start [ajint] Undocumented
** @param [r] pos [ajint] Undocumented
** @param [r] str [const AjPStr] Undocumented
** @param [u] seqout [AjPSeqout] Undocumented
** @param [r] around [ajint] Undocumented
** @param [r] stop [AjBool] If true report stop codons

** @@
******************************************************************************/

static void getorf_WriteORF(const AjPSeq seq,
                            ajint len, ajint seqlen, AjBool sense,
                            ajint find, ajint *orf_no, ajint start, ajint pos,
                            const AjPStr str, AjPSeqout seqout, ajint around, AjBool stop)
{
	AjPSeq new;
	AjPStr name  = NULL;       		/* name of the ORF */
	AjPStr value = NULL;       		/* string value of the ORF number */
	ajint s;
	ajint e;				/* start and end positions */
	AjPStr aroundstr = NULL;		/* holds sequence string around the
					   codon of interest */
	AjPStr stopstr = NULL;		/* holds sequence string of the STOP codon */
	ajint codonpos = 0;			/* holds position of start of codon
					   of interest */

	s = start + 1;
	e = pos + 1;

	/*
	** it is possible for an ORF in a circular genome to appear to start
	** past the end of the genome.
	** Move the reported positions back to start in the range 1..len
	** for readability.
	*/
	while (s > len) {
		s -= len;
		e -= len;
	}

	new = ajSeqNew();
	if (find == N_STOP2STOP || find == N_START2STOP ||
	        find == AROUND_INIT_STOP || find == AROUND_END_STOP ||
	        find == AROUND_START) {
		ajSeqSetNuc(new);
	} else {
		ajSeqSetProt(new);
	}


	/*
	** Set the start and end positions to report and get the sequence for
	** the AROUND* sequences
	*/
	if (find == AROUND_INIT_STOP) {
		codonpos = s - 3;
		s = codonpos - around;		/* 50 before the initial STOP */
		e = codonpos + around + 2;	/* 50 after the end of the STOP */
		if (s < 1) {
			return;
		}

		if (e > seqlen) {
			return;
		}
		ajStrAssignSubS(&aroundstr, ajSeqGetSeqS(seq), s - 1, e - 1);

	} else if (find == AROUND_START) {
		codonpos = s;
		s = codonpos - around;		/* 50 before the initial STOP */
		e = codonpos + around + 2;	/* 50 after the end of the STOP */
		if (s < 1) {
			return;
		}

		if (e > seqlen) {
			return;
		}
		ajStrAssignSubS(&aroundstr, ajSeqGetSeqS(seq), s - 1, e - 1);

	} else if (find == AROUND_END_STOP) {
		codonpos = e + 1;
		s = codonpos - around;		/* 50 before the initial STOP */
		e = codonpos + around + 2;	/* 50 after the end of the STOP */
		if (s < 1) {
			return;
		}

		if (e > seqlen) {
			return;
		}
		ajStrAssignSubS(&aroundstr, ajSeqGetSeqS(seq), s - 1, e - 1);
	}

 	if (stop) {
		s = s;		    /* start of the region */
		e = e + 3;	    /* end of the STOP */
		if (s < 1) {
			return;
		}

		if (e > seqlen) {
			return;
		}

		ajStrAssignSubS(&stopstr, ajSeqGetSeqS(seq), s - 1, e - 1);
	}

	/* set the name and description */
	ajStrAssignS(&name, ajSeqGetNameS(seq));
	ajStrAppendC(&name, "_");

	/* post-increment the ORF number for the next ORF */
	ajStrFromInt(&value, (*orf_no)++);

	ajStrAppendS(&name, value);
	ajSeqAssignNameS(new, name);

	/* set the description of the translation */
	ajStrAssignC(&name, "[");

	/* Reverse the reported positions if this is the reverse sense */
	if (!sense) {
		s = len - s + 1;
		e = len - e + 1;

		/*
		** shift the positions back into the range 1..len as far as possible
		       ** without going into negative numbers
		*/
		while (e > len) {
			s -= len;
			e -= len;
		}

		while (e < 0 || s < 0) {
			s += len;
			e += len;
		}
	}

	/* the base before the stop codon (numbering bases from 1) */
	ajStrFromInt(&value, s + ajSeqGetOffset(seq));

	ajStrAppendS(&name, value);
	ajStrAppendC(&name, " - ");

	/* the base before the stop codon (numbering bases from 1) */
	ajStrFromInt(&value, e + ajSeqGetOffset(seq));


	ajStrAppendS(&name, value);
	ajStrAppendC(&name, "] ");

	/* make it clear if this is the reverse sense */
	if (!sense) {
		ajStrAppendC(&name, "(REVERSE SENSE) ");
	}


	/*
	** make it clear if this is a circular genome and the ORF crosses
	** the breakpoint
	*/
	if (s > len || e > len) {
		ajStrAppendC(&name, "(ORF crosses the breakpoint) ");
	}


	if (find == AROUND_INIT_STOP || find == AROUND_START ||
	        find == AROUND_END_STOP) {
		ajStrAppendC(&name, "Around codon at ");
		ajStrFromInt(&value, codonpos);
		ajStrAppendS(&name, value);
		ajStrAppendC(&name, ". ");
	}

	ajStrAppendS(&name, ajSeqGetDescS(seq));
	AjPStr score;
	if (find == N_STOP2STOP || find == N_START2STOP || find == P_STOP2STOP || find == P_START2STOP) {
	    if (stop) {
	    	score = ajFmtStr(" sORFscore=%f", sORFscore(stopstr));
	    	ajStrAppendS(&name, score);
	    }
	}
	ajSeqAssignDescS(new, name);


	/* replace newstr in new */
	if (find == N_STOP2STOP || find == N_START2STOP || find == P_STOP2STOP || find == P_START2STOP) {
	    if (stop) {
			ajSeqAssignSeqS(new, stopstr);
		} else {
			ajSeqAssignSeqS(new, str);
		}
	} else {
		/* sequence to be 50 bases around the codon */
		ajSeqAssignSeqS(new, aroundstr);
	}

	ajSeqoutWriteSeq(seqout, new);

	ajSeqDel(&new);
	ajStrDel(&value);
	ajStrDel(&name);

	return;
}




/* @funcstatic getorf_AppORF **************************************************
**
** append aa to ORF sequence string
**
** @param [r] find [ajint] Find code
** @param [u] str [AjPStr*] Sequence string
** @param [r] chrseq [const char*] Undocumented
** @param [r] pos [ajint] Codon triplet position in chrseq
** @param [r] aa [char] Undocumented
** @@
******************************************************************************/

static void getorf_AppORF(ajint find, AjPStr *str,
                          const char *chrseq, ajint pos,
                          char aa)
{

	if (find == N_STOP2STOP || find == N_START2STOP ||
	        find == AROUND_INIT_STOP || find == AROUND_END_STOP) {
		ajStrAppendK(str, chrseq[pos]);
		ajStrAppendK(str, chrseq[pos + 1]);
		ajStrAppendK(str, chrseq[pos + 2]);
	} else if (find == P_STOP2STOP || find == P_START2STOP ||
	           find == AROUND_START) {
		ajStrAppendK(str, aa);
	}

	return;
}

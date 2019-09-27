/*
 * Copyright (C) 2014  Xcoo, Inc.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * This source implements public interfaces for using BWA functions from codes.
 * Some codes included in this source are based on the original BWA codes.
 */

#include "libbwa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bwa.h"
#include "bwase.h"
#include "utils.h"

extern void bwa_fprint_sam_hdr(FILE *stream, const bntseq_t *bns, const char *rg_line);

libbwa_samse_opt *libbwa_samse_opt_init(void)
{
    libbwa_samse_opt *o;
    o = calloc(1, sizeof(libbwa_samse_opt));
    o->n_occ = 3;
    o->rg_line = 0;
    return o;
}

void libbwa_samse_opt_destroy(libbwa_samse_opt *opt)
{
    free(opt);
}

// Copied from bwase.c
static int64_t pos_5(const bwa_seq_t *p)
{
    if (p->type != BWA_TYPE_NO_MATCH)
        return p->strand? pos_end(p) : p->pos;
    return -1;
}

// Based on bwa_print_sam1 in bwase.c
void libbwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2, FILE *out)
{
    int j;
    if (p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH)) {
        int seqid, nn, am = 0, flag = p->extra_flag;
        char XT;

        if (p->type == BWA_TYPE_NO_MATCH) {
            p->pos = mate->pos;
            p->strand = mate->strand;
            flag |= SAM_FSU;
            j = 1;
        } else j = pos_end(p) - p->pos; // j is the length of the reference in the alignment

        // get seqid
        nn = bns_cnt_ambi(bns, p->pos, j, &seqid);
        if (p->type != BWA_TYPE_NO_MATCH && p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len)
            flag |= SAM_FSU; // flag UNMAP as this alignment bridges two adjacent reference sequences

        // update flag and print it
        if (p->strand) flag |= SAM_FSR;
        if (mate) {
            if (mate->type != BWA_TYPE_NO_MATCH) {
                if (mate->strand) flag |= SAM_FMR;
            } else flag |= SAM_FMU;
        }
        err_fprintf(out, "%s\t%d\t%s\t", p->name, flag, bns->anns[seqid].name);
        err_fprintf(out, "%d\t%d\t", (int)(p->pos - bns->anns[seqid].offset + 1), p->mapQ);

        // print CIGAR
        if (p->cigar) {
            for (j = 0; j != p->n_cigar; ++j)
                err_fprintf(out, "%d%c", __cigar_len(p->cigar[j]), "MIDS"[__cigar_op(p->cigar[j])]);
        } else if (p->type == BWA_TYPE_NO_MATCH) err_fprintf(out, "*");
        else err_fprintf(out, "%dM", p->len);

        // print mate coordinate
        if (mate && mate->type != BWA_TYPE_NO_MATCH) {
            int m_seqid;
            long long isize;
            am = mate->seQ < p->seQ? mate->seQ : p->seQ; // smaller single-end mapping quality
            // redundant calculation here, but should not matter too much
            bns_cnt_ambi(bns, mate->pos, mate->len, &m_seqid);
            err_fprintf(out, "\t%s\t", (seqid == m_seqid)? "=" : bns->anns[m_seqid].name);
            isize = (seqid == m_seqid)? pos_5(mate) - pos_5(p) : 0;
            if (p->type == BWA_TYPE_NO_MATCH) isize = 0;
            err_fprintf(out, "%d\t%lld\t", (int)(mate->pos - bns->anns[m_seqid].offset + 1), isize);
        } else if (mate) err_fprintf(out, "\t=\t%d\t0\t", (int)(p->pos - bns->anns[seqid].offset + 1));
        else err_fprintf(out, "\t*\t0\t0\t");

        // print sequence and quality
        bwa_print_seq(out, p);
        err_putchar('\t');
        if (p->qual) {
            if (p->strand) seq_reverse(p->len, p->qual, 0); // reverse quality
            err_fprintf(out, "%s", p->qual);
        } else err_fprintf(out, "*");

        if (bwa_rg_id[0]) err_fprintf(out, "\tRG:Z:%s", bwa_rg_id);
        if (p->bc[0]) err_fprintf(out, "\tBC:Z:%s", p->bc);
        if (p->clip_len < p->full_len) err_fprintf(out, "\tXC:i:%d", p->clip_len);
        if (p->type != BWA_TYPE_NO_MATCH) {
            int i;
            // calculate XT tag
            XT = "NURM"[p->type];
            if (nn > 10) XT = 'N';
            // print tags
            err_fprintf(out, "\tXT:A:%c\t%s:i:%d", XT, (mode & BWA_MODE_COMPREAD)? "NM" : "CM", p->nm);
            if (nn) err_fprintf(out, "\tXN:i:%d", nn);
            if (mate) err_fprintf(out, "\tSM:i:%d\tAM:i:%d", p->seQ, am);
            if (p->type != BWA_TYPE_MATESW) { // X0 and X1 are not available for this type of alignment
                err_fprintf(out, "\tX0:i:%d", p->c1);
                if (p->c1 <= max_top2) err_fprintf(out, "\tX1:i:%d", p->c2);
            }
            err_fprintf(out, "\tXM:i:%d\tXO:i:%d\tXG:i:%d", p->n_mm, p->n_gapo, p->n_gapo+p->n_gape);
            if (p->md) err_fprintf(out, "\tMD:Z:%s", p->md);
            // print multiple hits
            if (p->n_multi) {
                err_fprintf(out, "\tXA:Z:");
                for (i = 0; i < p->n_multi; ++i) {
                    bwt_multi1_t *q = p->multi + i;
                    int k;
                    j = pos_end_multi(q, p->len) - q->pos;
                    nn = bns_cnt_ambi(bns, q->pos, j, &seqid);
                    err_fprintf(out, "%s,%c%d,", bns->anns[seqid].name, q->strand? '-' : '+',
                           (int)(q->pos - bns->anns[seqid].offset + 1));
                    if (q->cigar) {
                        for (k = 0; k < q->n_cigar; ++k)
                            err_fprintf(out, "%d%c", __cigar_len(q->cigar[k]), "MIDS"[__cigar_op(q->cigar[k])]);
                    } else err_fprintf(out, "%dM", p->len);
                    err_fprintf(out, ",%d;", q->gap + q->mm);
                }
            }
        }
        err_fputc('\n', out);
    } else { // this read has no match
        //ubyte_t *s = p->strand? p->rseq : p->seq;
        int flag = p->extra_flag | SAM_FSU;
        if (mate && mate->type == BWA_TYPE_NO_MATCH) flag |= SAM_FMU;
        err_fprintf(out, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", p->name, flag);
        //Why did this work differently to the version above??
        //for (j = 0; j != p->len; ++j) putchar("ACGTN"[(int)s[j]]);
        bwa_print_seq(out, p);
        err_fputc('\t', out);
        if (p->qual) {
            if (p->strand) seq_reverse(p->len, p->qual, 0); // reverse quality
            err_fprintf(out, "%s", p->qual);
        } else err_fprintf(out, "*");
        if (bwa_rg_id[0]) err_fprintf(out, "\tRG:Z:%s", bwa_rg_id);
        if (p->bc[0]) err_fprintf(out, "\tBC:Z:%s", p->bc);
        if (p->clip_len < p->full_len) err_fprintf(out, "\tXC:i:%d", p->clip_len);
        err_fputc('\n', out);
    }
}

// Based on bwa_sai2sam_se in bwase.c
void libbwa_sai2sam_se_core(const char *prefix, const char *fn_sa,
                            const char *fn_fa, int n_occ, const char *rg_line,
                            FILE *out)
{
    extern bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);
    int i, n_seqs, tot_seqs = 0, m_aln;
    bwt_aln1_t *aln = 0;
    bwa_seq_t *seqs;
    bwa_seqio_t *ks;
    bntseq_t *bns;
    FILE *fp_sa;
    gap_opt_t opt;
    char magic[4];

    // initialization
    bwase_initialize();
    bns = bns_restore(prefix);
    srand48(bns->seed);
    fp_sa = xopen(fn_sa, "r");

    m_aln = 0;
    err_fread_noeof(magic, 1, 4, fp_sa);
    if (strncmp(magic, SAI_MAGIC, 4) != 0) {
        fprintf(stderr, "[E::%s] Unmatched SAI magic. Please re-run `aln' with the same version of bwa.\n", __func__);
        exit(LIBBWA_E_UNMATCHED_SAI);
    }
    err_fread_noeof(&opt, sizeof(gap_opt_t), 1, fp_sa);
    bwa_fprint_sam_hdr(out, bns, rg_line);
    // set ks
    ks = bwa_open_reads(opt.mode, fn_fa);
    // core loop
    while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt.mode, opt.trim_qual)) != 0) {
        tot_seqs += n_seqs;

        // read alignment
        for (i = 0; i < n_seqs; ++i) {
            bwa_seq_t *p = seqs + i;
            int n_aln;
            err_fread_noeof(&n_aln, 4, 1, fp_sa);
            if (n_aln > m_aln) {
                m_aln = n_aln;
                aln = (bwt_aln1_t*)realloc(aln, sizeof(bwt_aln1_t) * m_aln);
            }
            err_fread_noeof(aln, sizeof(bwt_aln1_t), n_aln, fp_sa);
            bwa_aln2seq_core(n_aln, aln, p, 1, n_occ);
        }

        bwa_cal_pac_pos(bns, prefix, n_seqs, seqs, opt.max_diff, opt.fnr); // forward bwt will be destroyed here

        bwa_refine_gapped(bns, n_seqs, seqs, 0);

        for (i = 0; i < n_seqs; ++i)
            libbwa_print_sam1(bns, seqs + i, 0, opt.mode, opt.max_top2, out);

        bwa_free_read_seq(n_seqs, seqs);
    }

    // destroy
    bwa_seq_close(ks);
    bns_destroy(bns);
    err_fclose(fp_sa);
    free(aln);
}

// Based on bwa_sai2sam_se in bwase.c
int libbwa_samse(const char *db, const char *sai, const char *read, const char *out, const libbwa_samse_opt *opt)
{
    FILE *fpo;

    // Validate arguments
    if (!db || !sai || !read || !out || !opt) return LIBBWA_E_INVALID_ARGUMENT;

    char *prefix;
    if ((prefix = bwa_idx_infer_prefix(db)) == 0) {
        fprintf(stderr, "[%s] fail to locate the index\n", __func__);
        return LIBBWA_E_INDEX_ERROR;
    }
    fpo = xopen(out, "w");
    libbwa_sai2sam_se_core(prefix, sai, read, opt->n_occ, opt->rg_line, fpo);
    err_fclose(fpo);
    free(prefix);
    return LIBBWA_E_SUCCESS;
}

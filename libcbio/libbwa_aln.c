/*
 *  Copyright (C) 2014  Xcoo, Inc.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "libbwa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bwa.h"
#include "bwtaln.h"
#include "utils.h"

// Based on gap_init_opt in bwtaln.c
libbwa_aln_opt *libbwa_aln_opt_init(void)
{
    libbwa_aln_opt *o;
    o = (libbwa_aln_opt *)calloc(1, sizeof(libbwa_aln_opt));
    /* IMPORTANT: s_mm*10 should be about the average base error
       rate. Voilating this requirement will break pairing! */
    o->s_mm = 3; o->s_gapo = 11; o->s_gape = 4;
    o->max_diff = -1; o->max_gapo = 1; o->max_gape = 6;
    o->indel_end_skip = 5; o->max_del_occ = 10; o->max_entries = 2000000;
    o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
    o->seed_len = 32; o->max_seed_diff = 2;
    o->fnr = 0.04;
    o->n_threads = 1;
    o->max_top2 = 30;
    o->trim_qual = 0;
    return o;
}

void libbwa_aln_opt_destroy(libbwa_aln_opt *opt)
{
    free(opt);
}

void convert_aln_opt(const libbwa_aln_opt *src, gap_opt_t *dst)
{
    dst->s_mm = src->s_mm; dst->s_gapo = src->s_gapo; dst->s_gape = src->s_gape;
    dst->max_diff = src->max_diff; dst->max_gapo = src->max_gapo; dst->max_gape = src->max_gape;
    dst->indel_end_skip = src->indel_end_skip; dst->max_del_occ = src->max_del_occ; dst->max_entries = src->max_entries;
    dst->mode = src->mode;
    dst->seed_len = src->seed_len; dst->max_seed_diff = src->max_seed_diff;
    dst->fnr = src->fnr;
    dst->n_threads = src->n_threads;
    dst->max_top2 = src->max_top2;
    dst->trim_qual = src->trim_qual;
}

void libbwa_aln_core(const char *prefix, const char *fn_fa, const char *out ,const gap_opt_t *opt)
{
    int i, n_seqs, tot_seqs = 0;
    bwa_seq_t *seqs;
    bwa_seqio_t *ks;
    bwt_t *bwt;
    FILE *fpo;

    // initialization
    ks = bwa_open_reads(opt->mode, fn_fa);
    fpo = xopen(out, "w");

    { // load BWT
        char *str = (char*)calloc(strlen(prefix) + 10, 1);
        strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
        free(str);
    }

    // core loop
    err_fwrite(SAI_MAGIC, 1, 4, fpo);
    err_fwrite(opt, sizeof(gap_opt_t), 1, fpo);
    while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual)) != 0) {
        tot_seqs += n_seqs;

#ifdef HAVE_PTHREAD
        if (opt->n_threads <= 1) { // no multi-threading at all
            bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
        } else {
            pthread_t *tid;
            pthread_attr_t attr;
            thread_aux_t *data;
            int j;
            pthread_attr_init(&attr);
            pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
            tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
            for (j = 0; j < opt->n_threads; ++j) {
                data[j].tid = j; data[j].bwt = bwt;
                data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = opt;
                pthread_create(&tid[j], &attr, worker, data + j);
            }
            for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
            free(data); free(tid);
        }
#else
        bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
#endif

        for (i = 0; i < n_seqs; ++i) {
            bwa_seq_t *p = seqs + i;
            err_fwrite(&p->n_aln, 4, 1, fpo);
            if (p->n_aln) err_fwrite(p->aln, sizeof(bwt_aln1_t), p->n_aln, fpo);
        }

        bwa_free_read_seq(n_seqs, seqs);
    }

    // destroy
    bwt_destroy(bwt);
    bwa_seq_close(ks);
    err_fclose(fpo);
}

// Based on bwa_aln in bwtaln.c
int libbwa_aln(const char *db, const char *read, const char *out, const libbwa_aln_opt *opt_)
{
    gap_opt_t *opt;
    char *prefix;

    // Validate arguments
    if (!db || !read || !out || !opt_) return LIBBWA_E_INVALID_ARGUMENT;

    opt = gap_init_opt();
    convert_aln_opt(opt_, opt);

    if ((prefix = bwa_idx_infer_prefix(db)) == 0) {
        free(opt);
        return LIBBWA_E_INDEX_ERROR;
    }
    libbwa_aln_core(prefix, read, out, opt);
    free(opt); free(prefix);
    return LIBBWA_E_SUCCESS;
}

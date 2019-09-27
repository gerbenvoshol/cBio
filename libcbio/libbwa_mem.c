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
#include <math.h>
#include <time.h>

#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "kseq.h"
#include "utils.h"
KSEQ_DECLARE(gzFile)

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

extern void bwa_fprint_sam_hdr(FILE *stream, const bntseq_t *bns, const char *rg_line);

// Same as mem_opt_init in bwamem.c
libbwa_mem_opt *libbwa_mem_opt_init(void)
{
    libbwa_mem_opt *o;
    o = calloc(1, sizeof(libbwa_mem_opt));
    o->flag = 0;
    o->match_score = 1;
    o->mismatch_penalty = 4;
    o->o_del = o->o_ins = 6;
    o->e_del = o->e_ins = 1;
    o->band_width = 100;
    o->t = 30;
    o->zdrop = 100;
    o->pen_unpaired = 17;
    o->pen_clip5 = o->pen_clip3 = 5;
    o->min_seed_len = 19;
    o->split_width = 10;
    o->max_occ = 500;
    o->max_chain_gap = 10000;
    o->max_ins = 10000;
    o->mask_level = 0.50;
    o->drop_ratio = 0.50;
    o->xa_drop_ratio = 0.80;
    o->split_factor = 1.5;
    o->chunk_size = 10000000;
    o->n_threads = 1;
    o->max_hits = 5;
    o->max_matesw = 50;
    o->mask_level_redun = 0.95;
    o->min_chain_weight = 0;
    o->max_chain_extend = 1<<30;
    o->mapq_coef_len = 50;
    return o;
}

void libbwa_mem_opt_destroy(libbwa_mem_opt *opt)
{
    free(opt);
}

void convert_mem_opt(const libbwa_mem_opt *src, mem_opt_t *dst)
{
    dst->a = src->match_score;
    dst->b = src->mismatch_penalty;
    dst->o_del = src->o_del; dst->e_del = src->e_del;
    dst->o_ins = src->o_ins; dst->e_ins = src->e_ins;
    dst->pen_unpaired = src->pen_unpaired;
    dst->pen_clip5 = src->pen_clip5; dst->pen_clip3 = src->pen_clip3;
    dst->w = src->band_width;
    dst->zdrop = src->zdrop;
    dst->T = src->t;
    dst->flag = src->flag;
    dst->min_seed_len = src->min_seed_len;
    dst->min_chain_weight = src->min_chain_weight;
    dst->max_chain_extend = src->max_chain_extend;
    dst->split_factor = src->split_factor;
    dst->split_width = src->split_width;
    dst->max_occ = src->max_occ;
    dst->max_chain_gap = src->max_chain_gap;
    dst->n_threads = src->n_threads;
    dst->chunk_size = src->chunk_size;
    dst->mask_level = src->mask_level;
    dst->drop_ratio = src->drop_ratio;
    dst->XA_drop_ratio = src->xa_drop_ratio;
    dst->mask_level_redun = src->mask_level_redun;
    dst->mapQ_coef_len = src->mapq_coef_len;
    dst->mapQ_coef_fac = log(src->mapq_coef_len);
    dst->max_ins = src->max_ins;
    dst->max_matesw = src->max_matesw;
    dst->max_hits = src->max_hits;
}

// Modified based on main_mem in fastmap.c
int libbwa_mem(const char *db, const char *read, const char *mate, const char *out, const libbwa_mem_opt *opt_)
{
    mem_opt_t *opt;
    int fd, fd2, i, n, copy_comment = 0;
    gzFile fp, fp2 = 0;
    kseq_t *ks, *ks2 = 0;
    bseq1_t *seqs;
    bwaidx_t *idx;
    char *rg_line = 0;
    void *ko = 0, *ko2 = 0;
    int64_t n_processed = 0;
    mem_pestat_t *pes0 = 0;
    FILE *fpo;

    // Validate arguments
    if (!db || !read || !out || !opt_)
        return LIBBWA_E_INVALID_ARGUMENT;

    opt = mem_opt_init();
    convert_mem_opt(opt_, opt);

    bwa_fill_scmat(opt->a, opt->b, opt->mat);
    if ((idx = bwa_idx_load(db, BWA_IDX_ALL)) == 0) return LIBBWA_E_INDEX_ERROR; // FIXME: memory leak

    ko = kopen(read, &fd);
    if (ko == 0) {
        return LIBBWA_E_FILE_ERROR;
    }
    fp = gzdopen(fd, "r");
    ks = kseq_init(fp);
    if (mate != NULL) {
        if (opt->flag&MEM_F_PE) {
            if (bwa_verbose >= 2)
                fprintf(stderr, "[W::%s] when '-p' is in use, the second query file will be ignored.\n", __func__);
        } else {
            ko2 = kopen(mate, &fd2);
            if (ko2 == 0) {
                return LIBBWA_E_FILE_ERROR;
            }
            fp2 = gzdopen(fd2, "r");
            ks2 = kseq_init(fp2);
            opt->flag |= MEM_F_PE;
        }
    }

    fpo = xopen(out, "w");

    bwa_fprint_sam_hdr(fpo, idx->bns, rg_line);
    while ((seqs = bseq_read(opt->chunk_size * opt->n_threads, &n, ks, ks2)) != 0) {
        if ((opt->flag & MEM_F_PE) && (n&1) == 1) {
            if (bwa_verbose >= 2)
                fprintf(stderr, "[W::%s] odd number of reads in the PE mode; last read dropped\n", __func__);
            n = n>>1<<1;
        }
        if (!copy_comment)
            for (i = 0; i < n; ++i) {
                free(seqs[i].comment); seqs[i].comment = 0;
            }
        mem_process_seqs(opt, idx->bwt, idx->bns, idx->pac, n_processed, n, seqs, pes0);
        n_processed += n;
        for (i = 0; i < n; ++i) {
            err_fputs(seqs[i].sam, fpo);
            free(seqs[i].name); free(seqs[i].comment); free(seqs[i].seq); free(seqs[i].qual); free(seqs[i].sam);
        }
        free(seqs);
    }

    free(opt);
    bwa_idx_destroy(idx);
    kseq_destroy(ks);
    err_fclose(fpo);
    err_gzclose(fp); kclose(ko);
    if (ks2) {
        kseq_destroy(ks2);
        err_gzclose(fp2); kclose(ko2);
    }
    return LIBBWA_E_SUCCESS;
}

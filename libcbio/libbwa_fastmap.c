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

#include <stdlib.h>

#include "bwa.h"
#include "bwamem.h"
#include "kseq.h"
#include "utils.h"
KSEQ_DECLARE(gzFile)

libbwa_fastmap_opt *libbwa_fastmap_opt_init(void)
{
    libbwa_fastmap_opt *o;
    o = calloc(1, sizeof(libbwa_fastmap_opt));
    o->print_seq = 0;
    o->min_iwidth = 20;
    o->min_len = 17;
    return o;
}

void libbwa_fastmap_opt_destroy(libbwa_fastmap_opt *opt)
{
    free(opt);
}

// Based on main_fastmap in fastmap.c
int libbwa_fastmap(const char *db, const char *read, const char *out, const libbwa_fastmap_opt *opt)
{
    int i;
	kseq_t *seq;
	bwtint_t k;
	gzFile fp;
	smem_i *itr;
	const bwtintv_v *a;
	bwaidx_t *idx;
    FILE *fpo;

    // Validate arguments
    if (!db || !read || !out || !opt)
        return LIBBWA_E_INVALID_ARGUMENT;

	fp = xzopen(read, "r");
	seq = kseq_init(fp);
	if ((idx = bwa_idx_load(db, BWA_IDX_BWT|BWA_IDX_BNS)) == 0)
        return LIBBWA_E_INDEX_ERROR;

    fpo = xopen(out, "w");

	itr = smem_itr_init(idx->bwt);
	while (kseq_read(seq) >= 0) {
		err_fprintf(fpo, "SQ\t%s\t%ld", seq->name.s, seq->seq.l);
		if (opt->print_seq) {
			err_fputc('\t', fpo);
			err_fputs(seq->seq.s, fpo);
		} else err_fputc('\n', fpo);
		for (i = 0; i < seq->seq.l; ++i)
			seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
		smem_set_query(itr, seq->seq.l, (uint8_t*)seq->seq.s);
		while ((a = smem_next(itr)) != 0) {
			for (i = 0; i < a->n; ++i) {
				bwtintv_t *p = &a->a[i];
				if ((uint32_t)p->info - (p->info>>32) < opt->min_len) continue;
				err_fprintf(fpo, "EM\t%d\t%d\t%ld", (uint32_t)(p->info>>32), (uint32_t)p->info, (long)p->x[2]);
				if (p->x[2] <= opt->min_iwidth) {
					for (k = 0; k < p->x[2]; ++k) {
						bwtint_t pos;
						int len, is_rev, ref_id;
						len  = (uint32_t)p->info - (p->info>>32);
						pos = bns_depos(idx->bns, bwt_sa(idx->bwt, p->x[0] + k), &is_rev);
						if (is_rev) pos -= len - 1;
						bns_cnt_ambi(idx->bns, pos, len, &ref_id);
						err_fprintf(fpo, "\t%s:%c%ld", idx->bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - idx->bns->anns[ref_id].offset) + 1);
					}
				} else err_fputs("\t*", fpo);
				err_fputc('\n', fpo);
			}
		}
		err_fputs("//", fpo);
	}

	smem_itr_destroy(itr);
	bwa_idx_destroy(idx);
	kseq_destroy(seq);
    err_fclose(fpo);
	err_gzclose(fp);
	return LIBBWA_E_SUCCESS;
}

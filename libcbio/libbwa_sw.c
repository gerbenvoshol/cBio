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
#include <math.h>

#include "bwtsw2.h"
#include "bwa.h"
#include "utils.h"
#include "kstring.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

// Based on nt_comp_table in bwtsw2_aux.c
unsigned char libbwa_nt_comp_table[256] = {
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','T','V','G', 'H','N','N','C', 'D','N','N','M', 'N','K','N','N',
    'N','N','Y','S', 'A','N','B','W', 'X','R','N','N', 'N','N','N','N',
    'n','t','v','g', 'h','n','n','c', 'd','n','n','m', 'n','k','n','n',
    'n','n','y','s', 'a','n','b','w', 'x','r','n','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N'
};

extern int bsw2_resolve_duphits(const bntseq_t *bns, const bwt_t *bwt, bwtsw2_t *b, int IS);
extern int bsw2_resolve_query_overlaps(bwtsw2_t *b, float mask_level);

// Based on bsw2_init_opt in bwtsw2_aux.c
libbwa_sw_opt *libbwa_sw_opt_init(void)
{
    libbwa_sw_opt *o = (libbwa_sw_opt *)calloc(1, sizeof(libbwa_sw_opt));
    o->a = 1;
    o->b = 3;
    o->q = 5;
    o->r = 2;
    o->t = 30;

    o->bw = 50;

    o->max_ins = 20000;

    o->z = 1;
    o->is = 3;
    o->t_seeds = 5;
    o->hard_clip = 0;
    o->skip_sw = 0;

    o->mask_level = 0.50f;
    o->coef = 5.5f;

    o->qr = o->q + o->r;
    o->n_threads = 1;
    o->chunk_size = 10000000;

    o->max_chain_gap = 10000;

    o->cpy_cmt = 0;

    return o;
}

void libbwa_sw_opt_destroy(libbwa_sw_opt *opt)
{
    free(opt);
}

void _convert_sw_opt(const libbwa_sw_opt *src, bsw2opt_t *dst)
{
    dst->a = src->a;
    dst->b = src->b;
    dst->q = src->q;
    dst->r = src->r;
    dst->t = src->t;

    dst->bw = src->bw;

    dst->max_ins = src->max_ins;

    dst->z = src->z;
    dst->is = src->is;
    dst->t_seeds = src->t_seeds;
    dst->hard_clip = src->hard_clip;
    dst->skip_sw = src->skip_sw;

    dst->mask_level = src->mask_level;
    dst->coef = src->coef;

    dst->qr = src->qr;
    dst->n_threads = src->n_threads;
    dst->chunk_size = src->chunk_size;

    dst->max_chain_gap = src->max_chain_gap;

    dst->cpy_cmt = src->cpy_cmt;
}

// Copied from bwtsw2_aux.c
static void gen_cigar(const bsw2opt_t *opt, int lq, uint8_t *seq[2], int64_t l_pac, const uint8_t *pac, bwtsw2_t *b, const char *name)
{
    int i;
    int8_t mat[25];

    bwa_fill_scmat(opt->a, opt->b, mat);
    for (i = 0; i < b->n; ++i) {
        bsw2hit_t *p = b->hits + i;
        bsw2aux_t *q = b->aux + i;
        uint8_t *query;
        int beg, end, score;
        if (p->l) continue;
        beg = (p->flag & 0x10)? lq - p->end : p->beg;
        end = (p->flag & 0x10)? lq - p->beg : p->end;
        query = seq[(p->flag & 0x10)? 1 : 0] + beg;
        q->cigar = bwa_gen_cigar(mat, opt->q, opt->r, opt->bw, l_pac, pac, end - beg, query, p->k, p->k + p->len, &score, &q->n_cigar, &q->nm);
#if 0
        if (name && score != p->G) { // debugging only
            int j, glen = 0;
            for (j = 0; j < q->n_cigar; ++j)
                if ((q->cigar[j]&0xf) == 1 || (q->cigar[j]&0xf) == 2)
                    glen += q->cigar[j]>>4;
            fprintf(stderr, "[E::%s] %s - unequal score: %d != %d; (qlen, aqlen, arlen, glen, bw) = (%d, %d, %d, %d, %d)\n",
                    __func__, name, score, p->G, lq, end - beg, p->len, glen, opt->bw);
        }
#endif
        if (q->cigar && (beg != 0 || end < lq)) { // write soft clipping
            q->cigar = realloc(q->cigar, 4 * (q->n_cigar + 2));
            if (beg != 0) {
                memmove(q->cigar + 1, q->cigar, q->n_cigar * 4);
                q->cigar[0] = beg<<4 | 4;
                ++q->n_cigar;
            }
            if (end < lq) {
                q->cigar[q->n_cigar] = (lq - end)<<4 | 4;
                ++q->n_cigar;
            }
        }
    }
}

// Copied from bwtsw2_aux.c
static void merge_hits(bwtsw2_t *b[2], int l, int is_reverse)
{
    int i;
    if (b[0]->n + b[1]->n > b[0]->max) {
        b[0]->max = b[0]->n + b[1]->n;
        b[0]->hits = realloc(b[0]->hits, b[0]->max * sizeof(bsw2hit_t));
    }
    for (i = 0; i < b[1]->n; ++i) {
        bsw2hit_t *p = b[0]->hits + b[0]->n + i;
        *p = b[1]->hits[i];
        if (is_reverse) {
            int x = p->beg;
            p->beg = l - p->end;
            p->end = l - x;
            p->flag |= 0x10;
        }
    }
    b[0]->n += b[1]->n;
    bsw2_destroy(b[1]);
    b[1] = 0;
}

// Copied from bwtsw2_aux.c
static bwtsw2_t *bsw2_aln1_core(const bsw2opt_t *opt, const bntseq_t *bns, uint8_t *pac, const bwt_t *target,
                                int l, uint8_t *seq[2], bsw2global_t *pool)
{
    extern void bsw2_chain_filter(const bsw2opt_t *opt, int len, bwtsw2_t *b[2]);
    bwtsw2_t *b[2], **bb[2], **_b, *p;
    int k, j;
    bwtl_t *query;
    query = bwtl_seq2bwtl(l, seq[0]);
    _b = bsw2_core(bns, opt, query, target, pool);
    bwtl_destroy(query);
    for (k = 0; k < 2; ++k) {
        bb[k] = calloc(2, sizeof(void*));
        bb[k][0] = calloc(1, sizeof(bwtsw2_t));
        bb[k][1] = calloc(1, sizeof(bwtsw2_t));
    }
    for (k = 0; k < 2; ++k) { // separate _b into bb[2] based on the strand
        for (j = 0; j < _b[k]->n; ++j) {
            bsw2hit_t *q;
            p = bb[_b[k]->hits[j].is_rev][k];
            if (p->n == p->max) {
                p->max = p->max? p->max<<1 : 8;
                p->hits = realloc(p->hits, p->max * sizeof(bsw2hit_t));
            }
            q = &p->hits[p->n++];
            *q = _b[k]->hits[j];
            if (_b[k]->hits[j].is_rev) {
                int x = q->beg;
                q->beg = l - q->end;
                q->end = l - x;
            }
        }
    }
    b[0] = bb[0][1]; b[1] = bb[1][1]; // bb[*][1] are "narrow SA hits"
    bsw2_chain_filter(opt, l, b); // NB: only unique seeds are chained
    for (k = 0; k < 2; ++k) {
        bsw2_extend_left(opt, bb[k][1], seq[k], l, pac, bns->l_pac, pool->aln_mem);
        merge_hits(bb[k], l, 0); // bb[k][1] is merged to bb[k][0] here
        bsw2_resolve_duphits(0, 0, bb[k][0], 0);
        bsw2_extend_rght(opt, bb[k][0], seq[k], l, pac, bns->l_pac, pool->aln_mem);
        bsw2_resolve_duphits(0, 0, bb[k][0], 0);
        b[k] = bb[k][0];
        free(bb[k]);
    }
    merge_hits(b, l, 1); // again, b[1] is merged to b[0]
    bsw2_resolve_query_overlaps(b[0], opt->mask_level);
    bsw2_destroy(_b[0]); bsw2_destroy(_b[1]); free(_b);
    return b[0];
}

// Copied from bwtsw2_aux.c
static void flag_fr(bwtsw2_t *b[2])
{
    int i, j;
    for (i = 0; i < b[0]->n; ++i) {
        bsw2hit_t *p = b[0]->hits + i;
        p->flag |= 0x10000;
    }
    for (i = 0; i < b[1]->n; ++i) {
        bsw2hit_t *p = b[1]->hits + i;
        p->flag |= 0x20000;
    }
    for (i = 0; i < b[0]->n; ++i) {
        bsw2hit_t *p = b[0]->hits + i;
        for (j = 0; j < b[1]->n; ++j) {
            bsw2hit_t *q = b[1]->hits + j;
            if (q->beg == p->beg && q->end == p->end && q->k == p->k && q->len == p->len && q->G == p->G) {
                q->flag |= 0x30000; p->flag |= 0x30000;
                break;
            }
        }
    }
}

// Copied from bwtsw2_aux.c
static int fix_cigar(const bntseq_t *bns, bsw2hit_t *p, int n_cigar, uint32_t *cigar)
{
    // FIXME: this routine does not work if the query bridge three reference sequences
    int32_t coor, refl, lq;
    int x, y, i, seqid;
    bns_cnt_ambi(bns, p->k, p->len, &seqid);
    coor = p->k - bns->anns[seqid].offset;
    refl = bns->anns[seqid].len;
    x = coor; y = 0;
    // test if the alignment goes beyond the boundary
    for (i = 0; i < n_cigar; ++i) {
        int op = cigar[i]&0xf, ln = cigar[i]>>4;
        if (op == 1 || op == 4 || op == 5) y += ln;
        else if (op == 2) x += ln;
        else x += ln, y += ln;
    }
    lq = y; // length of the query sequence
    if (x > refl) { // then fix it
        int j, nc, mq[2], nlen[2];
        uint32_t *cn;
        bwtint_t kk = 0;
        nc = mq[0] = mq[1] = nlen[0] = nlen[1] = 0;
        cn = calloc(n_cigar + 3, 4);
        x = coor; y = 0;
        for (i = j = 0; i < n_cigar; ++i) {
            int op = cigar[i]&0xf, ln = cigar[i]>>4;
            if (op == 4 || op == 5 || op == 1) { // ins or clipping
                y += ln;
                cn[j++] = cigar[i];
            } else if (op == 2) { // del
                if (x + ln >= refl && nc == 0) {
                    cn[j++] = (uint32_t)(lq - y)<<4 | 4;
                    nc = j;
                    cn[j++] = (uint32_t)y<<4 | 4;
                    kk = p->k + (x + ln - refl);
                    nlen[0] = x - coor;
                    nlen[1] = p->len - nlen[0] - ln;
                } else cn[j++] = cigar[i];
                x += ln;
            } else if (op == 0) { // match
                if (x + ln >= refl && nc == 0) {
                    // FIXME: not consider a special case where a split right between M and I
                    cn[j++] = (uint32_t)(refl - x)<<4 | 0; // write M
                    cn[j++] = (uint32_t)(lq - y - (refl - x))<<4 | 4; // write S
                    nc = j;
                    mq[0] += refl - x;
                    cn[j++] = (uint32_t)(y + (refl - x))<<4 | 4;
                    if (x + ln - refl) cn[j++] = (uint32_t)(x + ln - refl)<<4 | 0;
                    mq[1] += x + ln - refl;
                    kk = bns->anns[seqid].offset + refl;
                    nlen[0] = refl - coor;
                    nlen[1] = p->len - nlen[0];
                } else {
                    cn[j++] = cigar[i];
                    mq[nc?1:0] += ln;
                }
                x += ln; y += ln;
            }
        }
        if (mq[0] > mq[1]) { // then take the first alignment
            n_cigar = nc;
            memcpy(cigar, cn, 4 * nc);
            p->len = nlen[0];
        } else {
            p->k = kk; p->len = nlen[1];
            n_cigar = j - nc;
            memcpy(cigar, cn + nc, 4 * (j - nc));
        }
        free(cn);
    }
    return n_cigar;
}

// Copied from bwtsw2_aux.c
static void write_aux(const bsw2opt_t *opt, const bntseq_t *bns, int qlen, uint8_t *seq[2], const uint8_t *pac, bwtsw2_t *b, const char *name)
{
    int i;
    // allocate for b->aux
    if (b->n<<1 < b->max) {
        b->max = b->n;
        kroundup32(b->max);
        b->hits = realloc(b->hits, b->max * sizeof(bsw2hit_t));
    }
    b->aux = calloc(b->n, sizeof(bsw2aux_t));
    // generate CIGAR
    gen_cigar(opt, qlen, seq, bns->l_pac, pac, b, name);
    // fix CIGAR, generate mapQ, and write chromosomal position
    for (i = 0; i < b->n; ++i) {
        bsw2hit_t *p = &b->hits[i];
        bsw2aux_t *q = &b->aux[i];
        q->flag = p->flag & 0xfe;
        q->isize = 0;
        if (p->l == 0) { // unique hit
            float c = 1.0;
            int subo;
            // fix out-of-boundary CIGAR
            q->n_cigar = fix_cigar(bns, p, q->n_cigar, q->cigar);
            // compute mapQ
            subo = p->G2 > opt->t? p->G2 : opt->t;
            if (p->flag>>16 == 1 || p->flag>>16 == 2) c *= .5;
            if (p->n_seeds < 2) c *= .2;
            q->qual = (int)(c * (p->G - subo) * (250.0 / p->G + 0.03 / opt->a) + .499);
            if (q->qual > 250) q->qual = 250;
            if (q->qual < 0) q->qual = 0;
            if (p->flag&1) q->qual = 0; // this is a random hit
            q->pqual = q->qual; // set the paired qual as qual
            // get the chromosomal position
            q->nn = bns_cnt_ambi(bns, p->k, p->len, &q->chr);
            q->pos = p->k - bns->anns[q->chr].offset;
        } else q->qual = 0, q->n_cigar = 0, q->chr = q->pos = -1, q->nn = 0;
    }
}

// Copied from bwtsw2_aux.c
static void update_mate_aux(bwtsw2_t *b, const bwtsw2_t *m)
{
    int i;
    if (m == 0) return;
    // update flag, mchr and mpos
    for (i = 0; i < b->n; ++i) {
        bsw2aux_t *q = &b->aux[i];
        q->flag |= 1; // paired
        if (m->n == 0) q->flag |= 8; // mate unmapped
        if (m->n == 1) {
            q->mchr = m->aux[0].chr;
            q->mpos = m->aux[0].pos;
            if (m->aux[0].flag&0x10) q->flag |= 0x20; // mate reverse strand
            if (q->chr == q->mchr) { // set insert size
                if (q->mpos + m->hits[0].len > q->pos)
                    q->isize = q->mpos + m->hits[0].len - q->pos;
                else q->isize = q->mpos - q->pos - b->hits[0].len;
            } else q->isize = 0;
        } else q->mchr = q->mpos = -1;
    }
    // update mapping quality
    if (b->n == 1 && m->n == 1) {
        bsw2hit_t *p = &b->hits[0];
        if (p->flag & BSW2_FLAG_MATESW) { // this alignment is found by Smith-Waterman
            if (!(p->flag & BSW2_FLAG_TANDEM) && b->aux[0].pqual < 20)
                b->aux[0].pqual = 20;
            if (b->aux[0].pqual >= m->aux[0].qual) b->aux[0].pqual = m->aux[0].qual;
        } else if ((p->flag & 2) && !(m->hits[0].flag & BSW2_FLAG_MATESW)) { // properly paired
            if (!(p->flag & BSW2_FLAG_TANDEM)) { // pqual is bounded by [b->aux[0].qual,m->aux[0].qual]
                b->aux[0].pqual += 20;
                if (b->aux[0].pqual > m->aux[0].qual) b->aux[0].pqual = m->aux[0].qual;
                if (b->aux[0].pqual < b->aux[0].qual) b->aux[0].pqual = b->aux[0].qual;
            }
        }
    }
}

// Copied from bwtsw2_aux.c
static void print_hits(const bntseq_t *bns, const bsw2opt_t *opt, bsw2seq1_t *ks, bwtsw2_t *b, int is_pe, bwtsw2_t *bmate)
{
    int i, k;
    kstring_t str;
    memset(&str, 0, sizeof(kstring_t));
    if (b == 0 || b->n == 0) { // no hits
        ksprintf(&str, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t", ks->name);
        for (i = 0; i < ks->l; ++i) kputc(ks->seq[i], &str);
        if (ks->qual) {
            kputc('\t', &str);
            for (i = 0; i < ks->l; ++i) kputc(ks->qual[i], &str);
        } else kputs("\t*", &str);
        kputc('\n', &str);
    }
    for (i = 0; b && i < b->n; ++i) {
        bsw2hit_t *p = b->hits + i;
        bsw2aux_t *q = b->aux + i;
        int j, beg, end, type = 0;
        // print mandatory fields before SEQ
        if (q->cigar == 0) q->flag |= 0x4;
        ksprintf(&str, "%s\t%d", ks->name, q->flag | (opt->multi_2nd && i? 0x100 : 0));
        ksprintf(&str, "\t%s\t%ld", q->chr>=0? bns->anns[q->chr].name : "*", (long)q->pos + 1);
        if (p->l == 0 && q->cigar) { // not a repetitive hit
            ksprintf(&str, "\t%d\t", q->pqual);
            for (k = 0; k < q->n_cigar; ++k)
                ksprintf(&str, "%d%c", q->cigar[k]>>4, (opt->hard_clip? "MIDNHHP" : "MIDNSHP")[q->cigar[k]&0xf]);
        } else ksprintf(&str, "\t0\t*");
        if (!is_pe) kputs("\t*\t0\t0\t", &str);
        else ksprintf(&str, "\t%s\t%d\t%d\t", q->mchr==q->chr? "=" : (q->mchr<0? "*" : bns->anns[q->mchr].name), q->mpos+1, q->isize);
        // get the sequence begin and end
        beg = 0; end = ks->l;
        if (opt->hard_clip && q->cigar) {
            if ((q->cigar[0]&0xf) == 4) beg += q->cigar[0]>>4;
            if ((q->cigar[q->n_cigar-1]&0xf) == 4) end -= q->cigar[q->n_cigar-1]>>4;
        }
        for (j = beg; j < end; ++j) {
            if (p->flag&0x10) kputc(libbwa_nt_comp_table[(int)ks->seq[ks->l - 1 - j]], &str);
            else kputc(ks->seq[j], &str);
        }
        // print base quality if present
        if (ks->qual) {
            kputc('\t', &str);
            for (j = beg; j < end; ++j) {
                if (p->flag&0x10) kputc(ks->qual[ks->l - 1 - j], &str);
                else kputc(ks->qual[j], &str);
            }
        } else kputs("\t*", &str);
        // print optional tags
        ksprintf(&str, "\tAS:i:%d\tXS:i:%d\tXF:i:%d\tXE:i:%d\tNM:i:%d", p->G, p->G2, p->flag>>16, p->n_seeds, q->nm);
        if (q->nn) ksprintf(&str, "\tXN:i:%d", q->nn);
        if (p->l) ksprintf(&str, "\tXI:i:%d", p->l - p->k + 1);
        if (p->flag&BSW2_FLAG_MATESW) type |= 1;
        if (p->flag&BSW2_FLAG_TANDEM) type |= 2;
        if (type) ksprintf(&str, "\tXT:i:%d", type);
        if (opt->cpy_cmt && ks->comment) {
            int l = strlen(ks->comment);
            if (l >= 6 && ks->comment[2] == ':' && ks->comment[4] == ':') {
                kputc('\t', &str); kputs(ks->comment, &str);
            }
        }
        kputc('\n', &str);
    }
    ks->sam = str.s;
    free(ks->seq); ks->seq = 0;
    free(ks->qual); ks->qual = 0;
    free(ks->name); ks->name = 0;
}

// Copied from bwtsw2_aux.c
static void update_opt(bsw2opt_t *dst, const bsw2opt_t *src, int qlen)
{
    double ll = log(qlen);
    int i, k;
    *dst = *src;
    if (dst->t < ll * dst->coef) dst->t = (int)(ll * dst->coef + .499);
    // set band width: the query length sets a boundary on the maximum band width
    k = (qlen * dst->a - 2 * dst->q) / (2 * dst->r + dst->a);
    i = (qlen * dst->a - dst->a - dst->t) / dst->r;
    if (k > i) k = i;
    if (k < 1) k = 1; // I do not know if k==0 causes troubles
    dst->bw = src->bw < k? src->bw : k;
}

// Copied from bwtsw2_aux.c
static void bsw2_aln_core(bsw2seq_t *_seq, const bsw2opt_t *_opt, const bntseq_t *bns, uint8_t *pac, const bwt_t *target, int is_pe)
{
    int x;
    bsw2opt_t opt;
    bsw2global_t *pool = bsw2_global_init();
    bwtsw2_t **buf;
    buf = calloc(_seq->n, sizeof(void*));
    for (x = 0; x < _seq->n; ++x) {
        bsw2seq1_t *p = _seq->seq + x;
        uint8_t *seq[2], *rseq[2];
        int i, l, k;
        bwtsw2_t *b[2];
        l = p->l;
        update_opt(&opt, _opt, p->l);
        if (pool->max_l < l) { // then enlarge working space for aln_extend_core()
            int tmp = ((l + 1) / 2 * opt.a + opt.r) / opt.r + l;
            pool->max_l = l;
            pool->aln_mem = realloc(pool->aln_mem, (tmp + 2) * 24);
        }
        // set seq[2] and rseq[2]
        seq[0] = calloc(l * 4, 1);
        seq[1] = seq[0] + l;
        rseq[0] = seq[1] + l; rseq[1] = rseq[0] + l;
        // convert sequences to 2-bit representation
        for (i = k = 0; i < l; ++i) {
            int c = nst_nt4_table[(int)p->seq[i]];
            if (c >= 4) { c = (int)(drand48() * 4); ++k; } // FIXME: ambiguous bases are not properly handled
            seq[0][i] = c;
            seq[1][l-1-i] = 3 - c;
            rseq[0][l-1-i] = 3 - c;
            rseq[1][i] = c;
        }
        if (l - k < opt.t) { // too few unambiguous bases
            buf[x] = calloc(1, sizeof(bwtsw2_t));
            free(seq[0]); continue;
        }
        // alignment
        b[0] = bsw2_aln1_core(&opt, bns, pac, target, l, seq, pool);
        for (k = 0; k < b[0]->n; ++k)
            if (b[0]->hits[k].n_seeds < opt.t_seeds) break;
        if (k < b[0]->n) {
            b[1] = bsw2_aln1_core(&opt, bns, pac, target, l, rseq, pool);
            for (i = 0; i < b[1]->n; ++i) {
                bsw2hit_t *p = &b[1]->hits[i];
                int x = p->beg;
                p->flag ^= 0x10, p->is_rev ^= 1; // flip the strand
                p->beg = l - p->end;
                p->end = l - x;
            }
            flag_fr(b);
            merge_hits(b, l, 0);
            bsw2_resolve_duphits(0, 0, b[0], 0);
            bsw2_resolve_query_overlaps(b[0], opt.mask_level);
        } else b[1] = 0;
        // generate CIGAR and print SAM
        buf[x] = bsw2_dup_no_cigar(b[0]);
        // free
        free(seq[0]);
        bsw2_destroy(b[0]);
    }
    if (is_pe) bsw2_pair(&opt, bns->l_pac, pac, _seq->n, _seq->seq, buf);
    for (x = 0; x < _seq->n; ++x) {
        bsw2seq1_t *p = _seq->seq + x;
        uint8_t *seq[2];
        int i;
        seq[0] = malloc(p->l * 2); seq[1] = seq[0] + p->l;
        for (i = 0; i < p->l; ++i) {
            int c = nst_nt4_table[(int)p->seq[i]];
            if (c >= 4) c = (int)(drand48() * 4);
            seq[0][i] = c;
            seq[1][p->l-1-i] = 3 - c;
        }
        update_opt(&opt, _opt, p->l);
        write_aux(&opt, bns, p->l, seq, pac, buf[x], _seq->seq[x].name);
        free(seq[0]);
    }
    for (x = 0; x < _seq->n; ++x) {
        if (is_pe) update_mate_aux(buf[x], buf[x^1]);
        print_hits(bns, &opt, &_seq->seq[x], buf[x], is_pe, buf[x^1]);
    }
    for (x = 0; x < _seq->n; ++x) bsw2_destroy(buf[x]);
    free(buf);
    bsw2_global_destroy(pool);
}

// Based on process_seqs in bwtsw2_aux.c
static void _process_seqs(bsw2seq_t *_seq, const bsw2opt_t *opt, const bntseq_t *bns, uint8_t *pac, const bwt_t *target, int is_pe, FILE *out)
{
    int i;
    is_pe = is_pe? 1 : 0;

#ifdef HAVE_PTHREAD
    if (opt->n_threads <= 1) {
        bsw2_aln_core(_seq, opt, bns, pac, target, is_pe);
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
            thread_aux_t *p = data + j;
            p->tid = j; p->_opt = opt; p->bns = bns; p->is_pe = is_pe;
            p->pac = pac; p->target = target;
            p->_seq = calloc(1, sizeof(bsw2seq_t));
            p->_seq->max = (_seq->n + opt->n_threads - 1) / opt->n_threads + 1;
            p->_seq->n = 0;
            p->_seq->seq = calloc(p->_seq->max, sizeof(bsw2seq1_t));
        }
        for (i = 0; i < _seq->n; ++i) { // assign sequences to each thread
            bsw2seq_t *p = data[(i>>is_pe)%opt->n_threads]._seq;
            p->seq[p->n++] = _seq->seq[i];
        }
        for (j = 0; j < opt->n_threads; ++j) pthread_create(&tid[j], &attr, worker, &data[j]);
        for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
        for (j = 0; j < opt->n_threads; ++j) data[j]._seq->n = 0;
        for (i = 0; i < _seq->n; ++i) { // copy the result from each thread back
            bsw2seq_t *p = data[(i>>is_pe)%opt->n_threads]._seq;
            _seq->seq[i] = p->seq[p->n++];
        }
        for (j = 0; j < opt->n_threads; ++j) {
            thread_aux_t *p = data + j;
            free(p->_seq->seq);
            free(p->_seq);
        }
        free(data); free(tid);
    }
#else
    bsw2_aln_core(_seq, opt, bns, pac, target, is_pe);
#endif

    // print and reset
    for (i = 0; i < _seq->n; ++i) {
        bsw2seq1_t *p = _seq->seq + i;
        if (p->sam) err_fprintf(out, "%s", p->sam);
        free(p->name); free(p->seq); free(p->qual); free(p->sam);
        p->tid = -1; p->l = 0;
        p->name = p->seq = p->qual = p->sam = 0;
    }
    err_fflush(out);
    _seq->n = 0;
}

// Based on bsw2_aln in bwtsw2_aux.c
void _sw_aln(const bsw2opt_t *opt, const bntseq_t *bns, bwt_t * const target, const char *fn, const char *fn2, FILE *out)
{
    gzFile fp, fp2;
    kseq_t *ks, *ks2;
    int l, is_pe = 0, i, n;
    uint8_t *pac;
    bsw2seq_t *_seq;
    bseq1_t *bseq;

    pac = calloc(bns->l_pac/4+1, 1);
    for (l = 0; l < bns->n_seqs; ++l)
        err_fprintf(out, "@SQ\tSN:%s\tLN:%d\n", bns->anns[l].name, bns->anns[l].len);
    err_fread_noeof(pac, 1, bns->l_pac/4+1, bns->fp_pac);
    fp = xzopen(fn, "r");
    ks = kseq_init(fp);
    _seq = calloc(1, sizeof(bsw2seq_t));
    if (fn2) {
        fp2 = xzopen(fn2, "r");
        ks2 = kseq_init(fp2);
        is_pe = 1;
    } else fp2 = 0, ks2 = 0, is_pe = 0;
    while ((bseq = bseq_read(opt->chunk_size * opt->n_threads, &n, ks, ks2)) != 0) {
        if (n > _seq->max) {
            _seq->max = n;
            kroundup32(_seq->max);
            _seq->seq = realloc(_seq->seq, _seq->max * sizeof(bsw2seq1_t));
        }
        _seq->n = n;
        for (i = 0; i < n; ++i) {
            bseq1_t *b = &bseq[i];
            bsw2seq1_t *p = &_seq->seq[i];
            p->tid = -1; p->l = b->l_seq;
            p->name = b->name; p->seq = b->seq; p->qual = b->qual; p->comment = b->comment; p->sam = 0;
        }
        free(bseq);
        _process_seqs(_seq, opt, bns, pac, target, is_pe, out);
    }
    // free
    free(pac);
    free(_seq->seq); free(_seq);
    kseq_destroy(ks);
    err_gzclose(fp);
    if (fn2) {
        kseq_destroy(ks2);
        err_gzclose(fp2);
    }
}

// Based on bwa_bwtsw2 in bwtsw2_main.c
int libbwa_sw(const char *db, const char *read, const char *mate,
              const char *out, const libbwa_sw_opt *opt_)
{
    bsw2opt_t *opt;
    bwaidx_t *idx;
    FILE *fpo;

    // Validate arguments
    if (!db || !read || !out || !opt_)
        return LIBBWA_E_INVALID_ARGUMENT;

    opt = bsw2_init_opt();
    _convert_sw_opt(opt_, opt);

    srand48(11);

    opt->t *= opt->a;
    opt->coef *= opt->a;

    fpo = xopen(out, "w");

    if ((idx = bwa_idx_load(db, BWA_IDX_BWT|BWA_IDX_BNS)) == 0) return LIBBWA_E_INDEX_ERROR;
    _sw_aln(opt, idx->bns, idx->bwt, read, mate != NULL ? mate : 0, fpo);
    bwa_idx_destroy(idx);
    free(opt);

    err_fclose(fpo);

    return LIBBWA_E_SUCCESS;
}

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
 */

#include "libbwa.h"

#include "utils.h"

char *bwa_pg;

void bwa_fprint_sam_hdr(FILE *stream, const bntseq_t *bns, const char *rg_line)
{
	int i;
	for (i = 0; i < bns->n_seqs; ++i)
		err_fprintf(stream, "@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
	if (rg_line) err_fprintf(stream, "%s\n", rg_line);
	err_fprintf(stream, "@PG\tID:%s\tPN:%s\tVN:%s\n", LIBBWA_PG_ID, LIBBWA_PG_PN, LIBBWA_PACKAGE_VERSION);
	err_fflush(stream);
}

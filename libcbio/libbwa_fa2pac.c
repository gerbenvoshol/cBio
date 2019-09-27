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

#include <zlib.h>

#include "bntseq.h"
#include "utils.h"

// Based on bwa_fa2pac in bntseq.c
int libbwa_fa2pac(const char *db, const char *prefix, int for_only)
{
    gzFile fp;

    // Validate arguments
    if (!db || !prefix || (for_only < 0 || for_only > 1))
        return LIBBWA_E_INVALID_ARGUMENT;

    fp = xzopen(db, "r");
    bns_fasta2bntseq(fp, prefix, for_only);
	err_gzclose(fp);
	return LIBBWA_E_SUCCESS;
}

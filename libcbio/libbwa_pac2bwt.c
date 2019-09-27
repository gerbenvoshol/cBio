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

#include "bwt.h"
#include "bwtindex.h"

// Based on bwa_pac2bwt in bwtindex.c
int libbwa_pac2bwt(const char *pac, const char *out, int use_is)
{
    bwt_t *bwt;

    // Validate arguments
    if (!pac || !out || (use_is < 0 || use_is > 1))
        return LIBBWA_E_INVALID_ARGUMENT;

    bwt = bwt_pac2bwt(pac, use_is);
    bwt_dump_bwt(out, bwt);
    bwt_destroy(bwt);
    return LIBBWA_E_SUCCESS;
}

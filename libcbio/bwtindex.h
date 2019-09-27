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
 * This header is added for using two functions, bwa_seq_len and bwt_pac2bwt, in
 * libbwa.
 */

#ifndef BWTINDEX_H
#define BWTINDEX_H

#include "bwt.h"

int64_t bwa_seq_len(const char *fn_pac);
bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is);

#endif

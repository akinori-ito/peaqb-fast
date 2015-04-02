/*
 * Copyright (c) 2003, Giuseppe Gottardi 'oveRet' <gottardi@ailinux.org>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <stdlib.h>
#include <math.h>
#include <common.h>
#include <reldistframes.h>


extern int bark;

double
reldistframes(double *Pnoise, double *M, double *reldisttmp, int n)
{
    int k;

    for(k=0;k<bark;k++) {
	if(10.0*log10((double)Pnoise[k]/M[k]) >= 1.5) {
	    *reldisttmp = *reldisttmp  + 1;
	    break;
	}
    }

    return ((double)*reldisttmp/n);
}

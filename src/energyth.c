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
#include <energyth.h>


int
energyth(signed int *test, signed int *ref, int hann)
{
    int k;
    double sum;

    sum = 0;
    for(k=0;k<hann/2;k++) {
	sum += p(test[hann/2 + k], 2.0);
	if(sum > ENERGYLIMIT)
	    return 1;
    }

    sum = 0;
    for(k=0;k<hann/2;k++) {
	sum += p(ref[hann/2 + k], 2.0);
	if(sum > ENERGYLIMIT)
	    return 1;
    }

    return 0;
}

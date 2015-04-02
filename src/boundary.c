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
#include <boundary.h>

int
boundary(signed int *ch1ref, signed int *ch1test, signed int *ch2ref,
         signed int *ch2test, int hann)
{
    int k, i, sum;
    int ch1t = 0, ch1r = 0, ch2t = 0, ch2r = 0;

    for(k=0;k<hann-BOUNDWIN+1;k++) {
	if(!ch1t) {
	    sum = 0;
	    for(i=0;i<BOUNDWIN;i++)
		sum += abs(ch1test[k+i]);
	    if(sum > BOUNDLIMIT)
		ch1t = 1;
	}

	if(!ch1r) {
	    sum = 0;
	    for(i=0;i<BOUNDWIN;i++)
		sum += abs(ch1ref[k+i]);
	    if(sum > BOUNDLIMIT)
		ch1r = 1;
	 }

	if(ch1t || ch1r) // || or &&
	    return 1;
    }

    if(ch2test == NULL && ch2ref == NULL)
	return 0;

    for(k=0;k<hann-BOUNDWIN+1;k++) {
	if(!ch2t) {
	    sum = 0;
	    for(i=0;i<BOUNDWIN;i++)
		sum += abs(ch2test[k+i]);
	    if(sum > BOUNDLIMIT)
		ch2t = 1;
	}

	if(!ch2r) {
	    sum = 0;
	    for(i=0;i<BOUNDWIN;i++)
		sum += abs(ch2ref[k+i]);
	    if(sum > BOUNDLIMIT)
		ch2r = 1;
	}

	if((ch1t || ch2t) || (ch1r || ch2r)) // || or &&
	    return 1;
    }

    return 0;
}

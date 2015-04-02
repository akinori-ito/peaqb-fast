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
#include <loudness.h>


extern double *fC;
extern int bark;

double
loudness(double *E)
{
    int k;
    double Ntot = 0;
    double s, N, Ethres;

    for(k=0;k<bark;k++) {
	s = p(10.0, (-2.0-2.05*atan((double)fC[k]/4000.0) - 0.75
	    *atan((double)p(fC[k]/1600.0, 2.0)))/10.0);
 	Ethres = p(10.0, 0.364*p(fC[k]/1000.0, -0.8));
	N = (double)CONST*p(Ethres/(s*10000.0), 0.23)
	    *(p(1.0-s+s*E[k]/Ethres, 0.23) - 1.0);
	if(N > 0)
	    Ntot += N;
    }

    Ntot *= (double)24.0/bark;

    return Ntot;
}

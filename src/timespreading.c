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
#include <timespreading.h>

#define T100 0.03
#define Tmin 0.008

extern int bark;
extern double *fC;

int
timespreading(double *E2, double *Etmp, int rate, double *E)
{
    int k;
    double a, T;

    for(k=0;k<bark;k++) {
	T = (double)Tmin + (double)(100.0/fC[k])*(double)(T100- Tmin);
	a = exp((double)-HANN/(double)(T*2.0*rate));
	Etmp[k] = Etmp[k]*a + (1.0-a)*E2[k];
	if(Etmp[k] >= E2[k])
	    E[k] = Etmp[k];
	else
	    E[k] = E2[k];
    }

    return 0;
}

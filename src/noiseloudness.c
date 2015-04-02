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
#include <levpatadapt.h>
#include <noiseloudness.h>


extern int bark;
extern double *fC;

double
noiseloudness(double *Modtest, double *Modref, struct levpatadaptout lev, 
	      double *nltmp, int n)
{
    int k;
    double Pthres, stest, sref, beta, num, denom;
    double nl = 0;

    for(k=0;k<bark;k++) {
	Pthres = p(10.0, 0.4*0.364*p(fC[k]/1000.0, -0.8));
	stest = (double)THRESFAC0*Modtest[k] + S0;
	sref = (double)THRESFAC0*Modref[k] + S0;
	if(lev.Eptest[k] == 0 && lev.Epref[k] == 0)
	    beta = 1.0;
	else
	if(lev.Epref[k] == 0)
	    beta = 0;
	else
	    beta = exp((double)-ALPHA*(lev.Eptest[k] - lev.Epref[k])
	           /lev.Epref[k]);
	
	num = stest*lev.Eptest[k] - sref*lev.Epref[k];
	denom = Pthres + sref*lev.Epref[k]*beta;
	if(num < 0)
	    num = 0;

	nl += p(Pthres/(E0*stest), 0.23)*(p(1.0 + num/denom, 0.23) - 1.0);
    }

    nl *= (double)24.0/bark;
    if(nl < 0)
	nl = 0;

    *nltmp += p(nl, 2.0);

    return sqrt((double)*nltmp/n);
}

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
#include <spreading.h>


extern double *fC;
extern int bark;

int
spreading(double *pp, double *e2)
{
    int k, j, u;
    double L;
    double denom, sum1, sum2, Eline, Su, Sl = 27.0;
    double pplog[bark],su1[bark],su2[bark],el[bark],denom1[bark],denom2[bark];

    for(j=0;j<bark;j++) {
	pplog[j] = 10.0*log10(pp[j]);
	su1[j] = -24.0-230.0/fC[j]+0.2*pplog[j];
	su2[j]= -24.0-230.0/fC[j]; // L = 0
	el[j] = p(10.0, pplog[j]/10.0);
	denom1[j] = 0;
	for(u=0;u<j;u++)
	    denom1[j] += p(10.0, -dz*(j-u)*Sl/10.0);
	for(u=j;u<bark;u++)
	    denom1[j] += p(10.0, dz*(u-j)*su1[j]/10.0);
	denom2[j] = 0;
	for(u=0;u<j;u++)
	    denom2[j] += p(10.0, -dz*(j-u)*Sl/10.0);
	for(u=j;u<bark;u++)
	    denom2[j] += p(10.0, dz*(u-j)*(-24.0-230.0/fC[j])/10.0);
    }

    for(k=0;k<bark;k++) {
	sum1 = 0;
	sum2 = 0;
	// Eline
	for(j=0;j<bark;j++) {
	    L = pplog[j];
	    Su = su1[j];
	    Eline = el[j];
	    if(k < j)
		Eline *= p(10.0, -dz*(j-k)*Sl/10.0);
	    else
		Eline *= p(10.0, dz*(k-j)*Su/10.0);
	    
	    Eline /= denom1[j];
	    sum1 += p(Eline, 0.4);
	}

	// Eline (tilde)
	for(j=0;j<bark;j++) {
	    Su = su2[j];
	    if(k < j)
		Eline = p(10.0, -dz*(j-k)*Sl/10.0);
	    else
		Eline = p(10.0, dz*(k-j)*Su/10.0);

	    Eline /= denom2[j];
	    sum2 += p(Eline, 0.4);
	}

	sum2 = p(sum2, 1.0/0.4);
	sum1 = p(sum1, 1.0/0.4);

	e2[k] = sum1/sum2;
    }

    return 0;
}

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
#include <detprob.h>


extern int bark;

double
detprob(double *Etestch1, double *Etestch2, double *Erefch1,
	double *Erefch2, double *Ptildetmp, double *PMtmp,
	double *Qsum, int *ndistorcedtmp, int hann)
{
    int k;
    double Etilderefch1, Etilderefch2, Etildetestch1, Etildetestch2;
    double L, s, e, b, a, pch1, pch2, pbin,
    	   qch1, qch2, qbin, P, c0, c1, ADBb;
    double prod = 1.0, Q = 0.0;

    for(k=0;k<bark;k++) {
	Etildetestch1 = 10.0*log10((double)Etestch1[k]);
	Etilderefch1 = 10.0*log10((double)Erefch1[k]);

	if(Etilderefch1 > Etildetestch1)
	    L = 0.3*Etilderefch1;
	else
	    L = 0.3*Etildetestch1;

	L += 0.7*Etildetestch1;

	if(L > 0)
	    s = 5.95072*p(6.39468/L, 1.71332)+9.01033*p(10.0, -11.0)
	        *p(L, 4.0)+5.05622*p(10.0, -6.0)*p(L, 3.0)-0.00102438
		*p(L, 2.0)+0.0550197*L-0.198719;
	else
	    s = p(10.0, 30.0);

	e = Etilderefch1 - Etildetestch1;

	if(Etilderefch1 > Etildetestch1)
	    b = 4.0;
	else
	    b = 6.0;

	a = (double)p(10.0, log10((double)log10((double)2.0))/b)/s;
	pch1 = 1.0 - p(10.0, -p(a*e, b));
	qch1 = abs((int)e)/s;  // don't touch this

	pbin = pch1;
	qbin = qch1;

	if(Etestch2 != NULL && Erefch2 != NULL) {
	    Etildetestch2 = 10.0*log10((double)Etestch2[k]);
	    Etilderefch2 = 10.0*log10((double)Erefch2[k]);

	    if(Etilderefch2 > Etildetestch2)
		L = 0.3*Etilderefch2;
	    else
		L = 0.3*Etildetestch2;

	    L += 0.7*Etildetestch2;

	    if(L > 0)
		s = 5.95072*p(6.39468/L, 1.71332)+9.01033*p(10.0, -11.0)
		    *p(L, 4.0)+5.05622*p(10.0, -6.0)*p(L, 3.0)
		    -0.00102438*p(L, 2.0)+0.0550197*L-0.198719;
	    else
		s = 1.0*p(10.0, 30.0);

	    e = Etilderefch2 - Etildetestch2;

	    if(e > 0)
		b = 4.0;
	    else
		b = 6.0;

	    a = (double)p(10.0, log10((double)log10((double)2.0))/b)/s;
	    pch2 = 1.0 - p(10.0, -p(a*e, b));
	    qch2 = abs((int)e)/s;  // don't touch this

	    if(pch2 > pch1)
		pbin = pch2;
	    if(qch2 > qch1)
		qbin = qch2;
	}

	prod *= (1.0 - pbin);
	Q += qbin;
    }

    P = 1.0 - prod;
    if(P > 0.5) {
	*Qsum += Q;
	(*ndistorcedtmp)++;
    }

    if(*ndistorcedtmp == 0)
	ADBb = 0;
    else
    if(*Qsum > 0)
	ADBb = log10((double)*Qsum / (*ndistorcedtmp));
    else
	ADBb = -0.5;

    c0 = p(0.9, hann/(2.0*1024.0));
    #if !defined(C1)
    c1 = p(0.99, hann/(2.0*1024.0));
    #else
    c1 = C1;
    #endif
    *Ptildetmp = (1.0 - c0)*P + (*Ptildetmp)*c0;
    if(*Ptildetmp > (*PMtmp)*c1)
	*PMtmp = *Ptildetmp;
    else
	*PMtmp = (*PMtmp)*c1;

    return ADBb;
}

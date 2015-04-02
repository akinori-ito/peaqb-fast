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

extern int bark;
extern double *fC;

struct levpatadaptout
levpatadapt(double *Etest, double *Eref, int rate,
	    struct levpatadaptin *tmp, int hann)
{
    int k, i, m1, m2;
    double T, levcorr, numlevcorr = 0, denomlevcorr = 0, R;
    double pattcoeffref, pattcoefftest;
    double Elref[BARK], Eltest[BARK], Rtest[BARK], Rref[BARK], a[BARK];
    struct levpatadaptout out;

    for(k=0;k<bark;k++) {
	T = (double)Tmin + (100.0/fC[k])*(T100 - Tmin);
	a[k] = exp((double)-hann/(2.0*rate * T));
	tmp->Ptest[k] = tmp->Ptest[k]*a[k] + (1.0-a[k])*Etest[k];
	tmp->Pref[k] = tmp->Pref[k]*a[k] + (1.0-a[k])*Eref[k];
	numlevcorr += sqrt(tmp->Ptest[k]*tmp->Pref[k]);
	denomlevcorr += tmp->Ptest[k];
    }

    levcorr = p(numlevcorr/denomlevcorr, 2.0);

    for(k=0;k<bark;k++) {
	if(levcorr > 1.0) {
	    Elref[k] = (double)Eref[k]/levcorr;
	    Eltest[k] = Etest[k];
	}
	else {
	    Eltest[k] = (double)Etest[k]*levcorr;
	    Elref[k] = Eref[k];
	}

	// Autocorrelation
	tmp->Rnum[k] *= a[k];
	tmp->Rdenom[k] *= a[k];
	tmp->Rnum[k]  += Elref[k]*Eltest[k];
	tmp->Rdenom[k] += Elref[k]*Elref[k];

	if(tmp->Rdenom[k] == 0 && tmp->Rnum[k] != 0) {
	    Rtest[k] = 0;
	    Rref[k] = 1.0;
	}
	else
	if(tmp->Rdenom[k] == 0 && tmp->Rnum[k] == 0) {
	    //copy from frequency band below
	    if(k) {
		Rtest[k] = Rtest[k-1];
		Rref[k] = Rref[k-1];
	    }
	    //if don't exist
	    else {
		Rtest[k] = 1.0;
		Rref[k] = 1.0;
	    }
	}
	else {
	    R = tmp->Rnum[k] / tmp->Rdenom[k];
	    if(R >= 1.0) {
		Rtest[k] = 1.0/R;
		Rref[k] = 1.0;
	    }
	    else {
		Rtest[k] = 1.0;
		Rref[k] = R;
	    }
	}
    }

    for(k=0;k<bark;k++) {
	m1 = M1;
	m2 = M2;
	pattcoefftest = 0;
	pattcoeffref = 0;

	if(m1 > k)
	    m1 = k;
	if(m2 > bark -k -1)
	    m2 = bark -k -1;

	for(i = -m1;i <= m2;i++) {
	    pattcoefftest += Rtest[k+i];
	    pattcoeffref += Rref[k+i];
	}

	tmp->PattCorrTest[k] = a[k]*tmp->PattCorrTest[k] + 
			       pattcoefftest*(1.0-a[k])/(m1+m2+1);
	tmp->PattCorrRef[k] = a[k]*tmp->PattCorrRef[k] + 
			      pattcoeffref*(1.0-a[k])/(m1+m2+1);

	out.Epref[k] = Elref[k] * tmp->PattCorrRef[k];
	out.Eptest[k] = Eltest[k] * tmp->PattCorrTest[k];
    }

    return out;
}

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
#include <string.h>
#include <fftw.h>
#include <common.h>
#include <harmstruct.h>


extern char *filetest;
extern int count;
extern fftw_plan plan2;

double Cfft[HANN/2];
int maxk;

double
harmstruct(double *ffttest, double *fftref, double *EHStmp, int rate,
	   double *Cffttmp, int p, int *n)
{
    int k, i;
    double F0[HANN/2], C[HANN/2], hannwin[HANN/2];
    double num, denoma, denomb, Csum = 0;
    fftw_complex in[HANN/2], out[HANN/2];
    double max;

    bzero(Cfft, 8 * HANN/2);

    for(k=0;k<p*2-1;k++) {
	if(!fftref[k] || !ffttest[k]) {  // skip log(0)
	    #if defined(SKIPFRAME) && !defined(ZERO)
	    (*n)--;
	    return 0;
	    #elif defined(ZERO)
	    if(!fftref[k]) {
		fftref[k] = ZERO;
		#ifdef DEBUG
		debug("Warning [%s:%d] in Harmstruct.c: "
		      "fftref[%d] is set around zero\n",
		      filetest, count, k, fftref[k]);
		#endif
	    }
	    if(!ffttest[k]) {
		ffttest[k] = ZERO;
		#ifdef DEBUG
		debug("Warning [%s:%d] in Harmstruct.c: "
		      "ffttest[%d] is set around zero\n",
		      filetest, count, k, ffttest[k]);
		#endif
	    }
	    F0[k] = log10(p(fftref[k], 2.0)) - log10(p(ffttest[k], 2.0));
	    #else
	    if(!fftref[k] && !ffttest[k])
		F0[k] = 0;
	    else {
		#ifdef DEBUG
		debug("Error [%s:%d] in Harmstruct.c: "
		      "log(fftref[%d] = %g  log(ffttest[%d] = %g\n",
		      filetest, count, k, fftref[k], k, ffttest[k]);
		#endif
		F0[k] = 0;

	    }
	    #endif
	}
	else
	    F0[k] = log10(p(fftref[k], 2.0)) - log10(p(ffttest[k], 2.0));
    }

    for(i=0;i<p;i++) {
	num = 0;
	denoma = 0;
	denomb = 0;
	for(k=0;k<p;k++) {
	    num += F0[k] * F0[i+k];
	    denoma += p(F0[k], 2.0);
	    denomb += p(F0[i+k], 2.0);
	}

	hannwin[i] = 0.5*sqrt((double)8.0/3.0)*(1.0 - cos((double)2.0
		     *M_PI*i/(p-1.0)));
	C[i] = num / (sqrt((double)denoma) * sqrt((double)denomb));

	#if !defined(AVGHANN)
	C[i] *= hannwin[i];
	#endif
	Csum += C[i];
    }

    for(i=0;i<p;i++) {
	C[i] -= (double)Csum/p;
	#if defined(AVGHANN)
	C[i] *= hannwin[i];
	#endif
	in[i].im = 0;
	in[i].re = C[i];
    }
    
    fftw_one(plan2, in, out);

    for(k=0;k<p/2;k++) {
	out[k].re *= (double)(1.0/p);
	out[k].im *= (double)(1.0/p);
	Cfft[k] = p(out[k].re, 2.0) + p(out[k].im, 2.0);
    }

    #ifdef EHSMODO2
    for(k=0;k<p/2;k++) {
        Cffttmp[k] += Cfft[k];
	Cfft[k] = Cffttmp[k]/(*n);
    }
    #endif

    #if !defined(GETMAX)
    i = 0+PATCH;
    while(1) {
	if(Cfft[i] >= Cfft[i+1]) {
	    while(i < p/2-1 && Cfft[i] >= Cfft[i+1])
	        i++;

	    if(i < p/2-1)
	        break;
	    else {
	        (*n)--;
	        return 0;
	    }
	}
	else {
	    while(i < p/2-1 && Cfft[i] <= Cfft[i+1])
	        i++;
	    while(i < p/2-1 && Cfft[i] >= Cfft[i+1])
	        i++;
	    if(i < p/2-1)
	        break;
	    else {
	        (*n)--;
	        return 0;
	    }
        }
    }
    #else
    i = 0;
    #endif

    max = 0;
    for(k=i+1;k<p/2;k++)
	if(Cfft[k] > max) {
	    max = Cfft[k];
	    maxk = k;
	}

    #ifdef EHSMODO2
    return max*1000.0;
    #endif

    (*EHStmp) += max;
    return ((*EHStmp)*1000.0/(*n));
}

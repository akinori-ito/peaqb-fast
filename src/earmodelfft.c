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
#include <fftw.h>
#include <common.h>
#include <earmodelfft.h>


extern double hannwindow[];
extern fftw_plan plan;

int
earmodelfft(signed int *ch, int lp, int hann, double *ffte, double *absfft)
{
    int k;
    double w, fac;
    fftw_complex in[HANN], out[HANN];
    
    fac = p(10.0, lp/20.0)/NORM;

    for(k=0;k<hann;k++) {
	in[k].re = hannwindow[k] * (double)ch[k];
	in[k].im = 0;
    }
    
    fftw_one(plan, in, out);
    
    for(k=0;k<hann/2;k++) {
	out[k].re *= (double)(fac/hann);
	out[k].im *= (double)(fac/hann);
	absfft[k] = sqrt((double)(p(out[k].re, 2.0) + p(out[k].im, 2.0)));

	w = -0.6*3.64*p(k * FREQADAP/1000.0, -0.8) +
	    6.5*exp((double)-0.6*p(k * FREQADAP/1000.0 - 3.3, 2.0)) -
	    0.001*p(k * FREQADAP/1000.0, 3.6);

	ffte[k] = absfft[k]*p(10.0, w/20.0);

    }

    return 0;
}

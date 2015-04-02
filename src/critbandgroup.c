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
#include <critbandgroup.h>

extern double *fL, *fC, *fU;
extern int bark;

int
critbandgroup(double *ffte, int rate, int hann, double *pe)
{
    double fres;
    int i, k;

    fres = (double)rate/hann;

    for(i=0;i<bark;i++) {
	pe[i] = 0;
	for(k=0;k<hann/2;k++) {
	    if(((double)(k-0.5)*fres) >= fL[i] 
	       && ((double)(k+0.5)*fres <= fU[i]))
		pe[i] += p(ffte[k], 2.0);
	    else
		if(((double)(k-0.5)*fres) < fL[i] 
		   && ((double)(k+0.5)*fres > fU[i]))
		    pe[i] += p(ffte[k], 2.0)*(fU[i]-fL[i])/fres;
		else
		if(((double)(k-0.5)*fres) < fL[i] 
		   && ((double)(k+0.5)*fres > fL[i]))
		    pe[i] += p(ffte[k], 2.0)*(double)((k+0.5)
		             *fres-fL[i])/fres;
		else
		if(((double)(k-0.5)*fres) < fU[i] 
		   && ((double)(k+0.5)*fres > fU[i]))
		    pe[i] += p(ffte[k], 2.0)*(fU[i]-(double)(k-0.5)
		             *fres)/fres;
	}

	if(pe[i] < p(10.0, -12.0))
	    pe[i] = p(10.0, -12.0);
    }

    return 0;
}

int
AddIntNoise(double *pe)
{
    int k;
    double Pthres;

    for(k=0;k<bark;k++) {
	Pthres = p(10.0, 0.4*0.364*p(fC[k]/1000.0, -0.8));
	pe[k] += Pthres;
    }

    return 0;
}

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
#include <bandwidth.h>

int
bandwidth(double *ffttest, double *fftref, int hann,
	  struct bandwidthout *out)
{
    int k, BwRef = 0, BwTest = 0;
    double Flevtest, Flevref;
    double ZeroThreshold;

    ZeroThreshold = 20.0 * log10((double)ffttest[ZEROTHRESHOLD]);

    for(k=ZEROTHRESHOLD;k<hann/2;k++) {
	Flevtest = 20.0 * log10((double)ffttest[k]);

	if(Flevtest > ZeroThreshold)
	    ZeroThreshold = Flevtest;
    }

    for(k=ZEROTHRESHOLD-1;k>=0;k--) {
	Flevref = 20.0 * log10((double)fftref[k]);

	if(Flevref >= 10.0+ZeroThreshold) {
	    BwRef = k + 1;
	    break;
	}
    }

    for(k=BwRef-1;k>=0;k--) {
	Flevtest = 20.0 * log10((double)ffttest[k]);

	if(Flevtest >= 5.0+ZeroThreshold) {
	    BwTest = k + 1;
	    break;
	}
    }

    if(BwRef > BwMAX) {
	out->sumBandwidthRefb += (double)BwRef;
	out->countref++;
    }

    if(BwTest > BwMAX) {
	out->sumBandwidthTestb += (double)BwTest;
	out->counttest++;
    }

    if(out->countref == 0)
	out->BandwidthRefb = 0;
    else
	out->BandwidthRefb = out->sumBandwidthRefb/(double)out->countref;

    if(out->counttest == 0)
	out->BandwidthTestb = 0;
    else
	out->BandwidthTestb = out->sumBandwidthTestb/(double)out->counttest;

    return 0;
}

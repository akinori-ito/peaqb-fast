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
#include <assert.h>
#include <math.h>
#include <common.h>
#include <moddiff.h>

extern double *fC;
extern int bark;

struct moddiffout
moddiff(double *Modtest, double *Modref, double *Etilderef)
{
    int k;
    struct moddiffout out;
    double Pthres;

    out.ModDiff1 = 0;
    out.ModDiff2 = 0;
    out.TempWt = 0;

    for(k=0;k<bark;k++) {
	// WinModDiff1B && AvgModDiff1B
	out.ModDiff1 += module(Modtest[k] - Modref[k])/(1.0 + Modref[k]);
	// AvgModDiff2B
	if(Modtest[k] > Modref[k])
	    out.ModDiff2 += module(Modtest[k] - Modref[k])
	    		    /(0.01 + Modref[k]);
	else
	    out.ModDiff2 += 0.1*module(Modtest[k] - Modref[k])
	    		    /(0.01 + Modref[k]);

	Pthres = p(10.0, 0.4*0.364*p(fC[k]/1000.0, -0.8));
	out.TempWt += Etilderef[k]/(Etilderef[k] + p(Pthres, 0.3)*100.0);
    }

    out.ModDiff1 *= (double)100.0/bark;
    out.ModDiff2 *= (double)100.0/bark;

    return out;
}

double
ModDiff1(struct moddiffout in, struct moddiffin *intmp, int n)
{
    int i;

    intmp->mod[intmp->Lcount] = in.ModDiff1;
    intmp->Lcount++;
    if(intmp->Lcount == L)
    intmp->Lcount = 0;

    if(n < L)
	return 0;

    intmp->modtmp = 0;
    for(i=0;i<L;i++)
	intmp->modtmp += sqrt((double)intmp->mod[i]);

    intmp->modtmp /= (double)L;

    intmp->win += p(intmp->modtmp, 4.0);

    return sqrt((double)intmp->win/(double)(n-L+1.0));
}

double
ModDiff2(struct moddiffout in, struct moddiffin *intmp)
{
    intmp->num2 += in.ModDiff1 * in.TempWt;
    intmp->denom2 += in.TempWt;

    return (intmp->num2/intmp->denom2);
}

double
ModDiff3(struct moddiffout in, struct moddiffin *intmp)
{
    intmp->num3 += in.ModDiff2 * in.TempWt;
    intmp->denom3 += in.TempWt;

    return (intmp->num3/intmp->denom3);
}

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
#include <neural.h>

double amin[11] = {393.916656, 361.965332, -24.045116, 1.110661, -0.206623,
		   0.074318, 1.113683, 0.950345, 0.029985, 0.000101, 0.0};
double amax[11] = {921.0, 881.131226, 16.212030, 107.137772, 2.886017,
		   13.933351, 63.257874, 1145.018555, 14.819740, 1.0, 1.0};
double wx[12][3] = {{-0.502657, 0.436333, 1.219602},
		    {4.307481, 3.246017, 1.123743},
		    {4.984241, -2.211189, -0.192096},
		    {0.051056, -1.762424, 4.331315},
		    {2.321580, 1.789971, -0.754560},
		    {-5.303901, -3.452257, -10.814982},
		    {2.730991, -6.111805, 1.519223},
		    {0.624950, -1.331523, -5.955151},
		    {3.102889, 0.871260, -5.922878},
		    {-1.051468, -0.939882, -0.142913},
		    {-1.804679, -0.503610, -0.620456},
		    {-2.518254, 0.654841, -2.207228}};
double wy[4] = {-3.817048, 4.107138, 4.629582, -0.307594};
double bmin = -3.98;
double bmax = 0.22;


struct out
neural(struct outframes processed)
{
    int j, i;
    struct out oveRet;
    double sum1, sum2;
    double x[11];

    x[0] = processed.BandwidthRefb;
    x[1] = processed.BandwidthTestb;
    x[2] = processed.TotalNMRb;
    x[3] = processed.WinModDiff1b;
    x[4] = processed.ADBb;
    x[5] = processed.EHSb;
    x[6] = processed.AvgModDiff1b;
    x[7] = processed.AvgModDiff2b;
    x[8] = processed.RmsNoiseLoudb;
    x[9] = processed.MFPDb;
    x[10] = processed.RelDistFramesb;

	// gcodcla.wav
	/*x[0] = 834.117;
	x[1] = 647.095;
	x[2] = -14.6048;
	x[3] = 6.89483;
	x[4] = 0.432969;
	x[5] = 0.503605;
	x[6] = 7.14863;
	x[7] = 24.9353;
	x[8] = 0.124738;
	x[9] = 0.968876;
	x[10] = 0.0485208;  // così è tutto ok*/

	// ccodsax.wav
	/*x[0] = 853.375;
	x[1] = 645.444;
	x[2] = -7.94882;
	x[3] = 11.4108;
	x[4] = 1.41971;
	x[5] = 0.491164;
	x[6] = 12.6383;
	x[7] = 44.7187;
	x[8] = 0.21807;
	x[9] = 1.15;//0.675505;
	x[10] = 0.556215;*/

    sum2 = 0;
    for(j=0;j<J;j++) {
	sum1 = 0;
	for(i=0;i<I;i++)
	    sum1 += wx[i][j]*((x[i] - amin[i])/(amax[i] - amin[i]));
	    sum2 += wy[j]*sig(wx[I][j] + sum1);
    }

    oveRet.DI = wy[J] + sum2;
    oveRet.ODG = bmin + (bmax - bmin)*sig(oveRet.DI);

    return oveRet;
}


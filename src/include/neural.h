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

#define sig(x) (double)(1.0/(1.0 + exp((double)-(x))))

#define I 11
#define J 3

struct outframes {
    double WinModDiff1b;
    double AvgModDiff1b;
    double AvgModDiff2b;
    double RmsNoiseLoudb;
    double BandwidthRefb;
    double BandwidthTestb;
    double TotalNMRb;
    double RelDistFramesb;
    double ADBb;
    double MFPDb;
    double EHSb;
};

struct out {
    double ODG;
    double DI;
};

/* Function prototypes */
struct out neural(struct outframes processed);
/* Prototypes end */

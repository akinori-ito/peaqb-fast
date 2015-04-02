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

#define T100 0.05
#define Tmin 0.008
#define M 8
// if M is odd
/*
#define M1 (M-1)/2
#define M2 M1
*/
// if M is even
#define M1 M/2 - 1
#define M2 M/2

struct levpatadaptout {
    double Epref[BARK];
    double Eptest[BARK];
};

struct levpatadaptin {
    double Ptest[BARK];
    double Pref[BARK];
    double PattCorrTest[BARK];
    double PattCorrRef[BARK];
    double Rnum[BARK];
    double Rdenom[BARK];
};

/* Function prototypes */
struct levpatadaptout levpatadapt(double *, double *, int,
				  struct levpatadaptin *, int);
/* Prototypes end */

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

#define L 4

struct moddiffin {
    double win;
    int Lcount;
    double modtmp;
    double mod[L];
    double num2;
    double denom2;
    double num3;
    double denom3;
};

struct moddiffout {
    double ModDiff1;
    double ModDiff2;
    double TempWt;
};

/* Function prototypes */
struct moddiffout moddiff(double *, double *, double *);
double ModDiff1(struct moddiffout, struct moddiffin *, int);
double ModDiff2(struct moddiffout, struct moddiffin *);
double ModDiff3(struct moddiffout, struct moddiffin *);
/* Prototypes end */

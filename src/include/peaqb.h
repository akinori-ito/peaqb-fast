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

#if HAVE_CONFIG
#include <config.h>
#endif

#define LOGRESULT "analized"
/*
 * Il file dove vengono scritti i risulatati.
 */

#ifdef DEBUG
#define LOGFILE "debugged.txt"
#endif
/*
 * Il file dove vengono loggate eventuali warning da logaritmi.
 */

#define OPT_REF      0x01
#define OPT_TEST        0x02

#define THRESHOLDDELAY 0.050
#define AVERAGINGDEALAY 0.5

#define B(f) 7 * asinh((double)f /650)
#define BI(z) 650 * sinh((double)z /7)


/* Function prototypes */
void fatalerr(char *,...);
void usage(char *);
void logvariable(const char *, double *, int);
void ProcessFrame(signed int *, signed int *, int, signed int *,
                  signed int *, int, int, int, int);
/* Prototypes end */

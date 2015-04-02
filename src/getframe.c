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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getframe.h>


#include <errno.h>
//extern int errno;

signed int
GetFrameValue(FILE *fp, int bytes)
{
    int intvalue;

    if (bytes <= 0)
	return 0;
    intvalue = ReadInt(fp, bytes);
    switch(bytes) {
	case 3:
	    if (intvalue & 0x00800000)
		intvalue |= 0xff000000;
	    break;
	case 2:
	    if (intvalue & 0x00008000)
		intvalue |= 0xffff0000;
	    break;
	case 1:
	    if (intvalue & 0x00000080)
		intvalue |= 0xffffff00;
	    break;
    }

    return (signed int) intvalue;
}

int
GetMonoFrame(FILE *fp, signed int *vect, int bytes, int hann)
{
    int i = 0;

    if (fp == NULL)
	return 0;

    if (fseek(fp, (-hann/2)*bytes, SEEK_CUR) == -1)
	fatalerr("err: %s", strerror(errno));

    while(!feof(fp) && i < hann) {
	vect[i] = GetFrameValue(fp, bytes);
	i++;
    }

    if(i < hann) {
	bzero(vect, hann*4);
	fseek(fp, -i*bytes, SEEK_END);
	return 0;
    }

    /* Number of samples wrote */
    return i ;
}

int
GetStereoFrame(FILE *fp, signed int *sx, signed int *dx, int bytes,
	       int hann)
{
    int i = 0, k = 0, count = 0;

    if (fp == NULL)
	return 0;

    if (fseek(fp, -hann*bytes, SEEK_CUR) == -1)
	fatalerr("err: %s", strerror(errno));

    while(!feof(fp) && count < hann*2) {
	if(!k) {
	    sx[i] = GetFrameValue(fp, bytes);
	    k = 1;
	    i--;
	    count++;
	}
	else {
	    dx[i] = GetFrameValue(fp, bytes);
	    k = 0;
	    count++;
	}
	i++;
    }
    if (count < hann*2) {
	bzero(sx, hann*4);
	bzero(dx, hann*4);
	fseek(fp, -count*bytes, SEEK_END);
	return 0;
    }

    /* Number of samples wrote */
    return count ;
}

int
LevelPression(char *f)
{
    int lp;

    while(*f != ':' && *f)
	f++;
    
    if(*f == ':') {
	lp = atoi(f + 1);
	*f = '\0';
	return lp;
    }
    else
	return LP;
}

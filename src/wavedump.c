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
#include <wavedump.h>


#include <errno.h>
//extern int errno;

int
HeaderDump(FILE *fp, const char *string)
{
    char buff[7];
    int k, offset = 0;

    while (!feof(fp)) {
	if(!fread(buff, 7, 1, fp))
	    fatalerr("err: error in WaveHeader");
	offset += 7;
	for(k=0;k<4;k++)
	    if (!strncmp((char *)(buff + k), string, 4)) {
		fseek(fp, k-3, SEEK_CUR);
		offset += k-3;
		return offset;
	    }
	fseek(fp, -3, SEEK_CUR);
	offset += -3;
    }
	
    return -1;
}

int
SampleRate(FILE *fp)
{
    int rate, offset;

    if ((offset = HeaderDump(fp,"fmt ")) == -1)
	return -1;
    if (fseek(fp, 8, SEEK_CUR) == -1)
	fatalerr("err: %s", strerror(errno));

    rate = ReadInt(fp,4);
    fseek(fp, -8 - offset, SEEK_CUR);

    return rate;
}

int
BitForSample(FILE *fp)
{
    int bit, offset;

    if ((offset = HeaderDump(fp,"fmt ")) == -1)
	return -1;
    if (fseek(fp, 18, SEEK_CUR) == -1)
	fatalerr("err: %s", strerror(errno));

    bit = ReadInt(fp,2);
    fseek(fp, -18 - offset, SEEK_CUR);

    return bit;
}

int
NumOfChan(FILE *fp)
{
    int chan, offset;

    if ((offset = HeaderDump(fp,"fmt ")) == -1)
	return -1;
    if (fseek(fp, 6, SEEK_CUR) == -1)
	fatalerr("err: %s", strerror(errno));

    chan = ReadInt(fp,2);
    fseek(fp, -6 - offset, SEEK_CUR);

    return chan;
}

int
FindData(FILE *fp)
{
    int offset;

    if (fp == NULL)
	return -1;

    if ((offset = HeaderDump(fp,"data")) == -1)
	return -1;
    if (fseek(fp, 4 + offset, SEEK_CUR) == -1)
	fatalerr("err: %s", strerror(errno));
    
    return 1;
}

#if defined(LITTLE) && !defined(BIG)
int
ReadInt(FILE *fp, int size)
{
    int l;
    unsigned char c;

    if (size <= 0)
	return 0;
    
    c = fgetc(fp);
    l = ((int) c) & 255;
    l |= (ReadInt(fp, size-1) << 8);
    
    return l;
}
#elif defined(BIG) && !defined(LITTLE)
int
ReadInt(FILE *fp, int size)
{
    int l;

    if (size <= 0)
	return 0;
        
    l = (ReadInt(fp, size-1) << 8);
    l |= ((int) fgetc(fp)) & 255;
    
    return l;
}
#endif

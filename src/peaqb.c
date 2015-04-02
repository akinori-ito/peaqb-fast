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
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <getopt.h>
#include <assert.h>
#include <math.h>
#include <fftw.h>
#include <common.h>
#include <wavedump.h>
#include <getframe.h>
#include <bandwidth.h>
#include <levpatadapt.h>
#include <moddiff.h>
#include <modulation.h>
#include <loudness.h>
#include <neural.h>
#include <nmr.h>
#include <detprob.h>
#include <energyth.h>
#include <harmstruct.h>
#include <boundary.h>
#include <critbandgroup.h>
#include <earmodelfft.h>
#include <noiseloudness.h>
#include <reldistframes.h>
#include <spreading.h>
#include <timespreading.h>
#include <threshold.h>
#include <peaqb.h>


#include <errno.h>
//extern int errno;
char *fileref, *filetest;
double hannwindow[HANN];
double Etesttmpch1[BARK], Etesttmpch2[BARK], Ereftmpch1[BARK],
       Ereftmpch2[BARK], Cffttmpch1[HANN/2], Cffttmpch2[HANN/2];
int delaytime1, delaytime2;
int count = 0;
int harmsamples = 1;
fftw_plan plan, plan2;

/* Bark Tables */
double *fL, *fC, *fU;
int bark;
int csv_output = 0;

struct levpatadaptin levinch1, levinch2;
struct modulationin modintestch1, modintestch2, modinrefch1, modinrefch2;
struct moddiffin moddiffinch1, moddiffinch2;
struct bandwidthout bandwidthch1, bandwidthch2;
struct outframes processed;


int
main(int argc, char *argv[])
{

    signed int ch1ref[HANN];
    signed int ch2ref[HANN];
    signed int ch1test[HANN];
    signed int ch2test[HANN];

    int opt_line = 0;
    int rateref, numchref, bitsampleref, lpref;
    int ratetest, numchtest, bitsampletest, lptest;
    int boundflag, totalframes = 0;
    FILE *fpref, *fptest;

    struct boundaryflag boundbe = {0, 0};
    struct out oveRet;

    /* Parse command line */
    if (argc < 3)
	usage(argv[0]);


    {
    int c = 0;

    while ((c = getopt(argc, argv, "r:t:h:c")) != EOF)

	switch (c) {
	    case 'h':
		usage(argv[0]);
		break;
	    case 'r':
		opt_line |= OPT_REF;
		fileref = optarg;
		break;
	    case 't':
		opt_line |= OPT_TEST;
		filetest = optarg;
		break;
            case 'c':
		csv_output = 1;
		break;
	}
    }

    /* Input control */
    if (!(opt_line & OPT_REF) || !*fileref)
	fatalerr ("err: -r/--reference <arg> required");

    if (!(opt_line & OPT_TEST) || !*filetest)
	fatalerr ("err: -t/--test <arg> required");

    lpref = LevelPression(fileref);
    lptest = LevelPression(filetest);

    /* Init routines */

    // make Hann Window  (2.1.3)
    {
    int k;

    for(k=0;k<HANN;k++)
	hannwindow[k] = 0.5*sqrt((double)8/3)*
			(1 - cos((double)2*M_PI*k/(HANN -1)));
    }

    // make Bark tables  (2.1.5)
    {
    int k;
    double zL, zU;
    double *zl, *zc, *zu;

    zL = B(Flow);
    zU = B(Fup);

    bark = ceil((zU - zL) / dz);

    fL = (double *)malloc(bark * sizeof(double));
    fC = (double *)malloc(bark * sizeof(double));
    fU = (double *)malloc(bark * sizeof(double));
    zl = (double *)malloc(bark * sizeof(double));
    zc = (double *)malloc(bark * sizeof(double));
    zu = (double *)malloc(bark * sizeof(double));
    assert(fL != NULL && fC != NULL && fU != NULL && zl != NULL 
    	   && zc != NULL && zu != NULL);

    for(k=0;k<bark;k++) {
	zl[k] = zL + k*dz;
	zu[k] = zL + (k+1)*dz;
	zc[k] = 0.5 * (zl[k] + zu[k]);
    }

    zu[bark-1] = zU;
    zc[bark -1] = 0.5 * (zl[bark-1] + zu[bark-1]);

    for(k=0;k<bark;k++) {
	fL[k] = BI(zl[k]);
	fU[k] = BI(zu[k]);
	fC[k] = BI(zc[k]);
    }

    free(zl);
    free(zu);
    free(zc);
    }

    // Initialize temp var
    memset(&levinch1, 0x00, sizeof(struct levpatadaptin));
    memset(&levinch2, 0x00, sizeof(struct levpatadaptin));

    memset(&modintestch1, 0x00, sizeof(struct modulationin));
    memset(&modintestch2, 0x00, sizeof(struct modulationin));
    memset(&modinrefch1, 0x00, sizeof(struct modulationin));
    memset(&modinrefch2, 0x00, sizeof(struct modulationin));

    memset(Etesttmpch1, 0x00, BARK * sizeof(double));
    memset(Etesttmpch2, 0x00, BARK * sizeof(double));
    memset(Ereftmpch1, 0x00, BARK * sizeof(double));
    memset(Ereftmpch2, 0x00, BARK * sizeof(double));
    memset(Cffttmpch1, 0x00, (HANN/2) * sizeof(double));
    memset(Cffttmpch2, 0x00, (HANN/2) * sizeof(double));

    memset(&moddiffinch1, 0x00, sizeof(struct moddiffin));
    memset(&moddiffinch2, 0x00, sizeof(struct moddiffin));

    memset(&bandwidthch1, 0x00, sizeof(struct bandwidthout));
    memset(&bandwidthch2, 0x00, sizeof(struct bandwidthout));

    // ref file
    if ((fpref = fopen(fileref,"r")) == NULL)
	fatalerr("err: %s", strerror(errno));
    if ((rateref = SampleRate(fpref)) == -1)
	fatalerr("err: error in WaveHeader");
    if ((numchref = NumOfChan(fpref)) == -1)
	fatalerr("err: error in WaveHeader");
    if ((bitsampleref = BitForSample(fpref)) == -1)
	fatalerr("err: error in WaveHeader");
    if(FindData(fpref) == -1)
	fatalerr("err: can't find Data Field");

    // test file
    if ((fptest = fopen(filetest,"r")) == NULL)
	fatalerr("err: %s", strerror(errno));
    if ((ratetest = SampleRate(fptest)) == -1)
	fatalerr("err: error in WaveHeader");
    if ((numchtest = NumOfChan(fptest)) == -1)
	fatalerr("err: error in WaveHeader");
    if ((bitsampletest = BitForSample(fptest)) == -1)
	fatalerr("err: error in WaveHeader");
    if(FindData(fptest) == -1)
	fatalerr("err: can't find Data Field");

    if (csv_output) {
        fprintf(stdout,"\n PEAQb Algorithm. Author Giuseppe Gottardi 'oveRet'"
	           " <gottardi@ailinux.org>\n");

        fprintf(stdout,"\nRef File, %s"
	               "\n  - Sample Rate, %d"
	               "\n  - Number Of Channel, %d"
	               "\n  - Bits for Sample, %d"
		       "\n  - Level Playback, %d\n\n", fileref, rateref,
		numchref, bitsampleref, lpref);

        fprintf(stdout,"\nTest File, %s"
	               "\n  - Sample Rate, %d"
	               "\n  - Number Of Channel, %d"
	               "\n  - Bits for Sample, %d"
		       "\n  - Level Playback, %d\n\n", filetest, ratetest,
		numchtest, bitsampletest, lptest);
	fprintf(stdout,"frame,"
#ifdef DATABOUND_BE
		       "TotalFrame,BoundaryFrom,BoundaryTo,"
#endif
		       "BandwidthRefb, BandwidthTestb, TotalNMRb, WinModDiff1b, ADBb, EHSb, AvgModDiff1b, AvgModDiff2b, RmsNoiseLoudb, MFPDb, RelDistFramesb, DI, ODG\n");
    }
    else {
        fprintf(stdout,"\n PEAQb Algorithm. Author Giuseppe Gottardi 'oveRet'"
	               " <gottardi@ailinux.org>\n");

        fprintf(stdout,"\nRef File %s"
	               "\n  - Sample Rate: %d"
	               "\n  - Number Of Channel: %d"
	               "\n  - Bits for Sample: %d"
		       "\n  - Level Playback: %d\n\n", fileref, rateref,
		numchref, bitsampleref, lpref);

        fprintf(stdout,"\nTest File %s"
	               "\n  - Sample Rate: %d"
	               "\n  - Number Of Channel: %d"
	               "\n  - Bits for Sample: %d"
		       "\n  - Level Playback: %d\n\n", filetest, ratetest,
		numchtest, bitsampletest, lptest);
    }

    // Processing
    if(ratetest != rateref)
	fatalerr("err: Can't process Wave Files with different Sample Rate");
    if(numchref != numchtest)
	fatalerr("err: Can't process Mono Wave with Stereo Wave");

    // Find delaytime1 for Loudness Threshold
    delaytime1 = ceilf((float)THRESHOLDDELAY*ratetest*2/HANN);
    // Find delaytime2 for Delayed Averaging
    delaytime2 = ceilf((float)AVERAGINGDEALAY*ratetest*2/HANN);

    // make fft plan
    plan = fftw_create_plan(HANN, FFTW_FORWARD, FFTW_MEASURE);

    while(harmsamples < (Fup/ratetest)*(HANN/2.0)/2.0)
	harmsamples *= 2;
    plan2 = fftw_create_plan(harmsamples, FFTW_FORWARD, FFTW_MEASURE);

    if(numchref == 1) {
	if (fseek(fpref, (HANN/2)*bitsampleref/8, SEEK_CUR) == -1)
	    fatalerr("err: %s", strerror(errno));
	if (fseek(fptest, (HANN/2)*bitsampletest/8, SEEK_CUR) == -1)
	    fatalerr("err: %s", strerror(errno));

	#ifdef DATABOUND_BE
	#undef DATABOUND_ONE
	{
	int i = 0, flag = 0, f1, f2;
	long dataref, datatest, br1, br2;

	dataref = ftell(fpref);
	datatest = ftell(fptest);

	while(1) {
	    br1 = ftell(fpref);
	    br2 = ftell(fptest);
	    f1 = GetMonoFrame(fpref, (signed int *)ch1ref,
		                bitsampleref/8, HANN);
	    f2 = GetMonoFrame(fptest, (signed int *)ch1test,
			        bitsampletest/8, HANN);
	    if(f1 && f2) {
	    totalframes++;
		if(boundary(ch1ref, ch1test, NULL, NULL, HANN) && !flag) {
		    boundbe.begin = totalframes;
		    flag = 1;
	 	}
	    }
	    else {
	        fseek(fptest, br1, SEEK_SET);
		fseek(fpref, br2, SEEK_SET);
		break;
	    }
	}

	fseek(fptest, -(HANN/2)*bitsampletest/8, SEEK_CUR);
	fseek(fpref, -(HANN/2)*bitsampleref/8, SEEK_CUR);
	while(i<totalframes) {
	    GetMonoFrame(fpref, (signed int *)ch1ref, bitsampleref/8, HANN);
	    GetMonoFrame(fptest, (signed int *)ch1test, bitsampletest/8, HANN);
	    fseek(fptest, -2*(HANN/2)*bitsampletest/8, SEEK_CUR);
	    fseek(fpref, -2*(HANN/2)*bitsampleref/8, SEEK_CUR);
	    i++;
	    if(boundary(ch1ref, ch1test, NULL, NULL, HANN)) {
		boundbe.end = totalframes-i;
		break;
	    }
	}

	fseek(fptest, datatest, SEEK_SET);
	fseek(fpref, dataref, SEEK_SET);
	}
	#endif

	while (GetMonoFrame(fpref, (signed int *)ch1ref,
			    bitsampleref/8, HANN)
	       && GetMonoFrame(fptest, (signed int *)ch1test,
	                    bitsampletest/8, HANN)) {

	    count++;
	    #ifdef DATABOUND_BE
	    if(count >= boundbe.begin && count <= boundbe.end)
		boundflag = 1;
	    else
		boundflag = 0;
	    #else
	    boundflag = boundary(ch1ref, ch1test, NULL, NULL, HANN);
	       #ifdef DATABOUND_ONE
	        {
	        static int flag1 = 0, flag2 = 0;

	        if(boundflag && !flag1)
		    flag1 = 1;
		if(!boundflag && flag1)
		    flag2 = 1;
		if(flag2)
		    boundflag = 0;
		}
	        #endif
	    #endif

	    ProcessFrame((signed int *)ch1ref,
	                 (signed int *)NULL, lpref,
			 (signed int *)ch1test,
			 (signed int *)NULL,
			 lptest, rateref, boundflag, HANN);

	    oveRet = neural(processed);

	    if (csv_output) {
	        fprintf(stdout,"%d,"
#ifdef DATABOUND_BE
			       "%d,%d,%d,"
#endif
			       "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
			   count,
			   #ifdef DATABOUND_BE
			   totalframes, boundbe.begin, boundbe.end,
			   #endif
			   processed.BandwidthRefb,
			   processed.BandwidthTestb, processed.TotalNMRb,
			   processed.WinModDiff1b, processed.ADBb,
			   processed.EHSb, processed.AvgModDiff1b,
			   processed.AvgModDiff2b,
			   processed.RmsNoiseLoudb, processed.MFPDb,
			   processed.RelDistFramesb,
			   oveRet.DI, oveRet.ODG);
	    }
	    else {
	        fprintf(stdout,"\nframe: %d"
			   #ifdef DATABOUND_BE
			   "/%d"
			   "\ndata boundary: %d -> %d"
			   #endif
			   "\nBandwidthRefb: %g"
	    		   "\nBandwidthTestb: %g"
			   "\nTotalNMRb %g"
			   "\nWinModDiff1b: %g"
			   "\nADBb: %g"
			   "\nEHSb: %g"
			   "\nAvgModDiff1b: %g"
			   "\nAvgModDiff2b %g"
			   "\nRmsNoiseLoudb: %g"
			   "\nMFPDb: %g"
			   "\nRelDistFramesb: %g"
			   "\nDI: %g"
			   "\nODG: %g\n",
			   count,
			   #ifdef DATABOUND_BE
			   totalframes, boundbe.begin, boundbe.end,
			   #endif
			   processed.BandwidthRefb,
			   processed.BandwidthTestb, processed.TotalNMRb,
			   processed.WinModDiff1b, processed.ADBb,
			   processed.EHSb, processed.AvgModDiff1b,
			   processed.AvgModDiff2b,
			   processed.RmsNoiseLoudb, processed.MFPDb,
			   processed.RelDistFramesb,
			   oveRet.DI, oveRet.ODG);
	    }
	}

	{
	FILE *res;

	res = fopen(LOGRESULT,"a+");
	fprintf(res,"\nFile: %s\n"
		    "\nframe: %d"
		    "\nBandwidthRefb: %g"
		    "\nBandwidthTestb: %g"
		    "\nTotalNMRb %g"
		    "\nWinModDiff1b: %g"
		    "\nADBb: %g"
		    "\nEHSb: %g"
		    "\nAvgModDiff1b: %g"
		    "\nAvgModDiff2b %g"
		    "\nRmsNoiseLoudb: %g"
		    "\nMFPDb: %g"
		    "\nRelDistFramesb: %g"
		    "\nDI: %g"
		    "\nODG: %g\n",
		    filetest, count, processed.BandwidthRefb,
		    processed.BandwidthTestb, processed.TotalNMRb,
		    processed.WinModDiff1b, processed.ADBb,
		    processed.EHSb, processed.AvgModDiff1b,
		    processed.AvgModDiff2b,
		    processed.RmsNoiseLoudb, processed.MFPDb,
		    processed.RelDistFramesb,
		    oveRet.DI, oveRet.ODG);

	fclose(res);
	}
    }

    if(numchref == 2) {
	if (fseek(fpref, HANN*bitsampleref/8, SEEK_CUR) == -1)
	    fatalerr("err: %s", strerror(errno));
	if (fseek(fptest, HANN*bitsampletest/8, SEEK_CUR) == -1)
	    fatalerr("err: %s", strerror(errno));

	#ifdef DATABOUND_BE
	#undef DATABOUND_ONE
	{
	int i = 0, flag = 0, f1, f2;
	long dataref, datatest, br1, br2;

	dataref = ftell(fpref);
	datatest = ftell(fptest);

	while(1) {
	    br1 = ftell(fpref);
	    br2 = ftell(fptest);
	    f1 = GetStereoFrame(fpref, (signed int *)ch1ref,
		                (signed int *)ch2ref, bitsampleref/8, HANN);
	    f2 = GetStereoFrame(fptest, (signed int *)ch1test,
			        (signed int *)ch2test, bitsampletest/8, HANN);
	    if(f1 && f2) {
	    totalframes++;
		if(boundary(ch1ref, ch1test, ch2ref, ch2test, HANN) && !flag) {
		    boundbe.begin = totalframes;
		    flag = 1;
	 	}
	    }
	    else {
	        fseek(fptest, br1, SEEK_SET);
		fseek(fpref, br2, SEEK_SET);
		break;
	    }
	}

	fseek(fptest, -HANN*bitsampletest/8, SEEK_CUR);
	fseek(fpref, -HANN*bitsampleref/8, SEEK_CUR);
	while(i<totalframes) {
	    GetStereoFrame(fpref, (signed int *)ch1ref,
	    			(signed int *)ch2ref, bitsampleref/8, HANN);
	    GetStereoFrame(fptest, (signed int *)ch1test,
	    			(signed int *)ch2test, bitsampletest/8, HANN);
	    fseek(fptest, -2*HANN*bitsampletest/8, SEEK_CUR);
	    fseek(fpref, -2*HANN*bitsampleref/8, SEEK_CUR);
	    i++;
	    if(boundary(ch1ref, ch1test, ch2ref, ch2test, HANN)) {
		boundbe.end = totalframes-i+1;
		break;
	    }
	}

	fseek(fptest, datatest, SEEK_SET);
	fseek(fpref, dataref, SEEK_SET);
	}
	#endif

	while (GetStereoFrame(fpref, (signed int *)ch1ref,
			      (signed int *)ch2ref, bitsampleref/8, HANN)
	       && GetStereoFrame(fptest, (signed int *)ch1test,
	       			 (signed int *)ch2test, bitsampletest/8, HANN)) {

	    count++;
	    #ifdef DATABOUND_BE
	    if(count >= boundbe.begin && count <= boundbe.end)
		boundflag = 1;
	    else
		boundflag = 0;
	    #else
	    boundflag = boundary(ch1ref, ch1test, ch2ref, ch2test, HANN);
	        #ifdef DATABOUND_ONE
	        {
	        static int flag1 = 0, flag2 = 0;

	        if(boundflag && !flag1)
		    flag1 = 1;
		if(!boundflag && flag1)
		    flag2 = 1;
		if(flag2)
		    boundflag = 0;
		}
	        #endif
	    #endif

	    ProcessFrame((signed int *)ch1ref,
	                 (signed int *)ch2ref, lpref,
		         (signed int *)ch1test,
			 (signed int *)ch2test,
			 lptest, rateref, boundflag, HANN);

	    oveRet = neural(processed);

	    fprintf(stdout,"\nframe: %d"
			   #ifdef DATABOUND_BE
			   "/%d"
			   "\ndata boundary: %d -> %d"
			   #endif
			   "\nBandwidthRefb: %g"
	    		   "\nBandwidthTestb: %g"
			   "\nTotalNMRb %g"
			   "\nWinModDiff1b: %g"
			   "\nADBb: %g"
			   "\nEHSb: %g"
			   "\nAvgModDiff1b: %g"
			   "\nAvgModDiff2b %g"
			   "\nRmsNoiseLoudb: %g"
			   "\nMFPDb: %g"
			   "\nRelDistFramesb: %g"
			   "\nDI: %g"
			   "\nODG: %g\n",
			   count,
			   #ifdef DATABOUND_BE
			   totalframes, boundbe.begin, boundbe.end,
			   #endif
			   processed.BandwidthRefb,
			   processed.BandwidthTestb, processed.TotalNMRb,
			   processed.WinModDiff1b, processed.ADBb,
			   processed.EHSb, processed.AvgModDiff1b,
			   processed.AvgModDiff2b,
			   processed.RmsNoiseLoudb, processed.MFPDb,
			   processed.RelDistFramesb,
			   oveRet.DI, oveRet.ODG);
	}

	{
	FILE *res;

	res = fopen(LOGRESULT,"a+");
	fprintf(res,"\nFile: %s\n"
		    "\nframe: %d"
		    "\nBandwidthRefb: %g"
		    "\nBandwidthTestb: %g"
		    "\nTotalNMRb %g"
		    "\nWinModDiff1b: %g"
		    "\nADBb: %g"
		    "\nEHSb: %g"
		    "\nAvgModDiff1b: %g"
		    "\nAvgModDiff2b %g"
		    "\nRmsNoiseLoudb: %g"
		    "\nMFPDb: %g"
		    "\nRelDistFramesb: %g"
		    "\nDI: %g"
		    "\nODG: %g\n",
		    filetest, count, processed.BandwidthRefb,
		    processed.BandwidthTestb, processed.TotalNMRb,
		    processed.WinModDiff1b, processed.ADBb,
		    processed.EHSb, processed.AvgModDiff1b,
		    processed.AvgModDiff2b,
		    processed.RmsNoiseLoudb, processed.MFPDb,
		    processed.RelDistFramesb,
		    oveRet.DI, oveRet.ODG);

	fclose(res);
	}  
    }

    fftw_destroy_plan(plan);
    fclose(fpref);
    fclose(fptest);
    return 0;
}


void
ProcessFrame(signed int *ch1ref, signed int *ch2ref, int lpref,
	     signed int *ch1test, signed int *ch2test, int lptest,
	     int rate, int boundflag, int hann)
{
    int k;
    static int ch = 1;
    double Ntotaltest, Ntotalref;
    struct levpatadaptout lev;
    struct moddiffout mod;
    struct processing processch1, processch2;

    /*if(!boundflag)
        printf("\n[%d] -> scartato", count);
    return;*/

    earmodelfft(ch1ref, lpref, hann, processch1.ffteref,
                processch1.fftref);
    earmodelfft(ch1test, lptest, hann, processch1.fftetest,
    		processch1.ffttest);

    critbandgroup(processch1.ffteref, rate, hann, processch1.ppref);
    AddIntNoise(processch1.ppref);

    critbandgroup(processch1.fftetest, rate, hann, processch1.pptest);
    AddIntNoise(processch1.pptest);

    for(k=0;k<hann/2;k++)
	processch1.fnoise[k] = module(processch1.ffteref[k])
	 		       - module(processch1.fftetest[k]);

    critbandgroup(processch1.fnoise, rate, hann, processch1.ppnoise);

    spreading(processch1.pptest, processch1.E2test);
    spreading(processch1.ppref, processch1.E2ref);

    timespreading(processch1.E2test, Etesttmpch1, rate, processch1.Etest);
    timespreading(processch1.E2ref, Ereftmpch1, rate, processch1.Eref);

    threshold(processch1.Eref, processch1.Mref);

    modulation(processch1.E2test, rate, &modintestch1, processch1.Modtest);
    modulation(processch1.E2ref, rate, &modinrefch1, processch1.Modref);

    // Data boundary
    if(boundflag) {
	static int countboundary = 1;
	static double RelDistFramesb = 0, nmrtmp = 0;

	bandwidth(processch1.ffttest, processch1.fftref, hann,
	          &bandwidthch1);
	processed.BandwidthRefb = bandwidthch1.BandwidthRefb;
	processed.BandwidthTestb = bandwidthch1.BandwidthTestb;

	processed.TotalNMRb = nmr(processch1.ppnoise, processch1.Mref,
				  &nmrtmp, countboundary);
	processed.RelDistFramesb = reldistframes(processch1.ppnoise,
						 processch1.Mref,
						 &RelDistFramesb,
						 countboundary);
	countboundary++;

	// Data boundary + Energy threshold
	if(energyth(ch1test, ch1ref, hann)) {
	    static int countenergy = 1;
	    static double EHStmp = 0;

	    processed.EHSb = harmstruct(processch1.ffttest,
	    				processch1.fftref,
					&EHStmp, rate, Cffttmpch1,
					harmsamples, &countenergy);
	    countenergy++;
	}
    }

    // Delayed Averaging
    if(count > delaytime2) {
	static double nltmp = 0;
	static int noise = 0, internal_count = 0, loudcounter = 0;

	mod = moddiff(processch1.Modtest, processch1.Modref,
	              (double *)&(modinrefch1.Etildetmp));
	processed.WinModDiff1b = ModDiff1(mod, &moddiffinch1,
					  count - delaytime2);
	processed.AvgModDiff1b = ModDiff2(mod, &moddiffinch1);
	processed.AvgModDiff2b = ModDiff3(mod, &moddiffinch1);

	Ntotaltest = loudness(processch1.Etest);
	Ntotalref = loudness(processch1.Eref);

	if(Ntotaltest > 0.1 || Ntotalref > 0.1) {
	    noise = 1;
	    #if defined(LOUDMODO2)
	    internal_count = 0;
	    #endif
	}

	// Delayed Averaging + loudness threshold
	if(noise && internal_count <= delaytime1) {
	    // skip 0.05 sec (about 3 frames)
	    internal_count++;
	    loudcounter++;
	}
	else {
	    lev = levpatadapt(processch1.Etest, processch1.Eref, rate,
			      &levinch1, hann);
	    processed.RmsNoiseLoudb = noiseloudness(processch1.Modtest,
	    					    processch1.Modref,
						    lev, &nltmp,
						    count - delaytime2
						    - loudcounter);
        }
    }

    /*{
    extern double Cfft[];
    extern int maxk;
    FILE *fp;

    logvariable("Cfftsx.txt", Cfft, 128);
    fp = fopen("Cfftsxmaxpos.txt", "a+");
    fprintf(fp,"%d\n",maxk);
    fclose(fp);
    }*/

    if(ch2ref && ch2test && *ch2ref && *ch2test) {
	ch = 2;

	earmodelfft(ch2ref, lpref, hann, processch2.ffteref,
		    processch2.fftref);
	earmodelfft(ch2test, lptest, hann, processch2.fftetest,
		    processch2.ffttest);

	critbandgroup(processch2.ffteref, rate, hann, processch2.ppref);
	AddIntNoise(processch2.ppref);

	critbandgroup(processch2.fftetest, rate, hann, processch2.pptest);
	AddIntNoise(processch2.pptest);

	for(k=0;k<hann/2;k++)
	    processch2.fnoise[k] = module(processch2.ffteref[k])
	    				  - module(processch2.fftetest[k]);

	critbandgroup(processch2.fnoise, rate, hann, processch2.ppnoise);

	spreading(processch2.pptest, processch2.E2test);
	spreading(processch2.ppref, processch2.E2ref);

	timespreading(processch2.E2test, Etesttmpch2, rate,
		      processch2.Etest);
	timespreading(processch2.E2ref, Ereftmpch2, rate,
		      processch2.Eref);

	threshold(processch2.Eref, processch2.Mref);

	modulation(processch2.E2test, rate, &modintestch2,
		   processch2.Modtest);
	modulation(processch2.E2ref, rate, &modinrefch2,
		   processch2.Modref);

	// Data boundary
	if(boundflag) {
	    static int countboundary = 1;
	    static double RelDistFramesb = 0, nmrtmp = 0;

	    bandwidth(processch2.ffttest, processch2.fftref, hann,
	    	      &bandwidthch2);
	    processed.BandwidthRefb += bandwidthch2.BandwidthRefb;
	    processed.BandwidthTestb += bandwidthch2.BandwidthTestb;
	    processed.BandwidthRefb /= 2.0;
	    processed.BandwidthTestb /= 2.0;

	    processed.TotalNMRb += nmr(processch2.ppnoise,
	    			       processch2.Mref, &nmrtmp,
				       countboundary);
	    processed.RelDistFramesb += reldistframes(processch2.ppnoise,
	    					      processch2.Mref,
						      &RelDistFramesb,
						      countboundary);
	    processed.TotalNMRb /= 2.0;
	    processed.RelDistFramesb /= 2.0;
	    countboundary++;

	    // Data boundary + Energy threshold
	    if(energyth(ch2test, ch2ref, hann)) {
		static int countenergy = 1;
		static double EHStmp = 0;

		processed.EHSb += harmstruct(processch2.ffttest,
					     processch2.fftref,
					     &EHStmp, rate, Cffttmpch2,
					     harmsamples, &countenergy);
		processed.EHSb /= 2.0;
		countenergy++;
	    }
	}

	// Delayed Averaging
	if(count > delaytime2) {
	    static double nltmp = 0;
	    static int noise = 0, internal_count = 0, loudcounter = 0;

	    mod = moddiff(processch2.Modtest, processch2.Modref,
	    		  (double *)&(modinrefch2.Etildetmp));
	    processed.WinModDiff1b += ModDiff1(mod, &moddiffinch2,
	    				       count - delaytime2);
	    processed.AvgModDiff1b += ModDiff2(mod, &moddiffinch2);
	    processed.AvgModDiff2b += ModDiff3(mod, &moddiffinch2);
	    processed.WinModDiff1b /= 2.0;
	    processed.AvgModDiff1b /= 2.0;
	    processed.AvgModDiff2b /= 2.0;

	    Ntotaltest = loudness(processch2.Etest);
	    Ntotalref = loudness(processch2.Eref);

	    if(Ntotaltest > 0.1 || Ntotalref > 0.1) {
		noise = 1;
		#if defined(LOUDMODO2)
		internal_count = 0;
		#endif
	    }

	    // Delayed Averaging + loudness threshold
	    if(noise && internal_count <= delaytime1) {
		// skip 0.05 sec (about 3 frames)
		internal_count++;
		loudcounter++;
	    }
	    else {
	        lev = levpatadapt(processch2.Etest, processch2.Eref, rate,
			          &levinch2, hann);
		processed.RmsNoiseLoudb += noiseloudness(processch2.Modtest,
							 processch2.Modref,
							 lev, &nltmp,
							 count - delaytime2
							 - loudcounter);
		processed.RmsNoiseLoudb /= 2.0;
	    }
	}
    }

    {
    static int ndistorcedtmp = 0;
    static double Ptildetmp = 0, PMtmp = 0, Qsum = 0;

    if(ch == 2)
	processed.ADBb = detprob(processch1.Etest, processch2.Etest,
				 processch1.Eref, processch2.Eref,
				 &Ptildetmp, &PMtmp, &Qsum,
				 &ndistorcedtmp, hann);
    else
	processed.ADBb = detprob(processch1.Etest, NULL, 
				 processch1.Eref, NULL, 
				 &Ptildetmp, &PMtmp, &Qsum,
				 &ndistorcedtmp, hann);
    processed.MFPDb = PMtmp;
    }
    /*
    #ifdef LOGVARIABLE
    logvariable("fftetestsx.txt", processch1.fftetest, hann/2);
    logvariable("ffterefsx.txt", processch1.ffteref, hann/2);
    logvariable("ffttestsx.txt", processch1.ffttest, hann/2);
    logvariable("fftrefsx.txt", processch1.fftref, hann/2);
    logvariable("Etestsx.txt", processch1.Etest, bark);
    logvariable("Erefsx.txt", processch1.Eref, bark);
    logvariable("E2testsx.txt", processch1.E2test, bark);
    logvariable("E2refsx.txt", processch1.E2ref, bark);
    logvariable("pptestsx.txt", processch1.pptest, bark);
    logvariable("pprefsx.txt", processch1.ppref, bark);
    logvariable("ppnoisesx.txt", processch1.ppnoise, bark);
    logvariable("Mrefsx.txt", processch1.Mref, bark);
    logvariable("Modtestsx.txt", processch1.Modtest, bark);
    logvariable("Modrefsx.txt", processch1.Modref, bark);

    logvariable("fftetestdx.txt", processch2.fftetest, hann/2);
    logvariable("ffterefdx.txt", processch2.ffteref, hann/2);
    logvariable("ffttestdx.txt", processch2.ffttest, hann/2);
    logvariable("fftrefdx.txt", processch2.fftref, hann/2);
    logvariable("Etestdx.txt", processch2.Etest, bark);
    logvariable("Erefdx.txt", processch2.Eref, bark);
    logvariable("E2testdx.txt", processch2.E2test, bark);
    logvariable("E2refdx.txt", processch2.E2ref, bark);
    logvariable("pptestdx.txt", processch2.pptest, bark);
    logvariable("pprefdx.txt", processch2.ppref, bark);
    logvariable("ppnoisedx.txt", processch2.ppnoise, bark);
    logvariable("Mrefdx.txt", processch2.Mref, bark);
    logvariable("Modtestdx.txt", processch2.Modtest, bark);
    logvariable("Modrefdx.txt", processch2.Modref, bark);
    #endif
    */
 
    /*{
    extern double Cfft[];
    extern int maxk;
    FILE *fp;

    logvariable("Cfftdx.txt", Cfft, 128);
    fp = fopen("Cfftdxmaxpos.txt", "a+");
    fprintf(fp,"%d\n",maxk);
    fclose(fp);
    }*/
    return;
}



void
fatalerr(char * pattern,...) /* Error handling routine */
{
    va_list ap;
    
    va_start(ap, pattern);

    fprintf(stderr,"PEAQ-");
    vfprintf(stderr,pattern,ap);
    fprintf(stderr," (exit forced).\n");

    va_end(ap);

    exit(-1);
}


#ifdef DEBUG
void
debug(char * pattern,...) /* Debug handling routine */
{
    FILE *log;
    va_list ap;

    va_start(ap, pattern);

    log = fopen(LOGFILE,"a+");
    vfprintf(log,pattern,ap);

    va_end(ap);
    fclose(log);

    return;
}
#endif


#ifdef LOGVARIABLE
void
logvariable(const char *filename, double *var, int len)
{
    FILE *fp;
    int k;

    #ifdef LOGALLFRAMES
    fp = fopen(filename,"a+");
    #else
    fp = fopen(filename,"w");
    #endif

    for(k=0;k<len;k++)
    	fprintf(fp,"%g\n",var[k]);

    fclose(fp);
    return;
}
#endif


void
usage(char * name) /* Print usage */
{
    fprintf(stderr, "PEAQ Algorithm. Giuseppe Gottardi 'oveRet'"
    		    "<gottardi@ailinux.org>\n\n");
    
    fprintf(stderr, "usage: %s <option>\n", name);
    fprintf(stderr, "     -r  reffile[:lp] (lp default = 92)\n"
    		    "     -t  testfile[:lp] (lp default = 92)\n"
                    "     -h  print this help\n");

    exit (0);
}

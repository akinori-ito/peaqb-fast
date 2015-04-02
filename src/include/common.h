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

#define DEBUG

// servono solo per l'allocazione della memoria
#define HANN 2048
#define BARK 109

#define DOUBLE

#if defined(DOUBLE)
    #define module(x) fabs((double) (x))
    #define p(x,y) pow((double)(x), (double)(y))
#elif defined(LDOUBLE)
    #define module(x) fabsl((long double) x)
    #define p(x,y) powl((long double)x, (long double)y)
#endif

/************************ detprob ************************/
#define C1 1.0   // Should be 1.0
/*
 * Se non è definito se lo calcola
 */
/************************** end **************************/

/*********************** harmstruct **********************/
#define AVGHANN   //(secondo i canadesi)
/*
 * Questa define modifica la sequenza delle operazioni
 * facendo si che prima della finestratura venga tolta
 * la componente continua di C.
 */
#define SKIPFRAME
/*
 * Questa define fa si che il frame contenente
 * un log(0) venga interamente scartato,
 * altrimenti il valore di F0[k] viene settato a 0.
 */
//#define ZERO 0.0001
/*
 * Questa define fa si che non venga scartato il log(0)
 * ma venga calcolato come log() di un numero prossimo a 0
 */
#define GETMAX
/*
 * Questa define fa si che venga saltata la procedura per
 * la determinazione della prima valle
 */
//#define EHSMODO2
/*
 * In questo modo vengono sommati tutti i vettori Cfft e solo
 * alla fine viene calcolata la prima valle e il picco sulla
 * distorsione armonica
 */

#define Fup 18000.0
#define Flow 80.0
#define PATCH 1
/****************************** end **************************/

/***************************** peaqb *************************/
#define LOGVARIABLE
/*
 * Pemrette di loggare alcune variabili frame per frame.
 */

#ifdef LOGVARIABLE
#define LOGALLFRAMES
#endif
/*
 * Logga tutti i frames scrivendoli uno dopo l'altro.
 */
//#define DATABOUND_BE
/*
 * Il data boundary in questo modo viene calcolato dall'inizio (Begin)
 * e dalla fine del file (End). --> [versione offline].
 */

//#define DATABOUND_ONE
/*
 * In questo modo viene preso solo il primo blocco di dati validi
 */

//#define LOUDMODO2
/* Loudness threshold:
 *
 * salta 50ms ogni volta che Ntotal > 0.1;
 * se è commentata salta solo una volta.
 */
/****************************** end **************************/

struct processing {
	double fftref[HANN/2];
	double ffttest[HANN/2];
	double ffteref[HANN/2];
	double fftetest[HANN/2];
	double fnoise[HANN/2];
	double pptest[BARK];
	double ppref[BARK];
	double ppnoise[BARK];
	double E2test[BARK];
	double E2ref[BARK];
	double Etest[BARK];
	double Eref[BARK];
	double Mref[BARK];
	double Modtest[BARK];
	double Modref[BARK];
};

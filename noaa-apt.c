#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <sys/stat.h>

#include <samplerate.h>
#include <sndfile.h>

/*                                            _   
 *  _ __   ___   __ _  __ _        __ _ _ __ | |_ 
 * | '_ \ / _ \ / _` |/ _` |_____ / _` | '_ \| __|
 * | | | | (_) | (_| | (_| |_____| (_| | |_) | |_ 
 * |_| |_|\___/ \__,_|\__,_|      \__,_| .__/ \__|
 *                                     |_|        
 * A simple, automatic program which takes an audio recording of one 
 * of the NOAA satellites, and does its best to automatically turn 
 * it into the best image it can, without hurting my brain too hard
 * or otherwise cluttering up the code with lots of useless stuff.
 *
 * Written by Mark VandeWettering
 */

/*--------------------------------------------------------------------*/
/* Hilbert Transformer... 					      */
/*--------------------------------------------------------------------*/

/* Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
   Command line: /www/usr/fisher/helpers/mkshape -h 63 -w -l */

#define NZEROS 62
#define GAIN   1.568293847e+00

static float xv[NZEROS+1];

static float xcoeffs[] =
  { +0.0025806452, +0.0000000000, +0.0030833182, +0.0000000000,
    +0.0043436817, +0.0000000000, +0.0064979527, +0.0000000000,
    +0.0096989225, +0.0000000000, +0.0141274550, +0.0000000000,
    +0.0200126097, +0.0000000000, +0.0276672484, +0.0000000000,
    +0.0375532412, +0.0000000000, +0.0504077050, +0.0000000000,
    +0.0675073918, +0.0000000000, +0.0912854169, +0.0000000000,
    +0.1270041052, +0.0000000000, +0.1884398887, +0.0000000000,
    +0.3263013526, +0.0000000000, +0.9976398888, +0.0000000000,
    -0.9976398888, -0.0000000000, -0.3263013526, -0.0000000000,
    -0.1884398887, -0.0000000000, -0.1270041052, -0.0000000000,
    -0.0912854169, -0.0000000000, -0.0675073918, -0.0000000000,
    -0.0504077050, -0.0000000000, -0.0375532412, -0.0000000000,
    -0.0276672484, -0.0000000000, -0.0200126097, -0.0000000000,
    -0.0141274550, -0.0000000000, -0.0096989225, -0.0000000000,
    -0.0064979527, -0.0000000000, -0.0043436817, -0.0000000000,
    -0.0030833182, -0.0000000000, -0.0025806452,
  };

#if 0
static void 
hilbertfilter()
{ 
    for (;;) { 
	float sum; int i;
        for (i = 0; i < NZEROS; i++) 
	    xv[i] = xv[i+1];
        xv[NZEROS] = next input value / GAIN;
        sum = 0.0;
        for (i = 0; i <= NZEROS; i++) sum += (xcoeffs[i] * xv[i]);
        next output value = sum;
    }
}
#endif

/*--------------------------------------------------------------------*/

#define SQR(x)		((x)*(x))
#define MIN(a,b)        ((a)<(b)?(a):(b))
#define MAX(a,b)        ((a)>(b)?(a):(b))
#define LERP(t, a, b)	((1-(t))*(a)+(t)*(b))

#define SR		(10*2080)
// #define SR		(9600)
#define NEW_WIDTH	(2080)

#define FILTER_LENGTH	(63)
float lpfilter[FILTER_LENGTH] = 
{
 8.86537E-04, 7.41036E-04, 3.35876E-04, -7.44832E-04, -2.27032E-03,
-3.60568E-03, -3.92266E-03, -2.60731E-03, 2.84200E-04, 3.81065E-03, 
6.35481E-03, 6.28896E-03, 2.86615E-03, -3.10756E-03, -9.23947E-03,
-1.23601E-02, -9.97371E-03, -1.74965E-03, 9.71073E-03, 1.94075E-02,
 2.18215E-02, 1.35359E-02, -4.49403E-03, -2.61691E-02, -4.16204E-02,
-4.05155E-02, -1.61846E-02, 3.10536E-02, 9.28420E-02, 1.54690E-01,
 2.00353E-01, 2.17165E-01, 2.00353E-01, 1.54690E-01, 9.28420E-02,
 3.10536E-02, -1.61846E-02, -4.05155E-02, -4.16204E-02, -2.61691E-02,
-4.49403E-03, 1.35359E-02, 2.18215E-02, 1.94075E-02, 9.71073E-03,
-1.74965E-03, -9.97371E-03, -1.23601E-02, -9.23947E-03, -3.10756E-03,
 2.86615E-03, 6.28896E-03, 6.35481E-03, 3.81065E-03, 2.84200E-04,
-2.60731E-03, -3.92266E-03, -3.60568E-03, -2.27032E-03, -7.44832E-04,
 3.35876E-04, 7.41036E-04, 8.86537E-04,
} ;

#define SYNC_FILTER_LENGTH 	(28)
float syncifilter[SYNC_FILTER_LENGTH] = { -1, -1, 1, 1, -1, -1, 1, 1,
-1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, } ;

float syncqfilter[SYNC_FILTER_LENGTH] = { 1, -1, -1, 1, 1, -1, -1, 1, 1,
-1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1 } ;


/* The waveforms of the syncA/syncB signals are both 39/4160 seconds long.
 * They are divided into segments which are of length T = 1.0 / 4160 seconds long.
 * The syncA signal has 7 pulses of 2T high, 2T low...
 * The syncB signal has 7 pulses of 3T high, 2T low...
 */

#define SYNC_LENGTH (39*SR/4160)

float syncA[SYNC_LENGTH] ;
float syncB[SYNC_LENGTH] ;

void
generate_sync()
{
    FILE *fp ;
    int i, p ;

    for (i=0; i<SYNC_LENGTH; i++)
        syncA[i] = syncB[i] = -1. ;

    for (p=0; p<7; p++) {
        for (i=(4+4*p)*SR/4160; i<(6+4*p)*SR/4160; i++)
            syncA[i] = 1.0 ;
    }

    for (p=0; p<7; p++) {
        for (i=(4+5*p)*SR/4160; i<(7+5*p)*SR/4160; i++)
            syncB[i] = 1.0 ;
    }

#if 0
    /* Print them out, so I can check them to make
     * sure I got them right.
     */
    fp = fopen("syncA.dat", "w") ;
    for (i=0; i<39*SR/4160; i++)
        fprintf(fp, "%f\n", syncA[i]) ;
    fclose(fp) ;

    fp = fopen("syncB.dat", "w") ;
    for (i=0; i<39*SR/4160; i++)
        fprintf(fp, "%f\n", syncB[i]) ;
    fclose(fp) ;
#endif

}

#define GAUSS_FILTER_WIDTH	(21)
#define GAUSS_SIGMA		(4)
float gaussfilter[GAUSS_FILTER_WIDTH] ;

void
initgaussfilter()
{
    int i ;
    float r2 ;
    for (i=0; i<GAUSS_FILTER_WIDTH; i++) {
	r2 = SQR((i-GAUSS_FILTER_WIDTH/2.) / GAUSS_SIGMA) ;
	gaussfilter[i] = exp(-r2) ;
    }
}

void
dofilter(float *filtertab, int ntaps, float *inp, int n, float *out)
{
    int i, j ;
    for (i=0; i<n-ntaps; i++) {
	out[i] = 0.0 ;
	for (j=0; j<ntaps; j++)
	    out[i] += filtertab[j] * inp[i+j] ;
    }
}

float
histfind(float f[], int n, float t) 
{
    int i ;
    float x ;

    for (i=0; i<n-1; i++) {
	if (f[i] <= t && f[i+1] > t) {
	    /* we want a value x such that...
	     * (1-x)*f[i] + x*f[i+1] = t
	     * solving for x, we get... 
	     * f[i] - x * f[i] + x*f[i+1] = t 
	     * x * (f[i+1] - f[i]) = t - f[i]
	     * x = (t-f[i]) / (f[i+1] - f[i]) ;
	     */
	    x = (t - f[i]) / (f[i+1] - f[i]) ;
	    return (float) (i+x) / n ;
	}
    }
    abort() ;
}

int
main(int argc, char *argv[])
{
    SNDFILE *sf ;
    SF_INFO sfinfo ;

    FILE *fp ;
    struct stat sbuf ;
    int width, height ;
    float *inp, *inpq, *inpi, *out, *out2, *outq, *outi, *synci, *syncq, imin, imax ;
    float *ip, *op ;
    int i, j, idx, acc ;
    int nsamples ;

    initgaussfilter() ;
    generate_sync() ;

    /* try opening our file of raw sound samples... */
    if ((fp = fopen(argv[1], "rb")) == NULL) {
	perror(argv[1]) ;
	exit(1) ;
    }

    if ((sf = sf_open(argv[1], SFM_READ, &sfinfo)) == NULL) {
	perror(argv[1]) ;
	exit(1) ;
    }

    fprintf(stderr, ":::: %s\n", argv[1]) ;
    fprintf(stderr, "     sample rate: %d hz \n", sfinfo.samplerate) ;
    fprintf(stderr, "     channels: %d\n", sfinfo.channels) ;
    fprintf(stderr, "     sections: %d\n", sfinfo.sections) ;
    fprintf(stderr, "     seekable: %s\n", sfinfo.seekable?"yes":"no") ; 
    fprintf(stderr, "     samples: %ld\n", (long) sfinfo.frames) ;

    if (sfinfo.channels != 1) {
        fprintf(stderr, "ERROR: sound file should only be one channel.\n") ;
        exit(1) ;
    }


    /* Without thinking too hard, resample... */

    float *tmpo =  (float *) calloc(sfinfo.frames, sizeof(float)) ;

    SRC_DATA srcdata ;
    srcdata.input_frames = sfinfo.frames ;
    srcdata.data_in = (float *) calloc(srcdata.input_frames, sizeof(float)) ;
    srcdata.output_frames = sfinfo.frames * SR / sfinfo.samplerate ;
    srcdata.data_out = (float *) calloc(srcdata.output_frames, sizeof(float)) ;
    srcdata.src_ratio = (float) SR / sfinfo.samplerate ;

    /* Read in all the samples... */
    sf_read_float(sf, srcdata.data_in, srcdata.input_frames) ;

    /* Do the resample */
    fprintf(stderr, "     Resampling to %d samples/second\n", SR) ;
    src_simple(&srcdata, SRC_SINC_FASTEST, 1) ;

    nsamples = srcdata.output_frames_gen ;
    fprintf(stderr, "     Input Frames Used: %ld\n", srcdata.input_frames_used) ;
    fprintf(stderr, "     Output Frames Generated: %ld\n", nsamples) ;



    /* There are two scanlines per second of recording */
    height = 2 * nsamples / SR ;

    /* Each scanline lasts 0.5 seconds */
    width = SR / 2 ;

    fprintf(stderr, "     Image is ~ %d lines, %d samples per line.\n", height, width) ;
    fprintf(stderr, "     check: %ld versus %d\n", nsamples, height*width) ;

    /* allocate space for input and output images */
    inp =  srcdata.data_out ;
    
    inpi = (float *) calloc(nsamples, sizeof(float)) ;
    inpq = (float *) calloc(nsamples, sizeof(float)) ;

    outi = (float *) calloc(nsamples, sizeof(float)) ;
    outq = (float *) calloc(nsamples, sizeof(float)) ;

    out  = (float *) calloc(nsamples, sizeof(float)) ;

    /* Use the Hilbert Transform to change inp to inpi/inpq */

    /* Mix with an oscillator... */
    float complex omega = 1. ;
    float complex domega = cexp(I*2.0*M_PI*2400/SR) ;

    for (i=0; i<nsamples; i++) {
	float complex x = inp[i] ;
	x *= omega ;
	inpi[i] = creal(x) ;
	inpq[i] = cimag(x) ;
	omega *= domega ;
    }

#if 0
    /* filter */
    for (i=0; i<nsamples; i++) {
        outi[i] = inpi[i] ;
        outq[i] = inpq[i] ;
    }
#else
    dofilter(lpfilter, FILTER_LENGTH, inpi, nsamples, outi) ;
    dofilter(lpfilter, FILTER_LENGTH, inpq, nsamples, outq) ;
#endif

    // Compute the amplitide, which is what we are interested in.
    for (i=0; i<nsamples; i++)
	out[i] = sqrt(SQR(outi[i])+SQR(outq[i])) ;

    fp = fopen("newsync.dat", "w") ;

    for (i=0; i<MIN(nsamples-SYNC_LENGTH,4*SR); i++) {
        float sumA = 0., sumB = 0., avg =0. ;
	
        for (j=0; j<SYNC_LENGTH; j++)
	    avg += out[i+j] ;
	avg /= SYNC_LENGTH ;

        for (j=0; j<SYNC_LENGTH; j++) {
            sumA += (out[i+j] - avg) * syncA[SYNC_LENGTH-1-j] ;
            sumB += (out[i+j] - avg) * syncB[SYNC_LENGTH-1-j] ;
        }
        fprintf(fp, "%d %f %f\n", i, sumA, sumB) ;
    }

    fclose(fp) ;

    /* create a different output buffer to hold 
     * the resized image... */
    out2  = (float *) calloc(height*NEW_WIDTH, sizeof(float)) ;
    synci = (float *) calloc(height*NEW_WIDTH, sizeof(float)) ;
    syncq = (float *) calloc(height*NEW_WIDTH, sizeof(float)) ;

    /* Two different ways of resampling..
     * They mostly are exactly the same..
     */
#if 0
    /* resample the _bad_ way... */
    for (i=0; i<height*NEW_WIDTH; i++) {
	float x = (float) i * width / NEW_WIDTH ;
	int ix = floor(x) ;
	float fx = x - ix ;
	out2[i] = LERP(fx, out[ix], out[ix+1]) ;
    }
#else
    /* resample a bad, but better way... */
    float sum = 0.0, w = 0.0 ;
    float r = (float) NEW_WIDTH / width ;
    for (i=0, j=0; i<nsamples && j<height*NEW_WIDTH; i++) {
        if (w + r >= 1.0) {
            /* add in the fractional bit */
            sum += (1.0 - w) * out[i] ;
            /* output the appropriate amount */
            out2[j++] = sum ;
            w = w + r - 1.0 ;
            sum = w * out[i] ;
        } else {
            sum += r * out[i] ;
            w +=r ;
        }
    }
#endif
        
#if 0
    /* now, try generating the sync.. */
    dofilter(syncifilter, SYNC_FILTER_LENGTH, out2, height*NEW_WIDTH, synci) ;
    dofilter(syncqfilter, SYNC_FILTER_LENGTH, out2, height*NEW_WIDTH, syncq) ;

    /* Find the amplitude of the sync signal */
    for (i=0; i<height*NEW_WIDTH; i++)
	syncq[i] = SQR(synci[i]) + SQR(syncq[i]) ;
    /* Low pass filter it... */
    dofilter(gaussfilter, GAUSS_FILTER_WIDTH, syncq, height * NEW_WIDTH, synci) ;

    fp = fopen("sync.dat", "w") ;
    for (i=0; i<height*NEW_WIDTH; i++) 
        fprintf(fp, "%f\n", synci[i]) ;
    fclose(fp) ;
#endif

#define HISTOGRAM_BINS 2000
    float hist[HISTOGRAM_BINS] ;
    for (i=0; i<HISTOGRAM_BINS; i++)
	hist[i] = 0 ;
    for (i=0; i<height * NEW_WIDTH; i++) {
	hist[(int)(out2[i]*(HISTOGRAM_BINS-1))] += 1.0 ;
    }

    for (i=0; i<HISTOGRAM_BINS; i++) 
	hist[i] /= height * NEW_WIDTH ;
    for (i=1; i<HISTOGRAM_BINS; i++) 
	hist[i] += hist[i-1] ;

    imin = histfind(hist, HISTOGRAM_BINS, 0.02) ;
    imax = histfind(hist, HISTOGRAM_BINS, 0.98) ;
    fprintf(stderr, ":::: histogram pass, %f %f...\n", imin, imax) ;

    for (i=0; i<height * NEW_WIDTH; i++) {
	if (out2[i] < imin)
	    out2[i] = 0.0 ;
	else if (out2[i] > imax) 
	    out2[i] = 1.0 ;
	else 
	    out2[i] = (out2[i] - imin) / (imax - imin) ;
    }

    /* float a = 5.12459e-05, b = 0.0534101, c = 1257.64 ; */
    /* float a = 6.0924e-05, b = 0.0343549, c = 1893.04 ; */
    /* float a = 4.54852e-05, b = 0.0976785, c = 1781.21 ; */
    float a = 4.32867e-05, b = 0.00770067, c = 2008.28 ;

    c += GAUSS_FILTER_WIDTH / 2. ;

    /* okay, try to do the rectified image... */
    fp = popen("cjpeg -quality 95 -gray -progressive", "w") ;
    fprintf(fp, "P5\n%d %d\n%d\n", NEW_WIDTH, height-1, 255) ;
    for (i=0, op = out2; i<height-1; i++, op+= NEW_WIDTH) {
	float sx = (c + i * (b + i * a)) ;
	float ex = (c + (i+1) * (b + (i+1) * a)) ;
	float lx = (ex-sx) + NEW_WIDTH ;
	/* now, we start the scanline at sx, and go to the next 
	 * scanline at ex...
	 */
	for (j=0; j<NEW_WIDTH; j++) {
            /* BAD RESAMPLE */
            float x = sx + lx * j / (NEW_WIDTH) ;
            int ix = (int) floor(x) ;
            float fx = x - ix ;
            float val = LERP(fx, op[ix], op[ix+1]) ;
            fputc((int) (255.*val), fp) ;
	}
    }

    pclose(fp) ;

    return 0 ;
}

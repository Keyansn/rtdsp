/*************************************************************************************
      DEPARTMENT OF ELECTRICAL AND ELECTRONIC ENGINEERING
      IMPERIAL COLLEGE LONDON 

       EE 3.19: Real Time Digital Signal Processing
      Dr Paul Mitcheson and Daniel Harvey

        PROJECT: Frame Processing

             ********* ENHANCE. C **********
Shell for speech enhancement 

  Demonstrates overlap-add frame processing (interrupt driven) on the DSK. 

 *************************************************************************************
              By Danny Harvey: 21 July 2006
Updated for use on CCS v4 Sept 2010
 ************************************************************************************/
/*
 * You should modify the code so that a speech enhancement project is built 
 *  on top of this template.
 */
/**************************** Pre-processor statements ******************************/
//  library required when using calloc
#include <stdlib.h>
//  Included so program can make use of DSP/BIOS configuration tool.  
#include "dsp_bios_cfg.h"

/* The file dsk6713.h must be included in every program that uses the BSL.  This 
   example also includes dsk6713_aic23.h because it uses the 
   AIC23 codec module (audio interface). */
#include "dsk6713.h"
#include "dsk6713_aic23.h"

// math library (trig functions)
#include <math.h>
#include <float.h>

/* Some functions to help with Complex algebra and FFT. */
#include "cmplx.h"      
#include "fft_functions.h"  

// Some functions to help with writing/reading the audio ports when using interrupts.
#include <helper_functions_ISR.h>

#define WINCONST 0.85185 /* 0.46/0.54 for Hamming window */
#define FSAMP 8000.0 /* sample frequency, ensure this matches Config for AIC */
#define FFTLEN 256 /* fft length = frame length 256/8000 = 32 ms*/
#define NFREQ (1+FFTLEN/2) /* number of frequency bins from a real FFT */
#define OVERSAMP 4 /* oversampling ratio (2 or 4) */  
#define FRAMEINC (FFTLEN/OVERSAMP) /* Frame increment */
#define CIRCBUF (FFTLEN+FRAMEINC) /* length of I/O buffers */

#define OUTGAIN 16000.0 /* Output gain for DAC */
#define INGAIN  (1.0/16000.0) /* Input gain for ADC  */
// PI defined here for use in your code 
#define PI 3.141592653589793
#define TFRAME FRAMEINC/FSAMP
#define time_counter 312

      

/******************************* Global declarations ********************************/

/* Audio port configuration settings: these values set registers in the AIC23 audio 
   interface to configure it. See TI doc SLWS106D 3-3 to 3-10 for more info. */
DSK6713_AIC23_Config Config = { \
/**********************************************************************/
/*   REGISTER            FUNCTION      SETTINGS         */ 
/**********************************************************************/\
    0x0017,  /* 0 LEFTINVOL  Left line input channel volume  0dB                   */\
    0x0017,  /* 1 RIGHTINVOL Right line input channel volume 0dB                   */\
    0x01f9,  /* 2 LEFTHPVOL  Left channel headphone volume   0dB                   */\
    0x01f9,  /* 3 RIGHTHPVOL Right channel headphone volume  0dB                   */\
    0x0011,  /* 4 ANAPATH    Analog audio path control       DAC on, Mic boost 20dB*/\
    0x0000,  /* 5 DIGPATH    Digital audio path control      All Filters off       */\
    0x0000,  /* 6 DPOWERDOWN Power down control              All Hardware on       */\
    0x0043,  /* 7 DIGIF      Digital audio interface format  16 bit                */\
    0x008d,  /* 8 SAMPLERATE Sample rate control        8 KHZ-ensure matches FSAMP */\
    0x0001   /* 9 DIGACT     Digital interface activation    On                    */\
/**********************************************************************/
};

// Codec handle:- a variable used to identify audio interface  
DSK6713_AIC23_CodecHandle H_Codec;

float *inbuffer, *outbuffer;   /* Input/output circular buffers */
float *inframe, *outframe;          /* Input and output frames */
float *inwin, *outwin;              /* Input and output windows */
float ingain, outgain; /* ADC and DAC gains */ 
float cpufrac; /* Fraction of CPU time used */
volatile int io_ptr=0;              /* Input/ouput pointer for circular buffers */
volatile int frame_ptr=0;
int count = 0;           /* Frame pointer */
complex *X;
float *absX;
float *M1;
float *M2;
float *M3;
float *M4;
float *temp;
float *noise;
float *noise_prev;
float *g;
complex *Y;
float K_factor = 0.80;
float P[FFTLEN];
float P_prev[FFTLEN];
float absXsq;
float noisesq;
float Psq;
float SNR=1;
float alphascale=1;
float alpha= 4;
float lambda= 0.015;


/* 		absX[k] = min(FLT_MAX, cabs(X[i]));
		absXsq = min(FLT_MAX, absX[k]*absX[k]);
		noisesq = min(FLT_MAX, noise[i]*noise[i]);
		Psq = min(FLT_MAX, P[i]*P[i]);	 */

/**********************************Switches**********************************/
int Optimisation1=1;
int Optimisation2=0;
int Optimisation3=1;
int Optimisation4=4;		
int Optimisation5=0;
int Optimisation6=1;
//int Optimisation7=0;
//int Optimisation8=0;
int unfiltered =0;
/**********************************************************************/

 /******************************* Function prototypes *******************************/
void init_hardware(void);     /* Initialize codec */ 
void init_HWI(void);            /* Initialize hardware interrupts */
void ISR_AIC(void);             /* Interrupt service routine for codec */
void process_frame(void);
float min (float a, float b);
float max (float a, float b);     /* Frame processing routine */
float lowpass(float x, float prev[], int k);
           
/********************************** Main routine ************************************/
void main()
{      

  int k=0; // used in various for loops
  
/*  Initialize and zero fill arrays */  



inbuffer = (float *) calloc(CIRCBUF, sizeof(float)); /* Input array */
outbuffer = (float *) calloc(CIRCBUF, sizeof(float)); /* Output array */
inframe = (float *) calloc(FFTLEN, sizeof(float)); /* Array for processing*/
outframe = (float *) calloc(FFTLEN, sizeof(float)); /* Array for processing*/
inwin = (float *) calloc(FFTLEN, sizeof(float)); /* Input window */
outwin = (float *) calloc(FFTLEN, sizeof(float)); /* Output window */
X = (complex *) calloc(FFTLEN, sizeof(complex));
absX = (float *) calloc(FFTLEN, sizeof(float));
M1 = (float *) calloc(FFTLEN, sizeof(float));
M2 = (float *) calloc(FFTLEN, sizeof(float));
M3 = (float *) calloc(FFTLEN, sizeof(float));
M4 = (float *) calloc(FFTLEN, sizeof(float));
noise = (float *) calloc(FFTLEN, sizeof(float));
noise_prev = (float *) calloc(FFTLEN, sizeof(float));
g = (float *) calloc(FFTLEN, sizeof(float));
Y = (complex *) calloc(FFTLEN, sizeof(complex));
temp = (float *) calloc(FFTLEN, sizeof(float));




/* initialize board and the audio port */
  init_hardware();
  
  /* initialize hardware interrupts */
  init_HWI();    
  
/* initialize algorithm constants */  
                       
  for (k=0;k<FFTLEN;k++)
{                           
inwin[k] = sqrt((1.0-WINCONST*cos(PI*(2*k+1)/FFTLEN))/OVERSAMP);
outwin[k] = inwin[k]; 
P[k] = 0;
P_prev[k] =0;
} 
  ingain=INGAIN;
  outgain=OUTGAIN;        


//	for(k = 0; k < FFTLEN; k++)
//	{
//		 nmb[k] = m1[k] = m2[k] = m3[k] = m4[k] = FLT_MAX;
//		 nmb_prev[k] = P_prev[k] = 0;					
//	}

  
  /* main loop, wait for interrupt */  
  while(1) process_frame();
}
    
/********************************** init_hardware() *********************************/  
void init_hardware()
{
    // Initialize the board support library, must be called first 
    DSK6713_init();
    
    // Start the AIC23 codec using the settings defined above in config 
    H_Codec = DSK6713_AIC23_openCodec(0, &Config);

/* Function below sets the number of bits in word used by MSBSP (serial port) for 
receives from AIC23 (audio port). We are using a 32 bit packet containing two 
16 bit numbers hence 32BIT is set for  receive */
MCBSP_FSETS(RCR1, RWDLEN1, 32BIT); 

/* Configures interrupt to activate on each consecutive available 32 bits 
from Audio port hence an interrupt is generated for each L & R sample pair */ 
MCBSP_FSETS(SPCR1, RINTM, FRM);

/* These commands do the same thing as above but applied to data transfers to the 
audio port */
MCBSP_FSETS(XCR1, XWDLEN1, 32BIT); 
MCBSP_FSETS(SPCR1, XINTM, FRM); 

}
/********************************** init_HWI() **************************************/ 
void init_HWI(void)
{
IRQ_globalDisable(); // Globally disables interrupts
IRQ_nmiEnable(); // Enables the NMI interrupt (used by the debugger)
IRQ_map(IRQ_EVT_RINT1,4); // Maps an event to a physical interrupt
IRQ_enable(IRQ_EVT_RINT1); // Enables the event
IRQ_globalEnable(); // Globally enables interrupts

}
        
/******************************** process_frame() ***********************************/  
void process_frame(void)
{
int k, m; 
int io_ptr0;   

/* work out fraction of available CPU time used by algorithm */    
cpufrac = ((float) (io_ptr & (FRAMEINC - 1)))/FRAMEINC;  
/* wait until io_ptr is at the start of the current frame */
while((io_ptr/FRAMEINC) != frame_ptr); 
/* then increment the framecount (wrapping if required) */ 
if (++frame_ptr >= (CIRCBUF/FRAMEINC)) frame_ptr=0;
  
  /* save a pointer to the position in the I/O buffers (inbuffer/outbuffer) where the 
  data should be read (inbuffer) and saved (outbuffer) for the purpose of processing */
  io_ptr0=frame_ptr * FRAMEINC;
/* copy input data from inbuffer into inframe (starting from the pointer position) */ 
 
m=io_ptr0;
for (k=0;k<FFTLEN;k++){                           
	inframe[k] = inbuffer[m] * inwin[k]; 
	if (++m >= CIRCBUF) m=0; /* wrap if required */ }
	 
/************************* DO PROCESSING OF FRAME  HERE **************************/
/* please add your code, at the moment the code simply copies the input to the 
ouptut with no processing */  


/***************** Applying fft ****************/ 
for( k = 0; k < FFTLEN; k++){                           
	X[k] = cmplx(inframe[k],0);} /* copy input straight into output */ 
	
	      
fft(FFTLEN,X);	  
for(k =0 ; k<FFTLEN; k++){
absX[k] = cabs(X[k]);
}

for( k = 0; k < FFTLEN; k++){ 
	
	P[k] = (1-K_factor)*absX[k] + (K_factor)*P_prev[k];
	P_prev[k]=P[k];
}
	  
if(count >= time_counter/4){
	count = 0;
	temp = M4;
	M4 = M3;
	M3 = M2;
	M2 = M1;
	M1 = temp;

	for ( k = 0 ; k < FFTLEN ; k++)
	{
		M1[k] = P[k]; 
		}
	}	  
count++;  
	  
	 
/***********************************************/
 
 
/*
for( k = 0; k < FFTLEN; k++){                           
	absX[k] = cabs(X[k]);
	
	if(Optimisation1)
	{
		// Square absX to put in power domain if o2 picked
		if(Optimisation2)
		{
			absX[k] *= absX[k];
		}
		
		// Perform low pass filter	
		P[k] = (1-K_factor)*absX[k] + (K_factor)*P_prev[k];
				
		// Store value for next lpf equation
		P_prev[k] = P[k];
	}
		
	if(Optimisation2)
	{
		P[k] = sqrt(P[k]);
	}
		
	else
		P[k] = absX[k];
		
		M1[k] = min(M1[k],P[k]);
	}
	*/
/***********************************************/
	
for( k = 0; k < FFTLEN; k++){ 
	M1[k] = min(M1[k],P[k]);
	noise[k] = alpha*(min(min(M1[k],M2[k]),min(M3[k],M4[k])));
	

	
	if(Optimisation6 && k<20){
		absX[k] = cabs(X[k]);
		SNR = (absX[k]/noise[k] -1)*(absX[k]/noise[k] -1);
		alphascale = min(FLT_MAX, -log(SNR));
		noise[k] = min(FLT_MAX, (noise[k]*alphascale));
	}
	
	
	if(Optimisation3){
		noise[k] = (1-K_factor)*noise[k] + (K_factor)*noise_prev[k];
		
		noise_prev[k] = noise[k];	
	}
		
		
		/***********************	OPTIMISATION 6 goes here	************************/
	//}	
//}	
	

/***********************************************/
	

//for( k = 0; k < FFTLEN; k++){

	/***********************	OPTIMISATION 5 goes here	************************/
	//if(Optimisation5){
		
		//absX[k] = min(cabs(X[k]), FLT_MAX);
		//absXsq = min(absX[k]*absX[k], FLT_MAX);
		//noisesq = min(noise[k]*noise[k], FLT_MAX);
		//Psq = min(P[k]*P[k], FLT_MAX);		
		
	//}
 	/*******************************************************************************/
	/*
	switch(Optimisation4){
		
		//no change to default
		case 0:
			if (Optimisation5)
				g = max(lambda, sqrt(1-(noisesq/absXsq)));
			else
				g = max(lambda, 1-(noise[k]/absX[k]));
			
			break;
	
		//(N/X,N/X)	
		case 1:
			if (Optimisation5)
				g = max(lambda*sqrt(noisesq/absXsq), sqrt(1-(noisesq/absXsq)));
			else
				g = max(lambda*noise[k]/absX[k], 1-(noise[k]/absX[k]));
			
			break;
	
		//(P/X,N/X)		
		case 2:
			if (Optimisation5)
			g = max(lambda*sqrt(Psq/absXsq), sqrt(1-(noisesq/absXsq)));
			else
				g = max(lambda*P[k]/absX[k], 1-(noise[k]/absX[k]));
			
			break;
	
		//(N/P,N/P)		
		case 3:
			if (Optimisation5)
				g = max(lambda*sqrt(noisesq/Psq), sqrt(1-(noisesq/Psq)));
			else
				g = max(lambda*noise[k]/P[k], 1-(noise[k]/P[k]));			
			
			break;
	
		//(L,N/P)	
		case 4:
			if (Optimisation5)
				g = max(lambda, sqrt(1-(noisesq/Psq)));
			else
				g[k] = max(lambda, 1-(noise[k]/P[k]));
				
			break;
	}*/

g[k] = max(lambda, 1-(noise[k]/P[k]));

Y[k] = rmul(g[k],X[k]);
		}	//end of loop

						
ifft(FFTLEN,Y);

for( k = 0; k < FFTLEN; k++){
	outframe[k] = Y[k].r;}	//check if it needs to be real or .r $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




 
/********************************************************************************/
    /* multiply outframe by output window and overlap-add into output buffer */  
                           
m=io_ptr0;
    
for (k=0;k<(FFTLEN-FRAMEINC);k++) { 
    /* this loop adds into outbuffer */                       
	outbuffer[m] = outbuffer[m]+outframe[k]*outwin[k];   
	if (++m >= CIRCBUF) m=0; /* wrap if required */
}         
    
for (;k<FFTLEN;k++) {                           
outbuffer[m] = outframe[k]*outwin[k];   /* this loop over-writes outbuffer */        
m++;
} 


}        
/*************************** INTERRUPT SERVICE ROUTINE  *****************************/

// Map this to the appropriate interrupt in the CDB file
   
void ISR_AIC(void)
{       
short sample;
/* Read and write the ADC and DAC using inbuffer and outbuffer */
sample = mono_read_16Bit();
inbuffer[io_ptr] = ((float)sample)*ingain;
/* write new output data */
mono_write_16Bit((int)(outbuffer[io_ptr]*outgain)); 
/* update io_ptr and check for buffer wraparound */    
if (++io_ptr >= CIRCBUF) io_ptr=0;
}

/******************************Custom functions******************************/

/*Minimum of two floats*/
float min (float a, float b){
	if	( a < b ){return a;}
	else	{return b;}}
	
/*Minimum of two floats*/
float max (float a, float b){
	if 	( a < b ){return b;}
	else	{return a;}}
	
///*Low Pass Filter*/
float lowpass(float x, float prev[], int k){
	return (1-K_factor)*x + (K_factor)*prev[k];}
	

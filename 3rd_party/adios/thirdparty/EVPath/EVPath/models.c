// TODO: include checks to make sure fftw (and/or imkl in the future) are installed before allowing calls to it

#include "models.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <fftw3.h>
#include <stdlib.h>
#include <string.h>

// TODO: add more possible fits, like ARIMA, or LSF, or whatever your heart desires
// also, inject these constants into EVPath

// possible models
const int MODEL_NO_FIT = 0;
const int MODEL_FFT = 1;
// const int MODEL_ARIMA = 2;
// ...

// set to 1 after model framework is added to EVPath
static int initialized = 0;

FMField model_field_list[] = {
	{"type", "integer", sizeof(int), FMOffset(model_ptr, type)},
	{"num_terms", "integer", sizeof(int), FMOffset(model_ptr, num_terms)},
	{"terms", "double[num_terms]", sizeof(double), FMOffset(model_ptr, terms)},
	{"err_term", "double", sizeof(double), FMOffset(model_ptr, err_term)},
	{NULL, NULL, 0, 0}
};

FMStructDescRec model_format_list[] = {
	{"model_rec", model_field_list, sizeof(model), NULL},
	{NULL, NULL}
};

// simple sin wave generator
// change hardcoded values for different sin wave
double sin_wav(double sample_rate) {
	double frequency = 1023.63; // hz -- cycles/second
	double amplitude = 5.4; // percentage peak-to-peak
	double phase = -M_PI/8; // amount to shift the sin wave
	double confidence = 1.; // defines the confidence interval
	static int i = 0;
	srand(5);
	double next;
	next = amplitude * (sin(2*M_PI*(i++) / sample_rate * frequency + phase)); // sin wave of certain frequency and amplitude at sample rate
	next += next * (1.-confidence) * ((double)rand()/RAND_MAX-.5) * 2; // confidence interval
	return next;
}


// TODO: not sure if this is the best way to organize
// you can take the fft, returning the entire thing, or fit a model to the fft, only returning a few important variables
// not sure how extensible this concept is to other models, or even if this is useful
// also, you can play around with what variables get passed in, and what gets passed out, eg array size

// takes the fft on the data set
// stores results in model struct
// type: MODEL_FFT
// num_terms: (size of data set * 16 (padded zeros) / 2) + 1
// terms: the results of the fft
// err_term: null
// notes:
//        modifies the length of data to have 16X padded values for increased resolution on FFT
void EVmodel_take_fft(double* data, model_ptr m) {
	int num_samples = (sizeof(data)/sizeof(data[0]));
	// scale the sample size up by 16
	int buffer_size = num_samples*16;
	double* buffered_data = fftw_malloc(sizeof(double)*buffer_size);
	memcpy(buffered_data, data, sizeof(double)*num_samples);
	// output of FFT is complex
	// range of answers is buffer_size/2
	m->num_terms = buffer_size/2+1;
	fftw_complex* fft = fftw_malloc(sizeof(fftw_complex)*m->num_terms);
	// prepare the fftw plan, giving it all necessary parameters
	fftw_plan plan_forward = fftw_plan_dft_r2c_1d(buffer_size, buffered_data, fft, FFTW_ESTIMATE);
	// execute the plan once
	fftw_execute(plan_forward);
	// and then discard the plan
	fftw_destroy_plan(plan_forward);
	fftw_free(buffered_data);
	
	m->terms = (double *)malloc(sizeof(double)*m->num_terms*2);
	// translate fftw_complex into scalar array of doubles
	for(int i=0; i<m->num_terms; ++i) {
		m->terms[i*2] = fft[i][0];
		m->terms[i*2+1] = fft[i][1];
	}
	
	fftw_free(fft);
}

// takes the fft on the data set then finds a model for that fft
// stores results in model struct
// type: MODEL_FFT
// num_terms: 7
// terms: [0] range of FFT
//        [1] interval size for FFT
//        [2] largest bin number
//        [3] largest contributing frequency
//        [4] magnitude of the frequency
//        [5][6] phase of the frequency
// err_term: null
void EVmodel_fit_fft(double* data, model_ptr m, double sample_rate) {
	int num_samples = sizeof(data)/sizeof(data[0]);
	EVmodel_take_fft(data, m);
	int buffer_size = (m->num_terms-1)*2;
	double* fft = m->terms;
	
	// find the largest contributing signal and its amplitude
	int largest_idx = 0;
	double largest_amp = 0.;
	for(int i=0; i<m->num_terms; i+=2) {
		double i_amp = sqrt(fft[i]*fft[i] + fft[i+1]*fft[i+1]);
		if(i_amp > largest_amp) {
			largest_idx = i;
			largest_amp = i_amp;
		}
	}
	free(m->terms);
	
	// the all-real input causes a symmetrical low and high peak in the fft
	// both the low and high contribute to amplitude, so multiply high by 2 for same effect
	// also divide by num of data points for scaling effect
	largest_amp *= 2./num_samples;
	
	double FFT_interval = (double)sample_rate/buffer_size;
	double FFT_range = sample_rate/2;
	double largest_interval = FFT_interval*largest_idx;
	
	m->terms = (double *)malloc(sizeof(double)*7);
	// range of FFT
	m->terms[0] = FFT_range;
	// interval size of FFT
	m->terms[1] = FFT_interval;
	// largest contributing bin number
	m->terms[2] = largest_idx;
	// largest contributing frequency
	m->terms[3] = largest_interval;
	// magnitude of this frequency
	m->terms[4] = largest_amp;
	// phase of this frequency
	m->terms[5] = fft[largest_idx];
	m->terms[6] = fft[largest_idx+1];
	// you can combine the phase with atan2 to get an appropriate single value phase
}


// TODO: the entire automatic fitting thing
// this would require diagnostics on the data to determine what kind of model might fit it best
// potentially possible, but probably not very scalable. it may not be wise to pursue this right away...

// attempts to find a good model to fit data
// currently defaults to fitting FFT unconditionally
void EVmodel_fit(double* data, model_ptr m, void *client_data) {
	m->type = MODEL_FFT;
	
	switch(m->type) {
		case MODEL_FFT:
		double sample_rate = *((double *)client_data);
		EVmodel_fit_fft(data, m, sample_rate);
		break;
		default: break;
	}
}

// TODO: the (fake) enum injection does not work. not sure why, although a true enum is incoming anyway so ignore this problem for now
// additionally, probably adding more direct access to the modeling functions in addition to generic EVmodel_fit

// initialize the EVmodel framework
// adds functionality to EVPath for modeling data
int EVmodel_init(CManager cm) {
	if(initialized == 0){
		static char* extern_string = "\
			/* fake enum */ \n\
			/* int MODEL_NO_FIT = 0; \n\
			int MODEL_FFT = 1; */\n\
			/* actual functions */ \n\
			double sin_wav(double sample_rate); \n\
			void EVmodel_fit(double* data, model_rec m);";
		
		static cod_extern_entry model_funcs[] = {
			{"EVmodel_fit", NULL}, // 0
			{"sin_wav", NULL},     // 1
			{NULL, NULL}
		};
		model_funcs[0].extern_value = (void *) (long) EVmodel_fit;
		model_funcs[1].extern_value = (void *) (long) sin_wav;
		
		FMStructDescList formats_to_add[] = {model_format_list, NULL};
		EVadd_standard_structs(cm, formats_to_add);
		EVadd_standard_routines(cm, extern_string, model_funcs);
		printf("added model components now...\n");
		initialized = 1;
	}
	return initialized;
}
#include <stdio.h>
#include <stdbool.h>
#include <sndfile.h>
#include <math.h>
#include <string.h>
#include <complex.h> //Complex must be included before fftw3 for fftw3 to use the c standard complex type
#include <fftw3.h>
#include "filter.h"
#include <volk/volk.h>

void Hilbert(float *in, lv_32fc_t *out, int num_elements, int N)
{
	double complex *InputArray = malloc(sizeof(double complex) * num_elements);
	for (int i = 0; i < num_elements; i++)
	{
		InputArray[i] = in[i] + 0 * I;
	}

	double complex *ProcessArray = malloc(sizeof(double complex) * num_elements);

	fftw_plan plan_forward = fftw_plan_dft_1d(N, ProcessArray, ProcessArray, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_backward = fftw_plan_dft_1d(N, ProcessArray, ProcessArray, FFTW_BACKWARD, FFTW_ESTIMATE);

	for (int ix = 0; ix < num_elements; ix += N)
	{
		//printf("I: %d num_elements: %d\n", i, num_elements);
		if (ix + N > num_elements)
			break;
		memcpy(ProcessArray, InputArray + ix, sizeof(double complex) * N);
		fftw_execute(plan_forward);
		int hN = N >> 1; // half of the length (N/2)
		int numRem = hN; // the number of remaining elements
		for (int i = 1; i < hN; ++i)
		{
			ProcessArray[i] = (creal(ProcessArray[i]) * 2) + (cimag(ProcessArray[i]) * 2) * I;
		}
		if (N % 2 == 0)
			numRem--;
		else if (N > 1)
		{
			ProcessArray[hN] = (creal(ProcessArray[hN]) * 2) + (cimag(ProcessArray[hN]) * 2) * I;
		}
		memset(&ProcessArray[hN + 1], 0, numRem * sizeof(double complex));
		fftw_execute(plan_backward);
		for (int i = 0; i < N; ++i)
		{
			ProcessArray[i] = (creal(ProcessArray[i]) / N) + (cimag(ProcessArray[i]) / N) * I;
		}
		for (int i = 0; i < N; i++)
		{
			//out[i] = (complex float) ((float)creal(ProcessArray[i]) / N) + ((float)cimag(ProcessArray[i]) / N)*I;
			out[ix + i] = creal(ProcessArray[i]) + cimag(ProcessArray[i]) * I;
		}
	}

	fftw_destroy_plan(plan_forward);
	fftw_destroy_plan(plan_backward);
	free(ProcessArray);
	free(InputArray);
}

int main()
{
	char *inFileName;
	SNDFILE *inFile;
	SF_INFO inFileInfo;

	inFileName = "out.wav";
	inFile = sf_open(inFileName, SFM_READ, &inFileInfo);
	long int samp_count = inFileInfo.frames;
	int samp_rate = inFileInfo.samplerate;
	float *samples = calloc(samp_count, sizeof(float));
	sf_readf_float(inFile, samples, samp_count); //??? fix this pointer thing pls why it produce warning error idk
	sf_close(inFile);
	printf("Be careful! This code will break with WAV files at 48000 Hz that are longer than about 12h. Actually it might be 6h but idk\n");
	printf("Sample Rate = %d Hz\n", samp_rate);
	printf("Sample Count = %ld\n", samp_count);
	sf_close(inFile);

	const int MixFrequency = 10; //user
	//const int Bandwidth = 1000;	   //user
	const int speedDivider = 2000;   //user

	float Bandwidth = (float)(samp_rate / 2) / speedDivider;
	if(speedDivider > samp_rate)
	{
		printf("Minimum Bandwidth is 0.5Hz!\n");
		return 0;
	}

	printf("input Bandwidth: %fHz\n", Bandwidth);




	//Apply filters to I and Q	
	printf("Applying filters\n");
	BWBandPass *filterR = create_bw_band_pass_filter(50, samp_rate, MixFrequency, MixFrequency + Bandwidth);
	for (int i = 0; i < samp_count; i++)
	{
		//samples[i] = bw_band_pass(filterR, samples[i]);
	}
	free_bw_band_pass(filterR);
	



	//Convert real array into complex, set Imaginary to zero
	printf("Converting to complex\n");
	lv_32fc_t *inputComplex = (lv_32fc_t *)volk_malloc(sizeof(lv_32fc_t) * samp_count, volk_get_alignment());
	//Hilbert(samples, inputComplex, samp_count, samp_count);
	for(int i = 0; i < samp_count; i++)
	{
		inputComplex[i] = samples[i] + 0 * I;
	}
	free(samples);

	printf("Mixing up\n");
	float sinAngle = 2.0 * 3.14159265359 * MixFrequency / samp_rate;
	lv_32fc_t phase_increment = lv_cmake(cos(sinAngle), sin(sinAngle));
	lv_32fc_t phase = lv_cmake(1.f, 0.0f);
	//volk_32fc_s32fc_x2_rotator_32fc(inputComplex, inputComplex, phase_increment, &phase, samp_count);

	

	char *outFileName = "demod_out_complex.wav";
	SNDFILE *outFile;
	SF_INFO outFileInfo = inFileInfo;
	outFileInfo.channels = 2;
	outFile = sf_open(outFileName, SFM_WRITE, &outFileInfo);
	sf_writef_float(outFile, (float *)inputComplex, samp_count);
	sf_close(outFile);

	printf("Converting back to real\n");
	float *outputReal = calloc(samp_count, sizeof(float));
	for (int i = 0; i < samp_count; i++)
	{
		outputReal[i] = creal(inputComplex[i]);
	}
	volk_free(inputComplex);

	
	printf("Writing to file\n");
	char *outFileName2 = "demod_out.wav";
	SNDFILE *outFile2;
	SF_INFO outFileInfo2 = inFileInfo;
	outFile2 = sf_open(outFileName2, SFM_WRITE, &outFileInfo2);
	sf_writef_float(outFile2, outputReal, samp_count);
	sf_close(outFile2);
	

	return 0;
}

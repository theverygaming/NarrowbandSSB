#include <stdio.h>
#include <stdbool.h>
#include <sndfile.h>
#include <math.h>
#include <string.h>
#include <complex>
#include <fftw3.h>
#include "filter.h"
#include <volk/volk.h>
#include "wavwriter.h"


void fft_bandpass_filter(float *in, long int samplecount, int samplerate, int fftN, float lowCutoff, float highCutoff)
{
	double *ProcessR = (double*)malloc(sizeof(double)*fftN);
	std::complex<double> *ProcessC = (std::complex<double>*)malloc(sizeof(std::complex<double>)*(fftN / 2 + 1));
	fftw_plan plan_forward = fftw_plan_dft_r2c_1d(fftN, ProcessR, (fftw_complex*)ProcessC, FFTW_ESTIMATE);
	fftw_plan plan_backward = fftw_plan_dft_c2r_1d(fftN, (fftw_complex*)ProcessC, ProcessR, FFTW_ESTIMATE);

	for(long int ix = 0; ix < samplecount; ix += fftN)
	{
		if (ix + fftN > samplecount)
			break;
		for(int i = 0; i < fftN; i++)
		{
			ProcessR[i] = (double)in[ix+i];
		}

		fftw_execute(plan_forward);
		for(int i = 0; i < (fftN / 2 + 1); i++)
		{
			double real = ProcessC[i].real();
			double imag = ProcessC[i].imag();
			float binFrequency = (float)(i+1) * (float)samplerate / (float)fftN;
			if(binFrequency < lowCutoff || binFrequency > highCutoff)
			{
				real = 0;
				imag = 0;
			}
			ProcessC[i] = std::complex<double>(real,imag);
			
		}
		fftw_execute(plan_backward);

		for(int i = 0; i < fftN; i++)
		{
			in[ix+i] = (float)ProcessR[i] / fftN;
		}
	}

	fftw_destroy_plan(plan_forward);
	fftw_destroy_plan(plan_backward);
	free(ProcessR);
	free(ProcessC);
}


void Hilbert(float *in, lv_32fc_t *out, int num_elements, int N)
{
	std::complex<double> *InputArray = (std::complex<double>*)malloc(sizeof(std::complex<double>) * num_elements);
	for (int i = 0; i < num_elements; i++)
	{
		InputArray[i] = std::complex<double>(in[i], 0);
	}

	std::complex<double> *ProcessArray = (std::complex<double>*)malloc(sizeof(std::complex<double>) * num_elements);

	fftw_plan plan_forward = fftw_plan_dft_1d(N, (fftw_complex*)ProcessArray, (fftw_complex*)ProcessArray, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_backward = fftw_plan_dft_1d(N, (fftw_complex*)ProcessArray, (fftw_complex*)ProcessArray, FFTW_BACKWARD, FFTW_ESTIMATE);

	for (int ix = 0; ix < num_elements; ix += N)
	{
		//printf("I: %d num_elements: %d\n", i, num_elements);
		if (ix + N > num_elements)
			break;
		memcpy(ProcessArray, InputArray + ix, sizeof(std::complex<double>) * N);
		fftw_execute(plan_forward);
		int hN = N >> 1; // half of the length (N/2)
		int numRem = hN; // the number of remaining elements
		for (int i = 1; i < hN; ++i)
		{
			ProcessArray[i] = std::complex<double>(((ProcessArray[i]).real() * 2), ((ProcessArray[i]).imag() * 2));
		}
		if (N % 2 == 0)
			numRem--;
		else if (N > 1)
		{
			ProcessArray[hN] = std::complex<double>(((ProcessArray[hN]).real() * 2), ((ProcessArray[hN]).imag() * 2));
		}
		memset(&ProcessArray[hN + 1], 0, numRem * sizeof(std::complex<double>));
		fftw_execute(plan_backward);
		for (int i = 0; i < N; ++i)
		{
			ProcessArray[i] = std::complex<double>(((ProcessArray[i]).real() / N), ((ProcessArray[i]).imag() / N));
		}
		for (int i = 0; i < N; i++)
		{
			out[ix + i] = std::complex<float>(ProcessArray[i].real(), ProcessArray[i].imag());
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
	float *samples = (float*)calloc(samp_count, sizeof(float));
	sf_readf_float(inFile, samples, samp_count); //??? fix this pointer thing pls why it produce warning error idk
	sf_close(inFile);
	printf("Be careful! This code will break with WAV files at 48000 Hz that are longer than about 12h. Actually it might be 6h but idk\n");
	printf("Sample Rate = %d Hz\n", samp_rate);
	printf("Sample Count = %ld\n", samp_count);
	sf_close(inFile);

	const int MixFrequency = 500; //user
	//const int Bandwidth = 1000;	   //user
	const int speedDivider = 10;   //user

	float Bandwidth = (float)(samp_rate / 2) / speedDivider;
	if(speedDivider > samp_rate)
	{
		printf("Minimum Bandwidth is 0.5Hz!\n");
		return 0;
	}

	printf("input Bandwidth: %fHz\n", Bandwidth);




	//Apply filters to I and Q	
	printf("Applying filters\n");
	//BWBandPass *filterR = create_bw_band_pass_filter(50, samp_rate, MixFrequency, (float)MixFrequency + 10);
	//BWLowPass *filterL = create_bw_low_pass_filter(50, samp_rate, MixFrequency + Bandwidth);
	//BWHighPass *filterH = create_bw_high_pass_filter(50, samp_rate, MixFrequency);
	fft_bandpass_filter(samples, samp_count, samp_rate, samp_count, MixFrequency, MixFrequency + Bandwidth);
	



	//Convert real array into complex, set Imaginary to zero
	printf("Converting to complex\n");
	lv_32fc_t *inputComplex = (lv_32fc_t *)volk_malloc(sizeof(lv_32fc_t) * samp_count, volk_get_alignment());
	Hilbert(samples, inputComplex, samp_count, samp_count);
	free(samples);

	printf("Mixing up\n");
	float sinAngle = 2.0 * 3.14159265359 * -MixFrequency / samp_rate;
	lv_32fc_t phase_increment = lv_cmake(cos(sinAngle), sin(sinAngle));
	lv_32fc_t phase = lv_cmake(1.f, 0.0f);
	volk_32fc_s32fc_x2_rotator_32fc(inputComplex, inputComplex, phase_increment, &phase, samp_count);



	
	
	

	/*char *outFileName = "demod_out_complex.wav";
	SNDFILE *outFile;
	SF_INFO outFileInfo = inFileInfo;
	outFileInfo.channels = 2;
	outFile = sf_open(outFileName, SFM_WRITE, &outFileInfo);
	sf_writef_float(outFile, (float *)inputComplex, samp_count);
	sf_close(outFile);*/

	//Downsample
	printf("Converting back to real + Resampling\n");
	float *outputReal = (float*)malloc((samp_count / speedDivider) * sizeof(float));
	int samplecounter = 0;
	for(int i = 0; i < samp_count; i++)
	{
		samplecounter++;
		if(samplecounter == speedDivider)
		{
			outputReal[i/speedDivider] = inputComplex[i].real();
			samplecounter = 0;
		}
	}
	volk_free(inputComplex);

	
	printf("Writing to file\n");
	char *outFileName2 = "demod_out.wav";
	SNDFILE *outFile2;
	SF_INFO outFileInfo2 = inFileInfo;
	outFile2 = sf_open(outFileName2, SFM_WRITE, &outFileInfo2);
	sf_writef_float(outFile2, outputReal, samp_count / speedDivider);
	sf_close(outFile2);
	

	return 0;
}

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
#include "dsp.h"

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

	inFileName = "in.wav";
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
	const int speedDivider = 10;   //user

	if(speedDivider > samp_rate)
	{
		printf("Minimum Bandwidth is 0.5Hz!\n");
		return 0;
	}

	

	//Convert real array into complex, set Imaginary to zero
	printf("Converting to complex\n");
	std::complex<float> *inputComplex1 = (std::complex<float>*)calloc(samp_count, sizeof(std::complex<float>));
	//Hilbert(samples, inputComplex, samp_count, samp_count);
	DSP::filters::fftbrickwallhilbert *hilbert = new DSP::filters::fftbrickwallhilbert(100, samp_count);
	hilbert->processSamples(samp_count, samples, inputComplex1);
	free(samples);


	printf("Upsampling\n");
	if (speedDivider > 1)
	{
		printf("Output Bandwidth: %fHz\n", (float)(samp_rate / 2) / speedDivider);
	}
	std::complex<float> *inputComplex = (std::complex<float>*)calloc(samp_count * speedDivider, sizeof(std::complex<float>));
	DSP::upsamplers::complexUpsampler upsampler(samp_count, speedDivider);
	upsampler.processSamples(inputComplex1, inputComplex);
	samp_count = samp_count * speedDivider;
	free(inputComplex1);


	printf("Mixing up\n");
	unsigned int alignment = volk_get_alignment();
	lv_32fc_t *outputComplex = (lv_32fc_t *)volk_malloc(sizeof(lv_32fc_t) * samp_count, alignment);
	float sinAngle = 2.0 * 3.14159265359 * MixFrequency / samp_rate;
	lv_32fc_t phase_increment = lv_cmake(cos(sinAngle), sin(sinAngle));
	lv_32fc_t phase = lv_cmake(1.f, 0.0f);
	volk_32fc_s32fc_x2_rotator_32fc(outputComplex, inputComplex, phase_increment, &phase, samp_count);
	free(inputComplex);

	/*char *outFileName = "out_complex.wav";
	SNDFILE *outFile;
	SF_INFO outFileInfo = inFileInfo;
	outFileInfo.channels = 2;
	outFile = sf_open(outFileName, SFM_WRITE, &outFileInfo);
	sf_writef_float(outFile, (float *)outputComplex, samp_count);
	sf_close(outFile);*/

	printf("Converting back to real\n");
	float *outputReal = (float*)calloc(samp_count, sizeof(float));
	for (int i = 0; i < samp_count; i++)
	{
		outputReal[i] = outputComplex[i].real();
	}
	volk_free(outputComplex);

	
	printf("Writing to output file\n");
	WavWriter writer("out.wav", 32, 1, samp_rate);
	if(writer.isOpen())
	{
		writer.writeData(outputReal, samp_count * sizeof(float));
	}
	else {
		printf("could not open output file");
	}
	writer.finish();

	return 0;
}

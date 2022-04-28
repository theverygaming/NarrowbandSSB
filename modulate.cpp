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
	DSP::filters::fftbrickwallhilbert *hilbert = new DSP::filters::fftbrickwallhilbert(300, samp_count);
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

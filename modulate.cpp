#include <stdio.h>
#include <stdbool.h>
#include <sndfile.h>
#include <math.h>
#include <string.h>
#include <complex>
#include <fftw3.h>
#include <volk/volk.h>
#include "wavwriter.h"
#include "dsp/dsp.h"
#include "dsp/gain.h"
#include "dsp/firfilters.h"
#include "dsp/fftfilters.h"
#include "dsp/resamplers.h"

int main()
{
	char *inFileName;
	SNDFILE *inFile;
	SF_INFO inFileInfo;

	inFileName = "in.wav";
	inFile = sf_open(inFileName, SFM_READ, &inFileInfo);
	sf_count_t samp_count = inFileInfo.frames;
	int samp_rate = inFileInfo.samplerate;
	float *samplesIn = (float*)calloc(samp_count, sizeof(float));
	sf_readf_float(inFile, samplesIn, samp_count);
	sf_close(inFile);

	const int MixFrequency = 500;
	const int speedDivider = 100;

	const int chunkSize = 1000;

	if(samp_count < chunkSize)
	{
		printf("Your input file does not have enough samples!, expect this to crash\n");
	}

	if(speedDivider > samp_rate / 2)
	{
		printf("wow that's slow, don't expect this to work well!\n");
	}
	if (speedDivider > 1)
	{
		printf("Output Bandwidth: %fHz\n", (float)(samp_rate / 2) / speedDivider);
	}


	WavWriter writer("out.wav", 32, 1, samp_rate);
	std::complex<float> *inputComplex1 = (std::complex<float>*)malloc(chunkSize * sizeof(std::complex<float>));
	std::complex<float> *inputComplex = (std::complex<float>*)malloc(chunkSize * speedDivider * sizeof(std::complex<float>));
	lv_32fc_t *outputComplex = (lv_32fc_t *)volk_malloc(sizeof(lv_32fc_t) * chunkSize * speedDivider, volk_get_alignment());
	float *outputReal = (float*)malloc(chunkSize * speedDivider * sizeof(float));

	dsp::gain::agc agc(20, samp_rate, 0.95);
	dsp::filters::fftbrickwallhilbert *hilbert = new dsp::filters::fftbrickwallhilbert(300, samp_count);
	dsp::resamplers::complexUpsampler upsampler(chunkSize, speedDivider, speedDivider * 10);
	
	float* bpfCoeffs = (float*)malloc(speedDivider * 10 * sizeof(float));
	dsp::filters::FIRcoeffcalc::calcCoeffs_band(dsp::filters::FIRcoeffcalc::bandpass, bpfCoeffs, speedDivider * 10, samp_rate, MixFrequency, MixFrequency + (samp_rate / 2 / speedDivider));
	dsp::filters::FIRfilter finalBpf(speedDivider * 10, bpfCoeffs, chunkSize * speedDivider);

	float sinAngle = 2.0 * 3.14159265359 * MixFrequency / samp_rate;
	lv_32fc_t phase_increment = lv_cmake(cos(sinAngle), sin(sinAngle));
	lv_32fc_t phase = lv_cmake(1.f, 0.0f);

	for(sf_count_t x = 0; x < samp_count; x += chunkSize)
	{
		float *samples = &samplesIn[x];

		agc.process(samples, samples, chunkSize);

		hilbert->processSamples(chunkSize, samples, inputComplex1);
		
		upsampler.processSamples(inputComplex1, inputComplex);


		volk_32fc_s32fc_x2_rotator_32fc(outputComplex, inputComplex, phase_increment, &phase, chunkSize * speedDivider);

		for (int i = 0; i < chunkSize * speedDivider; i++)
		{
			outputReal[i] = outputComplex[i].real();
		}
		
		finalBpf.filter(outputReal, outputReal, chunkSize * speedDivider);

		if (writer.isOpen())
		{
			writer.writeData(outputReal, chunkSize * speedDivider * sizeof(float));
		}
		else
		{
			printf("could not open output file");
			return 1;
		}
		if(x % 4000 == 0)
		{
			float percentage = ((float)(x + chunkSize) / samp_count) * 100;
			printf("%.1f%%\n", percentage);
		}
	}
	writer.finish();
	free(inputComplex1);
	delete hilbert;
	free(inputComplex);
	volk_free(outputComplex);
	free(outputReal);
	free(bpfCoeffs);
	free(samplesIn);


	return 0;
}

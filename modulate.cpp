#include <complex>
#include <dsp/fftfilters.h>
#include <dsp/firfilters.h>
#include <dsp/gain.h>
#include <dsp/mixer.h>
#include <dsp/resamplers.h>
#include <dsp/wav.h>
#include <fftw3.h>
#include <math.h>
#include <sndfile.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <volk/volk.h>

int main() {
    char *inFileName;
    SNDFILE *inFile;
    SF_INFO inFileInfo;

    inFileName = "in.wav";
    inFile = sf_open(inFileName, SFM_READ, &inFileInfo);
    sf_count_t samp_count = inFileInfo.frames;
    int samp_rate = inFileInfo.samplerate;
    float *samplesIn = (float *)calloc(samp_count, sizeof(float));
    sf_readf_float(inFile, samplesIn, samp_count);
    sf_close(inFile);

    const int MixFrequency = 500;
    const int speedDivider = 100;

    const int chunkSize = 1000;

    if (samp_count < chunkSize) {
        printf("Your input file does not have enough samples!, expect this to crash\n");
    }

    if (speedDivider > samp_rate / 2) {
        printf("wow that's slow, don't expect this to work well!\n");
    }
    if (speedDivider > 1) {
        printf("Output Bandwidth: %fHz\n", (float)(samp_rate / 2) / speedDivider);
    }

    dsp::wav::wavWriter writer("out.wav", 32, 1, samp_rate);

    if (!writer.isOpen()) {
        printf("could not open output file");
        return 1;
    }

    float *input = (float *)malloc(chunkSize * sizeof(float));
    float *output = (float *)malloc(chunkSize * speedDivider * sizeof(float));

    dsp::gain::agc agc(20, samp_rate, 0.95);
    dsp::resamplers::realUpsampler upsampler(chunkSize, speedDivider, speedDivider * 10);

    dsp::mixer::real_mixer mixer(MixFrequency, samp_rate);

    std::vector<float> bpfCoeffs =
        dsp::filters::FIRcoeffcalc::calcCoeffs_band(dsp::filters::FIRcoeffcalc::bandpass, speedDivider * 100, samp_rate, MixFrequency, MixFrequency + (samp_rate / 2 / speedDivider));
    dsp::filters::FIRfilter finalBpf(bpfCoeffs, chunkSize * speedDivider);

    for (sf_count_t x = 0; x < samp_count; x += chunkSize) {
        float *samples = &samplesIn[x];
        int chunkSize2 = chunkSize;
        if (x + chunkSize > samp_count - 1) {
            chunkSize2 = chunkSize - (x + chunkSize - (samp_count - 1));
        }

        agc.process(samples, samples, chunkSize2);

        upsampler.upsample(chunkSize2, samples, output);

        mixer.run(output, output, chunkSize2 * speedDivider);

        finalBpf.run(output, output, chunkSize2 * speedDivider);

        writer.writeData(output, chunkSize2 * speedDivider * sizeof(float));

        if (x % 4000 == 0) {
            float percentage = ((float)(x + chunkSize) / samp_count) * 100;
            printf("%.1f%%\n", percentage);
        }
    }

    writer.finish();

    return 0;
}

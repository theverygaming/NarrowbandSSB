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

    inFileName = "out.wav";
    inFile = sf_open(inFileName, SFM_READ, &inFileInfo);
    sf_count_t samp_count = inFileInfo.frames;
    int samp_rate = inFileInfo.samplerate;
    float *samplesIn = (float *)malloc(samp_count * sizeof(float));
    sf_readf_float(inFile, samplesIn, samp_count);
    sf_close(inFile);

    const int MixFrequency = 500;
    const int speedDivider = 50;

    const int chunkSize = 1000;

    float Bandwidth = (float)(samp_rate / 2) / speedDivider;
    if (samp_count < chunkSize) {
        printf("Your input file does not have enough samples!, expect this to crash\n");
    }

    if (speedDivider > samp_rate / 2) {
        printf("wow that's slow, don't expect this to work well!\n");
    }

    printf("input Bandwidth: %fHz\n", Bandwidth);

    std::vector<float> bpf_coeffs = dsp::filters::FIRcoeffcalc::calcCoeffs_band(dsp::filters::FIRcoeffcalc::bandpass, speedDivider * 10, samp_rate, MixFrequency, MixFrequency + Bandwidth);
    dsp::filters::FIRfilter bpf(bpf_coeffs, chunkSize);

    std::vector<float> lpf_coeffs = dsp::filters::FIRcoeffcalc::calcCoeffs_band(dsp::filters::FIRcoeffcalc::bandpass, speedDivider * 10, samp_rate, 0, Bandwidth);
    dsp::filters::FIRfilter lpf(lpf_coeffs, chunkSize);

    dsp::mixer::real_mixer mixer(MixFrequency, samp_rate);

    dsp::resamplers::realDownsampler downsampler(speedDivider);
    float *output = (float *)malloc((chunkSize / speedDivider) * sizeof(float));

    dsp::wav::wavWriter writer("demod_out.wav", 32, 1, samp_rate);

    if (!writer.isOpen()) {
        printf("could not open output file");
        return 1;
    }

    for (sf_count_t x = 0; x < samp_count; x += chunkSize) {
        float *samples = &samplesIn[x];
        int chunkSize2 = chunkSize;
        if (x + chunkSize > samp_count - 1) {
            chunkSize2 = chunkSize - (x + chunkSize - (samp_count - 1));
        }
        bpf.run(samples, samples, chunkSize2);

        mixer.run(samples, samples, chunkSize2);

        lpf.run(samples, samples, chunkSize2);

        downsampler.downsample(chunkSize2, samples, output);

        writer.writeData(output, (chunkSize2 / speedDivider) * sizeof(float));
        if (x % 100000 == 0) {
            float percentage = ((float)(x + chunkSize) / samp_count) * 100;
            printf("%.1f%%\n", percentage);
        }
    }
    writer.finish();
    free(samplesIn);
    free(output);

    return 0;
}

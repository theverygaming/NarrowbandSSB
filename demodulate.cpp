#include <complex>
#include <dsp/fftfilters.h>
#include <dsp/firfilters.h>
#include <dsp/gain.h>
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
    const int speedDivider = 10;
    const int chunkSize = 1000;

    float Bandwidth = (float)(samp_rate / 2) / speedDivider;
    if (samp_count < chunkSize) {
        printf("Your input file does not have enough samples!, expect this to crash\n");
    }

    if (speedDivider > samp_rate / 2) {
        printf("wow that's slow, don't expect this to work well!\n");
    }

    printf("input Bandwidth: %fHz\n", Bandwidth);

    float *bpfCoeffs = (float *)malloc(speedDivider * 10 * sizeof(float));
    dsp::filters::FIRcoeffcalc::calcCoeffs_band(dsp::filters::FIRcoeffcalc::bandpass, bpfCoeffs, speedDivider * 10, samp_rate, MixFrequency, MixFrequency + Bandwidth);
    dsp::filters::FIRfilter firstBpf(speedDivider * 10, bpfCoeffs, chunkSize);

    dsp::filters::fftbrickwallhilbert hilbert(300, chunkSize);

    float sinAngle = 2.0 * 3.14159265359 * -MixFrequency / samp_rate;
    lv_32fc_t phase_increment = lv_cmake(cos(sinAngle), sin(sinAngle));
    lv_32fc_t phase = lv_cmake(1.f, 0.0f);
    lv_32fc_t *inputComplex = (lv_32fc_t *)volk_malloc(sizeof(lv_32fc_t) * samp_count, volk_get_alignment());

    dsp::resamplers::realDownsampler downsampler(speedDivider);
    float *processReal = (float *)malloc(chunkSize * sizeof(float));
    float *outputReal = (float *)malloc((chunkSize / speedDivider) * sizeof(float));

    dsp::wav::wavWriter writer("demod_out.wav", 32, 1, samp_rate);
    for (sf_count_t x = 0; x < samp_count; x += chunkSize) {
        float *samples = &samplesIn[x];
        int chunkSize2 = chunkSize;
        if (x + chunkSize > samp_count - 1) {
            chunkSize2 = chunkSize - (x + chunkSize - (samp_count - 1));
        }
        firstBpf.filter(samples, samples, chunkSize2);

        hilbert.processSamples(chunkSize2, samples, inputComplex);

        volk_32fc_s32fc_x2_rotator_32fc(inputComplex, inputComplex, phase_increment, &phase, chunkSize2);

        for (int i = 0; i < chunkSize2; i++) {
            processReal[i] = inputComplex[i].real();
        }

        downsampler.downsample(chunkSize2, processReal, outputReal);

        if (writer.isOpen()) {
            writer.writeData(outputReal, chunkSize2 / speedDivider * sizeof(float));
        } else {
            printf("could not open output file");
            return 1;
        }
        if (x % 100000 == 0) {
            float percentage = ((float)(x + chunkSize) / samp_count) * 100;
            printf("%.1f%%\n", percentage);
        }
    }
    writer.finish();
    free(samplesIn);
    free(bpfCoeffs);
    free(outputReal);
    free(processReal);
    volk_free(inputComplex);

    return 0;
}

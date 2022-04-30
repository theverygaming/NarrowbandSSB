#pragma once
#include "firfilters.h"

namespace dsp
{
    namespace resamplers
    {
        class realUpsampler
        {
        public:
            realUpsampler(int multiplier, int chunkSize, int taps)
            {
                _multiplier = multiplier;
                _chunkSize = chunkSize;
                coeffs = (float *)malloc(taps * sizeof(float));
                filters::FIRcoeffcalc::calcCoeffs(dsp::filters::FIRcoeffcalc::lowpass, coeffs, taps, 48000, 24000 * ((float)1 / multiplier));
                lpf = new filters::FIRfilter(taps, coeffs, _chunkSize * multiplier);
                processarr = (float *)malloc(chunkSize * _multiplier * sizeof(float));
            }

            ~realUpsampler()
            {
                free(coeffs);
                free(processarr);
                delete lpf;
            }

            void upsample(int incount, float *in, float *out)
            {
                for (long int i = 0; i < incount; i++)
                {
                    processarr[i * _multiplier] = in[i];
                }
                lpf->filter(processarr, out, incount * _multiplier);
            }

        private:
            filters::FIRfilter *lpf;
            float *coeffs;
            float *processarr;
            int _multiplier;
            int _chunkSize;
        };
        class complexUpsampler
        {
        public:
            complexUpsampler(int chunkSize, int multiplier, int taps)
            {
                _chunkSize = chunkSize;
                _multiplier = multiplier;

                realInArr = (float *)calloc(_chunkSize, sizeof(float));
                imagInArr = (float *)calloc(_chunkSize, sizeof(float));

                realOutArr = (float *)malloc(_chunkSize * _multiplier * sizeof(float));
                imagOutArr = (float *)malloc(_chunkSize * _multiplier * sizeof(float));

                upReal = new resamplers::realUpsampler(_multiplier, _chunkSize, taps);
                upImag = new resamplers::realUpsampler(_multiplier, _chunkSize, taps);
            }

            ~complexUpsampler()
            {
                free(realInArr);
                free(imagInArr);

                free(realOutArr);
                free(imagOutArr);

                delete upReal;
                delete upImag;
            }

            void processSamples(std::complex<float> *in, std::complex<float> *out)
            {
                for (int i = 0; i < _chunkSize; i++)
                {
                    realInArr[i] = in[i].real();
                    imagInArr[i] = in[i].imag();
                }
                upReal->upsample(_chunkSize, realInArr, realOutArr);
                upImag->upsample(_chunkSize, imagInArr, imagOutArr);
                for (int i = 0; i < _chunkSize * _multiplier; i++)
                {
                    out[i] = {realOutArr[i], imagOutArr[i]};
                }
            }

        private:
            int _chunkSize;
            int _multiplier;
            float *realInArr;
            float *imagInArr;
            float *realOutArr;
            float *imagOutArr;
            resamplers::realUpsampler *upReal;
            resamplers::realUpsampler *upImag;
        };
        class realDownsampler
        {
        public:
            realDownsampler(int divider)
            {
                _divider = divider;
            }

            ~realDownsampler()
            {
            }

            void downsample(int incount, float *in, float *out)
            {
                for (int i = 0; i < incount; i += _divider)
                {
                    out[i / _divider] = in[i];
                }
            }

        private:
            int _divider;
        };
    }
}
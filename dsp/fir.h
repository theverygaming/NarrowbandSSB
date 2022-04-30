#pragma once
#include <volk/volk.h>
#include <string.h>
#include "math.h"
#include "window.h"

namespace dsp {
    namespace filters
    {
        class FIRcoeffcalc
        {
        public:
            enum filter_type
            {
                lowpass,
                highpass,
                bandpass,
                bandstop
            };

            static void calcCoeffs(filter_type type, float* impulseresponse, int taps, int samplerate, float frequency)
            {
                float sampletime = 1 / (float)samplerate;
                for (int i = 0; i < taps; i++)
                {
                    int x = i - (taps / 2);
                    if(x == 0) {
                        switch (type)
                        {
                        case lowpass:
                            impulseresponse[i] = 2 * frequency;
                            break;
                        case highpass:
                            break;
                        }
                        continue;
                    }
                    switch (type) {
                    case lowpass:
                        impulseresponse[i] = sin(2 * FL_M_PI * frequency * sampletime * x) / (FL_M_PI * sampletime * x);
                        break;
                    case highpass:
                        break;
                    }
                }
                for(int i = 0; i < taps; i++)
                {
                    //impulseresponse[i] *= sampletime;
                    impulseresponse[i] /= frequency * 2;
                    impulseresponse[i] *= windowfunctions::blackman(i, taps - 1);
                }
            }

            static void calcCoeffs_band(filter_type type, float* impulseresponse, int taps, int samplerate, float freq1, float freq2)
            {
                float sampletime = 1 / (float)samplerate;
                for (int i = 0; i < taps; i++)
                {
                    int x = i - (taps / 2);
                    if(x == 0) {
                        switch (type)
                        {
                        case bandpass:
                            impulseresponse[i] = 2 * freq2 - 2 * freq1;
                            break;
                        case bandstop:
                            break;
                        }
                        continue;
                    }
                    switch (type) {
                    case bandpass:
                        impulseresponse[i] = (sin(2 * FL_M_PI * freq2 * sampletime * x) - sin(2 * FL_M_PI * freq1 * sampletime * x)) / (FL_M_PI * sampletime * x);
                        break;
                    case bandstop:
                        break;
                    }
                }
                for(int i = 0; i < taps; i++)
                {
                    impulseresponse[i] *= sampletime;
                    //impulseresponse[i] /= 2 * freq2 - 2 * freq1;
                    impulseresponse[i] *= windowfunctions::blackman(i, taps - 1);
                }
            }

        private:
        };
        
        class FIRfilter {
        public:
            FIRfilter(int taps, float* coeffs, int chunkSize)
            {
                _coeffs = coeffs;
                _taps = taps;
                buffer = (float*)malloc(chunkSize * 2 * sizeof(float));
                bufferStart = &buffer[_taps];
            }

            ~FIRfilter()
            {
                free(buffer);
            }

            void filter(float* in, float* out, long int count) {
                memcpy(bufferStart, in, count * sizeof(float));
                for(long int i = 0; i < count; i++)
                {
                    volk_32f_x2_dot_prod_32f(&out[i], &buffer[i + 1], _coeffs, _taps);
                }
                memmove(buffer, &buffer[count], _taps * sizeof(float));
            }

        private:
            float* _coeffs;
            int _taps;
            float* buffer;
            float* bufferStart;
        };
    }
}
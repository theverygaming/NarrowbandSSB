#pragma once
#include <math.h>
#include <cstring>
#include <volk/volk.h>
#include <fftw3.h>
#include "filter.h"

#define FL_M_PI 3.1415926535f

namespace DSP {
    // https://github.com/AlexandreRouma/SDRPlusPlus/blob/master/core/src/dsp/utils/window_functions.h
    namespace windowfunctions {
        inline double blackman(double n, double N, double alpha = 0.16f) {
            double a0 = (1.0f - alpha) / 2.0f;
            double a2 = alpha / 2.0f;
            return a0 - (0.5f * cos(2.0f * FL_M_PI * (n / N))) + (a2 * cos(4.0f * FL_M_PI * (n / N)));
        }
    }


    namespace upsamplers {
        class realUpsampler {
        public:
            realUpsampler(int multiplier) {
                _multiplier = multiplier;
                sampFilter = create_bw_low_pass_filter(30, 100 * _multiplier, 100 / 2);
            }

            ~realUpsampler() {
                free_bw_low_pass(sampFilter);
            }

            void upsample(int incount, float *in, float *out)
            {
                for (long int i = 0; i < incount - 1; i++)
		        {
			        out[i * _multiplier] = in[i];
		        }
                for(long int i = 0; i < _multiplier - 1; i++)
		        {
			        out[i] = bw_low_pass(sampFilter, out[i]) * (_multiplier * 0.5);
		        }
            }
        
        private:
            BWLowPass *sampFilter;
            int _multiplier;
        };
        class complexUpsampler {
        public:
            complexUpsampler(int chunkSize, int multiplier) {
                _chunkSize = chunkSize;
                _multiplier = multiplier;
                
                realInArr = (float*)malloc(_chunkSize * sizeof(float));
                imagInArr = (float*)malloc(_chunkSize * sizeof(float));

                realOutArr = (float*)malloc(_chunkSize * _multiplier * sizeof(float));
                imagOutArr = (float*)malloc(_chunkSize * _multiplier * sizeof(float));

                upReal = new upsamplers::realUpsampler(_multiplier);
                upImag = new upsamplers::realUpsampler(_multiplier);
            }

            ~complexUpsampler() {
                free(realInArr);
                free(imagInArr);

                free(realOutArr);
                free(imagOutArr);
                
                delete upReal;
                delete upImag;
            }


            void processSamples(std::complex<float> *in, std::complex<float> *out) {
                for(int i = 0; i < _chunkSize; i++)
                {
                    realInArr[i] = in[i].real();
                    imagInArr[i] = in[i].imag();
                }
                upReal->upsample(_chunkSize, realInArr, realOutArr);
                upImag->upsample(_chunkSize, imagInArr, imagOutArr);
                for(int i = 0; i < _chunkSize * _multiplier; i++)
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
            upsamplers::realUpsampler *upReal;
            upsamplers::realUpsampler *upImag;
        };
    }
    

    namespace filters {
        class fftbrickwallbandpass {
        };

        // https://github.com/AlexandreRouma/SDRPlusPlus/blob/master/core/src/dsp/noise_reduction.h
        class fftbrickwallhilbert {
        public:
            fftbrickwallhilbert(int tapcount, int chunksize) {
                _chunksize = chunksize;
                _tapcount = tapcount;

                delay = (std::complex<float>*)fftwf_malloc(sizeof(std::complex<float>) * _chunksize);
                delay_start = &delay[_tapcount];
                memset(delay, 0, sizeof(std::complex<float>) * _chunksize);

                fft_window = (float*)fftwf_malloc(sizeof(float) * _tapcount);
                for (int i = 0; i < _tapcount; i++) {
                    fft_window[i] = windowfunctions::blackman(i, _tapcount - 1);
                }

                fft_in = (std::complex<float>*)fftwf_malloc(sizeof(std::complex<float>) * _tapcount);
                fft_cout = (std::complex<float>*)fftwf_malloc(sizeof(std::complex<float>) * _tapcount);
                fft_cin = (std::complex<float>*)fftwf_malloc(sizeof(std::complex<float>) * _tapcount);
                fft_fcout = (std::complex<float>*)fftwf_malloc(sizeof(std::complex<float>) * _tapcount);
                memset(fft_in, 0, sizeof(std::complex<float>) * _tapcount);
                memset(fft_cout, 0, sizeof(std::complex<float>) * _tapcount);
                memset(fft_cin, 0, sizeof(std::complex<float>) * _tapcount);
                memset(fft_fcout, 0, sizeof(std::complex<float>) * _tapcount);

                forwardPlan = fftwf_plan_dft_1d(_tapcount, (fftwf_complex*)fft_in, (fftwf_complex*)fft_cout, FFTW_FORWARD, FFTW_ESTIMATE);
                backwardPlan = fftwf_plan_dft_1d(_tapcount, (fftwf_complex*)fft_cin, (fftwf_complex*)fft_fcout, FFTW_BACKWARD, FFTW_ESTIMATE);

            }

            ~fftbrickwallhilbert() {
                fftwf_free(delay);
                fftwf_free(fft_window);

                fftwf_free(fft_in);
                fftwf_free(fft_cout);
                fftwf_free(fft_cin);
                fftwf_free(fft_fcout);
                
                fftwf_destroy_plan(forwardPlan);
                fftwf_destroy_plan(backwardPlan);
            }


            void processSamples(int count, float* inf, std::complex<float>* out) {
                for (int i = 0; i < count - _tapcount; i++)
	            {
		            delay_start[i] = std::complex<float>(inf[i], 0);
                    //delay[i] = std::complex<float>(inf[i], 0);
	            }
                //memcpy(delay_start, in, count * sizeof(std::complex<float>));
                for(int i = 0; i < count; i++)
                {
                    volk_32fc_32f_multiply_32fc((lv_32fc_t*)fft_in, (lv_32fc_t*)&delay[i], fft_window, _tapcount);

                    // Forward FFT
                    fftwf_execute(forwardPlan);

                    // Copy forward FFT out to Backward FFT input
                    //memcpy(fft_cin, fft_cout, sizeof(std::complex<float>) * _tapcount);

                    // "Filter"
                    for(int j = 0; j < _tapcount; j++)
                    {
                        if(j < (_tapcount / 2) )
                        {
                            fft_cin[j] = {fft_cout[j].real() * 2, fft_cout[j].imag() * 2};
                        }
                        else{
                            fft_cin[j] = {0,0};
                        }
                    }

                    // Backward FFT and write middle element to output buffer
                    fftwf_execute(backwardPlan);
                    out[i] = fft_fcout[_tapcount / 2];
                    if(i % 100000 == 0)
                    {
                        printf("%d/%d\n", i, count);
                    }
                }
                volk_32f_s32f_multiply_32f((float*)out, (float*)out, 1.0f / (float)_tapcount, count * 2);

                // Copy last values to delay
                memmove(delay, &delay[count], _tapcount * sizeof(std::complex<float>));
            }

        private:
            float* fft_window;

            int _chunksize;
            int _tapcount;

            fftwf_plan forwardPlan;
            fftwf_plan backwardPlan;

            std::complex<float> *delay_start;
            std::complex<float> *delay;
            
            std::complex<float>* fft_in;
            std::complex<float>* fft_cout;
            std::complex<float>* fft_cin;
            std::complex<float>* fft_fcout;
        };
    }
}
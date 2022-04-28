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

    namespace filters {
        class fftbrickwallbandpass {
        public:
            fftbrickwallbandpass(int tapcount, int chunksize, float lowCutoff, float highCutoff) {
                _chunksize = chunksize;
                _tapcount = tapcount;
                _lowCutoff = lowCutoff;
                _highCutoff = highCutoff;
                _tapcountC = _tapcount / 2 + 1;

                delay = (float*)malloc(sizeof(float) * _chunksize);
                delay_start = &delay[_tapcount];
                memset(delay, 0, sizeof(float) * _chunksize);

                fft_window = (float*)fftwf_malloc(sizeof(float) * _tapcount);
                for (int i = 0; i < _tapcount; i++) {
                    fft_window[i] = windowfunctions::blackman(i, _tapcount - 1);
                }

                fft_rin = (float*)malloc(sizeof(float) * _tapcount);
                fft_cout = (std::complex<float>*)fftwf_malloc(sizeof(std::complex<float>) * _tapcountC);
                fft_cin = (std::complex<float>*)fftwf_malloc(sizeof(std::complex<float>) * _tapcountC);
                fft_rout = (float*)malloc(sizeof(float) * _tapcount);
                memset(fft_rin, 0, sizeof(float) * _tapcount);
                memset(fft_cout, 0, sizeof(std::complex<float>) * _tapcountC);
                memset(fft_cin, 0, sizeof(std::complex<float>) * _tapcountC);
                memset(fft_rout, 0, sizeof(float) * _tapcount);

                forwardPlan = fftwf_plan_dft_r2c_1d(_tapcount, fft_rin, (fftwf_complex*)fft_cout, FFTW_ESTIMATE);
                backwardPlan = fftwf_plan_dft_c2r_1d(_tapcount, (fftwf_complex*)fft_cin, fft_rout, FFTW_ESTIMATE);

            }

            ~fftbrickwallbandpass() {
                free(delay);
                fftwf_free(fft_window);

                free(fft_rin);
                fftwf_free(fft_cout);
                fftwf_free(fft_cin);
                free(fft_rout);
                
                fftwf_destroy_plan(forwardPlan);
                fftwf_destroy_plan(backwardPlan);
            }


            void processSamples(int count, float* inf, float* out) {
                for (int i = 0; i < count - _tapcount; i++)
	            {
		            delay_start[i] = inf[i];
	            }
                //memcpy(delay_start, in, count * sizeof(std::complex<float>));
                for(int i = 0; i < count; i++)
                {
                    //volk_32fc_32f_multiply_32fc((lv_32fc_t*)fft_in, (lv_32fc_t*)&delay[i], fft_window, _tapcount);
                    for(int j = 0; j < _tapcount; j++)
                    {
                        fft_rin[j] = (&delay[i])[j] * fft_window[j];
                    }
                    // Forward FFT
                    fftwf_execute(forwardPlan);

                    // Copy forward FFT out to Backward FFT input
                    //memcpy(fft_cin, fft_cout, sizeof(std::complex<float>) * _tapcountC);

                    // "Filter"
                    for(int j = 0; j < _tapcountC; j++)
                    {
                        float binF = (float)(j+1) / _tapcountC;
                        if(binF > _lowCutoff && binF < _highCutoff) {
                            fft_cin[j] = fft_cout[j];
                        }
                        else {
                            fft_cin[j] = {0, 0};
                        }
                    }


                    // Backward FFT and write middle element to output buffer
                    fftwf_execute(backwardPlan);
                    //out[i] = inf[i];
                    out[i] = fft_rout[_tapcount / 2];
                    if(i % 100 == 0)
                    {
                        //printf("%d/%d in:%f out:%f\n", i, count, inf[i], out[i]);
                    }
                }
                volk_32f_s32f_multiply_32f((float*)out, (float*)out, 1.0f / (float)_tapcount, count);

                // Copy last values to delay
                memmove(delay, &delay[count], _tapcount * sizeof(float));
            }

        private:
            float* fft_window;

            int _chunksize;
            int _tapcount;
            int _tapcountC;

            float _lowCutoff;
            float _highCutoff;

            fftwf_plan forwardPlan;
            fftwf_plan backwardPlan;

            float* delay_start;
            float* delay;
            
            float* fft_rin;
            std::complex<float>* fft_cout;
            std::complex<float>* fft_cin;
            float* fft_rout;
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

    namespace upsamplers {
        class realUpsampler {
        public:
            realUpsampler(int multiplier, int chunkSize) {
                _multiplier = multiplier;
                bpf = new filters::fftbrickwallbandpass(100, chunkSize * _multiplier, 0, (float)1 / multiplier);
            }

            ~realUpsampler() {
                delete bpf;
            }

            void upsample(int incount, float *in, float *out)
            {
                for (long int i = 0; i < incount - 1; i++)
		        {
			        out[i * _multiplier] = in[i];
		        }
                bpf->processSamples(incount * _multiplier, out, out);
            }
        
        private:
            filters::fftbrickwallbandpass *bpf;
            int _multiplier;
        };
        class complexUpsampler {
        public:
            complexUpsampler(int chunkSize, int multiplier) {
                _chunkSize = chunkSize;
                _multiplier = multiplier;
                
                realInArr = (float*)calloc(_chunkSize, sizeof(float));
                imagInArr = (float*)calloc(_chunkSize, sizeof(float));

                realOutArr = (float*)malloc(_chunkSize * _multiplier * sizeof(float));
                imagOutArr = (float*)malloc(_chunkSize * _multiplier * sizeof(float));

                upReal = new upsamplers::realUpsampler(_multiplier, _chunkSize);
                upImag = new upsamplers::realUpsampler(_multiplier, _chunkSize);
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
}
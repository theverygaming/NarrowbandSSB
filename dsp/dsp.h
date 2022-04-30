#pragma once
#include <math.h>
#include <cstring>
#include <volk/volk.h>
#include <fftw3.h>
#include "window.h"
#include "fir.h"

namespace dsp {
    namespace filters {
        // https://github.com/AlexandreRouma/SDRPlusPlus/blob/master/core/src/dsp/noise_reduction.h
        class fftbrickwallhilbert {
        public:
            fftbrickwallhilbert(int tapcount, int chunksize) {
                _chunksize = chunksize;
                _tapcount = tapcount;

                delay = (std::complex<float>*)fftwf_malloc(sizeof(std::complex<float>) * _chunksize * 2);
                delay_start = &delay[_tapcount];
                memset(delay, 0, sizeof(std::complex<float>) * _chunksize * 2);

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
                for (int i = 0; i < count; i++)
	            {
		            delay_start[i] = {inf[i], 0};
	            }
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
                    //out[i] = delay_start[i];
                    out[i] = fft_fcout[_tapcount / 2];
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
            realUpsampler(int multiplier, int chunkSize, int taps) {
                _multiplier = multiplier;
                _chunkSize = chunkSize;
                //bpf = new filters::fftbrickwallbandpass(chunkSize * _multiplier, chunkSize * _multiplier, 0, (float)1 / multiplier);
                coeffs = (float*)malloc(taps * sizeof(float));
                filters::FIRcoeffcalc::calcCoeffs(dsp::filters::FIRcoeffcalc::lowpass, coeffs, taps, 48000, 24000 * ((float)1 / multiplier));
                lpf = new filters::FIRfilter(taps, coeffs, _chunkSize * multiplier);
                processarr = (float*)malloc(chunkSize * _multiplier * sizeof(float));
            }

            ~realUpsampler() {
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
            //filters::fftbrickwallbandpass *bpf;
            filters::FIRfilter *lpf;
            float* coeffs;
            float* processarr;
            int _multiplier;
            int _chunkSize;
        };
        class complexUpsampler {
        public:
            complexUpsampler(int chunkSize, int multiplier, int taps) {
                _chunkSize = chunkSize;
                _multiplier = multiplier;
                
                realInArr = (float*)calloc(_chunkSize, sizeof(float));
                imagInArr = (float*)calloc(_chunkSize, sizeof(float));

                realOutArr = (float*)malloc(_chunkSize * _multiplier * sizeof(float));
                imagOutArr = (float*)malloc(_chunkSize * _multiplier * sizeof(float));

                upReal = new upsamplers::realUpsampler(_multiplier, _chunkSize, taps);
                upImag = new upsamplers::realUpsampler(_multiplier, _chunkSize, taps);
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
#pragma once

#include <volk/volk.h>
#include <fftw3.h>
#include <cstring>
#include "window.h"


namespace dsp {
    namespace filters {
        // https://github.com/AlexandreRouma/SDRPlusPlus/blob/master/core/src/dsp/noise_reduction.h

        // This filter is very slow, replacement soon:tm:
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

    
}
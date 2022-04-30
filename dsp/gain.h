#pragma once
#include <volk/volk.h>
#include "math.h"

namespace dsp
{
    namespace gain
    {
        class agc
        {
        public:
            agc(float fallrate, int samplerate, float maxLevel)
            {
                _correctedFallRate = fallrate / samplerate;
                _maxLevel = maxLevel;
            }

            ~agc()
            {
            }

            void process(float *in, float *out, int count)
            {
                // https://github.com/AlexandreRouma/SDRPlusPlus/blob/master/core/src/dsp/processing.h
                level = pow(10, ((10.0f * log10f(level)) - (_correctedFallRate * count)) / 10.0f);

                if (level < 10e-14)
                {
                    level = 10e-14;
                }

                float absVal;
                for (int i = 0; i < count; i++)
                {
                    absVal = fabsf(in[i]);
                    if (absVal > level)
                    {
                        level = absVal;
                    }
                }

                volk_32f_s32f_multiply_32f(out, in, _maxLevel / level, count);
            }

        private:
            float _maxLevel;
            float _correctedFallRate;
            float level = 0;
        };
    }
}
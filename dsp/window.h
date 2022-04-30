#pragma once
#include "math.h"

namespace dsp {
    // https://github.com/AlexandreRouma/SDRPlusPlus/blob/master/core/src/dsp/utils/window_functions.h
    namespace windowfunctions {
        inline double blackman(double n, double N, double alpha = 0.16f) {
            double a0 = (1.0f - alpha) / 2.0f;
            double a2 = alpha / 2.0f;
            return a0 - (0.5f * cos(2.0f * FL_M_PI * (n / N))) + (a2 * cos(4.0f * FL_M_PI * (n / N)));
        }
    }
}
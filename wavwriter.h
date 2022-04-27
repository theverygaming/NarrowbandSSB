#pragma once

#include <fstream>
#include <cstring>

//Epicly copied from https://github.com/AlexandreRouma/SDRPlusPlus/blob/master/misc_modules/recorder/src/wav.h

class WavWriter {
public:
    WavWriter(std::string path, uint16_t bitDepth, uint16_t channelCount, uint32_t sampleRate) {
        outstream = std::ofstream(path.c_str(), std::ios::binary);
        header.formatHeaderLength = 16;
        header.sampleType = 1;
        if(bitDepth == 32) {header.sampleType = 3; } //32-bit is usually float
        header.channelCount = channelCount;
        header.sampleRate = sampleRate;
        header.bytesPerSecond = sampleRate * channelCount * (bitDepth / 8);
        header.bytesPerSample = (bitDepth / 8) * channelCount;
        header.bitDepth = bitDepth;
        outstream.write((char*)&header, sizeof(waveheader));
    }

    void finish() {
        header.fileSize = written + sizeof(waveheader) - 8;
        header.dataSize = written;
        outstream.seekp(0);
        outstream.write((char*)&header, sizeof(waveheader));
        outstream.close();
    }

    bool isOpen() {
        return outstream.is_open();
    }

    void writeData(void* data, size_t size) {
        outstream.write((char*)data, size);
        written += size;
    }

private:
    struct waveheader {
        char signature[4] = {'R','I','F','F'};
        uint32_t fileSize;
        char fileType[4] = {'W', 'A', 'V', 'E'};
        char formatMarker[4] = {'f', 'm', 't', ' '};
        uint32_t formatHeaderLength;
        uint16_t sampleType; 
        uint16_t channelCount;
        uint32_t sampleRate;
        uint32_t bytesPerSecond;
        uint16_t bytesPerSample;
        uint16_t bitDepth;
        char dataMarker[4] = {'d', 'a', 't', 'a'};
        uint32_t dataSize;
    };

    waveheader header;
    std::ofstream outstream;
    size_t written = 0;
};
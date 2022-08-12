#include "Image.hpp"

SImage::SImage(const char * filename) {
    _shape = {1,1,1};
    _pos = {0,0,0};
    _nFrames = 1;
    _nRepeats = 1;
}

void SImage::loadData() {
    data = new int16_t[1];//std::shared_ptr<int16_t[]>(new int16_t[1]);
}

// Array SImage::get(const char * property) {
//     return factory.createArray<double>({0});
// }
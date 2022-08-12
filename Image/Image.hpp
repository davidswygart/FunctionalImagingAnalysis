#include <cstdint>
#include <memory>

template <typename T> 
struct __declspec( dllexport ) Vec3 {
    T x;
    T y;
    T z;
};

template <typename T>
class __declspec( dllexport ) Image {
    protected:
        Vec3<int> _shape;
        Vec3<double> _pos;
        int _nFrames;
        int _nRepeats;
        T* data = nullptr;
        virtual void loadData() = 0;

    public:
        Vec3<int> shape() { return _shape; };
        int nx() { return _shape.x; };
        int ny() { return _shape.y; };
        int nz() { return _shape.z; };
        Vec3<double> pos() { return _pos; };
        int nFrames() { return _nFrames; };
        int nRepeats() { return _nRepeats; };
        int nBytes() { return _shape.x * _shape.y * _shape.z * _nRepeats * _nFrames * sizeof(T); };

        const T* get() {
            if (data == nullptr) {
                loadData();
            }
            return data;
        };
};

class __declspec( dllexport ) SImage : public Image<int16_t> { //implementation of Image for scanimage images
    private:
        void loadData() override;
    protected:
        char * metadata;

    public:
        SImage(const char * filename);
    
};
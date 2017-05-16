#ifndef JFFT_H
#define JFFT_H

#include <QDebug>
#include <complex>
#include <cmath>
#include <assert.h>
#include <QVector>

//Radix 2. Complex 1 dimentional in place FFT/IFFT
//included is also a SFT (Slow Fourier Transform) on my desktop a 16384 point transform took 1.2ms for my FFT and
//3.4s for my SFT. Thats about 3000 times faster
class JFFT
{
public:
    typedef std::complex<double> cpx_type;//this is the complex number definition. people really should use complex numbers more they are much better than real numbers
    typedef enum fft_direction_t{FORWARD,INVERSE}fft_direction_t;
    JFFT();
    bool init(int &fft_size);
    bool fft(cpx_type *x,int size,fft_direction_t fft_direction=FORWARD);//this is the main fft/ifft function and is the fastest
    bool fft_easy_to_understand(cpx_type *x,int size,fft_direction_t fft_direction=FORWARD);//this is the same but is slower and eisier to undrstand
    bool sft(cpx_type *x,int size,fft_direction_t fft_direction=FORWARD);
    int get_nfft(){return nfft;}

    //convenience functions and others if you want to add them. these are a few of the ones I use
    bool fft(std::vector<cpx_type> &x)
    {
        if(((int)x.size())!=nfft)
        {
            int tnfft=x.size();
            init(tnfft);
            x.resize(tnfft,0);
        }
        return fft(x.data(),x.size());
    }
    bool ifft(std::vector<cpx_type> &x)
    {
        if(((int)x.size())!=nfft)
        {
            int tnfft=x.size();
            init(tnfft);
            x.resize(tnfft,0);
        }
        return fft(x.data(),x.size(),INVERSE);
    }
    bool fft(QVector<cpx_type> &x)
    {
        if(x.size()!=nfft)
        {
            int tnfft=x.size();
            init(tnfft);
            x.resize(tnfft);
        }
        return fft(x.data(),x.size());
    }
    bool ifft(QVector<cpx_type> &x)
    {
        if(x.size()!=nfft)
        {
            int tnfft=x.size();
            init(tnfft);
            x.resize(tnfft);
        }
        return fft(x.data(),x.size(),INVERSE);
    }

private:

    //size in both power of 2 and number
    int nfft_2power=0;
    int nfft=0;

    //memory
    std::vector<cpx_type> TWIDDLE_mem;
    std::vector<cpx_type> TWIDDLE_INV_mem;
    std::vector<cpx_type> F;//used if the slow DFT is done

    void inline swap(cpx_type &a,cpx_type &b){cpx_type c_tmp=b;b=a;a=c_tmp;}


};


#endif // JFFT_H

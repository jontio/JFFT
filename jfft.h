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

//----------------


//an example of 1D FastFir (1D Fast convolution)
class JFastFir
{
public:
    JFastFir();
    void SetKernel(const JFFT::cpx_type *kernel,int size);
    void update_block(JFFT::cpx_type *buffer,int size);//process a block at a time this does not seem to be any faster than the convenience function
    JFFT::cpx_type update(JFFT::cpx_type in_val);//process one sample at a time
    JFFT::cpx_type update_easy_to_understand(JFFT::cpx_type in_val);//process one sample at a time. easy to understand

    //convenience functions
    void update(JFFT::cpx_type *buffer,int size)//process a block at a time version 1
    {
        for(int i=0;i<size;++i)
        {
            buffer[i]=update(buffer[i]);
        }
    }
    void update(std::vector<JFFT::cpx_type> &buffer)//process a block at a time
    {
        update(buffer.data(),buffer.size());
    }
    void update(QVector<JFFT::cpx_type> &buffer)//process a block at a time
    {
        update(buffer.data(),buffer.size());
    }
    void SetKernel(const std::vector<JFFT::cpx_type> &_kernel)//for a complex kernel as a vector
    {
        SetKernel(_kernel.data(),_kernel.size());
    }
    void SetKernel(const std::vector<double> &_kernel)//for a real kernel
    {
        std::vector<JFFT::cpx_type> tmp_kernel;
        tmp_kernel.resize(_kernel.size());
        for(int i=0;i<((int)_kernel.size());++i)tmp_kernel[i]=_kernel[i];
        SetKernel(tmp_kernel);
    }
    double update(double in_val)//process one sample at a time for a real signal
    {
        JFFT::cpx_type tmp_in_val=in_val;
        return update(tmp_in_val).real();
    }

private:
    JFFT fft;
    std::vector<JFFT::cpx_type> kernel;
    std::vector<JFFT::cpx_type> sigspace;//in out buffer
    std::vector<JFFT::cpx_type> remainder;//used for overlap

    JFFT::cpx_type *pkernel;
    JFFT::cpx_type *psigspace;
    JFFT::cpx_type *premainder;
    JFFT::cpx_type *psigspace_overlap;

    int kernel_non_zero_size=0;
    int remainder_size=0;

    //this is for in and out buffer
    int signal_non_zero_size=0;
    int sigspace_ptr=0;

    int nfft=0;//fft size

    //for block prosessing using version 2
    std::vector<JFFT::cpx_type> tmp_space;


};

//------------------------


//filter design
//all designs are using the window method and derived from the low pass filter

class JFilterDesign
{
public:
    JFilterDesign(){}
    static std::vector<JFFT::cpx_type> LowPassHanning(double FrequencyCutOff, double SampleRate, int Length);
    static std::vector<JFFT::cpx_type> HighPassHanning(double FrequencyCutOff, double SampleRate, int Length);
    static std::vector<JFFT::cpx_type> BandPassHanning(double LowFrequencyCutOff,double HighFrequencyCutOff, double SampleRate, int Length);
    static std::vector<JFFT::cpx_type> BandStopHanning(double LowFrequencyCutOff,double HighFrequencyCutOff, double SampleRate, int Length);
private:
    static double sinc_normalized(double val);
};

#endif // JFFT_H

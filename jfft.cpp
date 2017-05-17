#include "jfft.h"


using namespace std;

JFFT::JFFT()
{

}

bool JFFT::init(int &fft_size)
{

    //use a bigger FFT size if its not a power of 2
    nfft = 1;
    nfft_2power=0;
    while(nfft < fft_size)
    {
        nfft<<=1;
        nfft_2power++;
    }
    fft_size=nfft;
    //qDebug()<<nfft<<nfft_2power;

    //alocate mem
    TWIDDLE_mem.resize(nfft);//this is more memory than we need. what is the exact amount we need?
    TWIDDLE_INV_mem.resize(nfft);

    //load twiddles (these are just roots of unity, like cutting a cake)
    //looking at http://www.alwayslearn.com/DFT%20and%20FFT%20Tutorial/DFTandFFT_FFT_Butterfly_8_Input.html
    //they are inserted as W^0_2, W^0_4,W^1_4, W^0_8,W^1_8,W^2_8,W^3_8 ...
    //the pattern is 2 4 8 16 32 .. for the bottom number (cake cut into this many piecies N) and the top number
    //increases till not quite half way around the cake (ie less than N/2)
    cpx_type imag=cpx_type(0,1);
    cpx_type *TWIDDLE=TWIDDLE_mem.data();
    cpx_type *TWIDDLE_INV=TWIDDLE_INV_mem.data();
    int w=0;
    for(int N=2;N<=nfft;N<<=1)
    {
        for(int i=0;i<N/2;i++)
        {
            cpx_type twiddle=exp(-2.0*imag*M_PI*((double)i)/((double)N));
            cpx_type twiddle_inv=exp(2.0*imag*M_PI*((double)i)/((double)N));
            assert(w<nfft);
            TWIDDLE[w]=twiddle;
            TWIDDLE_INV[w]=twiddle_inv;
            w++;
            //qDebug()<<i<<N<<twiddle.real()<<twiddle.imag();
        }
    }

    return true;


}

bool JFFT::fft(cpx_type *x,int size,fft_direction_t fft_direction)
{
    assert(size==nfft);
    cpx_type *TWIDDLE;
    if(fft_direction==FORWARD)TWIDDLE=TWIDDLE_mem.data();
     else TWIDDLE=TWIDDLE_INV_mem.data();

    //for the ifft an alternitive tick is given at http://www.embedded.com/design/configurable-systems/4210789/DSP-Tricks--Computing-inverse-FFTs-using-the-forward-FFT
    //it would mean taking the conjugate of x before the forward fft is done then taking the conjugate after the fft is done
    //this would be a good idea for memory limited devices

    //bit reversal. fast (lets say 2x as fast as the slow one)
    //this does the bit reversal in an iterative way.
    // This is a bit hard to draw in ascii art
    //   1 2 3 4
    //    X   X
    //   2 1 4 3
    //    \ X /
    //     X X
    //    / X \
    //   4 3 2 1
    //
    // or maybe this is a better to describe the process
    //
    //  1 2   3 4   5 6   7 8
    //   X     X     X     X
    //  2 1   4 3   6 5   8 7
    //
    //  21 43   65 87
    //    X       X
    //  43 21   87 65
    //
    //  4321 8765
    //      X
    //  8765 4321
    //
    for (quint32 i=0;i<((quint32)nfft);++i)
    {

        quint32 y=i;
        y = (((y & 0xaaaaaaaa) >> 1) | ((y & 0x55555555) << 1));
        y = (((y & 0xcccccccc) >> 2) | ((y & 0x33333333) << 2));
        y = (((y & 0xf0f0f0f0) >> 4) | ((y & 0x0f0f0f0f) << 4));
        y = (((y & 0xff00ff00) >> 8) | ((y & 0x00ff00ff) << 8));
        y = ((y >> 16) | (y << 16)) >> (32-nfft_2power);

        if(y>i)
        {
            swap(x[i], x[y]);
            //qDebug()<<y<<i;
        }

    }

    //fft. using pointers harder to understand but is identical to the one in fft_easy_to_understand
    int nfill=0;
    cpx_type y;
    cpx_type *xkp;
    cpx_type *xlp;
    cpx_type *wp;
    for(int n=1;n<nfft;n<<=1)//run through each stage
    {
        int k=0;
        wp=TWIDDLE+n-1;
        xkp=x;
        xlp=x+k+n;
        while(k<nfft)
        {

            //buterfly
            y=(*wp)*(*xlp);
            (*xlp)=(*xkp)-y;
            (*xkp)+=y;

            //find next pointers for the buterfly
            xkp++;xlp++;k++;
            if(k&nfill)wp++;
             else
             {
                k+=n;
                xkp+=n;xlp+=n;
                wp=TWIDDLE+n-1;
             }

        }
        nfill<<=1;
        nfill|=1;
    }

    //scale if we are doing an inverse
    //this only scalling on the ifft matchs what MATLAB does
    if(fft_direction==INVERSE)
    {
        for(int i=0;i<nfft;++i)x[i]*=(1.0/((double)nfft));
    }

    return true;
}

//this is eaiser to understand
bool JFFT::fft_easy_to_understand(cpx_type *x,int size,fft_direction_t fft_direction)
{
    assert(size==nfft);
    cpx_type *TWIDDLE;
    if(fft_direction==FORWARD)TWIDDLE=TWIDDLE_mem.data();
     else TWIDDLE=TWIDDLE_INV_mem.data();

    //for the ifft an alternitive tick is given at http://www.embedded.com/design/configurable-systems/4210789/DSP-Tricks--Computing-inverse-FFTs-using-the-forward-FFT
    //it would mean taking the conjugate of x before the forward fft is done then taking the conjugate after the fft is done
    //this would be a good idea for memory limited devices

    //bit reversal. slow (say 10% of the CPU use for this function is used here)
    for (int i=0;i<nfft;++i)
    {
        int ti=i;
        int ti2=0;
        for(int j=0;j<nfft_2power;++j)
        {
            ti2<<=1;
            ti2|=(ti&1);
            ti>>=1;
        }
        if(ti2>i)
        {
            swap(x[i], x[ti2]);
            //qDebug()<<ti2<<i;
        }
    }

    //fft. most clear
    //look at the image at http://www.alwayslearn.com/DFT%20and%20FFT%20Tutorial/DFTandFFT_FFT_Butterfly_8_Input.html
    //k points to top part of a buterfly and l points to the bottom part of a buterfly. w points to the twiddle that
    //is currently needed. n is for how many buterfly are in a set, stage 1 has 1 stage 2 has 2 stage 3 has 4 stage 4
    //has 8 and so on.
    for(int n=1;n<nfft;n<<=1)//for for each stage
    {
        int k=0;
        int l=k+n;
        int w=n-1;
        while(k<nfft)
        {

            //qDebug()<<k<<l<<w;

            cpx_type y=x[k]-TWIDDLE[w]*x[l];
            x[k]+=TWIDDLE[w]*x[l];
            x[l]=y;

            k++;l++;
            if(!(k%n))
            {
                k+=n;l+=n;
                w=n-1;
            }
             else w++;
        }
    }

    //scale if we are doing an inverse
    //this only scalling on the ifft matchs what MATLAB does
    if(fft_direction==INVERSE)
    {
        for(int i=0;i<nfft;++i)x[i]*=(1.0/((double)nfft));
    }

    return true;
}

//DFT from definition
//if you were really wanting the best from it you should move Wf to the init function
//but as this is just a rough comparison between a fft and a slow ft implimentation this should do
bool JFFT::sft(cpx_type *x,int size,fft_direction_t fft_direction)
{
    cpx_type imag=cpx_type(0,1);
    cpx_type W;
    F.assign(size,0);
    if(fft_direction==FORWARD)W=exp(-2.0*imag*M_PI/((double)size));
     else W=exp(2.0*imag*M_PI/((double)size));

    //this makes std::pow(W,n*k)==Wf[(n*k)%size] and Wf[(n*k)%size] is faster
    std::vector<cpx_type> Wf;
    Wf.resize(size);
    for(int i=0;i<size;++i)
    {
        Wf[i]=std::pow(W,i);
    }

    for(int n=0;n<size;n++)
    {
        for(int k=0;k<size;k++)
        {
            //F[n]+=x[k]*std::pow(W,n*k);//way way too slow
            F[n]+=x[k]*Wf[(n*k)%size];
        }
    }

    if(fft_direction==INVERSE)
    {
        for(int i=0;i<size;++i)x[i]=F[i]*(1.0/((double)size));
    }
     else
     {
        for(int i=0;i<size;++i)x[i]=F[i];
     }

    return true;
}

#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);


    //JFFT usage example


    //load time data in this case two real sine waves
    QVector<JFFT::cpx_type> x;
    x.resize(16384);//use a power of 2 else JFFT will change it to be so
    for(int i=0;i<x.size();++i)
    {
        x[i]=JFFT::cpx_type(sin(((double)i)/2.0),0);
        x[i]+=JFFT::cpx_type(sin(((double)i)/4.6),0);
    }

    //tranform it to the frequency domain
    fft.fft(x);

    //display frequency information
    //if you have matlab or octave you can copy and paste the resulting text to display the frequency spectrum
    QString str;
    str+="f=[";
    for(int i=0;i<x.size();++i)
    {
        str+=QString::number(x[i].real());
        str+="+";
        str+=QString::number(x[i].imag());
        str+="i";
        if(i+1<x.size())str+=",";
         else str+="];\nplot(abs(f));";

    }
    ui->plainTextEdit->setPlainText(str);

    //tranform it back to the time domain using an inverse FFT
    fft.ifft(x);




}

MainWindow::~MainWindow()
{
    delete ui;
}

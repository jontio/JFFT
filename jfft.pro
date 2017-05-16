#-------------------------------------------------
#
# Project created by QtCreator 2017-05-16T14:18:40
#
#-------------------------------------------------

QT       += core gui

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

QMAKE_LFLAGS_RELEASE -= -O1

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = jfft
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    jfft.cpp

HEADERS  += mainwindow.h \
    jfft.h

FORMS    += mainwindow.ui

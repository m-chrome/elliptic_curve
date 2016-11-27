TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    ec.cpp

HEADERS += \
    ec.hpp

LIBS += -lgcrypt -lgmp

QMAKE_CXXFLAGS += -std=c++1y

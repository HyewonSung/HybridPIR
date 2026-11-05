CC  ?= gcc
CXX ?= g++

ARCHFLAGS ?= -march=native
OMPFLAGS  ?= -fopenmp

CFLAGS   = -no-pie -g3 -Ofast -Wall -DUSE_FFT -DNDEBUG $(ARCHFLAGS)
CXXFLAGS = -std=gnu++11 $(CFLAGS) $(OMPFLAGS)
LDFLAGS  = $(OMPFLAGS)

HEADERS = $(wildcard *.h) $(wildcard spqlios/*.h)

OBJS = \
spqlios/spqlios-fft-fma.o \
spqlios/spqlios-ifft-fma.o \
spqlios/spqlios-fft-impl.o \
spqlios/fft_processor_spqlios.o \
spqlios/lagrangehalfc_impl_fma.o \
hybridpir.o

EXE = hybridpir

%.o: %.s ${HEADERS}
	$(CC) -c ${CFLAGS} -o $@ $<

%.o: %.cpp ${HEADERS}
	$(CXX) -c ${CXXFLAGS} -o $@ $<

all: ${OBJS}
	$(CXX) ${CXXFLAGS} ${OBJS} -o ${EXE} ${LDFLAGS}

clean:
	rm -f ${OBJS} ${EXE}

.PHONY: all clean

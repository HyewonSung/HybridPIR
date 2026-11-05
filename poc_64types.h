#ifndef _POC_64TYPESH_
#define _POC_64TYPESH_

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>

#include "generic_utils.h"

typedef int64_t Torus64;

static const int64_t _two32 = INT64_C(1) << 32;

struct Torus64Polynomial {
    Torus64* const coefs;

    Torus64Polynomial(int N): coefs(new Torus64[N]) {}

    ~Torus64Polynomial() { delete[] coefs; }
};

struct IntPolynomiala {
    int* const coefs;

    IntPolynomiala(int N): coefs(new int[N]) {}

    ~IntPolynomiala() { delete[] coefs; }
};

#ifdef USE_FFT
struct LagrangeHalfCPolynomiala {
    double* const values;
    LagrangeHalfCPolynomiala(int N): values(new double[N]) {}
    ~LagrangeHalfCPolynomiala() {
        delete[] values;
    }
};

#else

struct LagrangeHalfCPolynomiala {
    IntPolynomiala* intPoly;
    Torus64Polynomial* torus64Poly;
    LagrangeHalfCPolynomiala(int N) {
        intPoly=0;
        torus64Poly=0;
    }
    void clear() {
        if (intPoly) { delete intPoly; intPoly=0; }
        if (torus64Poly) { delete torus64Poly; torus64Poly=0;}
    }
    void setIntPoly(const IntPolynomiala* a, const int N) {
        clear();
        intPoly = new IntPolynomiala(N);
        for (int i=0; i<N; i++) intPoly->coefs[i]=a->coefs[i];
    }
    void setTorus64Poly(const Torus64Polynomial* a, const int N) {
        clear();
        torus64Poly = new Torus64Polynomial(N);
        for (int i=0; i<N; i++) torus64Poly->coefs[i]=a->coefs[i];
    }
    void setZeroTorus64Poly(const int N) {
        clear();
        torus64Poly = new Torus64Polynomial(N);
        for (int i=0; i<N; i++) torus64Poly->coefs[i]=0;
    }
    ~LagrangeHalfCPolynomiala() {
    }
};

#endif

struct TLweSample64 {
    Torus64Polynomial* const a;
    Torus64Polynomial* const b;

    TLweSample64(int N): a(new_array1<Torus64Polynomial>(2,N)), b(&a[1]) {}

    ~TLweSample64() { delete_array1<Torus64Polynomial>(a); }
};

struct TLweSampleFFTa {
    LagrangeHalfCPolynomiala* const a;
    LagrangeHalfCPolynomiala* const b;

    TLweSampleFFTa(int N): a(new_array1<LagrangeHalfCPolynomiala>(2,N)), b(&a[1]) {}

    ~TLweSampleFFTa() { delete_array1<LagrangeHalfCPolynomiala>(a); }
};

struct TGswSample64 {
    TLweSample64** const samples;
    TLweSample64* const allsamples;

    TGswSample64(int l, int N):
        samples(new_array2<TLweSample64>(2,l,N)),
        allsamples(samples[0]) {}

    ~TGswSample64() {
        delete_array2<TLweSample64>(samples);
    }
};

struct TGswSampleFFTa {
    TLweSampleFFTa** const samples;
    TLweSampleFFTa* const allsamples;

    TGswSampleFFTa(int l, int N):
        samples(new_array2<TLweSampleFFTa>((2),l,N)),
        allsamples(samples[0]) {}

    ~TGswSampleFFTa() {
        delete_array2<TLweSampleFFTa>(samples);
    }
};

struct GlweCipher64 : public TLweSample64 {
    GlweCipher64(int N) : TLweSample64(N) {}
};

struct GLevCipher64 {
    int l;
    GlweCipher64** cts;

    GLevCipher64(int l_in, int N) : l(l_in) {
        cts = new GlweCipher64*[l];
        for (int i = 0; i < l; ++i)
            cts[i] = new GlweCipher64(N);
    }

    ~GLevCipher64() {
        for (int i = 0; i < l; ++i) delete cts[i];
        delete[] cts;
    }

    GlweCipher64* operator[](int i)       { return cts[i]; }
    const GlweCipher64* operator[](int i) const { return cts[i]; }
};

struct AutoKsKey64 {
    int d;
    GLevCipher64* glev;

    AutoKsKey64(int d_in, int l, int N) : d(d_in) {
        glev = new GLevCipher64(l, N);
    }
    ~AutoKsKey64() {
        delete glev;
    }
};

class Globals {
    public:

        static const int k;
        static const int N;
        static const int t;
        static const int smalln;
        static const int basebit;
        int bgbit;
        int l;

    uint64_t torusDecompOffset;

    uint64_t* torusDecompBuf;
      int* lwekey;
      IntPolynomiala* tlwekey;
        int64_t* in_key;

        TLweSample64**** privKS;

        Globals(int bgbit_in = 9, int l_in = 6);
};

#endif

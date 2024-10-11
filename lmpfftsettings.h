#ifndef LMP_FFT_SETTINGS_H
#define LMP_FFT_SETTINGS_H 
#ifdef FFT_FFTW
#ifndef FFT_FFTW3
#define FFT_FFTW3 
#endif
#endif
#if defined(FFT_FFTW3)
#define LMP_FFT_LIB "FFTW3"
#elif defined(FFT_MKL)
#define LMP_FFT_LIB "MKL FFT"
#elif defined(FFT_CUFFT)
#define LMP_FFT_LIB "cuFFT"
#elif defined(FFT_HIPFFT)
#define LMP_FFT_LIB "hipFFT"
#else
#define LMP_FFT_LIB "KISS FFT"
#endif
#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define FFT_PRECISION 1
#define LMP_FFT_PREC "single"
#define MPI_FFT_SCALAR MPI_FLOAT
#else
typedef double FFT_SCALAR;
#define FFT_PRECISION 2
#define LMP_FFT_PREC "double"
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <raylib.h>
#include <complex.h>
#include <assert.h>

/*
TODO:
- Possibly change name of ft_create (seems vague)
*/

typedef struct {
    size_t amplitude_size;
    float* amplitude;
    float* frequency;
} FreqDomain;

// MEL DATA STRUCT

#define SAMPLE_RATE 44100 
#define BUFFER_SIZE 4096
#define PI 3.14159265358979323846f

void fd_create(FreqDomain* fd, Wave* wave);
void fd_destroy(FreqDomain* fd);
void ft_create(Wave* wave, FreqDomain* fd);
void fft(short in[], double complex out[], size_t n);
void _fft(short in[], double complex out[], size_t n, size_t stride);

int main(void) {
    Wave wave_original = LoadWave("songs/Methods.mp3");
    Wave wave_copy = WaveCopy(wave_original);
    FreqDomain fd;

    WaveFormat(&wave_copy, SAMPLE_RATE, 16, 1);

    fd_create(&fd, &wave_copy);

    ft_create(&wave_copy, &fd);

    fd_destroy(&fd);

    UnloadWave(wave_copy);
    UnloadWave(wave_original);
}


void fd_create(FreqDomain* fd, Wave* wave) {
    fd->amplitude_size = 0;
    fd->amplitude = malloc((wave->frameCount + 1) * sizeof(fd->amplitude));
    fd->frequency = malloc(BUFFER_SIZE * sizeof(fd->frequency));
}

void fd_destroy(FreqDomain* fd) {
    free(fd->amplitude);
    fd->amplitude = NULL;
    free(fd->frequency);
    fd->frequency = NULL;
}

void ft_create(Wave* wave, FreqDomain* fd) {
    short PCM_buffer[BUFFER_SIZE];
    double complex fft_output[BUFFER_SIZE];
    size_t wave_iterator = 0;
    short* wave_data_ptr = wave->data;
    
    for(wave_iterator = 0; wave_iterator < wave->frameCount - BUFFER_SIZE;) {
        for(size_t j = 0; j < BUFFER_SIZE; ++j) {
            PCM_buffer[j] = wave_data_ptr[wave_iterator++];
        }
        fft(PCM_buffer, fft_output, BUFFER_SIZE);
        for(size_t j = 0; j < BUFFER_SIZE; ++j) {
            fd->amplitude[fd->amplitude_size++] = fft_output[j];
        }
    }
    
    // Final buffer will exceed frame count size, so zeros must be put into after
    for(size_t i = 0; i < BUFFER_SIZE; ++i) {
        if (i < wave->frameCount % BUFFER_SIZE) {
            PCM_buffer[i] = wave_data_ptr[wave_iterator++];
        } 
        else  {
            PCM_buffer[i] = 0;
        }
    }
    fft(PCM_buffer, fft_output, BUFFER_SIZE);
    for(size_t i = 0; i < BUFFER_SIZE; ++i) {
        fd->amplitude[fd->amplitude_size++] = fft_output[i];
    }
}


void fft(short in[], double complex out[], size_t n) {
    assert(n % 2 == 0);
    _fft(in, out, BUFFER_SIZE, 1);
}

void _fft(short in[], double complex out[], size_t n, size_t stride) {
    assert(n > 0);

    if (n == 1) { 
        out[0] = in[0];
        return; 
    }

    _fft(in, out, n / 2, stride * 2); 
    _fft(in + stride, out + n/2, n / 2, stride * 2); 
    for (size_t i = 0; i < (n / 2); ++i) { 
        double t = (double) i / n;
        double complex v = cexp(-2* I * PI* t) * out[i + n/2];
        double complex e = out[i]; 
        out[i] = e + v; 
        out[i + (n/2)] = e - v; 
    } 
}
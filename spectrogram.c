#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <raylib.h>
#include <complex.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

#include "arena.c"
#include "double_array.c"

/*
TODO:
- Overlapping signals 
*/

#define SAMPLE_RATE 44100 
#define WINDOW_SIZE 4096
#define BIN_WIDTH WINDOW_SIZE / 2 
#define HOP_SIZE 2048
#define PI 3.14159265358979323846f // Might be redefinition from Raylib
#define PCM16_MAX 32768.0

typedef struct {
    size_t time_slice_index;
    DoubleArray magnitude;
    DoubleArray frequency;
} FreqDomain;


// MEL DATA STRUCT


FreqDomain* FreqDomain_create(Arena* arena, Wave* wave);
void compute_frequency_spectrum(FreqDomain* fd);
void compute_magnitude_spectrum(FreqDomain* fd, Wave* wave);
void FreqDomain_write_buffer_to_file(FreqDomain* fd, size_t bin_index);
void fft(double in[], double complex out[], size_t n);
void _fft(double in[], double complex out[], size_t n, size_t stride);

int main(void) {
    Arena* perm_arena = arena_create(MiB(100));
    Wave wave_original = LoadWave("songs/Methods.mp3");
    Wave wave_copy = WaveCopy(wave_original);
    FreqDomain* fd; 

    WaveFormat(&wave_copy, SAMPLE_RATE, 16, 1);

    fd = FreqDomain_create(perm_arena, &wave_copy);
    
    compute_magnitude_spectrum(fd, &wave_copy);
    compute_frequency_spectrum(fd);
    FreqDomain_write_buffer_to_file(fd, fd->time_slice_index / 2);

    arena_destroy(perm_arena);

    UnloadWave(wave_copy);
    UnloadWave(wave_original);
}

FreqDomain* FreqDomain_create(Arena* arena, Wave* wave) {
    FreqDomain *fd = arena_push(arena, sizeof(FreqDomain), true);
    fd->time_slice_index = 0;
    fd->magnitude.length = 0;
    fd->magnitude.capacity = wave->frameCount + BIN_WIDTH; 
    fd->magnitude.items = arena_push(arena, fd->magnitude.capacity * sizeof(*fd->magnitude.items), true);
    fd->frequency.length = 0;
    fd->frequency.capacity = BIN_WIDTH;
    fd->frequency.items = arena_push(arena, fd->frequency.capacity * sizeof(*fd->frequency.items), true);
    
    return fd;
}

void compute_frequency_spectrum(FreqDomain* fd) {
    for(size_t i = 0; i < BIN_WIDTH; ++i) {
        DoubleArray_push(&fd->frequency, i * ((double)SAMPLE_RATE / WINDOW_SIZE)); 
    }
}

void compute_magnitude_spectrum(FreqDomain* fd, Wave* wave) {
    double PCM_buffer_windowed[WINDOW_SIZE];
    double complex fft_output[WINDOW_SIZE];
    double complex fft_output_normalized[WINDOW_SIZE];
    double hamming_window[WINDOW_SIZE];
    size_t wave_iterator = 0;
    short* wave_data_ptr = wave->data;
    
    // Create hamming window data (https://www.sciencedirect.com/topics/engineering/hamming-window)
    for(size_t i = 0; i < WINDOW_SIZE; i++) {
        hamming_window[i] = 0.54 - 0.46 * cos((2 * PI * i) / (WINDOW_SIZE - 1));
    }
    
    // Grab buffer data from song file
    for(wave_iterator = 0; wave_iterator < wave->frameCount - WINDOW_SIZE; wave_iterator+=HOP_SIZE) {
        for(size_t i = 0; i < WINDOW_SIZE; ++i) {
            double t = wave_data_ptr[wave_iterator + i] / PCM16_MAX; // Normalize data before storing
            PCM_buffer_windowed[i] = t * hamming_window[i];
        }
        fft(PCM_buffer_windowed, fft_output, WINDOW_SIZE);
        // Normalize Window Output
        for(size_t i = 0; i < BIN_WIDTH; ++i) {
            fft_output_normalized[i] = fft_output[i] / sqrt(WINDOW_SIZE);
            DoubleArray_push(&fd->magnitude, cabs(fft_output_normalized[i]));
        }
        fd->time_slice_index++;
    }
}

// void compute_frame_magnitude

void FreqDomain_write_buffer_to_file(FreqDomain* fd, size_t bin_index) {
    FILE *file;
    
    file = fopen("debug/sample.txt", "w");
    
    if (file ==  NULL) {
       perror("Error Opening File");
       return;
    }
    
    for (size_t i = 0; i < BIN_WIDTH; ++i) {
        fprintf(file, "%zu: (%lfHz, %lf)\n", i, DoubleArray_get(fd->frequency, i), DoubleArray_get(fd->magnitude, (bin_index * BIN_WIDTH) + i));
    }

    fclose(file);
    return;
}

void fft(double in[], double complex out[], size_t n) {
    assert(n % 2 == 0);
    _fft(in, out, n, 1);
}

void _fft(double in[], double complex out[], size_t n, size_t stride) {
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
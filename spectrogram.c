#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <raylib.h>
#include <complex.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/*
TODO:
- Segmentation 
- Windowing Signal
*/

#define SAMPLE_RATE 44100 
#define BUFFER_SIZE 4096
#define HALF_BUFFER_SIZE (BUFFER_SIZE / 2)
#define PI 3.14159265358979323846f

typedef struct {
    size_t amplitude_size;
    size_t sample_size;
    double* amplitude;
    double* frequency;
} FreqDomain;

// MEL DATA STRUCT


void fd_create(FreqDomain* fd, Wave* wave);
void fd_destroy(FreqDomain* fd);
void fd_set_frequency(FreqDomain* fd);
void fd_set_amplitude(FreqDomain* fd, Wave* wave);
void fd_write_buffer_to_file(FreqDomain* fd, size_t sample_index);
void fft(short in[], double complex out[], size_t n);
void _fft(short in[], double complex out[], size_t n, size_t stride);

int main(void) {
    Wave wave_original = LoadWave("songs/Methods.mp3");
    Wave wave_copy = WaveCopy(wave_original);
    FreqDomain fd;

    WaveFormat(&wave_copy, SAMPLE_RATE, 16, 1);

    fd_create(&fd, &wave_copy);

    fd_set_amplitude(&fd, &wave_copy);
    fd_set_frequency(&fd);
    fd_write_buffer_to_file(&fd, fd.sample_size / 2);

    fd_destroy(&fd);

    UnloadWave(wave_copy);
    UnloadWave(wave_original);
}


void fd_create(FreqDomain* fd, Wave* wave) {
    fd->amplitude_size = 0;
    fd->sample_size = 0;
    fd->amplitude = malloc((wave->frameCount + HALF_BUFFER_SIZE) * sizeof(fd->amplitude));
    fd->frequency = malloc(HALF_BUFFER_SIZE * sizeof(fd->frequency));
}

void fd_destroy(FreqDomain* fd) {
    free(fd->amplitude);
    fd->amplitude = NULL;
    free(fd->frequency);
    fd->frequency = NULL;
}

void fd_set_frequency(FreqDomain* fd) {
    for(size_t i = 0; i < HALF_BUFFER_SIZE; i++) {
        fd->frequency[i] = ((double)i / HALF_BUFFER_SIZE) * SAMPLE_RATE;
    }
}

void fd_set_amplitude(FreqDomain* fd, Wave* wave) {
    short PCM_buffer[BUFFER_SIZE];
    double complex fft_output[BUFFER_SIZE];
    size_t wave_iterator = 0;
    short* wave_data_ptr = wave->data;
    double debug_amplitude[HALF_BUFFER_SIZE];
    
    for(wave_iterator = 0; wave_iterator < wave->frameCount - BUFFER_SIZE;) {
        for(size_t j = 0; j < BUFFER_SIZE; ++j) {
            PCM_buffer[j] = wave_data_ptr[wave_iterator++];
        }
        fft(PCM_buffer, fft_output, BUFFER_SIZE);
        for(size_t j = 0; j < HALF_BUFFER_SIZE; ++j) {
            fd->amplitude[fd->amplitude_size++] = cabs(fft_output[j]);
            debug_amplitude[j] = cabs(fft_output[j]);
        }
        fd->sample_size++;
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
    for(size_t i = HALF_BUFFER_SIZE - 1; i < BUFFER_SIZE; ++i) {
        fd->amplitude[fd->amplitude_size++] = cabs(fft_output[i]);
    }
    fd->sample_size++; // Increase sample by one for final buffer
    fd += 0;
}

void fd_write_buffer_to_file(FreqDomain* fd, size_t sample_index) {
    FILE *file;
    
    struct stat st = {0};
    
    file = fopen("debug/sample.txt", "w");
    
    if (file ==  NULL) {
       perror("Error Opening File");
       return;
    }
    
    for (size_t i = 0; i < HALF_BUFFER_SIZE; ++i) {
        fprintf(file, "%zu: (%lfHz, %lf)\n", i, fd->frequency[i], fd->amplitude[(sample_index * HALF_BUFFER_SIZE) + i]);
    }

    // for (size_t i = 0; i < HALF_BUFFER_SIZE; ++i) {
    //     fprintf(file, "(%zu, %lfHz): (%lf, %lf)\n", i, fd->frequency[HALF_BUFFER_SIZE - i - 1], fd->amplitude[(HALF_BUFFER_SIZE * sample_index) + i], fd->amplitude[(HALF_BUFFER_SIZE * (sample_index + 1)) - i]);
    // }
    
    fclose(file);
    return;
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
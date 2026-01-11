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
- Overlapping signals 
*/

#define SAMPLE_RATE 44100 
#define FRAME_SIZE 4096
#define HALF_FRAME_SIZE (FRAME_SIZE / 2)
#define HOP_LENGTH 441
#define PI 3.14159265358979323846f // Might be redefinition from Raylib
#define PCM16_MAX 32767.0

typedef struct {
    double* items; 
    size_t length;
    size_t capacity;
} DoubleArray;


typedef struct {
    size_t sample_size;
    DoubleArray magnitude;
    DoubleArray frequency;
} FreqDomain;

// MEL DATA STRUCT

double DoubleArray_get(DoubleArray array, size_t index);
void DoubleArray_push(DoubleArray* array, double value);

void FreqDomain_create(FreqDomain* fd, Wave* wave);
void FreqDomain_destroy(FreqDomain* fd);
void FreqDomain_frequency_set(FreqDomain* fd);
void FreqDomain_magnitude_set(FreqDomain* fd, Wave* wave);
void FreqDomain_write_buffer_to_file(FreqDomain* fd, size_t sample_index);
void fft(double in[], double complex out[], size_t n);
void _fft(double in[], double complex out[], size_t n, size_t stride);

int main(void) {
    Wave wave_original = LoadWave("songs/Methods.mp3");
    Wave wave_copy = WaveCopy(wave_original);
    FreqDomain fd;

    WaveFormat(&wave_copy, SAMPLE_RATE, 16, 1);

    FreqDomain_create(&fd, &wave_copy);

    FreqDomain_magnitude_set(&fd, &wave_copy);
    FreqDomain_frequency_set(&fd);
    FreqDomain_write_buffer_to_file(&fd, fd.sample_size / 2);

    FreqDomain_destroy(&fd);

    UnloadWave(wave_copy);
    UnloadWave(wave_original);
}


void FreqDomain_create(FreqDomain* fd, Wave* wave) {
    fd->sample_size = 0;
    fd->magnitude.length = 0;
    fd->magnitude.capacity = wave->frameCount + HALF_FRAME_SIZE;
    fd->magnitude.items = malloc(fd->magnitude.capacity * sizeof(*fd->magnitude.items));
    fd->frequency.capacity = HALF_FRAME_SIZE;
    fd->frequency.items = malloc(fd->frequency.capacity * sizeof(*fd->frequency.items));
}

void FreqDomain_destroy(FreqDomain* fd) {
    free(fd->magnitude.items);
    fd->magnitude.items = NULL;
    free(fd->frequency.items);
    fd->frequency.items = NULL;
}

void FreqDomain_frequency_set(FreqDomain* fd) {
    for(size_t i = 0; i < HALF_FRAME_SIZE; ++i) {
        DoubleArray_push(&fd->frequency, i * ((double)SAMPLE_RATE / FRAME_SIZE));
    }
}

void FreqDomain_magnitude_set(FreqDomain* fd, Wave* wave) {
    short PCM_buffer[FRAME_SIZE];
    double PCM_buffer_normalized[FRAME_SIZE];
    double PCM_buffer_normalized_windowed[FRAME_SIZE];
    double complex fft_output[FRAME_SIZE];
    double hamming_window[FRAME_SIZE];
    size_t wave_iterator = 0;
    short* wave_data_ptr = wave->data;
    
    // Create hamming window data (https://www.sciencedirect.com/topics/engineering/hamming-window)
    for(size_t i = 0; i < FRAME_SIZE; i++) {
        hamming_window[i] = 0.54 - 0.46 * cos((2 * PI * i) / (FRAME_SIZE - 1));
    }
    
    // Grab buffer data from song file
    for(wave_iterator = 0; wave_iterator < wave->frameCount - FRAME_SIZE;) {
        for(size_t j = 0; j < FRAME_SIZE; ++j) {
            PCM_buffer[j] = wave_data_ptr[wave_iterator++];
            PCM_buffer_normalized[j] = PCM_buffer[j] / PCM16_MAX;
            PCM_buffer_normalized_windowed[j] = PCM_buffer_normalized[j] * hamming_window[j];
        }
        fft(PCM_buffer_normalized_windowed, fft_output, FRAME_SIZE);
        for(size_t j = HALF_FRAME_SIZE - 1; j < FRAME_SIZE; ++j) {
            DoubleArray_push(&fd->magnitude, cabs(fft_output[j]));
        }
        fd->sample_size++;
    }
    
    // Final buffer will exceed frame count size, so zeros must be put into after
    for(size_t i = 0; i < FRAME_SIZE; ++i) {
        if (i < wave->frameCount % FRAME_SIZE) {
            PCM_buffer[i] = wave_data_ptr[wave_iterator++];
            PCM_buffer_normalized[i] = PCM_buffer[i] / PCM16_MAX;
            PCM_buffer_normalized_windowed[i] = PCM_buffer_normalized[i] * hamming_window[i];
        } 
        else  {
            PCM_buffer[i] = 0;
        }
    }
    fft(PCM_buffer_normalized_windowed, fft_output, FRAME_SIZE);
    for(size_t i = HALF_FRAME_SIZE - 1; i < FRAME_SIZE; ++i) {
        DoubleArray_push(&fd->magnitude, cabs(fft_output[i]));
    }
    fd->sample_size++; // Increase sample by one for final buffer
    fd += 0;
}

void FreqDomain_write_buffer_to_file(FreqDomain* fd, size_t sample_index) {
    FILE *file;
    
    file = fopen("debug/sample.txt", "w");
    
    if (file ==  NULL) {
       perror("Error Opening File");
       return;
    }
    
    for (size_t i = 0; i < HALF_FRAME_SIZE; ++i) {
        fprintf(file, "%zu: (%lfHz, %lf)\n", i, DoubleArray_get(fd->frequency, i), DoubleArray_get(fd->magnitude, (sample_index * HALF_FRAME_SIZE) + i));
    }

    // for (size_t i = 0; i < HALF_FRAME_SIZE; ++i) {
    //     fprintf(file, "(%zu, %lfHz): (%lf, %lf)\n", i, fd->frequency[HALF_FRAME_SIZE - i - 1], fd->magnitude[(HALF_FRAME_SIZE * sample_index) + i], fd->magnitude[(HALF_FRAME_SIZE * (sample_index + 1)) - i]);
    // }
    
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

double DoubleArray_get(DoubleArray array, size_t index) {
    if (index < array.length) {
        return array.items[index];
    }
    return 0;
}

void DoubleArray_push(DoubleArray* array, double value) {
    if (array->length < array->capacity) {
        array->items[array->length++] = value;
        return;
    }
    return;
}
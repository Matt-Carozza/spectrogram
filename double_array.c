#include <stdio.h>

typedef struct {
    double* items; 
    size_t length;
    size_t capacity;
} DoubleArray;

double DoubleArray_get(DoubleArray array, size_t index);
void DoubleArray_push(DoubleArray* array, double value);

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
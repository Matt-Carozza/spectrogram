#include <stdio.h>
#include <stdint.h>

typedef int8_t  i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef i8  b8;
typedef i32 b32;

typedef struct { 
   u64 capacity; // Implementation assumes we start after struct pointer
   u64 position;
   u8 data[];
} Arena;

#define ARENA_ALIGN (sizeof(void*))
static inline u64 KiB(u64 size);
static inline u64 MiB(u64 size);
static inline u64 GiB(u64 size);
static inline u64 min(u64 x, u64 y);
static inline u64 max(u64 x, u64 y);
Arena* arena_create(u64 capacity);
void arena_destroy(Arena* arena);
void* arena_push(Arena* arena, u64 size, b32 non_zero);
void arena_pop_to(Arena* arena, u64 pos);
void arena_clear(Arena* arena);
static inline u64 align_up_pow2(u64 n, u64 p);


static inline u64 KiB(u64 size) {
   return (size << 10);
}

static inline u64 MiB(u64 size) {
   return (size << 20);
}

static inline u64 GiB(u64 size) {
   return (size << 30);
}

/**
 * Aligns 'n' up to the next multiple of 'p'.
 * 'p' MUST be a power of two.
 */
static inline u64 align_up_pow2(u64 n, u64 p) {
    return (n + (p - 1)) & ~(p - 1);
}

static inline u64 min(u64 x, u64 y) {
    return (x < y) ? x : y;
}

static inline u64 max(u64 x, u64 y) {
    return (x > y) ? x : y;
}



Arena* arena_create(u64 capacity) {
    Arena* arena = malloc(sizeof(*arena) + sizeof(u8) * capacity);
    
    arena->capacity = capacity;
    arena->position = 0;
    
    return arena;
}

void arena_destroy(Arena* arena) {
    free(arena);
}

void* arena_push(Arena* arena, u64 size, b32 zero) {
    /*
     Ensure allignment is good for the types we store (multiple of two dependent on system)
    */
    u64 pos_alligned = align_up_pow2(arena->position, ARENA_ALIGN); 

    if (size > arena->capacity - pos_alligned) { 
        return NULL; // Maybe add error message
    } 
    u64 new_pos      = pos_alligned + size;
    arena->position = new_pos;
    
    u8* out = arena->data + pos_alligned;
    
    // zero out the data we are allocating
    if (zero) {
        memset(out, 0, size);
    }
    
    return out;
}

void arena_pop_to(Arena* arena, u64 pos) {
    if (pos > arena->position) {
        pos = arena->position;
    }
    arena->position = pos;
}

void arena_clear(Arena* arena) {
    arena_pop_to(arena, 0);
}

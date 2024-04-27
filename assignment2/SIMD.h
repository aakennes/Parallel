#include "params.h"

#include<immintrin.h>
#include<stdint.h>

#define Vec(sz, T) __attribute((vector_size(sz))) T
#define RC(T, x) reinterpret_cast<T>(x)

#define i32 int32_t
#define u32 uint32_t
#define i64 int64_t
#define u64 uint64_t

#define i32x8 Vec(32, i32)
#define u32x8 Vec(32, u32)
#define i64x4 Vec(32, i64)
#define u64x4 Vec(32, u64)

u32x8 setu32x8(u32 x);
// u32x8 shuffle(const u32x8 &x);
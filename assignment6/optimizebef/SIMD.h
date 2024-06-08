#include<cstdio>
#include<iostream>
#include<cstring>


#include "Math.h"


#define MMint u32

#define P 104857601

#define R 104857599
#define r2 45971250


#define Vec(sz, T) __attribute((vector_size(sz))) T

#define i32x8 Vec(32, i32)
#define u32x8 Vec(32, u32)
#define i64x4 Vec(32, i64)
#define u64x4 Vec(32, u64)

// #define Mintx8 

static  u32 get_r(MMint v);

static  u32 pow_mod(u32 x, u64 y);

static  u32 reduce(u64 x);

MMint In(u32 v);

u32 get(MMint v);

MMint subM(MMint v);

MMint add(MMint v,MMint rhs);

MMint sub(MMint v,MMint rhs);

MMint mul(MMint v,MMint rhs);

MMint Pow(u32 x,i64 y);

MMint inv(MMint v);

MMint div(MMint v,MMint rhs);

MMint neg(MMint v);

u32x8 setu32x8(u32 x);

// #undef P
// #undef R
// #undef r2
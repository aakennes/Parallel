#include "SIMD.h"

u32x8 setu32x8(u32 x){return (u32x8){x, x, x, x, x, x, x, x};}
// template<int typ>u32x8 shuffle(const u32x8 &x){return RC(u32x8, _mm256_shuffle_epi32(RC(__m256i, x), typ));}
// template<int typ>u32x8 blend(const u32x8 &x, const u32x8 &y){return RC(u32x8, _mm256_blend_epi32(RC(__m256i, x), RC(__m256i, y), typ));}
// __m256i swaplohi128(const __m256i &x){return _mm256_permute2x128_si256(x, x, 1);}	
// u32x8& x8(u32 *data){return *RC(u32x8* ,data);}
// const u32x8& x8(const u32 *data){return *RC(const u32x8*, data);}
// __m256i loadu(const void* data){return _mm256_loadu_si256(RC(const __m256i_u*, data));}
// void storeu(const __m256i &x, void* data){return _mm256_storeu_si256(RC(__m256i_u*, data), x);}
// u64x4 mulu32x8_fus(const u32x8 &x, const u32x8 &y){return RC(u64x4, _mm256_mul_epu32(RC(__m256i, x), RC(__m256i, y)));}
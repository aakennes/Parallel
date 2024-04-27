#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>

#define MOD 998244353

// Helper function to print __m256i
void print256_num(__m256i var) {
    uint32_t val[8];
    _mm256_storeu_si256((__m256i*)val, var);
    printf("Values: %u %u %u %u %u %u %u %u\n", val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]);
}

__m256i mod_mul(__m256i a, __m256i b) {
    __m256i result = _mm256_mul_epu32(a, b);

    // 计算结果的模
    __m256i mod = _mm256_set1_epi64x(MOD);
    __m256i quotient = _mm256_srli_epi64(result, 31); // 高 33 位移位到低 33 位
    result = _mm256_add_epi64(result, quotient);
    result = _mm256_sub_epi64(result, mod);
    result = _mm256_add_epi64(result, _mm256_and_si256(_mm256_cmpgt_epi64(result, mod), mod));
    return result;
}

int main() {
    // Example values
    __m256i a = _mm256_set1_epi32(123456789);
    __m256i b = _mm256_set1_epi32(987654321);
    print256_num(a);
    // Perform modulo multiplication
    __m256i result = mod_mul(a, b);
    
    // Print results
    print256_num(result);

    return 0;
}

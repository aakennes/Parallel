#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <map>
#include<immintrin.h>
#include<sys/time.h>
#include<chrono>
#include"../params.h"
#define i32 int32_t
#define u32 uint32_t
#define i64 int64_t
#define u64 uint64_t
using u128 = __uint128_t;
namespace hehub {

// Utility function for modular exponentiation
u32 __pow_mod(u32 modulus, u32 base, size_t index) {
    u32 power = 1;
    size_t mask = 1;
    while (mask <= index) {
        mask <<= 1;
    }
    mask >>= 1;
    while (mask) {
        power = (u64)power * power % modulus;
        if (mask & index) {
            power = (u64)power * base % modulus;
        }
        mask >>= 1;
    }
    return power;
}

// Utility function to get 2^n-th root of unity
u32 __get_2nth_unity_root(u32 modulus, u32 n) {
    if ((modulus - 1) % (2 * n) != 0) {
        throw std::invalid_argument("2N doesn't divide (modulus - 1)");
    }
    u32 candidate = 2;
    while (true) {
        if (__pow_mod(modulus, candidate, (modulus - 1) / 2) == modulus - 1) {
            break;
        }
        candidate++;
    }
    return __pow_mod(modulus, candidate, (modulus - 1) / (2 * n));
}
inline u32 __bit_rev_naive_16(u32 x, int bit_len) {
    x = ((x & 0xFF00FF00) >> 8) | ((x & 0x00FF00FF) << 8);
    x = ((x & 0xF0F0F0F0) >> 4) | ((x & 0x0F0F0F0F) << 4);
    x = ((x & 0xCCCCCCCC) >> 2) | ((x & 0x33333333) << 2);
    x = ((x & 0xAAAAAAAA) >> 1) | ((x & 0x55555555) << 1);
    return x >> (16 - bit_len);
}
// Struct for NTT factors
struct NTTFactors {
    NTTFactors(u32 modulus, size_t log_dimension, bool for_inverse = false) {
        const size_t log_modulus = (u32)(log2(modulus) + 0.5);
        if (log_modulus > 59) {
            throw std::invalid_argument("NTT not supporting primes with bit size > 59 currently.");
        }
        size_t dimension = 1 << log_dimension;
        const u32 root_of_2nth = __get_2nth_unity_root(modulus, dimension);
		// printf("%lu\n",root_of_2nth);
        if (!for_inverse) {
            seq.resize(dimension);
            for (size_t i = 0; i < dimension; i++) {
                seq[i] = __pow_mod(modulus, root_of_2nth, __bit_rev_naive_16(i, log_dimension));
            }
        } else {
            seq.resize(dimension * 2);
            auto root_of_2nth_inv = __pow_mod(modulus, root_of_2nth, 2 * dimension - 1);
			// printf("%lu\n",root_of_2nth_inv);
            for (size_t l = 0; l < log_dimension; l++) {
                auto start = (1 << l) - 1;
                auto power_index_factor = 1 << (log_dimension - l);
                for (size_t i = 0; i < (1 << l); i++) {
                    auto idx = start + i;
                    seq[idx] = __pow_mod(modulus, root_of_2nth_inv, __bit_rev_naive_16(i, l) * power_index_factor);
                }
            }
            const u32 dimension_inv = modulus - ((modulus - 1) >> log_dimension);
            for (size_t i = 0; i < dimension; i++) {
                auto temp = __pow_mod(modulus, root_of_2nth_inv, i);
                auto idx = i + dimension;
                seq[idx] = (u64)temp * dimension_inv % modulus;
                seq[idx] -= (seq[idx] >= modulus) ? modulus : 0;
            }
            shuffled_indices.resize(dimension);
            for (size_t i = 0; i < dimension; i++) {
                shuffled_indices[i] = __bit_rev_naive_16(i, log_dimension);
            }
        }
    }

    NTTFactors(const NTTFactors &copying) = default;
    NTTFactors(NTTFactors &&moving) = default;
    ~NTTFactors() {}

    std::vector<u32> seq;
    std::vector<size_t> shuffled_indices;
};

// Cache for NTT factors
std::map<std::pair<u32, u32>, NTTFactors>& ntt_factors_cache() {
    static std::map<std::pair<u32, u32>, NTTFactors> global_ntt_factors_cache;
    return global_ntt_factors_cache;
}

std::map<std::pair<u32, u32>, NTTFactors>& intt_factors_cache() {
    static std::map<std::pair<u32, u32>, NTTFactors> global_intt_factors_cache;
    return global_intt_factors_cache;
}

// Find or create NTT factors
inline const auto& __find_or_create_ntt_factors(const u32 modulus, const size_t log_dimension, const bool for_inverse = false) {
    if (!for_inverse) {
        auto it = ntt_factors_cache().find(std::make_pair(modulus, log_dimension));
        if (it == ntt_factors_cache().end()) {
            ntt_factors_cache().insert(std::make_pair(std::make_pair(modulus, log_dimension), NTTFactors(modulus, log_dimension)));
            it = ntt_factors_cache().find(std::make_pair(modulus, log_dimension));
        }
        return it->second;
    } else {
        auto it = intt_factors_cache().find(std::make_pair(modulus, log_dimension));
        if (it == intt_factors_cache().end()) {
            intt_factors_cache().insert(std::make_pair(std::make_pair(modulus, log_dimension), NTTFactors(modulus, log_dimension, true)));
            it = intt_factors_cache().find(std::make_pair(modulus, log_dimension));
        }
        return it->second;
    }
}

// NTT Negacyclic inplace lazy
void ntt_negacyclic_inplace_lazy(const size_t log_dimension, const u32 modulus, u32 coeffs[]) {
    const size_t dimension = 1ULL << log_dimension;
    const auto& ntt_factors = __find_or_create_ntt_factors(modulus, log_dimension);
    size_t level, start, data_step, h, l;
    size_t idx = 1;
    u32 temp, zeta;
    for (level = 1, data_step = dimension; level <= log_dimension; level++, data_step >>= 1) {
        auto gap = data_step / 2;
        for (start = 0; start < dimension; start += data_step, idx++) {
            zeta = ntt_factors.seq[idx];
            for (l = start; l < start + gap; l++) {
                h = l + gap;
                temp = (u64)coeffs[h] * zeta % modulus;
                coeffs[h] = coeffs[l] + 2 * modulus - temp;
                coeffs[l] = coeffs[l] + temp;
            }
        }
    }
    const u32 log_modulus = (u32)(log2(modulus) + 0.5);
    const u32 div_fix = (modulus >= (1ULL << log_modulus)) ? 1 : 0;
    for (size_t i = 0; i < dimension; i++) {
        coeffs[i] -= ((coeffs[i] >> log_modulus) - div_fix) * modulus;
    }
}

// INTT Negacyclic inplace lazy
void intt_negacyclic_inplace_lazy(const size_t log_dimension, const u32 modulus, u32 values[]) {
    const size_t dimension = 1ULL << log_dimension;
    const auto& intt_factors = __find_or_create_ntt_factors(modulus, log_dimension, true);
    u32 values_shuffled[dimension];
    const auto& shuffled_indices = intt_factors.shuffled_indices;
    for (size_t i = 0; i < dimension; i++) {
        values_shuffled[i] = values[shuffled_indices[i]];
    }
    size_t level, start, data_step, h, l;
    size_t idx = 0;
    u32 temp, zeta;
    for (level = 1, data_step = dimension; level <= log_dimension; level++, data_step >>= 1) {
        auto gap = data_step / 2;
        for (start = 0; start < dimension; start += data_step, idx++) {
            zeta = intt_factors.seq[idx];
            for (l = start; l < start + gap; l++) {
                h = l + gap;
                temp = (u64)values_shuffled[h] * zeta % modulus;
                values_shuffled[h] = values_shuffled[l] + 2 * modulus - temp;
                values_shuffled[l] = values_shuffled[l] + temp;
            }
        }
    }
    for (size_t i = 0; i < dimension; i++) {
        values[i] = values_shuffled[shuffled_indices[i]];
    }
    const u32 log_modulus = (u32)(log2(modulus) + 0.5);
    const u32 div_fix = (modulus >= (1ULL << log_modulus)) ? 1 : 0;
    idx++;
    for (size_t i = 0; i < dimension; i++, idx++) {
        values[i] -= ((values[i] >> log_modulus) - div_fix) * modulus;
        values[i] = (u64)values[i] * intt_factors.seq[idx] % modulus;
    }
}

// Function to perform polynomial multiplication
std::vector<u32> multiply_polynomials(const std::vector<u32>& poly1, const std::vector<u32>& poly2, u32 modulus) {
    size_t n = 1;
    size_t log_dimension = 0;
    while (n < poly1.size() + poly2.size() - 1) {
        n <<= 1;
        log_dimension++;
    }

    std::vector<u32> a(n, 0);
    std::vector<u32> b(n, 0);

    for (size_t i = 0; i < poly1.size(); i++) {
        a[i] = poly1[i];
    }
    for (size_t i = 0; i < poly2.size(); i++) {
        b[i] = poly2[i];
    }

    ntt_negacyclic_inplace_lazy(log_dimension, modulus, a.data());
    ntt_negacyclic_inplace_lazy(log_dimension, modulus, b.data());

    for (size_t i = 0; i < n; i++) {
        a[i] = (u64)a[i] * b[i] % modulus;
    }

    intt_negacyclic_inplace_lazy(log_dimension, modulus, a.data());

    return a;
}

} // namespace hehub
int nn[11]={256,512,1024,2048,4096,8192,16384,32768,65536,131072};
int qq[5]={1409,3329,7681,12289};
int f[131075],g[131075];
int main() {
    // Test polynomials
    freopen("1.out","w",stdout);
    u32 modulus = Original_Q;

    for(int i=0;i<=9;++i){
        for(int j=0;j<=3;++j){
			long double ans=0;
            int cnt=100;
			
			std::string str1="../data/datain/NTT";
			std::string str2=std::to_string(nn[i]);
			std::string str3=std::to_string(qq[i]);
			std::string strin=str1+str2+"_"+str3+".in";
			char charArrayin[strin.size() + 1];
			std::copy(strin.begin(), strin.end(), charArrayin);
			charArrayin[strin.size()] = '\0';
			freopen(charArrayin,"r",stdin);
			int n=nn[i],m=n;
            int limit=1;
            while(limit < m + n - 2){
                limit <<= 1;
            }
			for(int i = 0; i < n; ++i){std::cin >> f[i];}
			for(int i = 0; i < m; ++i){std::cin >> g[i];}
			for(int k=1;k<=cnt;++k){
				auto Start=std::chrono::high_resolution_clock::now();
                std::vector<u32> poly1;
                std::vector<u32> poly2;
                poly1.resize(limit);
                poly2.resize(limit);
                for(int i = 0; i < n; ++i){poly1[i]=f[i];}
                for(int i = 0; i < m; ++i){poly2[i]=g[i];}
				std::vector<u32> result = hehub::multiply_polynomials(poly1, poly2, modulus);
                auto End=std::chrono::high_resolution_clock::now();
				std::chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
				ans+=elapsed.count();
			}
            std::cout<<ans/cnt<<std::endl;
        }
    }
    return 0;
}

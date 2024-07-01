#include "ntt.h"
#include "mod_arith.h"
#include "permutation.h"
#include <cmath>
#include <map>
// #include <omp.h>

int qpow(int a,int x,int p){
    int ans=1;
    while(x){
        if(x&1)ans=1LL*ans*a%p;
        a=1LL*a*a%p;
        x>>=1;
    }
    return ans%p;
}

#define i32 int32_t
#define u32 uint32_t
#define i64 int64_t
#define u64 uint64_t

namespace hehub {

//common
u64 __pow_mod(u64 modulus, u64 base, size_t index) {
    u64 power = 1;
    size_t mask = 1;
    while (mask <= index) {
        mask <<= 1;
    }
    mask >>= 1;
    while (mask) {
        power = (u128)power * power % modulus;
        if (mask & index) {
            power = (u128)power * base % modulus;
        }
        mask >>= 1;
    }
    return power;
}

u64 __get_2nth_unity_root(u64 modulus, u64 n) {
    if ((modulus - 1) % (2 * n) != 0) {
        throw std::invalid_argument("2N doesn't divide (modulus - 1)");
    }
    u64 candidate = 2;
    while (true) {
        if (__pow_mod(modulus, candidate, (modulus - 1) / 2) == modulus - 1) {
            break;
        }
        candidate++;
    }
    auto root = __pow_mod(modulus, candidate, (modulus - 1) / (2 * n));
    return root;
}

u64 __get_nth_unity_root(u64 modulus, u64 n) {
    if ((modulus - 1) % n != 0) {
        throw std::invalid_argument("N doesn't divide (modulus - 1)");
    }
    u64 candidate = 2;
    while (true) {
        if (__pow_mod(modulus, candidate, (modulus - 1)/2) == modulus - 1) {
            break;
        }
        candidate++;
    }
    return candidate;
}

struct NTTFactors {
    NTTFactors(u64 modulus, size_t log_dimension, bool for_inverse = false) {
        const size_t log_modulus = (u64)(log2(modulus) + 0.5);
        if (log_modulus > 59) {
            throw std::invalid_argument(
                "NTT not supporting primes with bit size > 59 currently.");
        }
        size_t dimension = 1 << log_dimension;
        
        const u64 wnn = __get_nth_unity_root(modulus, dimension);
        auto invwnn =
                __pow_mod(modulus, wnn, modulus-2);
        auto temp = __get_2nth_unity_root(modulus,dimension);
        const u64 root_of_2nth = temp;
        if(!for_inverse){
            rt.resize(dimension);
            rt_harvey.resize(dimension);
            rt[0] = 1;
            rt[1] = temp;
            auto rt_harvey_temp = ((u128)temp << 64) / modulus;
            for (size_t i = 2; i < dimension; i++) rt[i] = mul_mod_harvey_lazy(modulus,rt[i-1],temp,rt_harvey_temp);
            // for (size_t i = 2; i < dimension; i++) rt[i] = (u128)rt[i-1] * rt[1] % modulus;
            for (size_t i = 0, j = 0; i < dimension; i++) {
                if (i > j) {
                    std::swap(rt[i], rt[j]);
                }
                for (size_t t = dimension >> 1; (j ^= t) < t; t >>= 1);
            }
            for (size_t i = 0; i < dimension; i++) rt_harvey[i] = ((u128)rt[i] << 64) / modulus;
        }else{
            irt.resize(dimension);
            irt_harvey.resize(dimension);
            irt[0] = 1;
            irt[1] = __pow_mod(modulus, temp, 2 * dimension - 1);
            auto irt_harvey_temp = ((u128)irt[1] << 64) / modulus;
            // for (size_t i = 2; i < dimension; i++) irt[i] = (u128)irt[i-1] * irt[1] % modulus;
            for (size_t i = 2; i < dimension; i++) irt[i] = mul_mod_harvey_lazy(modulus,irt[i-1],irt[1],irt_harvey_temp);;
            // for (int i = 2; i < limit; i++) printf("%d ", irt[i]);
            for (size_t i = 0, j = 0; i < dimension; i++) {
                if (i > j) {
                    std::swap(irt[i], irt[j]);
                }
                
                for (size_t t = dimension >> 1; (j ^= t) < t; t >>= 1);
            }
            for (size_t i = 0; i < dimension; i++) irt_harvey[i] = ((u128)irt[i] << 64) / modulus;
        }
        
        // printf("%lu\n",shuffled_indices[1]);
    }

    NTTFactors(const NTTFactors &copying) = default;

    NTTFactors(NTTFactors &&moving) = default;

    ~NTTFactors() {}
    std::vector<u64> rt;
    std::vector<u64> irt;

    std::vector<u64> rt_harvey;
    std::vector<u64> irt_harvey;

};

std::map<std::pair<u64, u64>, NTTFactors> &ntt_factors_cache() {
    static std::map<std::pair<u64, u64>, NTTFactors> global_ntt_factors_cache;
    return global_ntt_factors_cache;
}

std::map<std::pair<u64, u64>, NTTFactors> &intt_factors_cache() {
    static std::map<std::pair<u64, u64>, NTTFactors> global_intt_factors_cache;
    return global_intt_factors_cache;
}

inline const auto &
__find_or_create_ntt_factors(const u64 modulus, const size_t log_dimension,
                             const bool for_inverse = false) {
    if (!for_inverse) {
        auto it =
            ntt_factors_cache().find(std::make_pair(modulus, log_dimension));
        if (it == ntt_factors_cache().end()) {
            ntt_factors_cache().insert(
                std::make_pair(std::make_pair(modulus, log_dimension),
                               NTTFactors(modulus, log_dimension)));
            it =
                ntt_factors_cache().find(std::make_pair(modulus, log_dimension));
        }
        return it->second;
    } else {
        auto it =
            intt_factors_cache().find(std::make_pair(modulus, log_dimension));
        if (it == intt_factors_cache().end()) {
            intt_factors_cache().insert(
                std::make_pair(std::make_pair(modulus, log_dimension),
                               NTTFactors(modulus, log_dimension, true)));
            it = intt_factors_cache().find(
                std::make_pair(modulus, log_dimension));
        }
        return it->second;
    }
}
#define BARRETT_ADD(a, b, p) ({ \
    u64 aa = (u64)(a) + (b); \
    u64 x = aa - (((u128)(aa) * (((u128)1 << 64) / p)) >> 64) * (p); \
    if (x >= (p)) x -= (p); \
    x; \
})
void ntt_negacyclic_inplace_lazy(const size_t log_dimension, const u64 modulus,
                                 u64 coeffs[]) {
    const size_t dimension = 1ULL << log_dimension;
    // generate or read from cache
    const auto &ntt_factors =
        __find_or_create_ntt_factors(modulus, log_dimension);
    for (int i = 1, l = dimension >> 1; i < dimension; i <<= 1, l >>= 1) {
        // #pragma omp parallel for num_threads(thread_count) schedule(static)
        for (int j = 0; j < i; j++) {
            auto w = ntt_factors.rt[i + j];
            auto iw = ntt_factors.rt_harvey[i + j];
            for (int k = j * (l << 1); k < j * (l << 1) + l; k++) {
                auto x = coeffs[k] ;
                auto y = mul_mod_harvey_lazy(modulus, coeffs[k + l], w,iw);
                // auto y = (u128)coeffs[k + l]* w%modulus;
                // coeffs[k] = ((u128)x + y)%modulus,
                coeffs[k] =  BARRETT_ADD(x,y,modulus);
                coeffs[k + l] = (x + modulus - y)% modulus;
                // i128 temp = 1;
                // coeffs[k + l] = x - y + modulus;
                // if(coeffs[k + l] < 0)coeffs[k + l]+=modulus;
            }
        }
	}
}

void intt_negacyclic_inplace_lazy(const size_t log_dimension, const u64 modulus,
                                  u64 values[]) {
    const size_t dimension = 1ULL << log_dimension;
    // printf("%lu\n",modulus);
    // generate or read from cache
    const auto &intt_factors =
        __find_or_create_ntt_factors(modulus, log_dimension, true);
    for (int l = 1, i = dimension >> 1; i; l <<= 1, i >>= 1) {
        // #pragma omp parallel for num_threads(thread_count) schedule(static)
        for (int j = 0; j < i; j++) {
            auto w = intt_factors.irt[i + j];
            auto iw = intt_factors.irt_harvey[i + j];
            for (int k = j * (l << 1); k < j * (l << 1) + l; k++) {
                auto x = values[k] ;
                auto y = values[k + l];
                
                auto temp = x + modulus - y ;
                values[k + l] = mul_mod_harvey_lazy(modulus, temp, w, iw);
                // values[k + l] = (u128)temp*w % modulus;
                // values[k] = ((u128)x + y) % modulus;
                values[k] = BARRETT_ADD(x,y,modulus);
                // values[k] = x + y;
                // values[k + l] = (u128)(x - y + modulus) * w % modulus;
            }
        }
	}
    auto invn=__pow_mod(modulus,dimension,modulus-2);
    auto invn_harvey=((u128)invn << 64) / modulus;
    for(size_t i = 0; i < dimension; i++){
        // values[i] = ((u128)values[i] * invn) % modulus;
        values[i] = mul_mod_harvey_lazy(modulus, values[i], invn, invn_harvey);
    }

    // const u64 log_modulus = (u64)(log2(modulus) + 0.5);
    // const u64 div_fix = (modulus >= (1ULL << log_modulus)) ? 1 : 0;
    // for (size_t i = 0; i < dimension; i++, idx++) {
    //     values[i] -= ((values[i] >> log_modulus) - div_fix) * modulus;
    //     values[i] =
    //         mul_mod_harvey_lazy(modulus, values[i], intt_factors.seq[idx],
    //                             intt_factors.seq_harvey[idx]);
    // }
}
#undef BARRETT_ADD
#undef thread_count
void cache_ntt_factors_strict(const u64 log_dimension,
                              const std::vector<u64> &moduli) {
    for (auto modulus : moduli) {
        __find_or_create_ntt_factors(modulus, log_dimension);
        __find_or_create_ntt_factors(modulus, log_dimension, true);
    }
}

} // namespace hehub

#include "ntt.h"
#include "mod_arith.h"
#include "permutation.h"
#include <cmath>
#include <map>

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
        rt.resize(dimension);
        irt.resize(dimension);
        dw.resize(log_modulus);
        const u64 wnn = __get_2nth_unity_root(modulus, dimension);
        auto invwnn =
                __pow_mod(modulus, wnn, modulus-2);
        rt[0] = wnn;
        for (int i = 1; i < dimension; i++) rt[i] = (u128)rt[i - 1] * rt[i - 1] % modulus;
        irt[0] = invwnn;
        for (int i = 1; i < dimension; i++) irt[i] = (u128)irt[i - 1] * irt[i - 1] % modulus;
        if(!for_inverse){
            dw[0] = rt[log_modulus - 3];
            for (int i = 1; i < log_modulus - 2; i++)
                dw[i] = (u128)dw[i - 1] * irt[log_modulus - 1 - i]%modulus * rt[log_modulus - 3 - i]%modulus;//,printf("%d ",dw[i]);
            dw[log_modulus - 2] = (u128)dw[log_modulus - 3] * irt[1]%modulus;
        }else{
            dw[0] = irt[log_modulus - 3];
            for (int i = 1; i < log_modulus - 2; i++)
                dw[i] = (u128)dw[i - 1] * rt[log_modulus - 1 - i]%modulus * irt[log_modulus - 3 - i]%modulus;
            dw[log_modulus - 2] = (u128)dw[log_modulus - 3] * rt[1]%modulus;
        }
        // printf("%lu\n",shuffled_indices[1]);
    }

    NTTFactors(const NTTFactors &copying) = default;

    NTTFactors(NTTFactors &&moving) = default;

    ~NTTFactors() {}
    std::vector<u64> rt;
    std::vector<u64> irt;

    std::vector<u64> dw;
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

void ntt_negacyclic_inplace_lazy(const size_t log_dimension, const u64 modulus,
                                 u64 coeffs[]) {
    const size_t dimension = 1ULL << log_dimension;
    const size_t log_modulus = (u64)(log2(modulus) + 0.5);
    // generate or read from cache
    const auto &ntt_factors =
        __find_or_create_ntt_factors(modulus, log_dimension);
    u64 one = 1, imag = ntt_factors.rt[log_modulus-2];
    if (log_dimension & 1) {
        for (int i = 0; i < dimension/2; i++) {
            int x = coeffs[i], y = coeffs[i + dimension/2];
            coeffs[i] = (x + y )% modulus;
            coeffs[i + dimension/2] = (x + modulus - y)%modulus;
        }
    }
    for (int e = log_dimension & ~1; e >= 2; e -= 2) {
        const int m = 1 << e, m4 = m >> 2;
        u64 w2 = one;
        for (int i = 0; i < dimension; i += m) {
            const u64 w1 = (u128)w2 * w2 %modulus, w3 = (u128)w1 * w2 %modulus;
            for (int j = i; j < i + m4; ++j) {
                u64 a0 = (u128)coeffs[j + m4 * 0] * one % modulus, a1 = (u128)coeffs[j + m4 * 1] * w1 % modulus;
                u64 a2 = (u128)coeffs[j + m4 * 2] * w2 % modulus, a3 = (u128)coeffs[j + m4 * 3] * w3 % modulus;
                u64 t01p = (a0 + a1)%modulus, t23p = (a2 + a3)%modulus;
                u64 t01m = (a0 + modulus - a1)%modulus, t23m = ((u128)a2 + modulus - a3)%modulus * imag%modulus;
                coeffs[j + m4 * 0] = (t01p + t23p)%modulus;
                coeffs[j + m4 * 2] = (t01p + modulus - t23p)%modulus;
                coeffs[j + m4 * 1] = (t01m + t23m)%modulus;
                coeffs[j + m4 * 3] = (t01m + modulus - t23m)%modulus;
            }
            w2 = (u128)w2*ntt_factors.dw[__builtin_ctz(~(i >> e))]%modulus;
        }
        // printf("%d ",w2);
    }

}

void intt_negacyclic_inplace_lazy(const size_t log_dimension, const u64 modulus,
                                  u64 values[]) {
    const size_t dimension = 1ULL << log_dimension;
    const size_t log_modulus = (u64)(log2(modulus) + 0.5);
    // printf("%lu\n",modulus);
    // generate or read from cache
    const auto &intt_factors =
        __find_or_create_ntt_factors(modulus, log_dimension, true);
    u64 one = 1, imag = intt_factors.irt[log_modulus-2];
    for (int e = 2; e <= log_dimension; e += 2) {
        const u64 m = 1 << e, m4 = m >> 2;
        u64 w2 = one;
        for (int i = 0; i < dimension; i += m) {
            const u64 w1 = (u128)w2 * w2 %modulus, w3 = (u128)w1 * w2 %modulus;
            for (int j = i; j < i + m4; ++j) {
                u64 a0 = values[j + m4 * 0], a1 = values[j + m4 * 1];
                u64 a2 = values[j + m4 * 2], a3 = values[j + m4 * 3];
                u64 t01p = (a0 + a1)%modulus, t23p = (a2 + a3)%modulus;
                u64 t01m = (a0 + modulus - a1)%modulus, t23m = ((u128)a2 + modulus - a3) * imag%modulus;
                values[j + m4 * 0] = ((u128)t01p + t23p)%modulus*one%modulus;
                values[j + m4 * 2] = ((u128)t01p + modulus - t23p)*w1%modulus;
                values[j + m4 * 1] = ((u128)t01m + t23m)%modulus*w2%modulus;
                values[j + m4 * 3] = ((u128)t01m + modulus- t23m)*w3%modulus;
            }
            w2 = (u128)w2*intt_factors.dw[__builtin_ctz(~(i >> e))]%modulus;
        }
        // printf("%d ",w2);
    }
    if (log_dimension & 1) {
        for (int i = 0; i < dimension/2; i++) {
            int x = values[i], y = values[i + dimension/2];
            values[i] = (x + y )% modulus;
            values[i + dimension/2] = (x + modulus - y)%modulus;
        }
    }

    auto invn=__pow_mod(modulus,dimension,modulus-2);
    
    // std::cout<<invn<<'\n';
    for(size_t i = 0; i < dimension; i++){
        values[i] = ((u128)values[i] * invn) % modulus;
    } 
}

void cache_ntt_factors_strict(const u64 log_dimension,
                              const std::vector<u64> &moduli) {
    for (auto modulus : moduli) {
        __find_or_create_ntt_factors(modulus, log_dimension);
        __find_or_create_ntt_factors(modulus, log_dimension, true);
    }
}

} // namespace hehub

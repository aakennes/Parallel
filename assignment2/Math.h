#include<cstdio>
#include<iostream>
#include<cstring>

#include"SIMD.h"


int qpow(int a,int x,int p);

void findw(int countw,int wn[]);

template <std::uint32_t P> 
struct MontgomeryModInt32 {

private:
    u32 v;

    static constexpr u32 get_r() {
        u32 iv = P;

        for (u32 i = 0; i != 4; ++i)
            iv *= 2 - P * iv;

        return iv;
    }

    static constexpr u32 r = -get_r(), r2 = -u64(P) % P;

public:
    static constexpr u32 pow_mod(u32 x, u64 y) {
        if ((y %= P - 1) < 0)
            y += P - 1;

        u32 res = 1;

        for (; y != 0; y >>= 1, x = u64(x) * x % P)
            if (y & 1)
                res = u64(res) * x % P;

        return res;
    }

    static constexpr u32 get_pr() {
        u32 tmp[32] = {}, cnt = 0;
        const u64 phi = P - 1;
        u64 m = phi;

        for (u64 i = 2; i * i <= m; ++i) {
            if (m % i == 0) {
                tmp[cnt++] = i;

                while (m % i == 0)
                    m /= i;
            }
        }

        if (m > 1)
            tmp[cnt++] = m;

        for (u64 res = 2; res <= phi; ++res) {
            bool flag = true;

            for (u32 i = 0; i != cnt && flag; ++i)
                flag &= pow_mod(res, phi / tmp[i]) != 1;

            if (flag)
                return res;
        }

        return 0;
    }

    MontgomeryModInt32() = default;
    ~MontgomeryModInt32() = default;
    constexpr MontgomeryModInt32(u32 v) : v(reduce(u64(v) * r2)) {}
    constexpr MontgomeryModInt32(const MontgomeryModInt32 &rhs) : v(rhs.v) {}
    static constexpr u32 reduce(u64 x) {
        return x + (u64(u32(x) * r) * P) >> 32;
    }
    constexpr u32 get() const {
        u32 res = reduce(v);
        return res - (P & -(res >= P));
    }
    explicit constexpr operator u32() const {
        return get();
    }
    explicit constexpr operator i32() const {
        return i32(get());
    }
    constexpr MontgomeryModInt32 &operator=(const MontgomeryModInt32 &rhs) {
        return v = rhs.v, *this;
    }
    constexpr MontgomeryModInt32 operator-() const {
        MontgomeryModInt32 res;
        return res.v = (P << 1 & -(v != 0)) - v, res;
    }
    constexpr MontgomeryModInt32 inv() const {
        return pow(-1);
    }
    constexpr MontgomeryModInt32 &operator+=(const MontgomeryModInt32 &rhs) {
        return v += rhs.v - (P << 1), v += P << 1 & -(i32(v) < 0), *this;
    }
    constexpr MontgomeryModInt32 &operator-=(const MontgomeryModInt32 &rhs) {
        return v -= rhs.v, v += P << 1 & -(i32(v) < 0), *this;
    }
    constexpr MontgomeryModInt32 &operator*=(const MontgomeryModInt32 &rhs) {
        return v = reduce(u64(v) * rhs.v), *this;
    }
    constexpr MontgomeryModInt32 &operator/=(const MontgomeryModInt32 &rhs) {
        return this->operator*=(rhs.inv());
    }
    friend MontgomeryModInt32 operator+(const MontgomeryModInt32 &lhs,
                                        const MontgomeryModInt32 &rhs) {
        return MontgomeryModInt32(lhs) += rhs;
    }
    friend MontgomeryModInt32 operator-(const MontgomeryModInt32 &lhs,
                                        const MontgomeryModInt32 &rhs) {
        return MontgomeryModInt32(lhs) -= rhs;
    }
    friend MontgomeryModInt32 operator*(const MontgomeryModInt32 &lhs,
                                        const MontgomeryModInt32 &rhs) {
        return MontgomeryModInt32(lhs) *= rhs;
    }
    friend MontgomeryModInt32 operator/(const MontgomeryModInt32 &lhs,
                                        const MontgomeryModInt32 &rhs) {
        return MontgomeryModInt32(lhs) /= rhs;
    }
    friend std::istream &operator>>(std::istream &is, MontgomeryModInt32 &rhs) {
        return is >> rhs.v, rhs.v = reduce(u64(rhs.v) * r2), is;
    }
    friend std::ostream &operator<<(std::ostream &os, const MontgomeryModInt32 &rhs) {
        return os << rhs.get();
    }
    constexpr MontgomeryModInt32 pow(i64 y) const {
        if ((y %= P - 1) < 0)
            y += P - 1; // phi(P) = P - 1, assume P is a prime number

        MontgomeryModInt32 res(1), x(*this);

        for (; y != 0; y >>= 1, x *= x)
            if (y & 1)
                res *= x;

        return res;
    }
};

#include<cstdio>
#include<iostream>
#include<cstring>


#include "../params.h"


int qpow(int a,int x,int p);

void findw(int countw,int wn[]);

template <std::uint32_t P> 
struct MontgomeryModInt32 {

private:
    u32 v;
    const static u32 Montr = 104857599;
    const static u32 Montr2 = 45971250; 
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
 
    MontgomeryModInt32() = default;
    ~MontgomeryModInt32() = default;
    constexpr MontgomeryModInt32(u32 v) : v(reduce(u64(v) * Montr2)) {}
    constexpr MontgomeryModInt32(const MontgomeryModInt32 &rhs) : v(rhs.v) {}
    static constexpr u32 reduce(u64 x) {
        return x + (u64(u32(x) * Montr) * P) >> 32;
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
        return is >> rhs.v, rhs.v = reduce(u64(rhs.v) * Montr2), is;
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

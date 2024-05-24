#include"SIMD.h"


#define P 104857601

#define R 104857599
#define r2 45971250

static  u32 get_r(MMint v) {
    u32 iv = P;

    for (u32 i = 0; i != 4; ++i)
        iv *= 2 - P * iv;

    return iv;
}

static  u32 pow_mod(u32 x, u64 y) {
    if ((y %= P - 1) < 0)
        y += P - 1;

    u32 res = 1;

    for (; y != 0; y >>= 1, x = u64(x) * x % P)
        if (y & 1)
            res = u64(res) * x % P;

    return res;
}

u32 reduce(u64 x) {
    return x + (u64(u32(x) * R) * P) >> 32;
}

 MMint In(u32 v) {
    u32 res = reduce(u64(v)*r2);
    return res;
}

 u32 get(MMint v) {
    u32 res = reduce(v);
    return res - (P & -(res >= P));
}
    
 MMint subM(MMint v) {
    MMint res=0;
    return res = (P << 1 & -(v != 0)) - v;
}
 MMint add(MMint v,MMint rhs) {
    v += rhs - (P << 1), v += P << 1 & -(i32(v) < 0);
    return v;
}
 MMint sub(MMint v,MMint rhs) {
    v -= rhs, v += P << 1 & -(i32(v) < 0);
    return v;
}
 MMint mul(MMint v,MMint rhs) {
    return reduce(u64(v) * rhs);
}
 MMint Pow(MMint x,i64 y) {
    if ((y %= P - 1) < 0)
        y += P - 1; // phi(P) = P - 1, assume P is a prime number
    // if(y==204800)std::cout<<y<<'\n';
    MMint res=In(1);
    // while(y){
    //     if(y&1)res = mul(res,x);
    //     x = mul(x,x);
    //     y>>=1;
    // }
    for (; y != 0; y >>= 1, x = mul(x,x))
        if (y & 1)
            res = mul(res,x);
    // if(xx==1532914088)std::cout<<res<<'\n';
    return res;
}
 MMint inv(MMint v) {
    return Pow(v,-1);
}
 MMint div(MMint v,MMint rhs) {
    return mul(v,inv(rhs));
}

MMint neg(MMint v){
    return (!v - 1) & (P - v);
}

u32x8 setu32x8(u32 x){return (u32x8){x, x, x, x, x, x, x, x};}

// #undef P
// #undef R
// #undef r2
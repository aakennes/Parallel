#include<cstdio>
#include<iostream>
#include<cstring>
#include<semaphore.h>
#include<pthread.h>

#include"readwrite.h"
// #include"SIMD.h"

struct NttCommonData {
    int *a;
    int *r;
    int wnn;
    int invwnn;
    int id;
    int limit;
    int type;
};

struct NttMontgomeryData {
    int *a;
    int *r;
    Mint wnn;
    Mint invwnn;
    int id;
    int limit;
    int type;
};

struct NttMontgomeryMintData {
    Mint *a;
    int *r;
    Mint wnn;
    Mint invwnn;
    int id;
    int limit;
    int type;
};

struct DifMintData {
    Mint *a;
    Mint *r;
    Mint wnn;
    Mint invwnn;
    int id;
    int limit;
    int type;
};

struct DifX4Data {
    int *a;
    int *dw;
    int id;
    int imag;
    int logn;
    int limit;
};

struct DifX4MintData {
    Mint *a;
    Mint *dw;
    int id;
    Mint imag;
    int logn;
    int limit;
};

// #define MMint u32

void* ntt_common_func(void *arg);
void ntt_common(int *a,int *b,int *ab,int *r,int n);
void ntt_common(int *a,int *r,int limit,int type);

void* ntt_Montgomery_func(void *arg);
void ntt_Montgomery(int *a,int *b,int *ab,int *r,int n);
void ntt_Montgomery(int *a,int *r,int limit,int type);

void* ntt_Montgomery_Mint_func(void *arg);
void ntt_Montgomery_Mint(Mint *a,Mint *b,Mint *ab,int *r,int n);
void ntt_Montgomery_Mint(Mint *a,int *r,int limit,int type);

void* ntt_dif_func(void *arg);
void* ntt_dit_func(void *arg);
void ntt_dif(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dif(int *a,int *rt,int limit,int type);

void* ntt_dif_Mint_func(void *arg);
void* ntt_dit_Mint_func(void *arg);
void ntt_dif_Mint(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n);
void ntt_dif_Mint(Mint *a,Mint *rt,int limit,int type);

void* ntt_dif_x4_func(void *arg);
void* ntt_dit_x4_func(void *arg);
void ntt_dif_x4(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dif_x4(int *a,int *rt,int *irt,int limit,int level);
void ntt_dit_x4(int *a,int *rt,int *irt,int limit,int level);

void* ntt_dif_x4_Mint_func(void *arg);
void* ntt_dit_x4_Mint_func(void *arg);
void ntt_dif_x4_Mint(Mint *a, Mint *rt, Mint *irt, int limit, int level);
void ntt_dit_x4_Mint(Mint *a, Mint *rt, Mint *irt, int limit, int level);
void ntt_dif_x4_Mint(Mint *a, Mint *b, Mint *ab, Mint *rt,Mint *irt, int n);

// void ntt_dif_x4_avx2(MMint *a,MMint *b,MMint *ab,MMint *rt,MMint *irt,Mintx8 *RT1,Mintx8 *RT2,Mintx8 *IRT1,Mintx8 *IRT2);
// void ntt_dit_x4_avx2(MMint *a,Mintx8 *rt,Mintx8 *irt,int limit,int level);
// void ntt_dif_x4_avx2(MMint *a,Mintx8 *rt,Mintx8 *irt,int limit,int level);
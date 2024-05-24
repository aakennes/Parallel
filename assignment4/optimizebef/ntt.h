#include<cstdio>
#include<iostream>
#include<cstring>
#include <omp.h>
#include"readwrite.h"

void ntt_common(int *a,int *b,int *ab,int *r,int n);
void ntt_common(int *a,int *r,int limit,int type);
void ntt_common_openmp(int *a,int *b,int *ab,int *r,int n);
void ntt_common_openmp(int *a,int *r,int limit,int type);

void ntt_Montgomery(int *a,int *b,int *ab,int *r,int n);
void ntt_Montgomery(int *a,int *r,int limit,int type);
void ntt_Montgomery_openmp(int *a,int *b,int *ab,int *r,int n);
void ntt_Montgomery_openmp(int *a,int *r,int limit,int type);

void ntt_Montgomery_Mint(Mint *a,Mint *b,Mint *ab,int *r,int n);
void ntt_Montgomery_Mint(Mint *a,int *r,int limit,int type);
void ntt_Montgomery_Mint_openmp(Mint *a,Mint *b,Mint *ab,int *r,int n);
void ntt_Montgomery_Mint_openmp(Mint *a,int *r,int limit,int type);

void ntt_dif(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dit(int *a,int *rt,int limit);
void ntt_dif(int *a,int *rt,int limit);
void ntt_dif_openmp(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dit_openmp(int *a,int *rt,int limit);
void ntt_dif_openmp(int *a,int *rt,int limit);

void ntt_dif_Mint(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n);
void ntt_dit_Mint(Mint *a,Mint *rt,int limit);
void ntt_dif_Mint(Mint *a,Mint *rt,int limit);
void ntt_dif_Mint_openmp(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n);
void ntt_dit_Mint_openmp(Mint *a,Mint *rt,int limit);
void ntt_dif_Mint_openmp(Mint *a,Mint *rt,int limit);

void ntt_dif_x4(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dit_x4(int *a,int *rt,int *irt,int limit,int level);
void ntt_dif_x4(int *a,int *rt,int *irt,int limit,int level);
void ntt_dif_x4_openmp(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dit_x4_openmp(int *a,int *rt,int *irt,int limit,int level);
void ntt_dif_x4_openmp(int *a,int *rt,int *irt,int limit,int level);

void ntt_dif_x4_Mint(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n);
void ntt_dit_x4_Mint(Mint *a,Mint *dw,Mint imag,int limit,int level);
void ntt_dif_x4_Mint(Mint *a,Mint *dw,Mint imag,int limit,int level);
void ntt_dif_x4_Mint_openmp(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n);
void ntt_dit_x4_Mint_openmp(Mint *a,Mint *dw,Mint imag,int limit,int level);
void ntt_dif_x4_Mint_openmp(Mint *a,Mint *dw,Mint imag,int limit,int level);
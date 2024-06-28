#include<cstdio>
#include<iostream>
#include<cstring>

#include"readwrite.h"


void ntt_common(int *a,int *b,int *ab,int *r,int n);
void ntt_common(int *a,int *r,int limit,int type);

void ntt_dif(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dit(int *a,int *rt,int limit);
void ntt_dif(int *a,int *rt,int limit);

void ntt_dif_x4(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dit_x4(int *a,int *rt,int *irt,int limit,int level);
void ntt_dif_x4(int *a,int *rt,int *irt,int limit,int level);

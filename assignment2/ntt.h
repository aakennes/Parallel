#include"params.h"
#include"readwrite.h"
#include"Math.h"

#include<cstdio>
#include<iostream>
#include<cstring>

void ntt_common(int *a,int *b,int *ab,int *r);
void ntt_common(int *a,int *r,int limit,int type);
void ntt_Montgomery(int *a,int *b,int *ab,int *r);
void ntt_Montgomery(int *a,int *r,int limit,int type);
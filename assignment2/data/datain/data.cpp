#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../../params.h"
#include <cstring>

int main(){
    int q = Original_Q, n = Original_N;

    FILE* fp=fopen("NTT256_1409.in","w+");
    
    fprintf(fp,"%d %d\n",n,q);
    for (int i = 0; i < n; i++){
        srand(clock() * (i+1));
        fprintf(fp,"%d ",rand() % q);
    }
    fprintf(fp,"\n");
    for (int i = 0; i < n; i++){
        srand(clock() * (i+1));
        fprintf(fp,"%d ",rand() % q);
    }
    fclose(fp);
}

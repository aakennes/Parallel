#include"readwrite.h"

void fRead(FILE *fp,int a[],int b[],int n,int q){
    std::string str1="../data/datain/NTT";
    std::string str2=std::to_string(n);
    std::string str3=std::to_string(q);
    std::string strin=str1+str2+"_"+str3+".in";
    char charArrayin[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), charArrayin);
    charArrayin[strin.size()] = '\0';
    fp=fopen(charArrayin, "r+");
    for (int i = 0; i < n; i++){   
        fscanf(fp,"%d",&a[i]);
        // std::cout<<a[i]<<'\n';
    }
    for (int i = 0; i < n; i++){   
        fscanf(fp,"%d",&b[i]);
    }
}


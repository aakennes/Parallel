#include"readwrite.h"

void fRead(FILE *fp,int a[],int b[],int n){
    std::string str1="data/datain/NTT";
    std::string str2=std::to_string((int)Original_N);
    // std::string str3=std::to_string((int)Original_Q);
    std::string str3=std::to_string(1409);
    std::string strin=str1+str2+"_"+str3+".in";
    std::cout<<strin<<'\n';
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

void fWrite(FILE *fp,int ab[],int n){
    std::string str1="data/dataout/NTT";
    std::string str2=std::to_string((int)Original_N);
    std::string str3=std::to_string((int)Original_Q);
    std::string strout=str1+str2+"_"+str3+".out";
    char charArrayout[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), charArrayout);
    charArrayout[strout.size()] = '\0';
    fp=fopen(charArrayout, "w+");
    for (int i = 0; i < n*2 ; i++){
        // std::cout<<ab[i]<<'\n';
        fprintf(fp,"%d ",ab[i]);
    }
}

void fRead(FILE *fp,Mint a[],Mint b[],int n){
    std::string str1="data/datain/NTT";
    std::string str2=std::to_string((int)Original_N);
    // std::string str3=std::to_string((int)Original_Q);
    std::string str3=std::to_string(1409);
    std::string strin=str1+str2+"_"+str3+".in";
    std::cout<<strin<<'\n';
    char charArrayin[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), charArrayin);
    charArrayin[strin.size()] = '\0';
    fp=fopen(charArrayin, "r+");
    for (int i = 0; i < n; i++){   
        int temp;
        fscanf(fp,"%d",&temp);
        a[i]=temp;
        // std::cout<<a[i]<<'\n';
    }
    for (int i = 0; i < n; i++){   
        int temp;
        fscanf(fp,"%d",&temp);
        b[i]=temp;
    }
}

void fWrite(FILE *fp,Mint ab[],int n){
    std::string str1="data/dataout/NTT";
    std::string str2=std::to_string((int)Original_N);
    std::string str3=std::to_string((int)Original_Q);
    std::string strout=str1+str2+"_"+str3+".out";
    char charArrayout[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), charArrayout);
    charArrayout[strout.size()] = '\0';
    fp=fopen(charArrayout, "w+");
    for (int i = 0; i < n*2 ; i++){
        // std::cout<<ab[i]<<'\n';
        fprintf(fp,"%d ",ab[i].get());
        // printf("%d ",ab[i].get());
    }
}
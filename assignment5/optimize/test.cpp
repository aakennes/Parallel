#include<bits/stdc++.h>
using namespace std;
int main(){
    system("rm result.out");
    system("touch result.out");
    for(int i=0;i<=9;++i){
        for(int j=0;j<=3;++j){
            FILE* fpnm=fopen("result.in","w");
            fprintf(fpnm,"%d %d\n",i,j);
            fclose(fpnm);

            system("rm test.out");
            system("touch test.out");
            system("make clean");
            system("make");
            
            for(int i=0;i<20;++i){
                system("mpirun -np 4 main");
            }
            double ans=0;
            FILE* fileread=fopen("test.out","r+");
            for(int i=0;i<20;++i){
                double x;
                fscanf(fileread,"%lf",&x);
                ans+=x;
            }
            fclose(fileread);
            FILE* fpwrite=fopen("result.out","a+");
            fprintf(fpwrite,"%lf\n",ans/20);
            fclose(fpwrite);
        }
    }
}
#include<bits/stdc++.h>
using namespace std;
int main(){
    freopen("3.csv","w",stdout);
    printf("n,cnt,优化前,优化1,优化2,\n");
    const int maxn=1e7,maxs=26;
    for(int n=maxs;n<=maxs;n++){
        system("del 1.out");
        system("g++ -pthread -o II II.cpp");
        ofstream In1("1.in");
        int m=maxn*200/pow(2,n)+n;
        In1<<n<<" "<<m;
        In1.close();
        system("II");
        ifstream In2("1.out");
        double sum1,sum2,sum3;
        In2>>sum1>>sum2>>sum3;
        In2.close();
        printf("%d,%d,%lf,%lf,%lf,\n",n,m,sum1,sum2,sum3);
    }
}
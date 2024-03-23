#include<bits/stdc++.h>
using namespace std;
int main(){
    freopen("2.csv","w",stdout);
    printf("n,cnt,优化前,优化后,\n");
    const int maxn=2e5;
    for(int n=500;n<=500;n+=50){
        system("rm 1.out");
        system("aarch64-linux-gnu-g++ -static -o I -march=armv8.2-a I.cpp");
        int m=maxn/n;
        ofstream In1("1.in");
        In1<<n<<" "<<m;
        In1.close();
        system("./I");
        ifstream In2("1.out");
        double sum1,sum2;
        In2>>sum1>>sum2;
        In2.close();
        printf("%d,%d,%lf,%lf,\n",n,m,sum1,sum2);
    }
}
#include<bits/stdc++.h>
using namespace std;
int main(){
    freopen("3.csv","w",stdout);
    printf("n,cnt,优化前,优化后,\n");
    const int maxn=2e5;
    for(int n=50;n<=5000;n+=50){
        system("del 1.out");
        system("g++ -o I I.cpp");
        int m=maxn/n;
        ofstream In1("1.in");
        In1<<n<<" "<<m;
        In1.close();
        system("I.exe");
        ifstream In2("1.out");
        double sum1,sum2,sum3;
        In2>>sum1>>sum2>>sum3;
        In2.close();
        printf("%d,%d,%lf,%lf,\n",n,m,sum1,sum2);
    }
}
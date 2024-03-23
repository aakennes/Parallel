#include<bits/stdc++.h>
using namespace std;
int main(){
    freopen("4.csv","w",stdout);
    printf("n,cnt,2*1,3*1,3*3,4*4,5*5,\n");
    const int maxn=1e7,maxs=26;
    for(int n=26;n<=26;n++){
        system("rm 1.out");
        system("clang++ -pthread -o II II.cpp");
        ofstream In1("1.in");
        int m=maxn*200/pow(2,n)+n;
        In1<<n<<" "<<m;
        In1.close();
        system("./II");
        ifstream In2("1.out");
        double sum1,sum2,sum3,sum4,sum5;
        In2>>sum1>>sum2>>sum3>>sum4>>sum5;
        In2.close();
        printf("%d,%d,%lf,%lf,%lf,%lf,%lf,\n",n,m,sum1,sum2,sum3,sum4,sum5);
    }
}
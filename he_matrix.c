#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#define A0 0.000268   //Bohr radius in 1/electron-volt

float SimpsonDouble1(float (*f)(float x, float y), float b, int n){

int i, j, ki, kj;
float h, s, x, y, r;
h = (float) b/n;
i = j = 0;
s=0;
for(x=0; x<=b; x+=h){
        if(x==0 || x==b)
            ki=1;
        else
        {
            i=1-i; //if the index of x is even, i=0, if it's odd i=1
            ki=2*(i+1);
        }
	for(y=0; y<=x; y+=h){
        r=f(x,y);
        if(y==0 || y==x)
            kj=1;
        else
        {
            j=1-j;  //if the index of y is even, j=0, if it's odd j=1
            kj=2*(j+1);
        }
		s = s + (float)(ki*kj*r);
	}
}
return s*h*h/9;
}

float SimpsonDouble2(float (*f)(float x, float y), float b, int n){

int i, j, ki, kj;
float h, s, x, y,r;
h = (float) b/n;
i = j = 0;
s=0;
for(x=0; x<=b; x+=h){
        if(x==0 || x==b)
            ki=1;
        else
        {
            i=1-i;
            ki=2*(i+1);
        }
	for(y=x; y<=b; y+=h){
        r=f(x,y);
        if(y==x || y==b)
            kj=1;
        else
        {
            j=1-j;
            kj=2*(j+1);
        }
		s = s + (float)(ki*kj*r);
	}
}
return s*h*h/9;
}

float funct1(float x){

return x/exp(x);
}

float funct2(float x){

return x*x/exp(x);
}

float funct3(float x){

return x*x*x/exp(x);
}

float funct4(float x){

float s = x*x;
return s*s/exp(x);
}

float funct5(float x){

float s = x*x;
return s*s*x/exp(x);
}

float funct6(float x){

float s = x*x;
return s*s*s/exp(x);
}

float funct10(float x, float y){

return funct3(x)*funct4(y);
}

float funct11(float x, float y){

return funct2(x)*funct5(y);
}

float funct12(float x, float y){

    return funct1(x)*funct6(y);
}

float funct20(float x, float y){

    return funct4(x)*funct3(y);
}

float funct21(float x, float y){

    return funct5(x)*funct2(y);
}

float funct22(float x, float y){

    return funct6(x)*funct1(y);
}

//function that calculates the matrix element taking as input the truncation value b and the number of segments n
float element1s1sH3d3d(float b, int n){
    return (SimpsonDouble1(funct10,b,n)+SimpsonDouble2(funct20,b,n)+(SimpsonDouble1(funct11,b,n)+SimpsonDouble2(funct21,b,n))/3+(SimpsonDouble1(funct12,b,n)+SimpsonDouble2(funct22,b,n))/5)/(40960*A0);
}

//function that displays convergence rates
void ConvergenceTest(float b, int start, int step, int number){

    int v = start+step*(number-1);
    float ek1,ek2;
    float real = 1103/(8*40960*A0);
    ek1 = element1s1sH3d3d(b,start);
    for(int i=start+step;i<=v;i+=step){


            ek2 = element1s1sH3d3d(b,i);

            printf("\n%d: %.5f\n", i-step, (ek2-real)/(ek1-real));
            ek1=ek2;
    }

}

int main(){

    printf("************* Calculating the Matrix Element <1s1s|H|3d3d>***************\n");
    printf("   Program written by Mounir AFIFI for Atomic Physics Project \n   2nd Year Master Computational Physics\n\n");

    printf("Enter n (number of segmentations - if n<2, 2 will be taken as default - if n>1000000, 1000000 will be taken by default): ");
    int n;
    //printf("\n%d\n",boolean);
    while(scanf("%d",&n)==0){
        printf("Please enter a positive integer: ");
        getchar();
    }
    if(n<2)
        n=2;
    else
        if(n>1000000)
            n=1000000;

    float trunc;
    float max = -log(FLT_MIN);
    printf("Please enter the truncation threshold (if the value is larger than 87.33 or smaller than or equal to 0, it will be 87.33 by default): ");
    while(scanf("%f",&trunc)==0){
        printf("Please enter a positive real number (if the value is larger than 87.33 or smaller than or equal to 0, it will be 87.33 by default): ");
        getchar();
    }
    if(trunc<=0||trunc>max)
        trunc=max;

    printf("\nNumber of segments: %d - Truncation threshold for integration: %.2f\n\n",n,trunc);

    printf("I10 = %f\n",SimpsonDouble1(funct10,trunc,n));
    printf("I20 = %f\n\n",SimpsonDouble2(funct20,trunc,n));
    printf("I11 = %f\n",SimpsonDouble1(funct11,trunc,n));
    printf("I21 = %f\n\n",SimpsonDouble2(funct21,trunc,n));
    printf("I12 = %f\n",SimpsonDouble1(funct12,trunc,n));
    printf("I22 = %f\n\n",SimpsonDouble2(funct22,trunc,n));
    printf("<1s1s|H|3d3d> = %f\n",element1s1sH3d3d(trunc,n));

    ConvergenceTest(max, 10000, 10000, 30);

}

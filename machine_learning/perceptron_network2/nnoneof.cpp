#include <iostream>
#include <stdio.h>
#include <math.h>
#include "mat.h"
#include "rand.h"

//double f(double x) { return 2*(1.0/double (1.0 + exp(-2*x))-1); };
double f(double x) { return (1.0/double (1.0 + exp(-4.0*x))); };
double g(double x){ return (x <0.5 ? 0.0 : 1.0);  }
double  eat=0.1;

int main ()
{    
     initRand();
    int h, inputs;
    scanf("%d", &inputs);
    scanf("%d",&h);
    Matrix X;// matrix with all inputs and target answers
    X.read();//read data in to the matrix
    X.normalizeCols();
    int max_r= X.maxRows();
    int max_c= X.maxCols();
    //dimension of answer matrix
    int results_cols =max_c - inputs;
    // with bisa colums
    Matrix S (max_r,inputs+1,"tranning sample");
    Matrix V(inputs+1,h,"first wieght");
    Matrix H(max_r,h,"H without bisa");
    Matrix H_(max_r,h+1,"H with bisa -1");
    Matrix Y(max_r,results_cols,"my_results");
    Matrix T(max_r,results_cols,"target results");
    Matrix W(h+1,results_cols,"second weight");
    Matrix delat_W, delat_H;

    // initialize the matrix
     S.constant(-1.0);
     S.insert(X.extract(0,0,0,inputs),0,0);
     T.constant(0);
     T.insert(X.extract(0,inputs,0,0),0,0);
     W.rand(-0.2,0.2);
     V.rand(-0.2,0.2);
     //initialize withe bisa
     H_.constant(-1.0);
     Matrix temp1, temp2,temp3, temp4, temp5;
     for (int i =0; i<10000; i++){
        // H = f(XV)   
        H = (S.dot(V)).map(f);
        //H.map(f);
        H_.insert(H, 0 ,0);
        // Y = f(H_W)       
        Y= (H_.dot(W)).map(f);
        //Y.map(f);
        // delate_W = (Y-T)*Y*(1-Y)
        
        temp1= Y;
         temp1.sub(T);
        temp2= Y;
         temp2.scalarPreSub(1);
        delat_W = (temp1.mult(Y)).mult(temp2);
        // delate_H = H_ *(1-H_)*(delate_W Wt)
        temp3 = H_;
        temp3.scalarPreSub(1);
        temp4 =delat_W;
        temp4 = temp4.dotT(W);
        delat_H = (temp3.mult(H_)).mult(temp4);
        // W-= eat*H_t delat_W
        W.sub((H_.Tdot(delat_W)).scalarMult(eat));
          //V − = αX+T delat_H−
          //delat_H.print();
          //delat_H.extract(0,0,max_r,h).print();
         temp5=delat_H.extract(0,0,max_r,h);
        V.sub((S.Tdot(temp5)).scalarMult(eat));
     }
     // H = f(XV)
        Matrix NX, NXX;
        NX.read();
        NXX=NX;
        NX.normalizeCols();
        Matrix NH(NX.maxRows(),h,"NH without bisa");
        Matrix NH_(NX.maxRows(),h+1,"NH with bisa -1");
        Matrix NS(NX.maxRows(),inputs+1,"tranning sample");
        Matrix NY(NX.maxRows(),results_cols,"tranning answers");
        NY.constant(0.0);
        NS.constant(-1.0);
        NS.insert(NX,0,0);
        NH.constant(0.0);
        NH = NS.dot(V);
        NH.map(f);
        NH_.constant(-1.0);
        NH_.insert(NH, 0 ,0);
        // Y = f(H_W)
        NY= NH_.dot(W);
        NY.map(f);
        //NY.map(g);
        double m=0;
        int locat;
        for (int i=0; i < NX.maxRows(); i++)
        {
            NXX.writeLine(i);
            m = -100000;
            for(int j=0; j<results_cols; j++)
            {
                if (NY.get(i,j)>m)
                {
                    m=NY.get(i,j);
                    locat=j;
                }

            }  
            printf("  %d\n",locat);
        }
          
     return 0;

}







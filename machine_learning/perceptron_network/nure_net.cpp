#include <iostream>
#include <stdio.h>
#include "mat.h"
#include "rand.h"

int main ()
{
	int rows, total_cols, features, results_cols;
	scanf("%d",&features);
	Matrix x;
	x.read();
	//x.print();

	int max_r=x.maxRows();
	int max_c=x.maxCols();
	results_cols = max_c - features;

	Matrix s(x.maxRows(), features+1,"Traning_sample");
	s.constant(0);
	s.insert(x.extract(0,0,0,features),0,0);
	Matrix ze(max_r,1,"00");
	
	ze.constant(-1);	
	//ze.print();
	s.insert(ze,0,features);
	//s.print();

	Matrix T(max_r,results_cols,"results");
	T.constant(0);
	T.insert(x.extract(0,features,0,0),0,0);
	//T.print();

	Matrix Y(max_r,results_cols,"my_results");
	Matrix W(features+1,results_cols,"weights");
	initRand();
	W.constant(randUnit());
	int count=0;
	double eat=0.1;
	Matrix _T;
	// training 
	do{
		Y=s.dot(W);
		for (int i=0;i<max_r;i++)
		{
			for(int j=0; j<results_cols;j++)
			{
				if (Y.get(i,j)>0.5)
					Y.set(i,j,1);
				else
					Y.set(i,j,0);
			}
		}
		if (!Y.equal(T))
		{
			//W+ = αXT(T − Y )
			_T=T;
			W.add(s.Tdot(_T.sub(Y)).scalarMult(eat));
		}
		else
			break;
		count++;
	}while (count<30000);
	//Y.print();

	//testing 
	x.read();
	s.insert(x.extract(0,0,0,features),0,0);
	s.insert(ze,0,features);
	Y=s.dot(W);
		for (int i=0;i<max_r;i++)
		{
			for(int j=0; j<results_cols;j++)
			{
				if (Y.get(i,j)>0.5)
					Y.set(i,j,1);
				else
					Y.set(i,j,0);
			}
		}
	//Y.print();
	max_r=x.maxRows();
	max_c=x.maxCols();
		for (int i=0;i<max_r;i++)
		{
			for(int j=0; j<features;j++)
			{
				printf("%10.2lf", x.get(i,j));
			}
			for(int k=0;k<results_cols;k++)
			{
				printf("%10.2lf", Y.get(i,k));
			}
			
			printf("\n");
		}
		
		return 0;
}



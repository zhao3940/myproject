#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "mat.h"
#include "rand.h"

using namespace std;



int features;

void build_kdtrees(Matrix &trees,int feature ,int start, int end){
    int mid_point = (start+end)/2;
//    cout <<"start: "<<start<<endl;
//    cout <<"mid: "<<mid_point<<endl;
//    cout <<"end: "<<end<<endl;
    if (feature == features+1)
        feature = 1;
    trees.sortRowsByCol(feature,start,end);    
    feature++;
    if (end -start > 1){
        build_kdtrees(trees,feature,start,mid_point-1);
        build_kdtrees(trees,feature,mid_point+1,end);    	
    }
}

double dis(Matrix trees, Matrix others, int pick_point ){
	double sum=0;
	double tem=0;
	for (int i=0; i < others.maxCols();i++ )
	{
		tem = trees.get(pick_point , i+1) - others.get(0,i);
		sum = sum + tem * tem;
	}
	return sqrt(sum);
}

void nearest(Matrix trees, Matrix item,int rowstart, int rowend,int feature, double &best, int &bestex){
	// for single node
	if (rowstart == rowend){
		cout<<"RANGE: "<<rowstart<<"  to  "<<rowend<<endl;
		double x=dis(trees, item,rowstart);
		if (x<best){
			best = x;
			bestex = rowend;
			if (feature == features +1)
				feature = 1;
		//	cout<<setprecision(3)<<"BESTLEAF:  "<<best<<"  "<<rowstart<<endl;
		printf("BESTLEAF: %.3lf  %d\n",best,rowstart);
		}
	}
	else{
		cout<<"RANGE: "<<rowstart<<" to "<<rowend<<endl;
		int pick_point = (rowend +rowstart) / 2;
		if(features+1 == feature)
			feature = 1;
		// which side do I search
		if (trees.get(pick_point,feature) >= item.get(0,feature-1))
		{
			if (pick_point -rowstart> 0)
			{
				nearest(trees, item, rowstart, pick_point-1, feature+1, best, bestex);
				double temp1 = item.get(0,feature-1);
				double temp2=trees.get(pick_point,feature);
				if (fabs(temp2-temp1)>best)
				{
		    //        cout << fixed << setprecision(3)<<"BESTLEAF:  "<<best<<"  "<<bestex<<endl;
		            return;
				}
			}
			if (dis(trees,item,pick_point) < best)
			{
				best = dis(trees,item,pick_point);
	            cout << fixed << setprecision(3)<<"BESTPARENT1("<<feature<<"): "<<best<<" "<<pick_point<<endl;
				bestex = pick_point;
				if (best == 0)
				{
					return;
				}
			}
			if(rowend - pick_point >0)
			{
				nearest(trees, item, pick_point+1, rowend,feature+1, best, bestex);
			}
		}
			else{
				if (rowend - pick_point > 0)
				{
					nearest(trees, item, pick_point+1, rowend ,feature+1, best, bestex);
					double temp1 = item.get(0,feature-1);
					double temp2=trees.get(pick_point,feature);
					if (fabs(temp2-temp1)>best)
					{
		            	//cout << fixed << setprecision(3)<<"BESTLEAF:  "<<best<<"  "<<bestex<<endl;
		            	return;
					}	
				}
				if (dis(trees,item,pick_point) < best)
				{
					best = dis(trees,item,pick_point);
	            	cout << fixed << setprecision(3)<<"BESTPARENT("<<feature<<"): "<<best<<" "<<pick_point<<endl;
					bestex = pick_point;
					if (best == 0)
						{
							return;
						}
				}
				if(pick_point -rowstart > 0)
				{
					nearest(trees, item, rowstart ,pick_point-1,feature+1, best, bestex);
				}
		}
	}
	
}
void print(Matrix trees, int row,int l)
{
		for(int k = l; k<trees.maxCols();k++){
			//cout << setprecision(2)<<"   "<<trees.get(row,k);
			printf(" %.2lf",trees.get(row,k));
	}
}

int main(){
    Matrix trees("kd tree");
    char **label;
    label = trees.readLabeledRow();
    int end = trees.maxRows() - 1;
    features = trees.maxCols() - 1;
    //cout<<"features: "<<features<<endl;
    cout<<"KDTree version of matrix";
	build_kdtrees(trees,1,0,end);
    trees.printLabeledRow(label);
    int bestex;
    //double best;
    Matrix data;
	data.read();
	int pick_point;
	double best = DBL_MAX;
	for (int i = 0; i < data.maxRows(); i++){
    	//nearest(answer[3], trees, data.subMatrix(i,0,1,data.maxCols()), pick_point,1,0,end);
    	best = DBL_MAX;
		cout << "SOLVE:  ";
		print(data,i,0);
        cout<<endl;
		nearest(trees, data.subMatrix(i,0,1,data.maxCols()), 0, end, 1, best, bestex);
		cout<<"Ans:   ";
		print(trees,bestex,1);
		cout<< setprecision(0)<<"  "<<int (trees.get(bestex,0))<<" "<<label[int(trees.get(bestex,0))]<<endl;
		//cout<<"dis: "<<answer[1]<<endl;	
		cout<<endl<<endl;
	}

    return 0;
}

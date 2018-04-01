#include <iostream>
#include <iostream>
#include <stdio.h>
#include "mat.h"
#include "rand.h"
using namespace std;

int main()
{
  Matrix sample_org;
  int k;
  scanf("%d",&k);
  sample_org.readImagePpm("","sample_input");
  sample_org.printSize("Pic");
  
  // center the data
  Matrix sample_center;
  sample_center = sample_org;
  //double const mean = sample_org.mean();
  //printf("mean: %f\n",mean);
  Matrix mean;  
  mean = sample_org.meanVec();
  sample_center.subRowVector(mean);
  Matrix stddevs;
  stddevs = sample_org.stddevVec();
  sample_center.divRowVector(stddevs);
  //sample_center.print();
  
  // computer covariance matrix
  Matrix covariaance_center;
  covariaance_center = sample_center.cov();
  //covariaance_center.print();
  
  // computer eigenvalues and eigenvectors pf covariaance_center matrix
  Matrix eigenvalue_v, eigenvector_w;
  eigenvector_w=covariaance_center; 
  //eigenvetor_w destroys self by replacing self with eigenvectors in rows
  //also return a new matrix constant eigenvalues;
  eigenvalue_v = eigenvector_w.eigenSystem();
  eigenvalue_v.printSize("Eigenvalues");
  // normalize eigenvector_w
  //eigenvector_w.normalizeCols();
  
  //sort eigenvectors by eigenvalu, and pick k
  // eigenvector_w already be sorted by eigenSystem()
  Matrix K_vectors;
  K_vectors = eigenvector_w.extract(0,0,k,0);
  //new data matrix
  Matrix new_matrix;
  new_matrix= sample_center.dotT(K_vectors);
  //new_matrix= sample_center.dotT(eigenvector_w);
  new_matrix.printSize("Encoded");
  // recovering data from new matrix
  Matrix recover_matrix;
  //recover_matrix = new_matrix.dot(eigenvector_w);
  recover_matrix = new_matrix.dot(K_vectors);
  recover_matrix.multRowVector(stddevs).addRowVector(mean);
  //recover_matrix.print();
 
  double dis=1;
  dis =sample_org.dist2(recover_matrix);
  dis /= (sample_org.maxCols() * sample_org.maxRows());
  printf("DIST: %e\n",dis);

  // write in to a image file
  recover_matrix.writeImagePpm("z-after.ppm","");
  return 0;
}

#ifndef MATH
#define MATH

// // // // // // // // // // // // // // // // 
//
// Some simple matrix operators.
//
//
// These routines are self contained except for random number generation.
// Comment these out if you want completely stand-alone functionality.
// Uses some code from Numerical Recipes (see comments).  This code is NOT
// designed to be blindingly fast but rather just get the job done.  The
// variable Matrix::debug can be set to true if you want to debug memory
// allocation to look for overallocation of matrices and memory leaks.
// NOTE: most routines overwrite self with the answer.  For example: add
// adds to self.  See further in this comment block.
//
// Author: Robert B. Heckendorn, University of Idaho, 2018
// Version: 2.5
// Date: Apr 7, 2018

// IMPORTANT: If running on MICROSOFT WINDOWS uncomment this define statement!
// #define WINDOWS

// 
// Most matrix library routines replace the contents of the matrix
// object (overwrite self).  That is: X.sub(Y) will replace X with X - Y.
//
// Only the following routines allocate a new matrix leaving the original
// matrix untouched.  For these routines you need to assign the result
// to a new matrix: X = Y.dot(Z)
// or print them out: Y.dot(Z).print()
// 
// argMinRow()
// cartesianRow(double (*f)(int size, double *x, double *y), Matrix &other)
// cov()
// cov(Matrix &other)
// dot(const Matrix &other)
// dotT(const Matrix &other)
// Tdot(const Matrix &other)
// eigenSystem() // DANGER this replaces the object AND returns a new matrix
// extract(int minr, int minc, int sizer, int sizec)
// meanVec()
// minRow()
// pickRows(int match, const Matrix &list, int &num)
// stddevVec()
// transpose()
// 

#include <stdlib.h>
#include <stdio.h>
#include <vector>      // supports submatrices
#include <string>      // matrix names are strings
#include "rand.h"      // portable random number generator.  Include exactly
                       // ONE of the random number cpp files in your compile
class Matrix;

// // // // // // // // // // // // // // // // 
//
// class MatrixRowIter
//
// Iterator for incrementing through rows in a matrix
//
class MatrixRowIter {
private:
    Matrix *mat;
    int r;
    Matrix *arow;
    bool more;

public:
    MatrixRowIter(Matrix *mat);
    ~MatrixRowIter();

public:
    Matrix *rowBegin();
    Matrix *rowNext();
    bool rowNotEnd();
    int row();
};



// // // // // // // // // // // // // // // // 
//
// class Matrix
//
// A simple class for matrix operations.  It has some nice debugging
// features like trying hard to check that the proper dimensions are used.
// It is very draconian about this so there is an important different
// between row vectors and column vectors.  I find this helps students get
// the math to work out correctly if you pay attention to this difference.
// The routines allow you to name a matrix.  The name is then used in debug
// output.  Other things checked include referencing out of bounds.
//

class Matrix {
friend class MatrixRowIter;
public:
    static bool debug;      // debugging flag

private:
    bool defined;           // does it have rows and cols defined
    bool submatrix;         // if submatrix then it does NOT own the row content of m (see deallocate)!!
    int maxr, maxc;
    double **m;             // the data
    std::string name;       // the name of the matrix or ""

private:  // private methods
    void allocate(int r, int c, std::string namex, bool isSubMatrix=false);
    bool deallocate();
    void reallocate(int othermaxr, int othermaxc, std::string namex);

// constructors
public:
    Matrix(std::string namex="");
    Matrix(int r, std::string namex="");                            // create a subMatrix columns unallocated
    Matrix(int r, int c, std::string namex="");
    Matrix(int r, int c, double initValue, std::string namex="");   // create and initialize
    Matrix(int r, int c, double *data, std::string namex="");       // create and initialize from array
    Matrix(const Matrix &other, std::string namex="");              // copy constructor
    Matrix(Matrix *other);                                          // for convenience
    ~Matrix();
    Matrix &operator=(const Matrix &other);

// basic error checking support
public:
    void checkBounds(int r, int c, std::string msg) const;
    void assertColVector(std::string) const;
    void assertColsEqual(const Matrix &other, std::string msg) const;
    void assertDefined(std::string msg) const;
    void assertRowIndexOK(int r, std::string msg) const;
    void assertColIndexOK(int c, std::string msg) const;
    void assertIndexOK(int, int, std::string) const;
    void assertOtherLhs(const Matrix &other, std::string msg) const;
    void assertOtherSizeMatch(const Matrix &other, std::string msg) const;
    void assertRowVector(std::string) const;
    void assertRowsEqual(const Matrix &other, std::string msg) const;
    void assertSize(int r, int c, std::string msg) const;
    void assertUsableSize(std::string msg) const;
    void assertSquare(std::string msg) const;

public:  // auxillary routines but not private (for speed, they do not check self!!)
    void swapRows(int i, int j);                  // utility to swap two rows
    bool lessRows(int i, int j) const;            // utility to compare two rows

// accessors
public: 
    int maxRows() const { return maxr; }  // DEPRICATED
    int maxCols() const { return maxc; }  // DEPRICATED
    int numRows() const { return maxr; }
    int numCols() const { return maxc; }
    double get(int r, int c) const;      // get element value
    double inc(int r, int c);            // increment element
    double dec(int r, int c);            // decrement element
    double set(int r, int c, double v);  // set element
    void setDefined();                   // make defined when you *KNOW* the array has been defined by other means
    void setName(std::string newName);   // set matrix name
    const std::string &getName(const std::string &defaultName="") const;
    void narrow(int newMaxCol);          // remove trailing columns (without proper deallocation)
    void shorten(int newMaxRow);         // remove trailing rows (without proper deallocation)

// operations on matrices
public: 
    bool isRowVector() const { return defined && maxr==1; }  // exactly one row
    bool isColVector() const { return defined && maxc==1; }  // exactly one col
    bool equal(const Matrix &other) const;       // are the two matrices equal?
    bool nearEqual(double epsilon, const Matrix &other) const; // matrices nearly equal?
    int countGreater(const Matrix &other) const; // count number of elements >
    void argMax(int &r, int &c) const;           // what location is the largest in whole array
    void argMin(int &r, int &c) const;           // what location is the smallest in whole array
    Matrix argMinRow() const;                    // constructs a column vector of the argmin in each row
    Matrix minRow() const;                       // constructs a column vector of the min in each row
    double max() const;                          // minimum in whole array
    double min() const;                          // maximum in whole array
    double mean() const;                         // mean of whole array
    double maxCol(int c) const;                  // maximum value in a column
    double minCol(int c) const;                  // minimum value in a column
    double meanCol(int c) const;                 // mean in a column
    double stddevCol(int c) const;               // standard deviation in a column
    int countEqCol(int c, double value) const;   // count number of items in column c equal to value
    int countNeqCol(int c, double value) const;  // count number of items in column c not equal to value
    double dist2(const Matrix &other) const;     // *SQUARE* of distance between two matrices
    Matrix pickRows(int match, const Matrix &list, int &num);      // pick rows which have list value == match
    double dot(int r, int c, const Matrix &other) const;  // dot of row of this with col of other -> double 
    double dist2(int r, int c, const Matrix &other) const;  // *SQUARE* of distance between row of this with col of other

    // element by element operators (modifies self)
    Matrix &abs();
    Matrix &add(const Matrix &other);
    Matrix &sub(const Matrix &other);
    Matrix &mult(const Matrix &other);
    Matrix &div(const Matrix &other);

    Matrix &swap(Matrix &other);    // swaps two matrices so also modifies other
    Matrix &rowInc(int r);          // increment the values in a given row by 1

    // scalar operators
    Matrix &constant(double x);            // this can be used to zero a matrix
    Matrix &constantDiagonal(double x);    // this can be used to set the diagonal to a constant but does set rest of matrix
    Matrix &constantCol(int c, double x);  // this can be used to zero a column
    Matrix &constantColRange(int c, double start, double step);   // assign all elements in col range starting at start and going by step
    Matrix &identity();                    // convert to an identity matrix  (must be square)
    Matrix &scalarMult(double x);          // multiply all elements by x
    Matrix &scalarAdd(double x);           // add to all elements x
    Matrix &scalarPreSub(double x);        // NOTE: this is x - self   not   self - x, can be used to negate
    Matrix &scalarPostSub(double x);       // NOTE: this is self - x

    // Vector operations
    Matrix &divColVector(const Matrix &other);
    Matrix &multColVector(const Matrix &other);

    Matrix &divRowVector(const Matrix &other);  // self[r] / (row vector other) for each row
    Matrix &multRowVector(const Matrix &other); // self[r] * (row vector other) for each row
    Matrix &addRowVector(const Matrix &other);  // self[r] + (row vector other) for each row
    Matrix &subRowVector(const Matrix &other);  // self[r] - (row vector other) for each row
    Matrix &addRowVector(int r, const Matrix &other); // add row vector matrix in other to the given row of self
    
    // min/max normalization by columns
    Matrix normalizeCols();                       // normalize and return array of min and max of each col
    Matrix &normalizeCols(Matrix &minMax);        // normalize based on an array of min and max for each col

    // mapping functions
    Matrix &map(double (*f)(double x));              // apply given function to all elements
    Matrix &mapCol(int c, double (*f)(double x));    // apply given function to all elements in col c
    Matrix &mapIndex(double (*f)(int r, int c, double x)); // apply function to (index, element)
    Matrix cartesianRow(double (*)(int, double*, double*), Matrix&);  // apply given function to the cartesian product of two vectors of row vectors

    // random initialization (random number generator must be initialized with initRand() )
    Matrix &randCol(int c, double min, double max);  // random reals in a column
    Matrix &rand(double min, double max);            // random reals in range 
    Matrix &rand(int min, int max);                  // random integers in range

    // insertion and extraction
    Matrix &sample(Matrix &out);  // extract random rows with replacement into existing matrix out
    Matrix &extract(int minr, int minc, int sizer, int sizec, Matrix &out);  // extract into existing matrix out (see other versions of extract)
    Matrix &insert(const Matrix &other, int minr, int minc);    // insert the matrix at minr, minc.   Overflow is ignored.
    Matrix &insertRowVector(int row, const Matrix&);

    // input/output
    void print(std::string msg="") const;      // print matrix and its name
    void printInt(std::string msg="") const;   // print matrix as integers (it will error if not.)
    void printLabeledRow(char **labels, std::string msg="") const;   // print matrix with rows labeled
    void printSize(std::string msg="") const;  // print just the matrix size
    void write();                        // write out matrix in a form that can be read in
    void writeLine(int r);               // write out an unadorned row
    void read();                         // read in a matrix
    char **readLabeledRow();             // read in a matrix plus row labels

private:
    char **readAux(bool labeled);        // general read routine both labeled and un

public:
    // WARNING: The following CONSTRUCT TO NEW MATRIX for the answer  (BEWARE MEMORY LEAKS!)
    // NOTE: the result of these functions should be used somewhere like in an assignment
    //  e.g. a.dot(b) is probably wrong.   while x = a.dot(b) stores the result.
    Matrix extract(int minr, int minc, int sizer, int sizec);
    Matrix extractStride(int minr, int minc, int stepr, int stepc);
    Matrix transpose();                    // classic transpose into new matrix (see transposeSelf below)
    Matrix dot(const Matrix &other);       // classic matrix multiply, inner product
    Matrix dotT(const Matrix &other);      // classic matrix multiply self * Transpose(other)
    Matrix Tdot(const Matrix &other);      // classic matrix multiply Transpose(self) * other
    Matrix meanVec();                      // creates a row vector of means of columns
    Matrix stddevVec();                    // creates a row vector of standard deviations of columns
    Matrix cov();                          // covariance matrix (BIASED covariance)
    Matrix cov(Matrix &other);             // covariance matrix (BIASED covariance)

    // alternation versions of operaters that DO NOT create new matrices
    Matrix &transposeSelf();                // transpose in place of SQUARE MATRIX

    // special operators (destroys arguments)
    int *LU();                              // LU decomposition in place
    Matrix &solve(Matrix &B);               // solve Ax = B returns solutions and inverse
    Matrix &inverse();                      // replace with inverse

    // eigenSystem() destroys self by replacing self with eigenvectors in rows.
    // Returns a new matrix with the eigenvalues in it.
    // Eigenvalues and vectors returned sorted from largest magnitude to smallest
    // WARNING: allocates new matrix for answer
    void tridiagonalize(double *&d, double *&e);
    Matrix eigenSystem();

    // sorting support
public:
    void selectSort(int lower, int upper);
    void qs(int lower, int upper);
    void selectSortCol(int c, int lower, int upper);
    void qsCol(int c, int lower, int upper);

public: 
    void sortRows();                            // sort rows in place
    void sortRows(int startRow, int endRow);    // sort rows in place in a range of rows
    void sortRowsByCol(int c);                  // sort rows in place on given column
    void sortRowsByCol(int c, int startRow, int endRow);     // sort rows in place in a range of rows

public:
    // subMatrices are an efficiency for creating submatrices by pointing into
    // the parent matrix rather than copying all the contents of the matrix.
    // Read the warnings in the .cpp file
    // 
    Matrix subMatrix(int minr, int minc, int sizer, int sizec) const;  // create a submatrix whose corner is (minr, minc) and size given
    Matrix subMatrixEq(int c, double value) const;         // create submatrix with rows whose column c has the given value
    Matrix subMatrixNeq(int c, double value) const;        // create submatrix with rows whose column c does not have the given value

    // image (picture) support (currently only supports 8 bit pgm and ppm formats)
    // output is in ascii formats (zzz: fix someday to use more compressed output)
    // 8 bit gray is one integer in the range 0-255 for each pixel
    // 8 bit color is three integers in a row in the range 0-255 for RGB in each pixel.
    // That is an 8 bit color square 100x100 pixels gens a 100x300 dimensional array
private:
    int byteValue(double x);
    Matrix readImage(char *expectedType, char *caller, std::string filename, std::string namex);

public: 
    Matrix readImagePgm(std::string filename, std::string namex);   // read a P2 or P5 pgm  (8 bit gray scale) file into self
    Matrix readImagePpm(std::string filename, std::string namex);   // read a P3 or P6 ppm  (8 bit color)
    void writeImagePgm(std::string filename, std::string comment);  // write a P2 pgm file  (8 bit gray scale)
    void writeImagePpm(std::string filename, std::string comment);  // write a P3 pgm file (8 bit color)
};

#endif


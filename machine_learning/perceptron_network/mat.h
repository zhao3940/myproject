#ifndef MATH
#define MATH

// // // // // // // // // // // // // // // // 
//
// Some simple matrix operators.
//
//
// These routines are self contained except for random number
// generation.   Comment these out if you want completely stand-alone
// functionality.  Uses some code from Numerical Recipes (see comments).
// This code is not designed to be blindingly fast but rather just get the job done.
//
// Author: Robert B. Heckendorn, University of Idaho, 2017
// Version: 2.3
// Date: Mar 7, 2017

// IMPORTANT: If running on MICROSOFT WINDOWS uncomment this define statement!
// #define WINDOWS

#include <string>
#include "rand.h"

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
// A simple class for matrix operations.   It has some nice debugging features like
// trying hard to check that the proper dimensions are used.   It is very draconian
// about this so there is an important different between row vectors and column vectors.
// I find this helps students get the math to work out correctly if you pay attention to
// this difference.   The routines allow you to name a matrix.  The name is then used
// in debug output.   Other things checked include referencing out of bounds.
//
class Matrix {
friend class MatrixRowIter;
private:
    bool defined;           // does it have rows and cols defined
    int maxr, maxc;
    double **m;             // the data
    std::string name;       // the name of the matrix or ""

private:  // private methods
    void allocate(int r, int c, std::string namex);
    bool deallocate();
    void reallocate(int othermaxr, int othermaxc, std::string namex);

// constructors
public:
    Matrix(std::string namex="");
    Matrix(int r, int c, std::string namex="");
    Matrix(int r, int c, double *data, std::string namex="");
    Matrix(const Matrix &other, std::string namex="");         // copy constructor
    Matrix(Matrix *other);                                     // just for convenience
    ~Matrix();
    Matrix &operator=(const Matrix &other);

// basic error checking support
public:
    void checkBounds(int r, int c, std::string msg) const;
    void assertColVector(std::string) const;
    void assertColsEqual(const Matrix &other, std::string msg) const;
    void assertDefined(std::string msg) const;
    void assertIndexOK(int, int, std::string) const;
    void assertOtherLhs(const Matrix &other, std::string msg) const;
    void assertOtherSizeMatch(const Matrix &other, std::string msg) const;
    void assertRowVector(std::string) const;
    void assertRowsEqual(const Matrix &other, std::string msg) const;
    void assertSize(int r, int c, std::string msg) const;
    void assertSquare(std::string msg) const;

public:  // auxillary routines but not private (they do not check their arguments)
    void swapRows(int i, int j);            // utility to swap two rows
    bool lessRows(int i, int j) const;            // utility to compare two rows

// accessors
public: 
    int maxRows() const { return maxr; }  // DEPRICATED
    int maxCols() const { return maxc; }  // DEPRICATED
    int numRows() const { return maxr; }
    int numCols() const { return maxc; }
    double get(int r, int c);            // get element value
    double inc(int r, int c);            // increment element
    double dec(int r, int c);            // decrement element
    double set(int r, int c, double v);  // set element
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
    Matrix &rowInc(int r);
    Matrix &rowAdd(int r, const Matrix &other);

    // scalar operators
    Matrix &constant(double x);            // this can be used to zero a matrix
    Matrix &constantDiagonal(double x);    // this can be used to set the diagonal to a constant but does set rest of matrix
    Matrix &constantCol(int c, double x);  // this can be used to zero a column
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
    
    // min/max normalization by columns
    Matrix normalizeCols();                       // normalize and return array of min and max of each col
    Matrix &normalizeCols(Matrix &minMax);        // normalize based on an array of min and max for each col

    // mapping functions
    Matrix &map(double (*f)(double x));              // apply given function to all elements
    Matrix &mapCol(int c, double (*f)(double x));    // apply given function to all elements in col c
    Matrix &mapIndex(double (*f)(int r, int c, double x)); // apply function to (index, element)
    Matrix cartesianRow(double (*)(int, double*, double*), Matrix&);  // apply given function to the cartesian product of two vectors of row vectors

    // random initialization
    Matrix &randCol(int c, double min, double max);  // random reals in a column
    Matrix &rand(double min, double max);            // random reals in range 
    Matrix &rand(int min, int max);                  // random integers in range

    // insertion and extraction
    Matrix &sample(Matrix &out);  // extract random rows with replacement into existing matrix out
    Matrix &extract(int minr, int minc, int sizer, int sizec, Matrix &out);  // extract into existing matrix out
    Matrix &insert(const Matrix &other, int minr, int minc);
    Matrix &insertRowVector(int row, const Matrix&);

    // input/output
    void print(std::string msg="");
    void printSize(std::string msg="");
    void write();
    void writeLine(int r);
    void read();

    // WARNING: The following CONSTRUCT TO NEW MATRIX for the answer  (BEWARE MEMORY LEAKS!)
    // NOTE: the result of these functions should be used somewhere like in an assignment
    //  e.g. a.dot(b) is probably wrong.   while x = a.dot(b) stores the result.
    Matrix extract(int minr, int minc, int sizer, int sizec);
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

    // eigen systems (destroys arguments and creates a new matrix)
    void tridiagonalize(double *&d, double *&e);
    Matrix eigenSystem();

    // sort
private: 
    void selectSort(int lower, int upper);
    void qs(int lower, int upper);

public: 
    void sortRows();
    };

#endif

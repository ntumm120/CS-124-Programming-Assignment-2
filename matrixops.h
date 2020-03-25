#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <math.h>

using namespace std;

class MatrixOps{
    public:
        void make(vector<vector<int> > &first, vector<vector<int> > &second, int n);
        void makeidentity(vector<vector<int> > &first, vector<vector<int> > &second, int n);
        void add(vector<vector<int> > &first, vector<vector<int> > &second,
             vector<vector<int> > &result, int n);
        void sub(vector<vector<int> > &first, vector<vector<int> > &second,
             vector<vector<int> > &result, int n);
        int powerround(int n);
        void pad(vector<vector<int> > &padfirst, vector<vector<int> > &padsecond,
                       vector<vector<int> > &first, vector<vector<int> > &second, int n);
        MatrixOps();
        int findPad(int n, int crossover);
        void makepad(vector<vector<int> > &first, int newdim);
        void removepad(vector<vector<int> > &first, int newdim);
        void printDiagonal(vector<vector<int> > &first, int dimension);
        void combine (std::vector< std::vector<int> > &A,
                              std::vector< std::vector<int> > &B,
                              int row, int col, int d);
        void divide(std::vector< std::vector<int> > &A, std::vector< std::vector<int> > &B,int row, int col, int d);
};

MatrixOps::MatrixOps() {
}

void MatrixOps::make(vector<vector<int> > &first, vector<vector<int> > &second, int n){
    // initialize matrices with 0s and 1s, can change later to 0, 1, -1
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            first[i][j] = rand() % 2;
            second[i][j] = rand() % 2;
        }
    }
}

void MatrixOps::makeidentity(vector<vector<int> > &first, vector<vector<int> > &second, int n){

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j){
                first[i][j] = 1;
                second[i][j] = 1;
            }
            else{
                first[i][j] = 0;
                second[i][j] = 0;
            }
        }
    }
}

void MatrixOps::add(vector<vector<int> > &first, vector<vector<int> > &second,
        vector<vector<int> > &result, int n){
    
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            result[i][j] = first[i][j] + second[i][j];
        }
    }
}

void MatrixOps::sub(vector<vector<int> > &first, vector<vector<int> > &second,
              vector<vector<int> > &result, int n){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            result[i][j] = first[i][j] - second[i][j];
        }
    }
}

int MatrixOps::powerround(int n){
    return pow(2, int(ceil(log2(n))));
}

void MatrixOps::pad(vector<vector<int> > &padfirst, vector<vector<int> > &padsecond,
                          vector<vector<int> > &first, vector<vector<int> > &second, int n){

    int paddedn = powerround(n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            padfirst[i][j] = first[i][j];
            padsecond[i][j] = second[i][j];
        }
    }
    for (int i = n; i < paddedn; i++){
        for (int j = n; j < paddedn; j++){
            padfirst[i][j] = 0;
            padsecond[i][j] = 0;
        }
    }
}

int MatrixOps::findPad(int n, int crossover){
    int temp = 0;
    while (n > crossover) {
        if (n % 2 == 0) {
            n /= 2;
        }
        else {
            n = (n + 1) / 2;
            temp++;
        }
    }
    return (n*pow(2,temp));
}


void MatrixOps::makepad(vector<vector<int> > &first, int newdim){
    first.resize(newdim);
    for (int i = 0; i < newdim; i++){
        first[i].resize(newdim);
    }
}

void MatrixOps::removepad(vector<vector<int> > &first, int newdim){
    first.resize(newdim);
    for (int i = 0; i < newdim; i++){
        first.resize(newdim);
    }
}

void MatrixOps::printDiagonal(vector<vector<int> > &first, int dimension){
    for (int i = 0; i < dimension; ++i) {
        printf("%d\n", first[i][i]);
    }
}


void MatrixOps::divide (std::vector< std::vector<int> > &A,
            std::vector< std::vector<int> > &B,
            int row, int col, int d)
{
    for (int i1 = 0, i2 = row; i1 < d; i1++, i2++)
    {
        for (int j1 = 0, j2 = col; j1 < d; j1++, j2++)
        {
            B[i1][j1] = A[i2][j2];
        }
    }
}

// combine the 4 matrices in a resulting matrix
void MatrixOps::combine (	std::vector< std::vector<int> > &A,
               std::vector< std::vector<int> > &B,
               int row, int col, int d)
{
    for (int i1 = 0, i2 = row; i1 < d; i1++, i2++)
    {
        for (int j1 = 0, j2 = col; j1 < d; j1++, j2++)
        {
            B[i2][j2] = A[i1][j1];
        }
    }
}

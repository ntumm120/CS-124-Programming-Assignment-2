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
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            padfirst[i][j] = first[i][j];
            padsecond[i][j] = second[i][j];
        }
    }
}
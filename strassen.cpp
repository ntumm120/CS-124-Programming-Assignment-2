#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include "matrixops.h"

using namespace std;

void conventional(vector<vector<int> > &first,
                         vector<vector<int> > &second, vector<vector<int> > &result, int n){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int dot  = 0;
            for (int k = 0; k < n; k++) {
                dot += first[i][k] * second[k][j];
            }
            result[i][j] = dot;
        }
    }
    for (int i = 0; i < n; i++){
        cout << result[i][i] << '\n';
    }
}

//To execute normal Strassen let crossover = 0, currently have crossover in command-line

void strassen(vector<vector<int> > &first, vector<vector<int> > &second, vector<vector<int> > &result, int n, int crossover){

    MatrixOps helper = *(new MatrixOps());
    if(n <= crossover){
        conventional(first, second, result, n);
    }
    else{
        int paddedsize = helper.powerround(n);
        int blocksize = paddedsize/2;

        vector<int> pads(paddedsize); 
        vector<vector<int> > padfirst(paddedsize, pads);
        vector<vector<int> > padsecond(paddedsize, pads);
        vector<vector<int> > padresult(paddedsize, pads);

        if(n != paddedsize){
            helper.pad(padfirst, padsecond, first, second, n);
        }
        else{
            padfirst = first;
            padsecond = second;
        }
        
        vector<int> block(blocksize);

        vector<vector<int> > firstq1(blocksize, block), firstq2(blocksize, block);
        vector<vector<int> > firstq3(blocksize, block), firstq4(blocksize, block);
        vector<vector<int> > secondq1(blocksize, block), secondq2(blocksize, block);
        vector<vector<int> > secondq3(blocksize, block), secondq4(blocksize, block);
        vector<vector<int> > thirdq1(blocksize, block), thirdq2(blocksize, block);
        vector<vector<int> > thirdq3(blocksize, block), thirdq4(blocksize, block);

        vector<vector<int> > ABlock(blocksize, block);
        vector<vector<int> > BBlock(blocksize, block);

        vector<vector<int> > m1(blocksize, block);
        vector<vector<int> > m2(blocksize, block);
        vector<vector<int> > m3(blocksize, block);
        vector<vector<int> > m4(blocksize, block);
        vector<vector<int> > m5(blocksize, block);
        vector<vector<int> > m6(blocksize, block);
        vector<vector<int> > m7(blocksize, block);

        //partitioning matrix into its quadrants/blocks
        for (int i=0; i<blocksize; i++){
            for (int j=0; j<blocksize; j++){
                firstq1[i][j] = padfirst[i][j];
                firstq2[i][j] = padfirst[i][j+blocksize];
                firstq3[i][j] = padfirst[i+blocksize][j];
                firstq4[i][j] = padfirst[i+blocksize][j+blocksize];
                secondq1[i][j] = padsecond[i][j];
                secondq2[i][j] = padsecond[i][j+blocksize];
                secondq3[i][j] = padsecond[i+blocksize][j];
                secondq4[i][j] = padsecond[i+blocksize][j+blocksize];
            }
        }

        //Constructing submatrices
        helper.add(firstq1, firstq4, ABlock, blocksize);
        helper.add(secondq1, secondq4, BBlock, blocksize);
        strassen(ABlock, BBlock, m1, blocksize, crossover);

        helper.add(firstq3, firstq4, ABlock, blocksize);
        strassen(ABlock, secondq1, m2, blocksize, crossover);

        helper.subtract(secondq2, secondq4, BBlock, blocksize);
        strassen(firstq1, BBlock, m3, blocksize, crossover);

        helper.subtract(secondq3, secondq1, BBlock, blocksize);
        strassen(firstq4, BBlock, m4, blocksize, crossover);

        helper.add(firstq1, firstq2, ABlock, blocksize);
        strassen(ABlock, secondq4, m5, blocksize, crossover);

        helper.subtract(firstq3, firstq1, ABlock, blocksize);
        helper.add(secondq1, secondq2, BBlock, blocksize);
        strassen(ABlock, BBlock, m6, blocksize, crossover);

        helper.subtract(firstq2, firstq4, ABlock, blocksize);
        helper.add(secondq2, secondq4, BBlock, blocksize);
        strassen(ABlock, BBlock, m7, blocksize, crossover);

        helper.add(m1, m4, ABlock, blocksize);
        helper.subtract(ABlock, m5, BBlock, blocksize);
        helper.add(BBlock, m7, secondq1, blocksize);

        helper.add(m3, m5, thirdq2, blocksize);

        helper.add(m2, m4, thirdq3, blocksize);

        helper.subtract(m1, m2, ABlock, blocksize);
        helper.add(ABlock, m3, BBlock, blocksize);
        helper.add(BBlock, m6, thirdq4, blocksize);

        //calculate result matrix
        for(int i=0; i<blocksize; i++){
            for(int j=0; j<blocksize; j++){
                padresult[i][j] = thirdq1[i][j];
                padresult[i][blocksize+j] = thirdq2[i][j];
                padresult[blocksize+i][j] = thirdq3[i][j];
                padresult[blocksize+i][blocksize+j] = thirdq4[i][j];
            }
        }

        //remove additional values from expanded matrix
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                result[i][j] = padresult[i][j];
            }
        }
    }

    for (int i = 0; i < n; i++){
       cout << result[i][i] << '\n';
    }

}

int main(int argc, char *argv[])
{
    //Need to write code to read a file of numbers into a matrix, for now we will use
    // randomly generated numbers to test
    
    //Still need to figure out experimental crossover point, possibly using time library,
    // not sure how to calculate numbers of operations

    long double startTime;
    long double conventionaltime;
    long double normaltime;
    long double varianttime;

    if(argc!= 3){
         cout << "You must have 3 command-line arguments." << "\n";
         return - 1;
    }

    int crossover = atoi(argv[1]);
    int dim = atoi(argv[2]);

    vector<vector<int> > first(dim,vector<int>(dim));
    vector<vector<int> > second(dim,vector<int>(dim));
    vector<vector<int> > result(dim,vector<int>(dim));

    MatrixOps help = *(new MatrixOps());
    help.make(first,second,dim);

    //Have to fix timing stuff, not sure if it works

    startTime = time(0);
    conventional(first,second,result,dim);
    conventionaltime = time(0) - startTime;
    
    //No crossover strassen (DOESNT WORK NEED TO FIGURE OUT WHY)
    //startTime = time(0);
    //strassen(first, second, result ,dim, 0);
    //normaltime = time(0) - startTime;
    
    //Crossover strassen
    startTime = time(0);
    strassen(first, second, result ,dim, crossover);
    varianttime = time(0) - startTime;


    cout << "Standard"<< "\n";
    cout << "Size "
         << dim << " : " << conventionaltime << "\n";
    cout << "\n";

    cout << "Normal strassen"<< "\n";
    cout << "Size "
         << dim << " : " << normaltime << "\n";
    cout << "\n";

    cout << "Variant strassen"<< "\n";
    cout << "Crossover used: "<< crossover << "\n";
    cout << "Size "
         << dim << " : " << varianttime << "\n";
    cout << "\n";

    return 0;
}
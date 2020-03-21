#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include "matrixops.h"
#include <fstream>

using namespace std;

void conventional(vector<vector<int> > &first,
                         vector<vector<int> > &second, vector<vector<int> > &result, int n){
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int dot = 0;
            for (int k = 0; k < n; ++k) {
                dot += first[i][k] * second[k][j];
            }
            result[i][j] = dot;
        }
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

        vector<vector<int> > p1(blocksize, block);
        vector<vector<int> > p2(blocksize, block);
        vector<vector<int> > p3(blocksize, block);
        vector<vector<int> > p4(blocksize, block);
        vector<vector<int> > p5(blocksize, block);
        vector<vector<int> > p6(blocksize, block);
        vector<vector<int> > p7(blocksize, block);

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
        strassen(ABlock, BBlock, p1, blocksize, crossover);

        helper.add(firstq3, firstq4, ABlock, blocksize);
        strassen(ABlock, secondq1, p2, blocksize, crossover);

        helper.sub(secondq2, secondq4, BBlock, blocksize);
        strassen(firstq1, BBlock, p3, blocksize, crossover);

        helper.sub(secondq3, secondq1, BBlock, blocksize);
        strassen(firstq4, BBlock, p4, blocksize, crossover);

        helper.add(firstq1, firstq2, ABlock, blocksize);
        strassen(ABlock, secondq4, p5, blocksize, crossover);

        helper.sub(firstq3, firstq1, ABlock, blocksize);
        helper.add(secondq1, secondq2, BBlock, blocksize);
        strassen(ABlock, BBlock, p6, blocksize, crossover);

        helper.sub(firstq2, firstq4, ABlock, blocksize);
        helper.add(secondq2, secondq4, BBlock, blocksize);
        strassen(ABlock, BBlock, p7, blocksize, crossover);

        helper.add(p1, p4, ABlock, blocksize);
        helper.sub(ABlock, p5, BBlock, blocksize);
        helper.add(BBlock, p7, secondq1, blocksize);

        helper.add(p3, p5, thirdq2, blocksize);

        helper.add(p2, p4, thirdq3, blocksize);

        helper.sub(p1, p2, ABlock, blocksize);
        helper.add(ABlock, p3, BBlock, blocksize);
        helper.add(BBlock, p6, thirdq4, blocksize);

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

}

void findOptimalconvThreshold() {

    for (int convThreshold = 8; convThreshold <= 256; convThreshold*=2){
        //for (convThreshold = 2; convThreshold <= 64; convThreshold++){
        double total = 0;
        //cout << "multiplying matrices, n = " << i << endl;
        for (int j = 0; j < 5; j ++){
            vector<vector<int> > &first();
            vector<vector<int> > &second();
            MatrixOps help = *(new MatrixOps());
            help.makeidentity(first, 256);
            help.makeidentity(second, 256);

            clock_t start;
            start = clock();

            vector<vector<int> > &result() = multiply(m1, m2);
            total += (std::clock() - start) / (double)(CLOCKS_PER_SEC);
            delete(m1);
            delete(m2);
            delete(m3);
        }

        cout << convThreshold << "\t" << total / 5 << endl;

        //cout << "finished multiplying.\n" << endl;

    }
}

vector<vector<int> >& FileToMatrix(char* inputfile, int dimension, int order){
    vector<vector<int>> matrix(dimension, vector<int>(dimension));

    ifstream inFile(inputfile);
    string line;
    for (int i = 0; i < pow(dimension,2)*order; ++i)
    {
        getline(inFile, line);
    }

    for (int i = 0; i < dimension; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            getline(inFile, line);
            matrix[i][j] = stoi(line);
        }
    }
    inFile.close();
    return matrix;
}

int main(int argc, char *argv[])
{
    //Wrote function to read in a file of numbers, for now we will use
    // randomly generated numbers to test
    
    //Still need to figure out experimental crossover point, possibly using time library,
    // not sure how to calculate number of operations

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


    vector<vector<int> > test(dim,vector<int>(dim));
    help.makeidentity(test, test, dim);


    //Have to fix timing stuff, not sure if it works

    //Conventional fully works on test matrices
    startTime = time(0);
    //conventional(test,test,result, dim);
    conventionaltime = 1000 * (time(0) - startTime);
    
    //No crossover strassen (DOESNT WORK NEED TO FIGURE OUT WHY)
    startTime = time(0);
    //strassen(first, second, result, dim, 0);
    normaltime = 1000 * (time(0) - startTime);
    
    //Crossover strassen doesn't work either
    startTime = time(0);
    strassen(test, test, result, dim, crossover);
    varianttime = 1000 * (time(0) - startTime);


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
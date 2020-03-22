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

//To execute normal Strassen let crossover = 0, currently have crossover in command-line, this strassen
//doesn't split up the quadrants in signature of method
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

//void multiply(vector<vector<int> > &, vector<vector<int> > &, vector<vector<int> > &, int , int , int , int , int , int , int);

void strassenMult(vector<vector<int> > &A, vector<vector<int> > &B, vector<vector<int> > &C,
        int topA, int leftA, int topB, int leftB, int topC, int leftC, int dimension, int crossover) {

    if (dimension > crossover) {
        MatrixOps helper = *(new MatrixOps());
        // C12 = A21 - A11
        helper.sub(A, B, C, topA + dimension / 2, leftA, topA, leftA, topC, leftC + dimension / 2, dimension / 2);
        // C21 = B11 + B12
        helper.add(B, B, C, topB, leftB, topB, leftB + dimension / 2, topC + dimension / 2, leftC, dimension / 2);
        // C22 = C12 * C21

        strassenMult(C, C, C, topC, leftC + dimension / 2, topC + dimension / 2, leftC, topC + dimension / 2,
                 leftC + dimension / 2, dimension / 2, crossover);

        //C12 = A12 - A22
        helper.sub(A, A, C, topA, leftA + dimension / 2, topA + dimension / 2, leftA + dimension / 2, topC,
                   leftC + dimension / 2, dimension / 2);
        //C21 = B21 + B22
        helper.add(B, B, C, topB + dimension / 2, leftB, topB + dimension / 2, leftB + dimension / 2,
                   topC + dimension / 2, leftC, dimension / 2);
        //C11 = C12 * C21

        strassenMult(C, C, C, topC, leftC + dimension / 2, topC + dimension / 2, leftC, topC, leftC, dimension / 2, crossover);

        //C12 = A11 + A22
        helper.add(A, A, C, topA, leftA, topA + dimension / 2, leftA + dimension / 2, topC, leftC + dimension / 2,
                   dimension / 2);
        //C21 = B11 + B22
        helper.add(B, B, C, topB, leftB, topB + dimension / 2, leftB + dimension / 2, topC + dimension / 2, leftC,
                   dimension / 2);

        vector<vector<int> > newA(dimension / 2, vector<int>(dimension / 2));
        vector<vector<int> > newB(dimension / 2, vector<int>(dimension / 2));
        helper.make(newA, newB, dimension / 2); // TODO deal with non-power of 2 case


        //newA = C12*C21
        strassenMult(C, C, newA, topC, leftC + dimension / 2, topC + dimension / 2, leftC, 0, 0, dimension / 2, crossover);
        //C11 = newA + C11
        helper.add(newA, C, C, 0, 0, topC, leftC, topC, leftC, dimension / 2);
        //C22 = newA + C22
        helper.add(newA, C, C, 0, 0, topC + dimension / 2, leftC + dimension / 2, topC + dimension / 2,
                   leftC + dimension / 2, dimension / 2);

        //newB = A21 + A22
        helper.add(A, A, newB, topA + dimension / 2, leftA, topA + dimension / 2, leftA + dimension / 2, 0, 0,
                   dimension / 2);
        //C21 = newB * B11

        strassenMult(newB, B, C, 0, 0, topB, leftB, topC + dimension / 2, leftC, dimension / 2, crossover);

        //C22 = C22 - C21
        helper.sub(C, C, C, topC + dimension / 2, leftC + dimension / 2, topC + dimension / 2, leftC,
                   topC + dimension / 2, leftC + dimension / 2, dimension / 2);
        //newA = B21 - B11
        helper.sub(B, B, newA, topB + dimension / 2, leftB, topB, leftB, 0, 0, dimension / 2);
        //newB = A22 * newA

        strassenMult(A, newA, newB, topA + dimension / 2, leftA + dimension / 2, 0, 0, 0, 0, dimension / 2, crossover);

        //C21 = C21 + newB
        helper.add(C, newB, C, topC + dimension / 2, leftC, 0, 0, topC + dimension / 2, leftC, dimension / 2);
        //C11 = C11 + newB
        helper.add(C, newB, C, topC, leftC, 0, 0, topC, leftC, dimension / 2);
        //newA = B12 - B22
        helper.sub(B, B, newA, topB, leftB + dimension / 2, topB + dimension / 2, leftB + dimension / 2, 0, 0,
                   dimension / 2);
        //C12 = A11 * newA

        strassenMult(A, newA, C, topA, leftA, 0, 0, topC, leftC + dimension / 2, dimension / 2, crossover);

        //C22 = C22 + C12
        helper.add(C, C, C, topC + dimension / 2, leftC + dimension / 2, topC, leftC + dimension / 2,
                   topC + dimension / 2, leftC + dimension / 2, dimension / 2);
        //newB = A11 + A12
        helper.add(A, A, newB, topA, leftA, topA, leftA + dimension / 2, 0, 0, dimension / 2);
        //newA = newB * B22

        strassenMult(newB, B, newA, 0, 0, topB + dimension / 2, leftB + dimension / 2, 0, 0, dimension / 2, crossover);

        //C12 = C12 + newA
        helper.add(C, newA, C, topC, leftC + dimension / 2, 0, 0, topC, leftC + dimension / 2, dimension / 2);
        //C11 = C11 - newA
        helper.sub(C, newA, C, topC, leftC, 0, 0, topC, leftC, dimension / 2);
    }
    else{
        conventional(A, B, C, dimension);
    }
}

/**void multiply(vector<vector<int> > &A, vector<vector<int> > &B, vector<vector<int> > &C, int topA, int leftA, int topB, int leftB, int topC, int leftC, int dimension, int crossover) {
    if (dimension > crossover) {
        strassenMult(A, B, C, topA, leftA, topB, leftB, topC, leftC, dimension);
    }
    else {
        conventional(A, B, C, dimension);
    }
}**/

vector<vector<int> > & multiply(vector<vector<int> > &A, vector<vector<int> > &B, vector<vector<int> > &C, int crossover){

    MatrixOps helper = *(new MatrixOps());
    int dimension = A.size();

    int padding = helper.findPad(dimension, crossover);
    helper.makepad(A, padding);
    helper.makepad(B, padding);
    helper.makepad(C, padding);
    strassenMult(A,B,C,0,0,0,0,0,0,padding, crossover);

    helper.removepad(A, dimension);
    helper.removepad(B, dimension);
    helper.removepad(C, dimension);
    return C;
}

void findOptimalconvThreshold() {

    for (int convThreshold = 8; convThreshold <= 256; convThreshold*=2){
        double total = 0;
        for (int j = 0; j < 5; j ++){
            vector<vector<int> > first(256, vector<int>(256));
            vector<vector<int> > second(256, vector<int>(256));
            vector<vector<int> > third(256, vector<int>(256));
            MatrixOps help = *(new MatrixOps());
            help.makeidentity(first, second, 256);

            clock_t start;
            start = clock();

            vector<vector<int> > result(256, vector<int>(256));
            result = multiply(first, second, third, convThreshold);
            total += (std::clock() - start) / (double)(CLOCKS_PER_SEC);
        }

        cout << convThreshold << "\t" << total / 5 << endl;

    }
}

vector<vector<int> >& FileToMatrix(vector<vector<int> >& matrix, char* inputfile, int dimension, int order){

    matrix.resize(dimension, vector<int>(dimension));
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

    int flag = atoi(argv[3]);
    if (flag == 0){
        findOptimalconvThreshold();
    }

    if (flag == 1) {
        long double startTime;
        long double conventionaltime;
        long double normaltime;
        long double varianttime;

        if (argc != 4) {
            cout << "You must have 3 command-line arguments." << "\n";
            return -1;
        }

        int crossover = atoi(argv[1]);
        int dim = atoi(argv[2]);

        vector<vector<int> > first(dim, vector<int>(dim));
        vector<vector<int> > second(dim, vector<int>(dim));
        vector<vector<int> > result(dim, vector<int>(dim));

        MatrixOps help = *(new MatrixOps());
        help.make(first, second, dim);


        vector<vector<int> > test(dim, vector<int>(dim));
        help.makeidentity(test, test, dim);

        vector<vector<int> > matResult(dim, vector<int>(dim));
        matResult = multiply(test, test, result, crossover);

        help.printDiagonal(matResult, dim);

        //Have to fix timing stuff, not sure if it works
        /**
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
        //strassen(test, test, result, dim, crossover);
        varianttime = 1000 * (time(0) - startTime);


        cout << "Standard" << "\n";
        cout << "Size "
             << dim << " : " << conventionaltime << "\n";
        cout << "\n";

        cout << "Normal strassen" << "\n";
        cout << "Size "
             << dim << " : " << normaltime << "\n";
        cout << "\n";

        cout << "Variant strassen" << "\n";
        cout << "Crossover used: " << crossover << "\n";
        cout << "Size "
             << dim << " : " << varianttime << "\n";
        cout << "\n";**/

        return 0;
    }
}
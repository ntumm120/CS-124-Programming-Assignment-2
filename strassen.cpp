#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

#include "matrixops.h"

using namespace std;

// optimized for cache efficiency
void conventional(vector<vector<int> > &first, vector<vector<int> > &second,
                  vector<vector<int> > &result, int n) {
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            int r = first[i][k];
            for (int j = 0; j < n; j++) {
                result[i][j] += r * second[k][j];
            }
        }
    }
}

// To execute normal Strassen let crossover = 0, currently have crossover in
// command-line, has some bugs
void strassen(vector<vector<int> > &first, vector<vector<int> > &second,
              vector<vector<int> > &result, int n, int crossover) {
    MatrixOps helper = *(new MatrixOps());
    if (n <= crossover) {
        conventional(first, second, result, n);
    } else {
        int paddedsize = helper.powerround(n);
        int blocksize = paddedsize / 2;

        vector<int> pads(paddedsize);
        vector<vector<int> > padfirst(paddedsize, pads);
        vector<vector<int> > padsecond(paddedsize, pads);
        vector<vector<int> > padresult(paddedsize, pads);

        if (n != paddedsize) {
            helper.pad(padfirst, padsecond, first, second, n);
        } else {
            padfirst = first;
            padsecond = second;
        }

        vector<int> block(blocksize);

        vector<vector<int> > firstq1(blocksize, block),
            firstq2(blocksize, block);
        vector<vector<int> > firstq3(blocksize, block),
            firstq4(blocksize, block);
        vector<vector<int> > secondq1(blocksize, block),
            secondq2(blocksize, block);
        vector<vector<int> > secondq3(blocksize, block),
            secondq4(blocksize, block);
        vector<vector<int> > thirdq1(blocksize, block),
            thirdq2(blocksize, block);
        vector<vector<int> > thirdq3(blocksize, block),
            thirdq4(blocksize, block);

        vector<vector<int> > ablock(blocksize, block);
        vector<vector<int> > bblock(blocksize, block);

        vector<vector<int> > p1(blocksize, block);
        vector<vector<int> > p2(blocksize, block);
        vector<vector<int> > p3(blocksize, block);
        vector<vector<int> > p4(blocksize, block);
        vector<vector<int> > p5(blocksize, block);
        vector<vector<int> > p6(blocksize, block);
        vector<vector<int> > p7(blocksize, block);

        // partitioning matrix into its quadrants/blocks
        for (int i = 0; i < blocksize; i++) {
            for (int j = 0; j < blocksize; j++) {
                firstq1[i][j] = padfirst[i][j];
                firstq2[i][j] = padfirst[i][j + blocksize];
                firstq3[i][j] = padfirst[i + blocksize][j];
                firstq4[i][j] = padfirst[i + blocksize][j + blocksize];
                secondq1[i][j] = padsecond[i][j];
                secondq2[i][j] = padsecond[i][j + blocksize];
                secondq3[i][j] = padsecond[i + blocksize][j];
                secondq4[i][j] = padsecond[i + blocksize][j + blocksize];
            }
        }

        // constructing submatrices
        helper.add(firstq1, firstq4, ablock, blocksize);
        helper.add(secondq1, secondq4, bblock, blocksize);
        strassen(ablock, bblock, p1, blocksize, crossover);

        helper.add(firstq3, firstq4, ablock, blocksize);
        strassen(ablock, secondq1, p2, blocksize, crossover);

        helper.sub(secondq2, secondq4, bblock, blocksize);
        strassen(firstq1, bblock, p3, blocksize, crossover);

        helper.sub(secondq3, secondq1, bblock, blocksize);
        strassen(firstq4, bblock, p4, blocksize, crossover);

        helper.add(firstq1, firstq2, ablock, blocksize);
        strassen(ablock, secondq4, p5, blocksize, crossover);

        helper.sub(firstq3, firstq1, ablock, blocksize);
        helper.add(secondq1, secondq2, bblock, blocksize);
        strassen(ablock, bblock, p6, blocksize, crossover);

        helper.sub(firstq2, firstq4, ablock, blocksize);
        helper.add(secondq2, secondq4, bblock, blocksize);
        strassen(ablock, bblock, p7, blocksize, crossover);

        helper.add(p1, p4, ablock, blocksize);
        helper.sub(ablock, p5, bblock, blocksize);
        helper.add(bblock, p7, secondq1, blocksize);

        helper.add(p3, p5, thirdq2, blocksize);

        helper.add(p2, p4, thirdq3, blocksize);

        helper.sub(p1, p2, ablock, blocksize);
        helper.add(ablock, p3, bblock, blocksize);
        helper.add(bblock, p6, thirdq4, blocksize);

        // calculate result matrix
        for (int i = 0; i < blocksize; i++) {
            for (int j = 0; j < blocksize; j++) {
                padresult[i][j] = thirdq1[i][j];
                padresult[i][blocksize + j] = thirdq2[i][j];
                padresult[blocksize + i][j] = thirdq3[i][j];
                padresult[blocksize + i][blocksize + j] = thirdq4[i][j];
            }
        }

        // remove additional values from expanded matrix
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = padresult[i][j];
            }
        }
    }
}

// This one works, has efficient padding implementation
void newstrassen(vector<vector<int> > &a, vector<vector<int> > &b,
                 vector<vector<int> > &c, int dim, int crossover) {
    if (dim <= crossover) {
        conventional(a, b, c, dim);
    } else {
        MatrixOps helper = *(new MatrixOps());
        if (dim % 2 == 1) {
            dim = dim + 1;
            a.resize(dim);
            b.resize(dim);
            c.resize(dim);

            for (int i = 0; i < dim; i++) {
                a[i].resize(dim);
                b[i].resize(dim);
                c[i].resize(dim);
            }
        }
        int newdim = dim / 2;

        vector<int> quad(newdim);
        vector<vector<int> > a11(newdim, quad), a12(newdim, quad),
            a21(newdim, quad), a22(newdim, quad), b11(newdim, quad),
            b12(newdim, quad), b21(newdim, quad), b22(newdim, quad),
            c11(newdim, quad), c12(newdim, quad), c21(newdim, quad),
            c22(newdim, quad);

        helper.divide(a, a11, 0, 0, newdim);
        helper.divide(a, a12, 0, newdim, newdim);
        helper.divide(a, a21, newdim, 0, newdim);
        helper.divide(a, a22, newdim, newdim, newdim);
        helper.divide(b, b11, 0, 0, newdim);
        helper.divide(b, b12, 0, newdim, newdim);
        helper.divide(b, b21, newdim, 0, newdim);
        helper.divide(b, b22, newdim, newdim, newdim);

        vector<vector<int> > first(newdim, quad), second(newdim, quad),
            P1(newdim, quad), P2(newdim, quad), P3(newdim, quad),
            P4(newdim, quad), P5(newdim, quad), P6(newdim, quad),
            P7(newdim, quad);

        helper.add(a11, a22, first, newdim);
        helper.add(b11, b22, second, newdim);

        newstrassen(first, second, P1, newdim, crossover);
        helper.add(a21, a22, first, newdim);

        newstrassen(first, b11, P2, newdim, crossover);
        helper.sub(b12, b22, second, newdim);

        newstrassen(a11, second, P3, newdim, crossover);
        helper.sub(b21, b11, second, newdim);

        newstrassen(a22, second, P4, newdim, crossover);
        helper.add(a11, a12, first, newdim);
        newstrassen(first, b22, P5, newdim, crossover);

        helper.sub(a21, a11, first, newdim);
        helper.add(b11, b12, second, newdim);
        newstrassen(first, second, P6, newdim, crossover);

        helper.sub(a12, a22, first, newdim);
        helper.add(b21, b22, second, newdim);
        newstrassen(first, second, P7, newdim, crossover);

        helper.add(P1, P4, first, newdim);
        helper.add(first, P7, second, newdim);
        helper.sub(second, P5, c11, newdim);

        helper.add(P3, P5, c12, newdim);
        helper.add(P2, P4, c21, newdim);
        helper.sub(P1, P2, first, newdim);
        helper.add(P3, P6, second, newdim);
        helper.add(first, second, c22, newdim);

        helper.combine(c11, c, 0, 0, newdim);
        helper.combine(c12, c, 0, newdim, newdim);
        helper.combine(c21, c, newdim, 0, newdim);
        helper.combine(c22, c, newdim, newdim, newdim);
    }
}

void convertFile(vector<vector<int> > &first, vector<vector<int> > &second,
                 char *inputfile, int dim) {
    first.resize(dim, vector<int>(dim));
    second.resize(dim, vector<int>(dim));

    ifstream inFile(inputfile);
    string line;

    for (int i = 0; i < (2 * dim * dim); i++) {
        getline(inFile, line);
        if (i < (dim * dim)) {
            first[i / dim][i % dim] = stoi(line);
        } else {
            second[(i - (dim * dim)) / dim][i % dim] = stoi(line);
        }
    }

    inFile.close();
}

void findcrossover() {
    for (int crossover = 128; crossover <= 256; crossover *= 2) {
        double total = 0;
        for (int j = 0; j < 3; j++) {
            vector<vector<int> > first(800, vector<int>(800));
            vector<vector<int> > second(800, vector<int>(800));
            vector<vector<int> > third(800, vector<int>(800));
            MatrixOps help = *(new MatrixOps());
            help.makeidentity(first, second, 800);

            clock_t start;
            start = clock();

            newstrassen(first, second, third, 800, crossover);
            clock_t end = clock();
            total += (end - start) / (double)(CLOCKS_PER_SEC);
        }
        cout << crossover << "\t" << total / 3 << endl;
    }
}

float triangles(float p, int dim, int crossover) {
    vector<vector<int> > graph(dim, vector<int>(dim));
    vector<vector<int> > graph_squared(dim, vector<int>(dim));
    vector<vector<int> > graph_cubed(dim, vector<int>(dim));

    // Set up random number generator
    default_random_engine generator(time(0));
    std::bernoulli_distribution distribution(p);

    for (int i = 0; i < dim; i++) {
        for (int j = i; j < dim; j++) {
            bool r = distribution(generator);
            if (r) {
                graph[i][j] = 1;
                graph[j][i] = 1;
            } else {
                graph[i][j] = 0;
                graph[j][i] = 0;
            }
        }
    }

    newstrassen(graph, graph, graph_squared, dim, crossover);
    newstrassen(graph, graph_squared, graph_cubed, dim, crossover);
    float total = 0;
    for (int i = 0; i < dim; i++) {
        total += graph_cubed[i][i];
    }
    return total / 6;
}

int main(int argc, char *argv[]) {
    // Wrote function to read in a file of numbers, for now we will use
    // randomly generated numbers to test

    // Still need to figure out experimental crossover point, made function to
    // do it

    if (argc != 5) {
        cout << "You must have 5 command-line arguments."
             << "\n";
        return -1;
    }

    int flag = atoi(argv[1]);
    char *infile = argv[4];
    int crossover = atoi(argv[3]);
    int dim = atoi(argv[2]);

    // Have to read first d^2 numbers into one matrix and the rest into the
    // other
    if (flag == 0) {
        MatrixOps help = *(new MatrixOps());
        vector<vector<int> > first(dim, vector<int>(dim));
        vector<vector<int> > second(dim, vector<int>(dim));
        vector<vector<int> > result(dim, vector<int>(dim));

        convertFile(first, second, infile, dim);
        newstrassen(first, second, result, dim, crossover);
        help.printDiagonal(result, dim);

    } else if (flag == 1) {
        findcrossover();
    } else if (flag == 2) {
        long double startTime;
        long double conventionaltime;
        long double normaltime;
        long double varianttime;

        vector<vector<int> > first(dim, vector<int>(dim));
        vector<vector<int> > second(dim, vector<int>(dim));
        vector<vector<int> > result(dim, vector<int>(dim));

        MatrixOps help = *(new MatrixOps());
        help.make(first, second, dim);

        vector<vector<int> > test(dim, vector<int>(dim));
        help.makeidentity(test, test, dim);

        newstrassen(test, test, result, dim, crossover);

        help.printDiagonal(result, dim);

        /**
        //conventional fully works on test matrices
        startTime = time(0);
        //conventional(test,test,result, dim);
        conventionaltime = 1000 * (time(0) - startTime);

        //No crossover strassen (DOESNT WORK NEED TO FIGURE OUT WHY)
        startTime = time(0);
        //strassen(first, second, result, dim, 0);
        normaltime = 1000 * (time(0) - startTime);

        //crossover strassen doesn't work either
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
        cout << "crossover used: " << crossover << "\n";
        cout << "Size "
             << dim << " : " << varianttime << "\n";
        cout << "\n";**/

        return 0;
    } else if (flag == 4) {
        float p = atof(argv[4]);
        float t = triangles(p, dim, crossover);
        printf("%f\n", t);
        return 0;
    }
}
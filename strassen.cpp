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
        graph[i][i] = 0;
        for (int j = i + 1; j < dim; j++) {
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
        return 0;
    } else if (flag == 4) {
        float p = atof(argv[4]);
        float total = 0;
        for (int i = 0; i < 5; i++) {
            total += triangles(p, dim, crossover);
        }
        float average = total / 5;
        printf("%f\n", average);
        return 0;
    }
}
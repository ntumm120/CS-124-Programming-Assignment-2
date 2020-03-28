#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <random>
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
    void combine (vector< vector<int> > &A,
                  vector< vector<int> > &B,
                  int row, int col, int d);
    void divide(vector< vector<int> > &A, vector<vector<int> > &B,int row, int col, int d);
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

void MatrixOps::combine (vector< vector<int> > &first,
                         vector<vector<int> > &second, int row, int col, int dim)
{
    for (int i = 0, r = row; i < dim; i++, r++)
    {
        for (int j = 0, c = col; j < dim; j++, c++)
        {
            second[r][c] = first[i][j];
        }
    }
}

void MatrixOps::divide(vector< vector<int> > &first,
                        vector< vector<int> > &second, int row, int col, int dim)
{
    for (int i = 0, r = row; i < dim; i++, r++)
    {
        for (int j = 0, c = col; j < dim; j++, c++)
        {
            second[i][j] = first[r][c];
        }
    }
}

//optimized for cache efficiency
void conventional(vector<vector<int> > &first,
                         vector<vector<int> > &second, vector<vector<int> > &result, int n){
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            int r = first[i][k];
            for (int j = 0; j < n; j++) {
                result[i][j] += r * second[k][j];
            }
        }
    }
}

//This works, has efficient padding implementation
void newstrassen(vector< vector<int> > &firstmat, vector< vector<int> > &secondmat,
                   vector< vector<int> > &resultmat, int dim, int crossover){

    if (dim <= crossover)
    {
        conventional(firstmat, secondmat, resultmat, dim);
    }
    else
    {
        MatrixOps helper = *(new MatrixOps());
        if (dim % 2 == 1){
            dim = dim + 1;
            firstmat.resize( dim ,vector<int>(dim ,0));
            secondmat.resize( dim ,vector<int>(dim ,0));
            resultmat.resize( dim ,vector<int>(dim ,0));
        }
        int newdim = round(dim/2);

        vector<int> quad(newdim);
        vector< vector<int> > a11(newdim, quad), a12(newdim, quad), a21(newdim, quad), 
        a22(newdim, quad), b11(newdim, quad), b12(newdim, quad), b21(newdim, quad), b22(newdim, quad),
        c11(newdim, quad), c12(newdim, quad), c21(newdim, quad), c22(newdim, quad), first(newdim, quad),
        second(newdim, quad), P1(newdim, quad), P2(newdim, quad), P3(newdim, quad), P4(newdim, quad),
        P5(newdim, quad), P6(newdim, quad), P7(newdim, quad);
        
        helper.divide(firstmat, a11, 0 , 0, newdim);
        helper.divide(firstmat, a12, 0 , newdim, newdim);
        helper.divide(firstmat, a21, newdim, 0, newdim);
        helper.divide(firstmat, a22, newdim, newdim, newdim);
        helper.divide(secondmat, b11, 0 , 0, newdim);
        helper.divide(secondmat, b12, 0 , newdim, newdim);
        helper.divide(secondmat, b21, newdim, 0, newdim);
        helper.divide(secondmat, b22, newdim, newdim, newdim);

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

        helper.combine(c11, resultmat, 0 , 0, newdim);
        helper.combine(c12, resultmat, 0 , newdim, newdim);
        helper.combine(c21, resultmat, newdim, 0, newdim);
        helper.combine(c22, resultmat, newdim, newdim, newdim);
    }
}


void convertFile(vector<vector<int> >& first, vector<vector<int> >& second, char* inputfile, int dim){

    first.resize(dim, vector<int>(dim));
    second.resize(dim, vector<int>(dim));

    ifstream inFile(inputfile);
    string line;

    for (int i = 0; i < (2*dim*dim); i++) {
        getline(inFile, line);
        if (i < (dim * dim)) {
            first[i / dim][i % dim] = stoi(line);
        }
        else{
            second[(i - (dim * dim))/dim][i % dim] = stoi(line);
        }
    }

    inFile.close();
}


void findcrossover() {

    for (int crossover = 128; crossover <= 256; crossover+=2){
        double total = 0;
        for (int j = 0; j < 1; j++){
            vector<vector<int> > first(10000, vector<int>(10000));
            vector<vector<int> > second(10000, vector<int>(10000));
            vector<vector<int> > third(10000, vector<int>(10000));
            MatrixOps help = *(new MatrixOps());
            help.makeidentity(first, second, 10000);

            clock_t start;
            start = clock();

            newstrassen(first, second, third, 10000, crossover);
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
    bernoulli_distribution distribution(p);

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

    if (argc != 4) {
        cout << "You must have 4 command-line arguments."
             << "\n";
        return -1;
    }

    int flag = atoi(argv[1]);
    char *infile = argv[3];
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
    }
    else if (flag == 1) {
        findcrossover();
    }
    else if (flag == 2) {
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
    }
    else if (flag == 3) {
        float p1 = .05;
        float t1 = triangles(p1, dim, crossover);
        printf("%f\n", t1);
        return 0;
    }
}

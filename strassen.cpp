/**
*Programming Assignemnt 2, Source Code
*Chris Zhu, Rakesh Nori
**/
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <climits>
#include <time.h>
#include <chrono>
#include <thread>
#include <pthread.h>

using namespace std;

struct stra_mtx {
    //start indices inclusive, end indices exclusive
    int rowStart;
    int rowEnd;
    int colStart;
    int colEnd;

    stra_mtx(int a, int b, int c, int d){
        rowStart = a;
        rowEnd = b;
        colStart = c;
        colEnd = d;
    }

    int r_split(){
        return rowStart + (rowEnd - rowStart) / 2;
    }

    int c_split(){
        return colStart + (colEnd - colStart) / 2;
    }

    int size(){
        return rowEnd - rowStart;
    }

};

void genRandMats(int low, int high, int n, string fileName){

    ofstream file(fileName.c_str());
    for (int i = 0; i < n * n * 2; ++i){

        int curr = (rand() % (high - low + 1)) + low;
        file << curr << endl;
    }
    file.close(); 
}

// matrix_multiply(c, a, b)
//    NOTE: assumes `a`, `b`, and `c` are square matrices with the same dimensions.
//    Computes the matrix product `a * b` and stores it in `c` through the conventional method.
//    Inspired by Eddie Kohler's Network 8 Lecture from CS 61
vector<vector<int>> matrix_multiply(vector<vector<int>> &a, vector<vector<int>> &b) {

    int n = a.size();

    vector<vector<int>> c(n , vector<int> (n, 0)); 

    for (size_t i = 0; i < n; ++i){
        vector<int> currRowA = a[i];
        for (size_t k = 0; k < n; ++k){
            vector<int> currRowB  = b[k];
            double currElemA = currRowA[k];
            for (size_t j = 0; j < n; ++j){
                c[i][j] += currElemA * currRowB[j];   
            }
        }
    }

    return c;
}

vector<vector<int>> top_left(vector<vector<int>> m){
    vector<vector<int>> out;
    for(int row = 0; row < m.size() / 2; row++){
        vector<int> temp;
        for(int col = 0; col < m.size() / 2; col++){
            temp.push_back(m[row][col]);
        }
        out.push_back(temp);
    }
    return out;
}

vector<vector<int>> bottom_left(vector<vector<int>> m){
    vector<vector<int>> out;
    for(int row = m.size() / 2; row < m.size(); row++){
        vector<int> temp;
        for(int col = 0; col < m.size() / 2; col++){
            temp.push_back(m[row][col]);
        }
        out.push_back(temp);
    }
    return out;
}

vector<vector<int>> top_right(vector<vector<int>> m){
    vector<vector<int>> out;
    for(int row = 0; row < m.size() / 2; row++){
        vector<int> temp;
        for(int col = m.size() / 2; col < m.size(); col++){
            temp.push_back(m[row][col]);
        }
        out.push_back(temp);
    }
    return out;
}

vector<vector<int>> bottom_right(vector<vector<int>> m){
    vector<vector<int>> out;
    for(int row = m.size() / 2; row < m.size(); row++){
        vector<int> temp;
        for(int col = m.size() / 2; col < m.size(); col++){
            temp.push_back(m[row][col]);
        }
        out.push_back(temp);
    }
    return out;
}

void print_mtx(vector<vector<int>> m){
    for(vector<int> row : m){
        for(int col : row){
            cout << col << " ";
        }
        cout << endl;
    }
}

vector<vector<int>> add(vector<vector<int>> m, vector<vector<int>> n){
    vector<vector<int>> out;
    for(int row = 0; row < m.size(); row++){
        vector<int> temp;
        for(int col = 0; col < m.size(); col++){
            temp.push_back(m[row][col] + n[row][col]);
        }
        out.push_back(temp);
    }
    return out;
}

vector<vector<int>> sub(vector<vector<int>> m, vector<vector<int>> n){
    vector<vector<int>> out;
    for(int row = 0; row < m.size(); row++){
        vector<int> temp;
        for(int col = 0; col < m.size(); col++){
            temp.push_back(m[row][col] - n[row][col]);
        }
        out.push_back(temp);
    }
    return out;
}

vector<vector<int>> join(vector<vector<int>> topLeft, vector<vector<int>> topRight,
    vector<vector<int>> bottomLeft, vector<vector<int>> bottomRight){
        vector<vector<int>> out;
        for(int row = 0; row < topLeft.size(); row++){
            vector<int> tempRow(topLeft.size() * 2, 0);
            for(int col = 0; col < topLeft.size(); col++){
                tempRow[col] = topLeft[row][col];
                tempRow[col + topLeft.size()] = topRight[row][col];
            }
            out.push_back(tempRow);
        }
        for(int row = 0; row < topLeft.size(); row++){
            vector<int> tempRow(topLeft.size() * 2, 0);
            for (int col = 0; col < bottomLeft.size(); col++){
                tempRow[col] = bottomLeft[row][col];
                tempRow[col + topLeft.size()] = bottomRight[row][col];
            }
            out.push_back(tempRow);
        }
        return out;
    }

vector<vector<int>> pad(vector<vector<int>> m){
    int extraZeros = ((int) pow(2, ceil(log2(m.size())))) - m.size();
    vector<vector<int>> out;
    for(int row = 0; row < m.size(); row++){
        vector<int> tempRow;
        for(int i = 0; i < m.size(); i++){
            tempRow.push_back(m[row][i]);
        }
        for(int i = 0; i < extraZeros; i++){
            tempRow.push_back(0);
        }
        out.push_back(tempRow);
    }
    for(int i = 0; i < extraZeros; i++){
        out.push_back(vector<int> (m.size() + extraZeros, 0));
    }
    return out;
}

vector<vector<int>> strassen_multiply(vector<vector<int>> m, vector<vector<int>> n, int thresh){
    if(m.size() < thresh){
        
        return matrix_multiply(m, n);
        //conventional
    }
    else{
        vector<vector<int>> a = top_left(m);
        vector<vector<int>> b = top_right(m);
        vector<vector<int>> c = bottom_left(m);
        vector<vector<int>> d = bottom_right(m);
        vector<vector<int>> e = top_left(n);
        vector<vector<int>> f = top_right(n);
        vector<vector<int>> g = bottom_left(n);
        vector<vector<int>> h = bottom_right(n);

       // print_mtx(add(a, b));

        vector<vector<int>> p1 = strassen_multiply(a, sub(f, h), thresh);
        vector<vector<int>> p2 = strassen_multiply(add(a, b), h, thresh);
        vector<vector<int>> p3 = strassen_multiply(add(c, d), e, thresh);
        vector<vector<int>> p4 = strassen_multiply(d, sub(g, e), thresh);
        vector<vector<int>> p5 = strassen_multiply(add(a, d), add(e, h), thresh);
        vector<vector<int>> p6 = strassen_multiply(sub(b, d), add(g, h), thresh);
        vector<vector<int>> p7 = strassen_multiply(sub(a, c), add(e, f), thresh);
        
        vector<vector<int>> topLeft = add(sub(add(p5, p4), p2),p6);
        vector<vector<int>> topRight = add(p1, p2);
        vector<vector<int>> bottomLeft = add(p3, p4);
        vector<vector<int>> bottomRight = sub(add(p5, p1), add(p3, p7));

        return join(topLeft, topRight, bottomLeft, bottomRight);

        /**
        stra_mtx a = stra_mtx(m_ind.rowStart, m_ind.r_split(), m_ind.colStart, m_ind.c_split());
        stra_mtx b = stra_mtx(m_ind.rowStart, m_ind.r_split(), m_ind.c_split(), m_ind.colEnd);
        stra_mtx c = stra_mtx(m_ind.r_split(), m_ind.rowEnd, m_ind.colStart, m_ind.c_split());
        stra_mtx d = stra_mtx(m_ind.r_split(), m_ind.rowEnd, m_ind.c_split(), m_ind.colEnd);
        stra_mtx e = stra_mtx(n_ind.rowStart, n_ind.r_split(), n_ind.colStart, n_ind.c_split());
        stra_mtx f = stra_mtx(n_ind.rowStart, n_ind.r_split(), n_ind.c_split(), n_ind.colEnd);
        stra_mtx g = stra_mtx(n_ind.r_split(), n_ind.rowEnd, n_ind.colStart, n_ind.c_split());
        stra_mtx h = stra_mtx(n_ind.r_split(), n_ind.rowEnd, n_ind.c_split(), n_ind.colEnd);
        **/
    }

}


int main(int argc, char *argv[]){	
    srand((unsigned)time(NULL));
    int flag = stoi(argv[1]);
    int dim = stoi(argv[2]);
    string fileName = argv[3];

    genRandMats(0, 2, dim, fileName.c_str());

    fstream file; 
    file.open(fileName.c_str()); 

    vector<vector<int>> m(dim , vector<int> (dim));
    vector<vector<int>> n(dim , vector<int> (dim));
    int curr;
    
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            file >> curr;
            m[i][j] = curr;
        }
    }

    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            file >> curr;
            n[i][j] = curr;

        }
    }
    
    m = pad(m);
    n = pad(n);
    
    cout << "Matrix 1: " << endl;
    //print_mtx(m);
    cout << "Matrix 2: " << endl; 
    //print_mtx(n);

    int thresh = 5;

    /**
    while(thresh <= 1200){
        auto start = chrono::high_resolution_clock::now(); 
        print_mtx(strassen_multiply(m, n, thresh));
        auto stop = chrono::high_resolution_clock::now(); 
	    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start); 
        cout << "DURATION for strassen: " << duration.count() << " ms" << endl; 
        thresh += 100;
    }

    **/

    auto start = chrono::high_resolution_clock::now(); 
    strassen_multiply(m, n, thresh);
    //print_mtx(strassen_multiply(m, n, thresh));
    auto stop = chrono::high_resolution_clock::now(); 
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start); 
    cout << "DURATION for strassen: " << duration.count() << " ms" << endl;

}
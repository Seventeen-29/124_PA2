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
vector<vector<int> > matrix_multiply(vector<vector<int> > const &a, vector<vector<int> > const &b) {

    int n = a.size();

    vector<vector<int> > c(n , vector<int> (n, 0)); 

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

vector<vector<int> > top_left(vector<vector<int> > const &m){
    vector<vector<int> > out;
    for(int row = 0; row < m.size() / 2; row++){
        vector<int> temp;
        for(int col = 0; col < m.size() / 2; col++){
            temp.push_back(m[row][col]);
        }
        out.push_back(temp);
    }
    return out;
}

vector<vector<int> > bottom_left(vector<vector<int> > const &m){
    vector<vector<int> > out;
    for(int row = m.size() / 2; row < m.size(); row++){
        vector<int> temp;
        for(int col = 0; col < m.size() / 2; col++){
            temp.push_back(m[row][col]);
        }
        out.push_back(temp);
    }
    return out;
}

vector<vector<int> > top_right(vector<vector<int> > const &m){
    vector<vector<int> > out;
    for(int row = 0; row < m.size() / 2; row++){
        vector<int> temp;
        for(int col = m.size() / 2; col < m.size(); col++){
            temp.push_back(m[row][col]);
        }
        out.push_back(temp);
    }
    return out;
}

vector<vector<int> > bottom_right(vector<vector<int> > const &m){
    vector<vector<int> > out;
    for(int row = m.size() / 2; row < m.size(); row++){
        vector<int> temp;
        for(int col = m.size() / 2; col < m.size(); col++){
            temp.push_back(m[row][col]);
        }
        out.push_back(temp);
    }
    return out;
}

void print_mtx(vector<vector<int> > const &m){
    for(vector<int> row : m){
        for(int col : row){
            cout << col << " ";
        }
        cout << endl;
    }
}

void print_diag(vector<vector<int> > const &m, int dim){
    for(int i = 0; i < dim; ++i){
        cout << m[i][i] << endl;
    }
}

vector<vector<int> > add(vector<vector<int> > const &m, vector<vector<int> > const &n){
    vector<vector<int> > out;
    for(int row = 0; row < m.size(); row++){
        vector<int> temp;
        for(int col = 0; col < m.size(); col++){
            temp.push_back(m[row][col] + n[row][col]);
        }
        out.push_back(temp);
    }
    return out;
}

vector<vector<int> > sub(vector<vector<int> > const &m, vector<vector<int> > const &n){
    vector<vector<int> > out;
    for(int row = 0; row < m.size(); row++){
        vector<int> temp;
        for(int col = 0; col < m.size(); col++){
            temp.push_back(m[row][col] - n[row][col]);
        }
        out.push_back(temp);
    }
    return out;
}

vector<vector<int> > join(vector<vector<int> > const &topLeft, vector<vector<int> > const &topRight,
    vector<vector<int> > const &bottomLeft, vector<vector<int> > const &bottomRight){
        vector<vector<int> > out;
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

void pad(vector<vector<int> > &m){
    vector<int> newRow;
    for (int i = 0; i < m.size(); i++){
        m[i].push_back(0);
        newRow.push_back(0);
    }
    newRow.push_back(0);
    m.push_back(newRow);
}

vector<vector<int> > strassen_multiply(vector<vector<int> > &m, vector<vector<int> > &n, int thresh){
    if(m.size() < thresh){
        
        return matrix_multiply(m, n);
        //conventional
    }
    else{
        if (m.size() % 2 != 0)
            pad(m);
        if (n.size() % 2 != 0)
            pad(n);
        
        vector<vector<int> > a = top_left(m);
        vector<vector<int> > b = top_right(m);
        vector<vector<int> > c = bottom_left(m);
        vector<vector<int> > d = bottom_right(m);
        vector<vector<int> > e = top_left(n);
        vector<vector<int> > f = top_right(n);
        vector<vector<int> > g = bottom_left(n);
        vector<vector<int> > h = bottom_right(n);

       // print_mtx(add(a, b));
       vector<vector<int> > sub_fh = sub(f, h);
       vector<vector<int> > add_ab = add(a, b);
       vector<vector<int> > add_cd = add(c, d);
       vector<vector<int> > sub_ge = sub(g, e);
       vector<vector<int> > add_ad = add(a, d);
       vector<vector<int> > add_eh = add(e, h);
       vector<vector<int> > sub_bd = sub(b, d);
       vector<vector<int> > add_gh = add(g, h);
       vector<vector<int> > sub_ac = sub(a, c);
       vector<vector<int> > add_ef = add(e, f);



        vector<vector<int> > p1 = strassen_multiply(a, sub_fh, thresh);
        vector<vector<int> > p2 = strassen_multiply(add_ab, h, thresh);
        vector<vector<int> > p3 = strassen_multiply(add_cd, e, thresh);
        vector<vector<int> > p4 = strassen_multiply(d, sub_ge, thresh);
        vector<vector<int> > p5 = strassen_multiply(add_ad, add_eh, thresh);
        vector<vector<int> > p6 = strassen_multiply(sub_bd, add_gh, thresh);
        vector<vector<int> > p7 = strassen_multiply(sub_ac, add_ef, thresh);
        
        vector<vector<int> > topLeft = add(sub(add(p5, p4), p2),p6);
        vector<vector<int> > topRight = add(p1, p2);
        vector<vector<int> > bottomLeft = add(p3, p4);
        vector<vector<int> > bottomRight = sub(add(p5, p1), add(p3, p7));

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
    srand48((unsigned)time(NULL));
    int flag = stoi(argv[1]);
    int dim = stoi(argv[2]);
    string fileName = argv[3];

    //cout << "FLAG " << flag << endl;

    // option to generate random matrices on the fly
    if (flag == 1)
        genRandMats(0, 2, dim, fileName.c_str()); // generates random ints from 0 to 2 (inclusive)
    

    // task 3: vertices
    if (flag == 3){
        int vertices = 1024;
        cout << "type value of p: " << endl;
        double p;
        cin >> p;
        vector<vector<int> > A(vertices, vector<int> (vertices, 0));
        for (int i = 0; i < vertices; i++)
        {
            for (int j = 0; j < i; j++)
            {
                if (drand48() < p){
                    A[i][j] = 1;
                    A[j][i] = 1;
                }
            }
        }
        //print_mtx(A);
        vector<vector<int> > a_squared = strassen_multiply(A, A, 512);
        vector<vector<int> > res = strassen_multiply(a_squared, A, 512);
        int sum = 0;
        for (int i = 0; i < res.size(); i++){
            sum += res[i][i];
        }
        //print_diag(strassen_multiply(strassen_multiply(A, A, 512), A, 512));
        cout << "NUM TRIANGLES: " << (sum / 6) << endl;
        return 0;
    }
    fstream file; 
    file.open(fileName.c_str()); 

    vector<vector<int> > m(dim , vector<int> (dim));
    vector<vector<int> > n(dim , vector<int> (dim));
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
    
    //m = pad(m);
    //n = pad(n);
    
    //cout << "Matrix 1: " << endl;
    //print_mtx(m);
    //cout << "Matrix 2: " << endl; 
    //print_mtx(n);

    // for testing optimal n_0
    if (flag == 1){
        vector<int> thresh;
        //thresh.push_back(4);
        //thresh.push_back(8);
        thresh.push_back(250);
        thresh.push_back(260);
        thresh.push_back(270);
        thresh.push_back(280);
        thresh.push_back(290);
        thresh.push_back(300);
        //thresh = { 256, 512, 768, 896, 1024, 1152, 1280, 1536, 2048 };

        for (int elem : thresh){
            cout << elem << endl;
            genRandMats(0, 2, elem + 1, fileName.c_str());
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
            int sum = 0;
            int sum2 = 0;
            int trials = 10;
            vector<chrono::milliseconds> results();
            for (int i = 0; i < trials; i++){
                 auto start = chrono::high_resolution_clock::now(); 
                strassen_multiply(m, n, elem);
                //print_mtx(strassen_multiply(m, n, thresh));
                auto stop = chrono::high_resolution_clock::now(); 
                auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start); 
                sum += duration.count();

                auto start2 = chrono::high_resolution_clock::now(); 
                matrix_multiply(m, n);
                auto stop2 = chrono::high_resolution_clock::now(); 
                auto duration2 = chrono::duration_cast<chrono::milliseconds>(stop2 - start2); 
                sum2 += duration2.count();
            }
            sum /= trials;
            sum2 /= trials;

            cout << elem << " AVG DURATION for strassen: " << sum << " ms" << endl;

            cout << elem << " AVG DURATION for conventional: " << sum2 << " ms" << endl;
        }
    }
    else if (flag == 0){
        print_diag(strassen_multiply(m, n, 512), dim);
    }

/**
    auto start = chrono::high_resolution_clock::now(); 
    strassen_multiply(m, n, thresh0);
    //print_mtx(strassen_multiply(m, n, thresh));
    auto stop = chrono::high_resolution_clock::now(); 
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start); 
    cout << "DURATION for strassen: " << duration.count() << " ms" << endl;

    start = chrono::high_resolution_clock::now(); 
    strassen_multiply(m, n, thresh1);
    //print_mtx(strassen_multiply(m, n, thresh));
    stop = chrono::high_resolution_clock::now(); 
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start); 
    cout << "DURATION for strassen: " << duration.count() << " ms" << endl;

    start = chrono::high_resolution_clock::now(); 
    strassen_multiply(m, n, thresh2);
    //print_mtx(strassen_multiply(m, n, thresh));
    stop = chrono::high_resolution_clock::now(); 
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start); 
    cout << "DURATION for strassen: " << duration.count() << " ms" << endl;

    start = chrono::high_resolution_clock::now(); 
    strassen_multiply(m, n, thresh3);
    //print_mtx(strassen_multiply(m, n, thresh));
    stop = chrono::high_resolution_clock::now(); 
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start); 
    cout << "DURATION for strassen: " << duration.count() << " ms" << endl;

    **/

}
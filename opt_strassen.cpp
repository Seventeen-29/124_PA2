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

    stra_mtx(int rS, int rE, int cS, int cE){
        rowStart = rS;
        rowEnd = rE;
        colStart = cS;
        colEnd = cE;
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

vector<vector<int> > top_left(vector<vector<int> > &m){
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

vector<vector<int> > bottom_left(vector<vector<int> > &m){
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

vector<vector<int> > top_right(vector<vector<int> > &m){
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

vector<vector<int> > bottom_right(vector<vector<int> > &m){
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
    for(int i = 0; i < m.size(); i++){
        m[i].push_back(0);
    }
    m.push_back(vector<int> (m.size() + 1, 0));
}

void trim(vector<vector<int> > &m){
    m.pop_back();
    for(int i = 0; i < m.size(); i++){
        m[i].pop_back();
    }
}

//performs m [m_ind] - n[n_ind] = res, all are square matrices of same size
vector<vector<int> > sub_ind(vector<vector<int> > &m, vector<vector<int> > &n, stra_mtx m_ind, stra_mtx n_ind){
    vector<vector<int> > res(m_ind.size(), vector<int> (m_ind.size(), 0));
    for(int row = 0; row < m_ind.size(); row++){
        for(int col = 0; col < m_ind.size(); col++){
            res[row][col] = m[m_ind.rowStart + row][m_ind.colStart + col] -
                            n[n_ind.rowStart + row][n_ind.colStart + col];
        }
    }
    return res;
}

//performs m [m_ind] + n[n_ind] = res, all are square matrices of same size
vector<vector<int> > add_ind(vector<vector<int> > &m, vector<vector<int> > &n, stra_mtx m_ind, stra_mtx n_ind){
    vector<vector<int> > res(m_ind.size(), vector<int> (m_ind.size(), 0));
    for(int row = 0; row < m_ind.size(); row++){
        for(int col = 0; col < m_ind.size(); col++){
            res[row][col] = m[m_ind.rowStart + row][m_ind.colStart + col] +
                            n[n_ind.rowStart + row][n_ind.colStart + col];
        }
    }
    return res;
}

//adds n to m and stores result back in m
void add_inplace(vector<vector<int> > &m, vector<vector<int> > &n){
    for(int row = 0; row < m.size(); row++){
        for(int col = 0; col < m.size(); col++){
            m[row][col] += n[row][col];
        }
    }
}

//subtracts n from m and stores result back in m
void sub_inplace(vector<vector<int> > &m, vector<vector<int> > &n){
    for(int row = 0; row < m.size(); row++){
        for(int col = 0; col < m.size(); col++){
            m[row][col] -= n[row][col];
        }
    }
}

vector<vector<int> > strassen_multiply(vector<vector<int> > &m, vector<vector<int> > &n, int thresh){
    if(m.size() < thresh){
        return matrix_multiply(m, n);
    }
    else{
        bool isOdd = false;
        if (m.size() % 2 == 1){
            pad(m);
            pad(n);
            isOdd = true;
        }

        int mid = m.size() / 2;
        int end = m.size();

        stra_mtx tL = stra_mtx(0, mid, 0, mid);
        stra_mtx tR = stra_mtx(0, mid, mid, end);
        stra_mtx bL = stra_mtx(mid, end, 0, mid);
        stra_mtx bR = stra_mtx(mid, end, mid, end);

        vector<vector<int> > sub_fh = sub_ind(n, n, tR, bR);
        vector<vector<int> > add_ab = add_ind(m, m, tL, tR);
        vector<vector<int> > add_cd = add_ind(m, m, bL, bR);
        vector<vector<int> > sub_ge = sub_ind(n, n, bL, tL);
        vector<vector<int> > add_ad = add_ind(m, m, tL, bR);
        vector<vector<int> > add_eh = add_ind(n, n, tL, bR);
        vector<vector<int> > sub_bd = sub_ind(m, m, tR, bR);
        vector<vector<int> > add_gh = add_ind(n, n, bL, bR);
        vector<vector<int> > sub_ac = sub_ind(m, m, tL, bL);
        vector<vector<int> > add_ef = add_ind(n, n, tL, tR);

        vector<vector<int> > a = top_left(m);
        vector<vector<int> > h = bottom_right(n);
        vector<vector<int> > e = top_left(n);
        vector<vector<int> > d = bottom_right(m);

        vector<vector<int> > p1 = strassen_multiply(a, sub_fh, thresh);
        vector<vector<int> > p2 = strassen_multiply(add_ab, h, thresh);
        vector<vector<int> > p3 = strassen_multiply(add_cd, e, thresh);
        vector<vector<int> > p4 = strassen_multiply(d, sub_ge, thresh);
        vector<vector<int> > p5 = strassen_multiply(add_ad, add_eh, thresh);
        vector<vector<int> > p6 = strassen_multiply(sub_bd, add_gh, thresh);
        vector<vector<int> > p7 = strassen_multiply(sub_ac, add_ef, thresh);
        
        //topLeft is p6
        add_inplace(p6, p5);
        add_inplace(p6, p4);
        sub_inplace(p6, p2);
        //bottomRight is p5
        add_inplace(p5, p1);
        sub_inplace(p5, p3);
        sub_inplace(p5, p7);
        //topRight is p1
        add_inplace(p1, p2);
        //bottomLeft is p3
        add_inplace(p3, p4);

        vector<vector<int> > result = join(p6, p1, p3, p5);

        if(isOdd){
            trim(result);
        }
        return result;
    }

}

int main(int argc, char *argv[]){	
    srand((unsigned)time(NULL));
    srand48((unsigned)time(NULL));
    int flag = stoi(argv[1]);
    int dim = stoi(argv[2]);
    string fileName = argv[3];

    if (flag == 4){
        vector<vector<int> > mtx1 = {{2, 3, 1}, {7, 1, 2}, {8, 2, 4}};
        vector<vector<int> > mtx2 = {{1, 5, 1}, {6, 9, 5}, {8, 0, 6}};
        print_mtx(strassen_multiply(mtx1, mtx2, 2));
        return 69;
    }

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
 
    if (flag == 1){
        vector<int> thresh;
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

}
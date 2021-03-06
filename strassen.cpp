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
        for (size_t k = 0; k < n; ++k){
            for (size_t j = 0; j < n; ++j){
                c[i][j] += a[i][k] * b[k][j]; 
            }
        }
    }
    return c;
}

void print_diag(vector<vector<int> > const &m){
    for(int i = 0; i < m.size(); ++i){
        cout << m[i][i] << endl;
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

//recursive strassen matrix mult algo
vector<vector<int> > strassen_multiply(vector<vector<int> > &m, vector<vector<int> > &n, int thresh){
    if(m.size() < thresh){
        return matrix_multiply(m, n);
    }
    else{
        bool needsTrim = false;
        if (m.size() % 2 == 1){
            for(int i = 0; i < m.size(); i++){
                m[i].push_back(0);
            }
            m.push_back(vector<int> (m.size() + 1, 0));
            for(int i = 0; i < n.size(); i++){
                n[i].push_back(0);
            }   
            n.push_back(vector<int> (n.size() + 1, 0));
            needsTrim = true;
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

        vector<vector<int> > a, d, e, h;
        for(int row = 0; row < m.size() / 2; row++){
            vector<int> temp, temp2;
            for(int col = 0; col < m.size() / 2; col++){
                temp.push_back(m[row][col]);
                temp2.push_back(n[row][col]);
            }
            a.push_back(temp);
            e.push_back(temp2);
        }

        for(int row = m.size() / 2; row < m.size(); row++){
            vector<int> temp, temp2;
            for(int col = m.size() / 2; col < m.size(); col++){
                temp.push_back(n[row][col]);
                temp2.push_back(m[row][col]);
            }
            h.push_back(temp);
            d.push_back(temp2);
        }

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

        vector<vector<int> > result;
        for(int row = 0; row < p6.size(); row++){
            vector<int> tempRow(p6.size() * 2, 0);
            for(int col = 0; col < p6.size(); col++){
                tempRow[col] = p6[row][col];
                tempRow[col + p6.size()] = p1[row][col];
            }
            if(needsTrim){
                tempRow.pop_back();
            }
            result.push_back(tempRow);
        }
        for(int row = 0; row < p6.size(); row++){
            vector<int> tempRow(p6.size() * 2, 0);
            for (int col = 0; col < p3.size(); col++){
                tempRow[col] = p3[row][col];
                tempRow[col + p6.size()] = p5[row][col];
            }
            if(needsTrim){
                tempRow.pop_back();
            }
            result.push_back(tempRow);
        }
        if(needsTrim){
            result.pop_back();
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
    fstream file;
    file.open(fileName.c_str()); 
    int curr = 0;
    
    if (flag == 0){
        vector<vector<int> > m(dim, vector<int> (dim, 0));
        vector<vector<int> > n(dim, vector<int> (dim, 0));
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
        print_diag(strassen_multiply(m, n, 100));
    }
    
    if (flag == 1){
        cout << "N = EVEN" << endl;
        fstream file1;
        //change testing bounds for n here
        for(int n = 2; n < 70; n += 2){
            genRandMats(0, 2, n, fileName.c_str());
            file1.open(fileName.c_str());
            vector<vector<int> > mtx1(n, vector<int> (n, 0));
            vector<vector<int> > mtx2(n, vector<int> (n, 0));
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    file1 >> curr;
                    mtx1[i][j] = curr;
                }
            }
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    file1 >> curr;
                    mtx2[i][j] = curr;

                }       
            }
            file1.close();

            double sum = 0.0;
            double sum2 = 0.0;
            int trials = 100;

            for (int i = 0; i < trials; i++){
                auto start = chrono::high_resolution_clock::now(); 
                strassen_multiply(mtx1, mtx2, n);
                auto stop = chrono::high_resolution_clock::now(); 
                auto duration = chrono::duration_cast<chrono::microseconds>(stop - start); 
                sum += duration.count();
                auto start2 = chrono::high_resolution_clock::now(); 
                matrix_multiply(mtx1, mtx2);
                auto stop2 = chrono::high_resolution_clock::now(); 
                auto duration2 = chrono::duration_cast<chrono::microseconds>(stop2 - start2);
                sum2 += duration2.count();
            }
            sum /= trials;
            sum2 /= trials;

            cout << n << "," << sum - sum2 << "," << endl;
        }
    }

    if (flag == 2){
        cout << "N = ODD" << endl;;
        fstream file1;
        //change testing bounds for n here
        for(int n = 3; n < 70; n += 2){
            genRandMats(0, 2, n, fileName.c_str());
            file1.open(fileName.c_str());
            vector<vector<int> > mtx1(n, vector<int> (n, 0));
            vector<vector<int> > mtx2(n, vector<int> (n, 0));
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    file1 >> curr;
                    mtx1[i][j] = curr;
                }
            }
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    file1 >> curr;
                    mtx2[i][j] = curr;

                }       
            }
            file1.close();

            double sum = 0.0;
            double sum2 = 0.0;
            int trials = 100;

            for (int i = 0; i < trials; i++){
                auto start = chrono::high_resolution_clock::now(); 
                strassen_multiply(mtx1, mtx2, n);
                auto stop = chrono::high_resolution_clock::now(); 
                auto duration = chrono::duration_cast<chrono::microseconds>(stop - start); 
                sum += duration.count();
                auto start2 = chrono::high_resolution_clock::now(); 
                matrix_multiply(mtx1, mtx2);
                auto stop2 = chrono::high_resolution_clock::now(); 
                auto duration2 = chrono::duration_cast<chrono::microseconds>(stop2 - start2);
                sum2 += duration2.count();
            }
            sum /= trials;
            sum2 /= trials;

            cout << n << "," << sum - sum2 << "," << endl;
        }
    }

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
        vector<vector<int> > a_squared = strassen_multiply(A, A, 100);
        vector<vector<int> > res = strassen_multiply(a_squared, A, 100);
        int sum = 0;
        for (int i = 0; i < res.size(); i++){
            sum += res[i][i];
        }
        //print_diag(strassen_multiply(strassen_multiply(A, A, 512), A, 512));
        cout << "NUM TRIANGLES: " << (sum / 6) << endl;
        return 0;
    }

}
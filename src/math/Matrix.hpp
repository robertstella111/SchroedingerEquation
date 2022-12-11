#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <chrono>

class Mat{
    private:
        std::vector<std::vector<double>> mat;
        unsigned dim;
    public:


        Mat(unsigned i) {
            std::vector<double>  buff(i, 0);
            for(unsigned k = 0; k < i; k++) {
                mat.push_back(buff);
            }
            dim = i;
        }

        void createIdentity() {

            for(unsigned k = 0; k < mat.size(); k++) {
                for(unsigned i = 0; i < mat.size(); i++) {
                    mat[k][i] = 0;
                }
                mat[k][k] = 1;
            }
        }

        void print() {
            for(unsigned k = 0; k < mat.size(); k++) {
                for(unsigned i = 0; i < mat.size(); i++) {
                    std::cout << mat[k][i] << " ";
                }
                std::cout << std::endl;
            }
        }

        double& insert(unsigned row, unsigned column) {
            return mat[row][column];
        }

        void max(unsigned *row, unsigned *column) {
            double buff = 0;
            for(unsigned i = 0; i < dim; i++) {
                for(unsigned j = i+1; j < dim; j++) {
                    if(std::abs(mat[i][j])>buff) {
                        *row = i;
                        *column = j;
                        buff = std::abs(mat[i][j]);
                    }
                }
            }
        }

        double& getCoeff(unsigned i, unsigned j) {
            return mat[i][j];
        }

        Mat operator*(Mat right) {
            Mat res = Mat(dim);
            double val;
            for(unsigned i = 0; i < dim; i++) {
                for(unsigned j = 0; j < dim; j++) {
                    val = 0;
                    for(unsigned k = 0; k < dim; k++) {
                        val += mat[i][k]*right.mat[k][j];
                    } 
                    res.insert(i,j) = val;  
                }
            }
            return res;

        }

        unsigned getDim() {
            return dim;
        }
};

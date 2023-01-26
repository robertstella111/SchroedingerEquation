#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <chrono>



class EigenValue{
    private:
        Mat *result;
        Mat eigvec = Mat(3);
        std::vector<double> eigvals;
        double epsilon = 1e-2;
    public:
        void setEpsilon(double eps) {
            epsilon = eps;
        }
        void solveJacobi(unsigned number, Mat *mat) {
            unsigned row_index, column_index;
            result = mat;
            unsigned dim = mat->getDim() - 1;
            
            eigvec = Mat(mat->getDim());
            eigvec.createIdentity();
            result->initMaxIndex();
            double angle;
            unsigned counter1 = dim;
            for(unsigned k = 0; k < number; k++) {
            for(unsigned i = 0; i < counter1;) {
                for(unsigned j = 0; j < i +1;  j++){
                    row_index = j % (dim);
                    column_index = result->index(row_index);
                    if(std::abs(result->getCoeff(row_index, row_index) - result->getCoeff(column_index,column_index)) < 1e-7) {
                        if(result->getCoeff(row_index,column_index) > 0) angle = M_PI/4;
                        else angle = -M_PI/4;
                    }
                    else {
                        angle =  0.5*std::atan((2*result->getCoeff(row_index, column_index))/(result->getCoeff(column_index, column_index) - result->getCoeff(row_index, row_index)));
                        if((result->getCoeff(row_index, row_index) - result->getCoeff(column_index, column_index)) < 0) angle += M_PI/2;
                    }
                
                    result->Transform1D(angle, row_index, column_index);  //rotate matrix A
                    eigvec.EigVec(angle,row_index, column_index); //unit matrix
                    
                }
                
                if(i + 5 < dim) i = i+5;
                else if (i + 1 == dim) break;
                else i = dim -1;
            }
            if(result->maxEntry() <= epsilon) break;
            }
            
            for(unsigned i = 0; i < result->getDim(); i++) eigvals.push_back(result->getCoeff(i, i));           
    
        }

        std::vector<double> EigenVector(unsigned i) {
            return eigvec.getEigenvector(i);
        }

        double getEigenValue(unsigned i) {
            return eigvals[i];
        }
};
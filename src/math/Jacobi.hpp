#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <chrono>



class EigenValue{
    private:
        SparseMat result = SparseMat(3);
        SparseMat eigvec = SparseMat(3);
        std::vector<double> eigvals;
    public:
        
        void solveJacobi(unsigned number, SparseMat mat) {
            unsigned row_index, column_index;
            result = mat;
            eigvec = SparseMat(mat.getDim());
            eigvec.createIdentity();
            double angle;
            unsigned i;
            for(i = 0; i < number; i++) {
                row_index = 0;
                column_index = 0;
                result.max(&row_index, &column_index); 
                if(row_index == 0 && column_index == 0) break;
                if(std::abs(result.getCoeff(row_index, row_index) - result.getCoeff(column_index,column_index)) < 1e-7) {
                    if(result.getCoeff(row_index,column_index) > 0) angle = M_PI/4;
                    else angle = -M_PI/4;
                }
                else {
                    angle =  0.5*std::atan((2*result.getCoeff(row_index, column_index))/(result.getCoeff(column_index, column_index) - result.getCoeff(row_index, row_index)));
                    if((result.getCoeff(row_index, row_index) - result.getCoeff(column_index, column_index)) < 0) angle += M_PI/2;
                }
              
                result.Transform1D(angle, row_index, column_index); 
                eigvec.EigVec(angle,row_index, column_index);
            }
           
            for(unsigned i = 0; i < result.getDim(); i++) eigvals.push_back(result.getCoeff(i, i));
    
        }

        std::vector<double> EigenVector(unsigned i) {
            return eigvec.getEigenvector(i);
        }

        double getEigenValue(unsigned i) {
            return eigvals[i];
        }
};
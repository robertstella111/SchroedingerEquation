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
        double epsilon = 1e-2;
    public:
        void setEpsilon(double eps) {
            epsilon = eps;
        }
        void solveJacobi(unsigned number, SparseMat mat) {
            unsigned row_index, column_index;
            result = mat;
            unsigned dim = mat.getDim() - 1;
            
            eigvec = SparseMat(mat.getDim());
            eigvec.createIdentity();
            result.initMaxIndex();
            double angle;
            unsigned i;
            for(i = 0; i < number; i++) {
                row_index = i % (dim);
                column_index = result.index(row_index);
                if(i % (20*dim) == 0) {
                    if(result.maxEntry() <= epsilon) break;
                }
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
           // for(unsigned i = 0; i < result.getDim(); i++) std::cout << result.getCoeff(i,i) << " " << result.maxRowIndex(i) << std::endl;
            
    
        }

        std::vector<double> EigenVector(unsigned i) {
            return eigvec.getEigenvector(i);
        }

        double getEigenValue(unsigned i) {
            return eigvals[i];
        }
};
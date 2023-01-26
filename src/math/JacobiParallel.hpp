#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <chrono>



class EigenValueParallel{
    private:
        MatParallel *result;
        MatParallel eigvec = MatParallel(1,3);
        std::vector<double> eigvals;
        double epsilon = 1e-2;
        unsigned NUM_THREADS = 1;
    public:
        void setEpsilon(double eps) {
            epsilon = eps;
        }

        EigenValueParallel(unsigned n) {
            NUM_THREADS = n;
        }

        void solveJacobi(unsigned number, MatParallel *mat) {
            unsigned row_index, column_index;
            std::vector<unsigned> all_rows(NUM_THREADS), all_column(NUM_THREADS);
            result = mat;
            unsigned dim = mat->getDim();
            eigvec = MatParallel(NUM_THREADS,mat->getDim());
            eigvec.createIdentity();
            result->initMaxIndex();
            std::vector<double> angles(NUM_THREADS);
            double angle = 0;
        
            //first part of algorithm
            for(unsigned outerCounter = 0; outerCounter < number; outerCounter++) {
                for(unsigned counter1 = dim - dim/NUM_THREADS; counter1 < dim; counter1 = counter1 + NUM_THREADS) {
                    for(unsigned counter2 = 0; counter2 < counter1 + NUM_THREADS;  counter2 = counter2 + NUM_THREADS){
                        for(unsigned s = 0; s < NUM_THREADS; s++) {
                            all_column[s] = dim;
                            all_rows[s] = dim;
                        }
                        all_rows[0] = counter2;
                        all_column[0] = result->getMaxEntryRow(counter2, epsilon);
                        if(all_column[0] == counter2) {
                            continue;
                        } 
                        for(unsigned s =1; s < NUM_THREADS; s++) {
                                for(unsigned searchCounter = 0; searchCounter < dim; searchCounter = searchCounter + NUM_THREADS) {
                                unsigned currentRow = counter2 + searchCounter;
                                currentRow+= s;
                                
                                //
                                while(std::find(all_column.begin(), all_column.end(), currentRow) != all_column.end() || currentRow +NUM_THREADS > dim) {
                                    if(currentRow+2*NUM_THREADS - 1 >= dim) currentRow = s;
                                    else currentRow+=NUM_THREADS;
                                }

                                all_rows[s] = currentRow;
                                result->getMaxEntryRow2(s, currentRow,epsilon,&all_column, &all_rows);
                                if(all_column[s] != currentRow) break;
                            }

                        }
                       
                        #pragma omp parallel private(row_index, column_index,angle) shared(result, all_rows, all_column, angles) 
                            {
                            angle = 0;
                            unsigned c3 = omp_get_thread_num();                           

                            row_index = all_rows[c3];
                            column_index = all_column[c3];
                            if(row_index != dim) {
                                column_index = all_column[c3];
                                if(std::abs(result->getCoeff(row_index, row_index) - result->getCoeff(column_index,column_index)) < 1e-7) {
                                    if(result->getCoeff(row_index,column_index) > 0) angle = M_PI/4;
                                    else angle = -M_PI/4;
                                }
                                else {
                                    angle =  0.5*std::atan((2*result->getCoeff(row_index, column_index))/(result->getCoeff(column_index, column_index) - result->getCoeff(row_index, row_index)));
                                    if((result->getCoeff(row_index, row_index) - result->getCoeff(column_index, column_index)) < 0) angle += M_PI/2;
                                }
                            }
                            #pragma omp barrier
                            angles[c3] = angle;
                            if(row_index%NUM_THREADS == column_index%NUM_THREADS) result->left(row_index, column_index, angle);
                    }

                    for(unsigned leftCounter = 0; leftCounter < NUM_THREADS; leftCounter++) {
                        if(all_rows[leftCounter]%NUM_THREADS != all_column[leftCounter]%NUM_THREADS) result->left(all_rows[leftCounter], all_column[leftCounter], angles[leftCounter]);
                    }
                    for(unsigned rightCounter = 0; rightCounter < NUM_THREADS; rightCounter++) {
                        result->right(all_rows[rightCounter], all_column[rightCounter], angles[rightCounter]);
                    }
                    for(unsigned rightCounter = 0; rightCounter < NUM_THREADS; rightCounter++) {
                        eigvec.EigVec(angles[rightCounter],all_rows[rightCounter], all_column[rightCounter]);
                    }
                    }
                    
                }
            }
            //fine tuning in last lines of matrix
            for(unsigned k = 0; k < number; k++) {
                for(unsigned i = 0; i <dim;  i++) {
                    for(unsigned j = 0; j < i+1;  j++){
                            if(j+1 == dim) break;
                            angle = 0;
                            row_index = j;
                            
                            if(row_index != dim) {
                                column_index = result->getMaxEntryRow(j,epsilon);
                                if(std::abs(result->getCoeff(row_index, row_index) - result->getCoeff(column_index,column_index)) < 1e-7) {
                                    if(result->getCoeff(row_index,column_index) > 0) angle = M_PI/4;
                                    else angle = -M_PI/4;
                                }
                                else {
                                    angle =  0.5*std::atan((2*result->getCoeff(row_index, column_index))/(result->getCoeff(column_index, column_index) - result->getCoeff(row_index, row_index)));
                                    if((result->getCoeff(row_index, row_index) - result->getCoeff(column_index, column_index)) < 0) angle += M_PI/2;
                                }
                            }
                            else continue;
                            if(column_index == row_index)continue;
                            result->left(row_index, column_index, angle);
                            result->right(row_index, column_index, angle);
                            eigvec.EigVec(angle, row_index, column_index);
                    }  
                }
            }  

            for(unsigned i = 0; i < result->getDim(); i++) eigvals.push_back(result->getCoeff(i, i));
           
        }
        
        double getEigenValue(unsigned i) {
            return eigvals[i];
        }

        std::vector<double> EigenVector(unsigned i) {
            return eigvec.getEigenvector(i);
        }
};
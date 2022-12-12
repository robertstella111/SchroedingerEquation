#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <chrono>

class SparseMat {
    private:
        std::vector<std::vector<double>> list;
        std::vector<unsigned> max_index_row;    
        std::vector<double> max_value_row;  
        unsigned dim;
    public:

    SparseMat(unsigned i) {
            std::vector<double>  buff(i, 0);
            for(unsigned k = 0; k < i; k++) {
                list.push_back(buff);
            }
            dim = i;
    }
/*
    SparseMat(unsigned i, std::map<unsigned, std::map<unsigned, double>>* buff) {
            list = *buff;
    }
*/  
    double maxRowIndex(unsigned i) {
        return max_value_row[i];
    }
    void initMaxIndex() {
        for(unsigned k = 0; k < list.size() - 1; k++) {
            max_index_row.push_back(k);
            max_value_row.push_back(0);
            for(unsigned i = k+1; i < list.size(); i++) {
                if(std::abs(list[k][i]) > max_value_row[k]) {
                     max_index_row[k] = i;
                    max_value_row[k] = std::abs(list[k][i]);
                }
            }
        }
        max_index_row.push_back(dim);
        max_value_row.push_back(0);
    }

    void createIdentity() {
            for(unsigned k = 0; k < dim; k++) {
                list[k][k] = 1;
            }
        }

        void print() {
            for(unsigned k = 0; k < list.size(); k++) {
                for(unsigned i = 0; i < list.size(); i++) {
                    std::cout << list[k][i] << " ";
                }
                std::cout << std::endl;
            }
        }

        double& insert(unsigned row, unsigned column) {
            return list[row][column];
        }

        void max(unsigned *row, unsigned *column) {
            double buff = 0;
            for(unsigned i = 0; i < dim; i++) {
                for(unsigned j = i+1; j < dim; j++) {
                    if(std::abs(list[i][j])>buff) {
                        *row = i;
                        *column = j;
                        buff = std::abs(list[i][j]);
                    }
                }
            }
        }

        unsigned index(unsigned row) {
            return max_index_row[row];
        }
        double maxEntry() {
            double max = 0;
            for(unsigned i = 0; i < dim; i++) {
                if(max_value_row[i] > max) max = max_value_row[i];
            }
            return max;
        }



        double& getCoeff(unsigned i, unsigned j) {
           return list[i][j];
        }

        void deleteCoeff(unsigned i, unsigned j) {
             list[i][j] = 0;
            
        }

        void sortColumns() {
        }

        void sortRows() {
        }

        unsigned getDim() {
            return dim;
        }

        void MaxEntryRow(unsigned row) {
            double buffer = 0;
            unsigned index_buffer = row+1;
            for(unsigned i = row + 1; i < dim; i++) {
                if(std::abs(list[row][i]) > buffer) {
                    index_buffer = i;
                    buffer = std::abs(list[row][i]);
                }
            }
            max_index_row[row] = index_buffer;
            max_value_row[row]=buffer;
        }

        void Transform1D(double angle, unsigned p, unsigned q) {
            double  s,c, app, aqq, buff;
            s = std::sin(angle);
            c = std::cos(angle);
            app = pow(c,2)*list[p][p]+pow(s,2)*list[q][q]-2*s*c*list[p][q];
            aqq = pow(s,2)*list[p][p]+pow(c,2)*list[q][q]+2*s*c*list[p][q];
            
            list[p][q] = 0;
            MaxEntryRow(p);
            //update in api : i > p 

            for(unsigned i = 0; i < p; i++) {
                //p > i, q > i
                list[p][i] = list[i][p]; //store old
                buff = c*list[p][i] - s*list[i][q];
                list[i][p] = buff; //update new
                if (p == max_index_row[i])  MaxEntryRow(i);
                else if (std::abs(list[i][p])>max_value_row[i]) {
                    max_index_row[i] = p;
                    max_value_row[i] = std::abs(list[i][p]);
                }
            }
            for(unsigned i = p+1; i < q; i++) {
                //i > p, q > i
               
                list[i][p] = list[p][i]; //store old
                buff = c*list[i][p] - s*list[i][q];
                list[p][i] = buff; //update new
                if (i == max_index_row[p])  MaxEntryRow(p);
                else if (std::abs(list[p][i])>max_value_row[p]) {
                    max_index_row[p] = i;
                    max_value_row[p] = std::abs(list[p][i]);
                }
                
            }
            for(unsigned i = q+1; i < dim; i++) {
                //p < i, q < i
                list[i][p] = list[p][i]; //store old
                buff = c*list[i][p] - s*list[q][i];
                list[p][i] = buff; //update new
                if (i == max_index_row[p])  MaxEntryRow(p);
                else if (std::abs(list[p][i])>max_value_row[p]) {
                    max_index_row[p] = i;
                    max_value_row[p] = std::abs(list[p][i]);
                }
            }


            for(unsigned i = 0; i < p; i++) {
                //p > i, q > i
                buff =  s*list[p][i]+c*list[i][q];
                list[i][q] = buff;
                if (q == max_index_row[i])  MaxEntryRow(i);
                else if (std::abs(list[i][q])>max_value_row[i]) {
                    max_index_row[i] = q;
                    max_value_row[i] = std::abs(list[i][q]);
                }
            }
            for(unsigned i = p+1; i < q; i++) {
                //i > p, q > i
                buff =  s*list[i][p]+c*list[i][q];
                list[i][q] = buff;
                if (q == max_index_row[i])  MaxEntryRow(i);
                else if (std::abs(list[i][q])>max_value_row[i]) {
                    max_index_row[i] = q;
                    max_value_row[i] = std::abs(list[i][q]);
                }
            }
            for(unsigned i = q+1; i < dim; i++) {
                //p < i, q < i
                buff = s*list[i][p]+c*list[q][i];
                list[q][i] = buff;
                if (i == max_index_row[q])  MaxEntryRow(q);
                else if (std::abs(list[q][i])>max_value_row[q]) {
                    max_index_row[q] = i;
                    max_value_row[q] = std::abs(list[q][i]);
                }
            }


            list[p][p] = app;
            list[q][q] = aqq;
        }


        void EigVec(double angle, unsigned p, unsigned q) {
            double  s,c;
            s = std::sin(angle);
            c = std::cos(angle);
            for(unsigned v = 0; v < dim; v++) {
                double buff = list[p] [v];
                list[p][v] = c*list[p][v] - s*list[q][v];
                list[q][v] = s*buff    + c*list[q][v];
            }
        }

        std::vector<double> getEigenvector(unsigned k) {
            std::vector<double> res;
            for(unsigned i = 0; i < dim; i++) res.push_back(list[k][i]);
            return res;
        }

        
};

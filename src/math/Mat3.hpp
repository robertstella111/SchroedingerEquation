#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <list>
#include <algorithm>
#include <map>
#include <chrono>
#include <memory>


class SparseMatPointer {
    private:
        std::unique_ptr<std::unique_ptr<double[]>[]> list;
        std::vector<unsigned> max_index_row;    
        std::vector<double> max_value_row;  
        unsigned dim;
    public:

    SparseMatPointer(unsigned n) {
        dim = n;
        std::cout << "Constructor" << std::endl;
        list =  std::make_unique<std::unique_ptr<double[]>[]>(n);

        for(unsigned i = 0; i < dim; i++) {
          list.get()[i] = std::make_unique<double[]>(n);
        }
    }

    void print() {
            auto p1 = list.get();
            for(unsigned k = 0; k < dim; k++) {
                auto p2 = p1[k].get();
                for(unsigned i = 0; i < dim; i++) {
                    std::cout << p2[i] << " ";
                }
                std::cout << std::endl;
            }
    }

    void createIdentity() {
        auto p1 = list.get();
        for(unsigned k = 0; k < dim; k++) {
            auto p2 = p1[k].get();
            for(unsigned i = 0; i < dim; i++) {
                p2[i] = 0;
            }
            p2[k] = 1;
        }
    }
    double maxRowIndex(unsigned i) {
        return max_value_row[i];
    }


    void initMaxIndex() {
        auto p1 = list.get();
        for(unsigned k = 0; k < dim - 1; k++) {
            max_index_row.push_back(k);
            max_value_row.push_back(0);
            auto p2 = p1[k].get();
            for(unsigned i = k+1; i < dim; i++) {
                if(std::abs(p2[i]) > max_value_row[k]) {
                    max_index_row[k] = i;
                    max_value_row[k] = std::abs(p2[i]);
                }
            }
        }
        max_index_row.push_back(dim);
        max_value_row.push_back(0);
    }


    double& insert(unsigned row, unsigned column) {
        return list.get()[row].get()[column];
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
    double& getCoeff(unsigned row, unsigned column) {
        return list.get()[row].get()[column];
    }

    void deleteCoeff(unsigned row, unsigned column) {
        list.get()[row].get()[column] = 0;
            
    }

    void MaxEntryRow(unsigned row) {
        double buffer = 0;
        unsigned index_buffer = row+1;
        auto p = list.get()[row].get();
        for(unsigned i = row + 1; i < dim; i++) {
            if(std::abs(p[i]) > buffer) {
                index_buffer = i;
                buffer = std::abs(p[i]);
            }
        }
        max_index_row[row] = index_buffer;
        max_value_row[row]=buffer;
    }
    void Transform1D(double angle, unsigned p, unsigned q) {
            double  s,c, app, aqq, buff, buff2;
            s = std::sin(angle);
            c = std::cos(angle);
            app = pow(c,2)*getCoeff(p,p)+pow(s,2)*getCoeff(q,q)-2*s*c*getCoeff(p,q);
            aqq = pow(s,2)*getCoeff(p,p)+pow(c,2)*getCoeff(q,q)+2*s*c*getCoeff(p,q);
            
            getCoeff(p,q) = 0;
            MaxEntryRow(p);
            //update in api : i > p 

            for(unsigned i = 0; i < p; i++) {
                //p > i, q > i
                getCoeff(p,i) = getCoeff(i,p); //store old
                buff = c*getCoeff(p,i) - s*getCoeff(i,q);
                getCoeff(i,p) = buff; //update new
                buff2 = std::abs(buff);
                if (p == max_index_row[i])  MaxEntryRow(i);
                else if (buff2>max_value_row[i]) {
                    max_index_row[i] = p;
                    max_value_row[i] = buff2;
                }
            }
            for(unsigned i = p+1; i < q; i++) {
                //i > p, q > i
               
                getCoeff(i,p) = getCoeff(p,i); //store old
                buff = c*getCoeff(i,p) - s*getCoeff(i,q);
                getCoeff(p,i) = buff; //update new
                buff2 = std::abs(buff);
                if (i == max_index_row[p])  MaxEntryRow(p);
                else if (buff2>max_value_row[p]) {
                    max_index_row[p] = i;
                    max_value_row[p] = buff2;
                }
                
            }
            for(unsigned i = q+1; i < dim; i++) {
                //p < i, q < i
                getCoeff(i,p) = getCoeff(p,i); //store old
                buff = c*getCoeff(i,p) - s*getCoeff(q,i);
                getCoeff(p,i) = buff; //update new
                buff2 = std::abs(buff);
                if (i == max_index_row[p])  MaxEntryRow(p);
                else if (buff2>max_value_row[p]) {
                    max_index_row[p] = i;
                    max_value_row[p] = buff2;
                }
            }


            for(unsigned i = 0; i < p; i++) {
                //p > i, q > i
                buff =  s*list.get()[p].get()[i]+c*list.get()[i].get()[q];
                list.get()[i].get()[q] = buff;
                buff2 = std::abs(buff);
                if (q == max_index_row[i])  MaxEntryRow(i);
                else if (buff2>max_value_row[i]) {
                    max_index_row[i] = q;
                    max_value_row[i] =buff2;
                }
            }
            for(unsigned i = p+1; i < q; i++) {
                //i > p, q > i
                buff =  s*list.get()[i].get()[p]+c*list.get()[i].get()[q];
                list.get()[i].get()[q] = buff;
                buff2 = std::abs(buff);
                if (q == max_index_row[i])  MaxEntryRow(i);
                else if (buff2>max_value_row[i]) {
                    max_index_row[i] = q;
                    max_value_row[i] = buff2;
                }
            }
            for(unsigned i = q+1; i < dim; i++) {
                //p < i, q < i
                buff = s*list.get()[i].get()[p]+c*list.get()[q].get()[i];
                list.get()[q].get()[i] = buff;
                buff2 = std::abs(buff);
                if (i == max_index_row[q])  MaxEntryRow(q);
                else if (buff2>max_value_row[q]) {
                    max_index_row[q] = i;
                    max_value_row[q] = buff2;
                }
            }


            list.get()[p].get()[p] = app;
            list.get()[q].get()[q] = aqq;
        }


        unsigned getDim() {
            return dim;
        }
        
};

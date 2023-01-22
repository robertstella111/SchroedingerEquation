#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <list>
#include <algorithm>
#include <map>
#include <chrono>
#include <memory>
#include <algorithm>
#include <omp.h>


typedef std::pair<unsigned int, double> IndexValuePair;

IndexValuePair myMin(IndexValuePair a, IndexValuePair b){
    return a.second > b.second ? a : b;
}


    #pragma omp declare reduction \
        (minPair:IndexValuePair:omp_out=myMin(omp_out, omp_in)) \
        initializer(omp_priv = IndexValuePair(0, -1))



class SparseMatParallel {
    private:
        std::unique_ptr<std::unique_ptr<double[]>[]> list;
        std::vector<unsigned> max_index_row;    
        std::vector<double> max_value_row;  
        unsigned dim;
        unsigned NUM_THREADS = 1;
    public:

    SparseMatParallel(unsigned Threads, unsigned n) {
        dim = n;
        NUM_THREADS = Threads;
        list =  std::make_unique<std::unique_ptr<double[]>[]>(n);
        omp_set_num_threads(NUM_THREADS);
        std::cout << "dim " << dim << std::endl; 
       #pragma omp parallel 
            {
            for(unsigned i = 0; i < dim; i = i+NUM_THREADS) {
                list.get()[i+omp_get_thread_num()] = std::make_unique<double[]>(dim);
                }
            
        }
        
       
        
    }
    void printRow(unsigned row) {
            auto p1 = list.get();
                auto p2 = p1[row].get();
                for(unsigned i = 0; i < dim; i++) {
                    std::cout << p2[i] << std::endl;
                }
                std::cout << std::endl;
            
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

        void max(unsigned *row, unsigned *column) {
            double buff = 0;
            for(unsigned i = 0; i < dim; i++) {
                for(unsigned j = i+1; j < dim; j++) {
                    auto buff2 = getCoeff(i,j);
                    if(std::abs(buff2)>buff) {
                        *row = i;
                        *column = j;
                        buff = std::abs(buff2);
                    }
                }
            }
        }


        void maxParallel(unsigned start, unsigned *row, unsigned *column) {
            double buff = 0;
            #pragma omp parallel shared(buff, row, column)
                {
                    for(unsigned i = start; i < dim; i += NUM_THREADS) {
                        for(unsigned k = i + omp_get_thread_num() + 1; k < dim; k++) {
                            auto buff2 = getCoeff(i,k);
                            if(std::abs(buff2)>buff) {
                                *row = i;
                                *column = k;
                                buff = std::abs(buff2);
                            }
                        }
                    }
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
        double buffer = -1;
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

    unsigned getMaxEntryRow(unsigned row, double epsilon) {
        double buffer = epsilon;
        unsigned index_buffer = row;
        auto p = list.get()[row].get();
        for(unsigned i = row+1; i < dim; i++) {
            if(std::abs(p[i]) > buffer) {
                index_buffer = i;
                buffer = std::abs(p[i]);
            }
        }
        return index_buffer;
        
    }


    unsigned getMaxEntryRowParallel(unsigned row, double epsilon) {
        IndexValuePair minValueIndex(row, epsilon);
        #pragma omp parallel reduction(minPair:minValueIndex)
            {
                for(unsigned i = row; i < dim; i+=NUM_THREADS) {
                    if(i+omp_get_thread_num() > row && i+omp_get_thread_num() < dim) {
                        if (std::abs(getCoeff(i+omp_get_thread_num(), row)) > minValueIndex.second) {
                            minValueIndex.first = i+omp_get_thread_num();
                           minValueIndex.second= std::abs(getCoeff(i+omp_get_thread_num(), row));
                        }

                    }
                    
                }
            }

        return minValueIndex.first;
        
    }


    void getMaxEntryRow2(unsigned thread, unsigned row, double epsilon, std::vector<unsigned> *columns, std::vector<unsigned> *rows) {
        double buffer = epsilon;
        unsigned index_buffer = row;
        auto p = list.get()[row].get();
        for(unsigned i = row+1; i < dim; i++) {
            if(std::abs(p[i]) > buffer) {
                
                if(std::find(columns->begin(), columns->end(), i) == columns->end() && 
                std::find(rows->begin(), rows->end(), i) == rows->end()){
                    index_buffer = i;
                    buffer = std::abs(p[i]);
                }
               
            }
        } 
        (*columns)[thread] = index_buffer;
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


        void right(unsigned p, unsigned q, double angle) {
            //p < q
            double  s,c, buff, buff2;
            s = std::sin(angle);
            c = std::cos(angle);
            #pragma omp parallel
                {
                    for(unsigned i = 0; i < dim; i = i + NUM_THREADS) {
                        buff = getCoeff(i + omp_get_thread_num(),p)*c - getCoeff(i + omp_get_thread_num(),q)*s;
                        buff2 = getCoeff(i+ omp_get_thread_num(),p)*s + getCoeff(i + omp_get_thread_num(),q)*c;
                        getCoeff(i + omp_get_thread_num(),q) = buff2;
                        getCoeff(i + omp_get_thread_num(),p) = buff;
                    }
            }
            

        }

                    

        void left(unsigned p, unsigned q, double angle) {
            //p < q
            double  s,c, buff, buff2;
            s = std::sin(angle);
            c = std::cos(angle);
            if(p % NUM_THREADS != q % NUM_THREADS) {
               for(unsigned i = 0; i < dim; i ++) {
                //no parallel
                    buff = getCoeff(p,i)*c - getCoeff(q,i)*s;
                    buff2 = getCoeff(p,i)*s + getCoeff(q,i)*c;
                    getCoeff(q,i) = buff2;
                    getCoeff(p,i) = buff;
                }  
            
            }
            
            else {
               for(unsigned i = 0; i < dim; i ++) {
                //no parallel
                    buff = getCoeff(p,i)*c - getCoeff(q,i)*s;
                    buff2 = getCoeff(p,i)*s + getCoeff(q,i)*c;
                    getCoeff(q,i) = buff2;
                    getCoeff(p,i) = buff;
                }  
            }
            
        }

        void EigVec(double angle, unsigned p, unsigned q) {
            double  s,c;
            s = std::sin(angle);
            c = std::cos(angle);
            for(unsigned v = 0; v < dim; v++) {
                double buff = getCoeff(p, v);
                getCoeff(p, v) = c*getCoeff(p,v) - s*getCoeff(q,v);
                getCoeff(q,v) = s*buff    + c*getCoeff(q,v);
            }
        }

        std::vector<double> getEigenvector(unsigned k) {
            std::vector<double> res;
            for(unsigned i = 0; i < dim; i++) res.push_back(getCoeff(k,i));
            return res;
        }

};

                            

#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <chrono>

class SparseMat {
    private:
        std::map<unsigned, std::map<unsigned, double>> list;    
        unsigned dim;
    public:

    SparseMat(unsigned i) {
            dim = i;
    }

    SparseMat(unsigned i, std::map<unsigned, std::map<unsigned, double>>* buff) {
            list = *buff;
    }


    void createIdentity() {
            list.clear();
            for(unsigned k = 0; k < dim; k++) {
                list[k][k] = 1;
            }
        }

        void print() {
            for(auto it = list.begin(); it != list.end(); it++) {
                for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
                    std::cout << it->first << " " << it2->first << " " << it2->second << std::endl;
                }
            }
        }

        double& insert(unsigned row, unsigned column) {
            return list[row][column];
        }

        void max(unsigned *row, unsigned *column) {
            double max = 0;
            auto end1 = list.end();
            for(auto it = list.begin(); it !=end1; it++) {
                auto end2 = it->second.end();
                for(auto it2 = it->second.begin(); it2 != end2; it2++) {
                    if(it->first < it2->first &&std::abs(it2->second) > max) {
                        max = std::abs(it2->second);
                        *row = it->first;
                        *column = it2->first;
                    }
                }
            }
        }

        double& getCoeff(unsigned i, unsigned j) {
           return list[i][j];
        }

        void deleteCoeff(unsigned i, unsigned j) {
             list[i].erase(j);
             if(list[i].size()==0) list.erase(i);
            
        }

        void sortColumns() {
        }

        void sortRows() {
        }


        SparseMat mult(SparseMat *right) {
            auto res = SparseMat(dim);
            auto end1 = list.end();
            for(auto it1 = list.begin(); it1 != end1; it1++){
                auto end2 = it1->second.end();
               for(auto it2 = it1->second.begin(); it2 != end2; it2++){
                auto end3 = right->list.end();
                    for(auto it3 = right->list.begin(); it3 !=end3; it3++){
                        auto end4 = it3->second.end();
                        for(auto it4 = it3->second.begin(); it4 != end4; it4++){
                            if(it2->first == it3->first) res.insert(it1->first,it4->first) += it2->second*it4->second;
                        } 
                    }
                } 
            }
            return res;
        }


        unsigned getDim() {
            return dim;
        }


        void Transform1D(double angle, unsigned p, unsigned q) {
            double  s,c, app, aqq;
            s = std::sin(angle);
            c = std::cos(angle);
            app = pow(c,2)*list[p][p]+pow(s,2)*list[q][q]-2*s*c*list[p][q];
            aqq = pow(s,2)*list[p][p]+pow(c,2)*list[q][q]+2*s*c*list[p][q];
            std::map<unsigned,std::map<unsigned, double>> buff;

            for(auto it = list[p].begin(); it != list[p].end(); it++) {
               if(it->first != p && it->first != q && std::abs(it->second) > 1e-8) {
                    buff[p][it->first] += c*it->second;
                    buff[q][it->first] += s*it->second;
               }
            }

            for(auto it = list[q].begin(); it != list[q].end(); it++) {
               if(it->first != p && it->first != q ) {
                    buff[p][it->first] -= s*it->second;
                    buff[q][it->first] += c*it->second;
               }
            }
            //update elements
            for(auto it1 = buff.begin(); it1 != buff.end(); it1++) {
                for(auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {      
                        list[it1->first][it2->first] = it2->second;
                        list[it2->first][it1->first] = it2->second;
                   
                }
            }
            list[p][p] = app;
            list[q][q] = aqq;
            deleteCoeff(p,q);
            deleteCoeff(q,p);
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

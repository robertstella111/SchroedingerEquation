#include <iostream>

#include <vector>
#include <map>


class Schroedinger1D{
    private:
        double b = 1, a = 0,h;
        unsigned length = 20;
        std::map<double, std::vector<double>> function;
        EigenValue solver;
    
    public:
        Schroedinger1D() {
            h = (b-a)/(length-1);
        }
        
        void setNumPoints(unsigned l) {
            length = l;
            h = (b-a)/(length-1);
        }

         std::vector<double> getEigenvektor(unsigned k) {
            std::vector<double> res;
            auto    buff = solver.EigenVector(k);
            double betrag = 0;
            double a, b, lambda;

            for(unsigned i = 0; i < buff.size()-1; i++) {
                a = buff[i];
                b = buff[i+1];
                a = pow(a,2);
                b = pow(b,2);
                betrag += ((a+b)/2)*h;
             }
            lambda = std::sqrt(1/betrag);
            res.push_back(0);
            for(unsigned i = 0; i < buff.size(); i++) res.push_back(lambda*buff[i]);
            res.push_back(0);
            return res;
            
        }


        void solve(unsigned number) {
            //Eigen::MatrixXd mat(length-2, length-2);
            SparseMat mat = SparseMat(length-2);
            std::cout << "Initializing matrix" << std::endl;
            for(unsigned i = 0; i < length-2; i++) {
                if(i == 0) {
                    mat.insert(i,i)= 2/pow(h,2);
                    mat.insert(i,i+1) = -1/pow(h,2);

                }
                else if(i == length -3) {
                    mat.insert(i,i)= 2/pow(h,2);
                    mat.insert(i,i-1) = -1/pow(h,2);
                }
                else {
                    mat.insert(i,i)= 2/pow(h,2);
                    mat.insert(i,i+1) = -1/pow(h,2);
                    mat.insert(i,i-1) = -1/pow(h,2);
                }
            }
            std::cout << "Solve eigenvalue problem" << std::endl;
            
            solver.solveJacobi(number,mat);
            std::cout << "solved!" << std::endl;
            for(unsigned k = 0; k < length-2; k++)  {
                function.insert(std::make_pair(solver.getEigenValue(k), getEigenvektor(k)));
            }
            
        }

    
        std::vector<double> getXKoord() {
            std::vector<double> res;
            for(unsigned i = 0; i < length; i++) res.push_back(a+i*h);
            return res;
            
        }

        std::vector<double> getEigenValuesSort(){
            std::vector<double> res;
            for(auto it = function.begin(); it != function.end(); it ++) res.push_back(it->first);
            return res;
            
        }

        std::vector<double> getEigenvektorSort(unsigned k) {
            std::vector<double> res;
            unsigned counter = 0;

            for(auto it = function.begin(); it != function.end(); it ++) {
                if(counter == k) {
                    return it->second;
                }
                counter++;
            }
            
            return res;
            
        }


        void setDomain(double a, double b) {
            (*this).a = a;
            (*this).b = b;
            h = (b-a)/(length-1);
        }
        
    
};

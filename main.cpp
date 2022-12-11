#include "io.hpp"
#include <chrono>



int main() {
  
    unsigned length = 10;
    double a,b,h;
    SparseMat mat = SparseMat(length-2);
    a = 0;
    b = 1;
    h = (b-a)/(length-1);
    for(unsigned i = 0; i < length-2; i++) {
                if(i == 0) {
                    mat.insert(i,i)= 2/pow(h,2);
                    mat.insert(i,i+1) = -1/pow(h,2);
                }
                else if(i == length -3) {
                    mat.insert(i,i-1) = -1/pow(h,2);
                    mat.insert(i,i)= 2/pow(h,2);
                    
                }
                else {
                    mat.insert(i,i-1) = -1/pow(h,2);
                    mat.insert(i,i)= 2/pow(h,2);
                    mat.insert(i,i+1) = -1/pow(h,2);
                   
                }
            }

    EigenValue solver;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    solver.solveJacobi(100000,mat);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "EigenValues: error relativ" << std::endl;
    std::vector<double> buff;
    for(unsigned i = 0; i < length-2; i++) buff.push_back(solver.getEigenValue(i));
    std::sort(buff.begin(), buff.end());
    for(unsigned i = 0; i < length-2; i++) std::cout << 100*(buff[i]-pow((i+1)*M_PI,2))/pow((i+1)*M_PI,2)<< " %" << std::endl;
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " Milliseconds" << std::endl;
}
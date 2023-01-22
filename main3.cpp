#include <omp.h>
#include "src/io.hpp"
#include <chrono>

#define THREADS 16
int main() {
    unsigned length = 1010;
    double a,b,h;
    SparseMatParallel mat = SparseMatParallel(THREADS, length - 2);
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
    std::cout << mat.getDim() << std::endl; 
    EigenValueParallel solver(THREADS);
    
    solver.setEpsilon(1e-4);
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    solver.solveJacobi(2,&mat);
    
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "EigenValues: error relativ" << std::endl;
    std::vector<double> buff;
    for(unsigned i = 0; i < length-2; i++) buff.push_back(solver.getEigenValue(i));
    std::sort(buff.begin(), buff.end());
    for(unsigned i = 0; i < length-2; i++) std::cout << 100*(buff[i]-pow((i+1)*M_PI,2))/pow((i+1)*M_PI,2)<< " %" << std::endl;
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " Milliseconds" << std::endl;
    std::cout << "Number of available threads: " << omp_get_max_threads() << std::endl;
    
}


//5.42%
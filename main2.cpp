#include "src/io.hpp"
#include <chrono>



int main() {
      
    unsigned length = 1002;
    double a,b,h;
    SparseMatPointer mat = SparseMatPointer(length-2);
    a = 0;
    b = 1;
    h = (b-a)/(length-1);
   // mat.createIdentity();

    
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
    //mat.print();
    EigenValueP solver;
    solver.setEpsilon(1e2);
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    solver.solveJacobi(10,&mat);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "EigenValues: error relativ" << std::endl;
    std::vector<double> buff;
    for(unsigned i = 0; i < length-2; i++) buff.push_back(solver.getEigenValue(i));
    std::sort(buff.begin(), buff.end());
    for(unsigned i = 0; i < length-2; i++) std::cout << 100*(buff[i]-pow((i+1)*M_PI,2))/pow((i+1)*M_PI,2)<< " %" << std::endl;
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " Milliseconds" << std::endl;
}


//5.42%
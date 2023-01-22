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


int main() {
    unsigned dim = 10000;

        unsigned update = 2;
        auto list =  std::make_unique<std::unique_ptr<double[]>[]>(dim);
        omp_set_num_threads(update);
        unsigned  penis2;
        
        #pragma omp parallel 
            {
            for(unsigned i = 0; i < dim; i = i+update) {
                list.get()[i+omp_get_thread_num()] = std::make_unique<double[]>(dim);
                }
            
        }
       std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


        omp_set_num_threads(update);
        #pragma omp parallel  firstprivate(update) firstprivate(penis2)
            {
            unsigned penis = omp_get_thread_num();
            penis2 = omp_get_thread_num();
            for(unsigned i = 0; i < dim; i = i+update) {
                double buff = 0;
                auto buff2 = list.get();
                for(unsigned k = 0; k < dim; k = k + update) buff += buff2[k + penis].get()[i+penis];
                buff2[i+penis].get()[i+penis] = buff;
            }
             
            
            
        }
        #pragma omp parallel firstprivate(penis2)
            {
            std::cout << penis2 << std::endl; 
        }

/*

        #pragma omp parallel
        if(omp_get_thread_num() == 0) {
            for(unsigned i = 0; i < grenze; i++) {
                for(unsigned k = i+1; k < dim; k++) list.get()[i].get()[i] += list.get()[i].get()[k];
            }
        }
        else {
            for(unsigned i = grenze; i < dim; i++) {
                 for(unsigned k = i+1; k < dim; k++) list.get()[i].get()[i] += list.get()[i].get()[k];
            }
        }
        
        */
         std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " Milliseconds" << std::endl;

}
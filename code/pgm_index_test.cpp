
#include <fstream>
#include "pgm_index.hpp"
#include "greedyplrBasedPGM.hpp"
#include "optimalPLR(C)BasedPGM.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>

#define Filename "../data/fbdata.txt"
template<typename K>
void load_data(std::vector<K> &data) {

    // Load data from a file or generate it
    std::ifstream file(Filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << Filename << std::endl;
        return;
    }
    K value;
    
    while(file >> value) {
        data.push_back(value);
    }
    file.close();
}

template<size_t epsilon = 2, size_t epsilon_recursive = 4>
void experiment_linear(std::vector<double> data, int threads_num = 16) {
   // Build the whole PGM-Index
   // Get the time of constructing the whole PGM-Index
   auto start = std::chrono::high_resolution_clock::now(); 
   pgm::PGMIndex<double, epsilon, epsilon_recursive> pgm_index(data.begin(), data.end());
   auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Time taken to build the index(linear): " << duration.count() << " seconds" << std::endl; 
   // Search for the first 10 elements
    double time_sum = 0;
   for (size_t i = 0; i < 10; ++i) {
        auto start_search = std::chrono::high_resolution_clock::now();
        auto result = pgm_index.search(data[i]);
        auto end_search = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> search_duration = end_search - start_search;
        time_sum += search_duration.count();
    }
    auto avg_time = time_sum / 10;
    std::cout << "Average time taken to search: " << avg_time << " seconds" << std::endl;
    // Print the size of the index
    std::cout << "Size of the index: " << pgm_index.size_in_bytes() << " bytes" << std::endl;
    // Print the number of segments of each level
    auto levels_offsets = pgm_index.get_levels_offsets();
    std::cout << "Number of segments: ";
    for (size_t i = 0; i < levels_offsets.size() - 1; ++i) {
        std::cout << levels_offsets[i + 1] - levels_offsets[i] << " ";
    }
    std::cout << std::endl;
    // Print the height of the index
    std::cout << "Height of the index: " << pgm_index.height() << std::endl;
}

template<size_t epsilon = 2, size_t epsilon_recursive = 4>
void experiment_greedy(std::vector<double> data, int threads_num = 16) {
    // Build the whole PGM-Index
    // Get the time of constructing the whole PGM-Index
    auto start = std::chrono::high_resolution_clock::now(); 
    Greedy::PGMIndex<double,epsilon,epsilon_recursive> greedy_model(data.begin(), data.end());
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Time taken to build the index(greedy): " << duration.count() << " seconds" << std::endl;
    
    double time_sum = 0;
    // Search for the first 10 elements
    
    for (size_t i = 0; i < 10; ++i) {
        auto start_search = std::chrono::high_resolution_clock::now();
        auto result = greedy_model.search(data[i]);
        auto end_search = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> search_duration = end_search - start_search;
        time_sum += search_duration.count();
    }
    auto avg_time = time_sum / 10;
    std::cout << "Average time taken to search: " << avg_time << " seconds" << std::endl;
    // Print the size of the index
    std::cout << "Size of the index: " << greedy_model.size_in_bytes() << " bytes" << std::endl;
    // Print the number of segments of each level
    auto levels_offsets = greedy_model.get_levels_offsets();
    std::cout << "Number of segments: ";
    for(size_t i = 0; i < levels_offsets.size() - 1; ++i) {
        std::cout << levels_offsets[i + 1] - levels_offsets[i] << " ";
    }
    std::cout << std::endl;
    // Print the height of the index
    std::cout << "Height of the index: " << greedy_model.height() << std::endl;
    // Check for epsilon

}



template<size_t epsilon = 2, size_t epsilon_recursive = 4>
void experiment_optimal(std::vector<double> data, int threads_num = 16) {

    // Preprocess the data
    using Point = PGM_C::internal::Point;
    std::vector<Point> points;
    points.reserve(data.size());
    // Convert the data to points
    for (size_t i = 0; i < data.size(); ++i) {
        points.emplace_back(data[i], i);
    }

    // Build the whole PGM-Index
    // Get the time of constructing the whole PGM-Index
    auto start = std::chrono::high_resolution_clock::now(); 
    PGM_C::PGMIndex<double,epsilon,epsilon_recursive> opt_model(points);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Time taken to build the index(optimal): " << duration.count() << " seconds" << std::endl;

     
    double time_sum = 0;
    // // Search for the first 10 elements
    
    for (size_t i = 0; i < 10; ++i) {
        auto start_search = std::chrono::high_resolution_clock::now();
        auto result = opt_model.search(data[i]);
        auto end_search = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> search_duration = end_search - start_search;
        time_sum += search_duration.count();
        // Print the result
        // printf("The real position of the key %lf is %ld\n", data[i], i);
        // printf("The approximate position of the key %lf is %ld\n", data[i], result.pos);    
        // printf("The range of the key %lf is [%ld, %ld]\n", data[i], result.lo, result.hi);
    }
    auto avg_time = time_sum / 10;
    std::cout << "Average time taken to search: " << avg_time << " seconds" << std::endl;
    // Print the size of the index
    std::cout << "Size of the index: " << opt_model.size_in_bytes() << " bytes" << std::endl;
    // Print the number of segments of each level
    auto levels_offsets = opt_model.get_levels_offsets();
    std::cout << "Number of segments: ";
    for(size_t i = 0; i < levels_offsets.size() - 1; ++i) {
        std::cout << levels_offsets[i + 1] - levels_offsets[i] << " ";
    }
    std::cout << std::endl;
    // Print the height of the index
    std::cout << "Height of the index: " << opt_model.height() << std::endl;
    // // Check for epsilon

   
}

int main() {
    // Load the data
    std::vector<double> data;
    load_data(data);
    printf("The size of the data is %ld\n",data.size());
    constexpr size_t epi = 4;
    constexpr size_t epi_recursive = 2;
    // experiment_linear<epi, epi_recursive>(data);
    experiment_optimal<epi, epi_recursive>(data);
    experiment_greedy<epi, epi_recursive>(data);
    experiment_linear<epi, epi_recursive>(data);
    return 0;
}

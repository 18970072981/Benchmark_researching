#include "greedy.hpp"
#include <iostream>
#include <vector>
#include <deque>
#include <limits>
#include <algorithm>
#include <fstream>
#include <chrono>



int main(){
    
    // Read the data from data1.txt
    std::ifstream file("../data/data1.txt");
    std::vector<double> data(20000000);
    
    // double maxi = -std::numeric_limits<double>::max();
    // double mini = std::numeric_limits<double>::max();
    for(int i = 0;i<20000000;i++){
        file >> data[i];
    }

    double epsilon = 2.0;

    // //Check the origin data
    // for(int i = 0;i<=10;i++){
    //     std::cout << "Data " << i << " : " << data[i] << std::endl;
    // }

    using Segment = GreedyPiecewiseLinearModel<double,long unsigned int>::CanonicalSegment::Segment;
    auto result_segments_serial = std::vector<Segment>();
    auto result_segments_parallel = std::vector<Segment>();
    // Initialize the epsilon
    
    // Initialize the GreedyPiecewiseLinearModel
    auto in = [&](size_t i) { return data[i]; };
    auto out_serial = [&result_segments_serial](const auto &cs) { 
        auto segment = cs.get_Canonicalsegment();
        result_segments_serial.push_back(segment); 
    };
    auto out_parallel = [&result_segments_parallel](const auto &cs) { 
        auto segment = cs.get_Canonicalsegment();
        result_segments_parallel.push_back(segment); 
    };


    // Compare the time taken for the serial and parallel version
    auto start = std::chrono::high_resolution_clock::now();
    size_t num_segments = make_segmentation(data.size(), epsilon, in, out_serial);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> serial_duration = end - start;
    std::cout << "The size of the segements is " << num_segments << std::endl;
    std::cout << "Serial Segmentation Time: " << serial_duration.count() << " seconds" << std::endl;
    std::cout << "The size of the serial segments is " << result_segments_serial.size() << std::endl;

    // for(int i = 0;i<=10;i++){
    //     std::cout << "The slope is " << result_segments_serial[i].slope << " and the intercept is " << result_segments_serial[i].intercept << std::endl;
    // }

    start = std::chrono::high_resolution_clock::now();
    num_segments = make_segmentation_par(data.size(), epsilon, in, out_parallel);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_duration = end - start;
    std::cout << "The size of the segements is " << num_segments << std::endl;
    std::cout << "Parallel Segmentation Time: " << parallel_duration.count() << " seconds" << std::endl;
    std::cout << "The size of the parallel segments is " << result_segments_parallel.size() << std::endl;
    return 0;
}
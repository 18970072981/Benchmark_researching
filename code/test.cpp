#include "greedy.hpp"
#include <iostream>
#include <vector>
#include <deque>
#include <limits>
#include <algorithm>
#include <fstream>
#include <chrono>
#include "linear_model.hpp"
#include "optimalPLR_Foundation.hpp"

#define SAMPLE_SIZE 20000000
#define FILE_NAME "../data/osmdata.txt"
#include <string> 

/**
 * @brief The result of a certain experiment (augmented with the time taken)
 * 
 */
struct Result
{
    int seg_serial,seg_parallel;
    double time_serial,time_parallel;

    // Initialize the result
    Result(int seg_serial,int seg_parallel,double time_serial,double time_parallel):
        seg_serial(seg_serial),seg_parallel(seg_parallel),time_serial(time_serial),time_parallel(time_parallel){}
    Result():seg_serial(0),seg_parallel(0),time_serial(0),time_parallel(0){}
    Result(const Result &other):seg_serial(other.seg_serial),seg_parallel(other.seg_parallel),time_serial(other.time_serial),time_parallel(other.time_parallel){}
};

/**
 * @brief Print out the results of the experiment
 * 
 * 
 * @param results 
 */
void print(std::vector<Result> results){
    std::cout << "#Segments (Serial): "<<std::endl;
    for(int i = 0;i<results.size();i++){
        std::cout << results[i].seg_serial  << ",";
    }
    std::cout << std::endl;
    std::cout << "#Segments (Parallel): ";
    for(int i = 0;i<results.size();i++){
        std::cout << results[i].seg_parallel  << ",";
    }
    std::cout << std::endl;
    std::cout << "Time (Serial): ";
    for(int i = 0;i<results.size();i++){
        std::cout << results[i].time_serial  << ",";
    }
    std::cout << std::endl;
    std::cout << "Time (Parallel): ";
    for(int i = 0;i<results.size();i++){
        std::cout << results[i].time_parallel  << ",";
    }
}

/**
 * @brief printout the #segments and the time taken for the serial and parallel version 
 *  
 */
void printout(std::vector<Result> results_greedy,
                std::vector<Result>results_linear,
                std::vector<Result>results_customized){
    
    // Printout the result of the GreedyPLR
    std::cout << "GreedyPLR :" << std::endl;
    print(results_greedy);
    std::cout << "OptimalPLR (PGM-Index):" << std::endl;
    print(results_linear);
    std::cout << "OptimalPLR (Customized):" << std::endl;
    print(results_customized);    
}


/**
 * @brief load the data from a certain file
 * 
 * @param data 
 */
void load_data( std::vector<double> &data){
    // Load the data from the file
    std::ifstream file(FILE_NAME);
    for(int i = 0;i<SAMPLE_SIZE;i++){
        file >> data[i];
    }
}
/**
 * @brief Check the data for the first 10 elements
 * 
 * @param data 
 */
void checkData(std::vector<double> data){
    for(int i = 0;i<=10;i++){
        std::cout << "Data " << i << " : " << data[i] << std::endl;
    }
}


/**
 * @brief experiment with the GreedyPLR
 * 
 * @param data 
 * @param epsilon 
 * @param threads_num 
 * @return Result of the GreedyPLR
 */
Result experiment_GreedyPLR(std::vector<double> data,double epsilon = 32,int threads_num = 16){

    
    // Intialize the results
    Result result;

    // Initialize the GreedyPiecewiseLinearModel and OptimalPiecewiseLinearModel and OptimalPLR  
    using Segment= greedy::internal::GreedyPiecewiseLinearModel<double,int>::CanonicalSegment::Segment;
    auto result_segments_serial = std::vector<Segment>();
    auto result_segments_parallel = std::vector<Segment>();
    
    // Initialize the epsilon
    // intialize the input and output function
    auto in = [&](size_t i) { return data[i]; };
    auto out_serial = [&result_segments_serial](const auto &cs) { 
        auto segment = cs.get_Canonicalsegment();
        result_segments_serial.push_back(segment); 
    };
    auto out_parallel = [&result_segments_parallel](const auto &cs) { 
        auto segment = cs.get_Canonicalsegment();
        result_segments_parallel.push_back(segment); 
    };

    // Start the experiment of Greedyplr
    // The time taken for the serial version and the number of the segments
    auto start = std::chrono::high_resolution_clock::now();
    size_t num_segments = greedy::internal::make_segmentation(data.size(), epsilon, in, out_serial);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> serial_duration = end - start;
    result.time_serial = serial_duration.count();
    result.seg_serial = num_segments;
    
   
    start = std::chrono::high_resolution_clock::now();
    num_segments = greedy::internal::make_segmentation_par(data.size(), epsilon, in, out_parallel,threads_num);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_duration = end - start;
    // The time taken for the parallel version and the number of the segments
    result.time_parallel = parallel_duration.count();
    result.seg_parallel = num_segments;
    // checkforepi_greplr(data,result_segments_serial);
    return result; 
           
}


/**
 * @brief Experiment with the OptimalPLR (PGM-Index)
 * 
 * @param data 
 * @param epsilon 
 * @param threads_num 
 * @return Result of the OptimalPLR
 *
 * 
 */
Result experiment_linear(std::vector<double> data,double epsilon = 32,int threads_num = 16){
    // Intialize the results
    Result result;

    // Initialize the OptimalPiecewiseLinearModel  
    using Segment = pgm::internal::OptimalPiecewiseLinearModel<double,int>::CanonicalSegment;
   std::vector<std::pair<double, double>> serial_segments;
    std::vector<std::pair<double, double>> parallel_segments;
    
    // Initialize the epsilon
    // intialize the input and output function
    auto in = [&data](size_t i) { return data[i]; };
    auto out_serial = [&serial_segments](const auto &segment) { 
        auto result = segment.get_floating_point_segment(0);
        serial_segments.push_back(result);
    };
    auto out_parallel = [&parallel_segments](const auto &segment) { 
        auto result = segment.get_floating_point_segment(0);
        parallel_segments.push_back(result);
    };

    // Start the experiment of OptimalPLR(PGM-Index)
    // The time taken for the serial version and the number of the segments
    auto start = std::chrono::high_resolution_clock::now();
    size_t num_segments = pgm::internal::make_segmentation(data.size(), epsilon, in, out_serial);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> serial_duration = end - start;
    result.time_serial = serial_duration.count();
    result.seg_serial = num_segments;
    
    // The time taken for the parallel version and the number of the segments
    start = std::chrono::high_resolution_clock::now();
    num_segments = pgm::internal::make_segmentation_par(data.size(), epsilon, in, out_parallel,threads_num);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_duration = end - start;
    // The time taken for the parallel version and the number of the segments
    result.time_parallel = parallel_duration.count();
    result.seg_parallel = num_segments;
    
    return result; 
}

/**
 * @brief Experiment with the OptimalPLR(Customized)
 * 
 * @param data 
 * @param epsilon 
 * @param threads_num 
 * @return Result of the OptimalPLR (Customized)  
 */
Result experiment_OptimalPLR(std::vector<double> data,double epsilon = 32,int threads_num = 16){
    // Intialize the results
    Result result;

    // Initialize the OptimalPLR  
    using Segment = PGM::internal::Segment;
    using Point = PGM::internal::Point;
    using Model = PGM::internal::OptimalPLR;

    // Preprocess the data
    std::vector<Point> points(data.size());
    for(int i = 0;i<data.size();i++){
        points[i].x = data[i];
        points[i].y = i;
    }    

    // Initialize the epsilon
    auto start1 = std::chrono::high_resolution_clock::now();

    Model opt(epsilon);
    std::vector<Segment> out_segments = opt.segmentData(points);
    // std::cout << "The size of the segements is " << out_segments.size() << std::endl;
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = end1 - start1;
    result.time_serial =duration1.count();
    result.seg_serial = out_segments.size();
    // The time taken for the parallel version
    out_segments.clear();
    auto start = std::chrono::high_resolution_clock::now();

    size_t num_segments = make_segmentation_par(points.size(), epsilon, points, out_segments,threads_num);
    // std::cout << "The size of the segements is " << num_segments << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    result.time_parallel = duration.count();
    result.seg_parallel = num_segments;


    return result; 
}


/**
 * 
 * @brief Experiment with the greedyPLR and OptimalPLR (PGM-Index & Customized) 
 * 
 * @param data 
 * @param epsilon 
 * @param threads_num 
 */
template <typename X, typename Y>
void experiment(std::vector<double> data,int start, int end,int threads_num = 16)
{
    // Experiment with the greedyPLR and OptimalPLR (PGM-Index & Customized)
    // Firstly, we need to check the data
    checkData(data);
       
    // Initialize the results
    std::vector<Result> results_greedy;
    std::vector<Result> results_linear;
    std::vector<Result> results_optimal_customized;

    // Run the experiment for the greedyPLR
    for(int i = start;i<end;i++){
        double epsilon = 1<<i;
        printf("The epsilon is %f\n",epsilon);
        
        // Run the experiment for the greedyPLR
        Result result_greedy = experiment_GreedyPLR(data,epsilon,threads_num);
        results_greedy.push_back(result_greedy); 
       
        // Run the experiment for the OptimalPLR (Customized)
        Result result_optimal = experiment_OptimalPLR(data,epsilon,threads_num);
        results_optimal_customized.push_back(result_optimal);
        // Run the experiment for the OptimalPLR (PGM-Index)
        Result result_linear = experiment_linear(data,epsilon,threads_num);
        results_linear.push_back(result_linear);
    
    }

    // Print the results
    printout(results_greedy,results_linear,results_optimal_customized);

}

int main(){
    
    // Read the data from data1.txt
    std::vector<double> data(SAMPLE_SIZE);
    load_data(data);


    // Experiment with the greedyPLR and OptimalPLR (PGM-Index & Customized)
    // experiment(data,start,end,threads_num);
    experiment<double, int>(data,1,13);
    // experiment(data,8,13,1);
    // experiment(data,8,13,4); 
    return 0;
}


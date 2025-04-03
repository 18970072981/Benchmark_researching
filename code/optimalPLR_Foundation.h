/* Write your docstring here*/
/**
 * @file OptimalPLR.h
 * @author Jiayong Qin (2057729401@qq.com)
 * @brief implement optimalPLR
 * @version 1.0
 * @date 2025-2-28
 *
 * @copyright Copyright (c) 2025
 *
 */
#include <iostream>
#include <vector>
#include <deque>
#include <limits>
#include <algorithm>
#include <fstream>
#include <omp.h>
#include <thread>

#include <chrono>

typedef int K;

struct Point
{
    /**
     * @brief Data used to represent a point in 2D space
     *
     */
    double x, y;
};

struct Segment
{
    /**
     * @brief Segements
     *
     */
    double slope, intercept;
    int start_idx, end_idx;
    // [start_idx, end_idx]
};

class OptimalPLR
{

public:
    /**
     * @brief Construct a new OptimalPLR object
     *
     * @param error error bounds of this implementation
     */
    OptimalPLR(double error) : epsilon(error) {}

    std::vector<Segment> segmentData(const std::vector<Point> &points)
    {
        std::vector<Segment> segments;
        if (points.empty())
            return segments;

        // Initialization: the first two points acting as the initial hull and segement
        int start = 0;
        Point sa = {points[0].x, points[0].y + epsilon};
        Point sc = {points[1].x, points[1].y - epsilon};
        Point sb = {points[0].x, points[0].y - epsilon};
        Point sd = {points[1].x, points[1].y + epsilon};

        double rho_max = computeSlope(sb, sd);
        double rho_min = computeSlope(sa, sc);

        std::deque<Point> upperHull = {sa, sd};
        std::deque<Point> lowerHull = {sb, sc};

        for (size_t i = 2; i < points.size(); ++i)
        {
            Point s = {points[i].x, points[i].y};

            if (!isWithinBounds(s, rho_max, sc, sd, rho_min))
            {
                // Check if this point is within the bounds , if not , we need to update the segement and initialize again
                Point so = findIntersection(sa, sc, sb, sd);
                double finalSlope = (rho_min + rho_max) / 2;
                segments.push_back({finalSlope, computeIntercept(so, finalSlope), start, int(i)});

                start = i;
                sa = {points[start].x, points[start].y + epsilon};
                sc = {points[i + 1].x, points[i + 1].y - epsilon};
                sb = {points[start].x, points[start].y - epsilon};
                sd = {points[i + 1].x, points[i + 1].y + epsilon};

                rho_max = computeSlope(sb, sd);
                rho_min = computeSlope(sa, sc);
                upperHull.clear();
                lowerHull.clear();
                upperHull = {sa, sd};
                lowerHull = {sb, sc};
                i++;

                continue;
            }
            else
            {
                Point s_upper;
                Point s_lower;
                bool upper = 0;
                bool lower = 0;
                

                if (s.y - epsilon - sa.y> rho_min * (s.x - sa.x) )
                {
                    // Find the max slope point and delete the points before it so that we can update the min slope
                    lower = 1;
                    s_lower.x = s.x;
                    s_lower.y = s.y - epsilon;
                    sa = findMaxSlopePoint(upperHull, s_lower);
                    sc = s_lower;
                    deletePointsBefore(upperHull, sa);
                    rho_min = computeSlope(sa, sc);
                    
                    
                }

                if (s.y + epsilon - sb.y < rho_max * (s.x - sb.x))
                {
                    // Find the min slope point and delete the points before it so that we can update the segement and the hull
                    upper = 1;
                    s_upper.x = s.x;
                    s_upper.y = s.y + epsilon;
                    sb = findMinSlopePoint(lowerHull, s_upper);
                    sd = s_upper;
                    deletePointsBefore(lowerHull, sb);
                    rho_max = computeSlope(sb, sd);
                    
                   
                }
                if(upper){
                    updateConvexHull(upperHull, s_upper, true);
                    sd = upperHull.back();
                }
                if(lower){
                    updateConvexHull(lowerHull, s_lower, false);
                    sc = lowerHull.back();
                }
            }
            

        }

        Point so = findIntersection(sa, sc, sb, sd);
        double finalSlope = (rho_min + rho_max) / 2;
        segments.push_back({finalSlope, computeIntercept(so, finalSlope), start, int(points.size() - 1)});

        return segments;
    }

private:
    double epsilon;

    /**
     * @brief compute the slope between two points
     *
     * @param p1
     * @param p2
     * @return double
     */
    double computeSlope(const Point &p1, const Point &p2)
    {
        if (p2.x == p1.x)
            return std::numeric_limits<double>::infinity();
        return (p2.y - p1.y) / (p2.x - p1.x);
    }

    /**
     * @brief compute the intercept of a line
     *
     * @param p
     * @param slope
     * @return double
     */
    double computeIntercept(const Point &p, double slope)
    {
        return p.y - slope * p.x;
    }

    /**
     * @brief To Check if the point is within the bounds
     *
     * @param s
     * @param rho_max
     * @param sa
     * @param sc
     * @param rho_min
     * @param sb
     * @param sd
     * @return true
     * @return false
     */
    bool isWithinBounds(const Point &s, double rho_max, Point sc, Point sd, double rho_min)
    {
        return (s.y + epsilon >= rho_min * (s.x - sc.x) + sc.y) && (s.y - epsilon <= rho_max * (s.x - sd.x) + sd.y);
    }

    /**
     * @brief update the convex hull
     *
     * @param hull
     * @param s
     * @param upper if the hull is upper or not
     */
    void updateConvexHull(std::deque<Point> &hull, const Point &s, bool upper)
    {
        while (hull.size() >= 2)
        {
            auto &p1 = hull[hull.size() - 2];
            auto &p2 = hull[hull.size() - 1];
            if (isRedundant(p1, p2, s, upper))
            {
                hull.pop_back();
            }
            else
            {
                break;
            }
        }
        hull.push_back(s);
    }

    /**
     * @brief delete the redundant points using triangle area method
     *
     * @param p1
     * @param p2
     * @param p3
     * @param upper
     * @return true
     * @return false
     */
    bool isRedundant(const Point &p1, const Point &p2, const Point &p3, bool upper)
    {
        double slope1 = (p2.y - p1.y) / (p2.x - p1.x);
        double slope2 = (p3.y - p2.y) / (p3.x - p2.x);
        return upper ? (slope1 >= slope2) : (slope1 <= slope2);
    }

    /**
     * @brief Find the min slope point
     *
     * @param hull
     * @param s
     * @return Point
     */
    Point findMinSlopePoint(std::deque<Point> &hull, const Point &s)
    {

        return *std::min_element(hull.begin(), hull.end(), [&](const Point &a, const Point &b)
                                 { return computeSlope(a, s) < computeSlope(b, s); });
    }

    /**
     * @brief find the max slope point
     *
     * @param hull
     * @param s
     * @return Point
     */
    Point findMaxSlopePoint(std::deque<Point> &hull, const Point &s)
    {
        return *std::max_element(hull.begin(), hull.end(), [&](const Point &a, const Point &b)
                                 { return computeSlope(a, s) < computeSlope(b, s); });
    }

    /**
     * @brief delete point before the reference point
     *
     * @param hull
     * @param reference
     */
    void deletePointsBefore(std::deque<Point> &hull, const Point &reference)
    {
        while (!hull.empty() && hull.front().x < reference.x)
        {
            hull.pop_front();
        }
    }

    /**
     * @brief Find the intersection point of two lines
     *
     * @param sa
     * @param sc
     * @param sb
     * @param sd
     * @return Point
     */
    Point findIntersection(Point sa, Point sc, Point sb, Point sd)
    {
        double a1 = computeSlope(sa, sc);
        double b1 = computeIntercept(sa, a1);

        double a2 = computeSlope(sb, sd);
        double b2 = computeIntercept(sb, a2);

        if (a1 == a2) // Parallel lines
            return {0, 0}; // or handle as needed
        double x_intersect = (b2 - b1) / (a1 - a2);
        double y_intersect = a1 * x_intersect + b1;

        return {x_intersect, y_intersect};
    }
};



template <typename T>
T calculate_residual(const T &a, const T &b)
{
    return a > b ? a - b : b - a;
}

/**
 * @brief Make segmentation of the data parallelly
 *
 * @tparam Fin the input data type
 * @tparam Fout the output data type
 * @param n the size of the data
 * @param epsilon error boundary
 * @param in input data
 * @param out output data
 * @return size_t the number of the segments
 */
template <typename Fin, typename Fout>
size_t make_segmentation_par(size_t n, size_t epsilon, Fin in, Fout &out)
{
    auto parallelism = std::min(std::min(omp_get_num_procs(), omp_get_max_threads()), 20);
    std::cout << parallelism << std::endl;
    auto chunk_size = n / parallelism;
    auto c = 0;
    OptimalPLR opt(epsilon);

#pragma omp parallel
    {
        std::vector<Segment> local_segments;
#pragma omp for reduction(+ : c) schedule(dynamic)
        for (int i = 0; i < parallelism; ++i)
        {
            size_t start = i * chunk_size;
            size_t end = (i == parallelism - 1) ? n : (i + 1) * chunk_size;
            std::vector<Point> segment_points(in.begin() + start, in.begin() + end);
            auto segments = opt.segmentData(segment_points);
            c += segments.size();
            local_segments.insert(local_segments.end(), segments.begin(), segments.end());
        }
#pragma omp critical
        out.insert(out.end(), local_segments.begin(), local_segments.end());
    }

    return c;
}

/**
 * @brief check the correctness for each segments
 */
void checkForEpsilon(const std::vector<Point> &points, const std::vector<Segment> &segments)
{
    int num = 0;
    for (int i = 0; i <= 10; i++)
    {

        int start = segments[i].start_idx;
        int end = segments[i].end_idx;
        double rho = segments[i].slope; // Ensure slope is double
        double b = segments[i].intercept;

        // Check index bounds
        if (start < 0 || end >= points.size() || start > end)
        {
            std::cerr << "Error: Segment " << i << " has invalid indices (" << start << ", " << end << ")\n";
            continue; // Skip invalid segments
        }

        double max_residual = std::numeric_limits<double>::lowest();
        for (size_t j = start; j < end; j++)
        {
            double residual = std::abs(points[j].y - (rho * points[j].x + b));
            max_residual = std::max(max_residual, residual);
        }

        std::cout << "Maximum residual for Segment " << i+1 << " is: " << max_residual << std::endl;
    }
}


#include <iostream>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Random_points_in_square_2<Point_2> Random_points_iterator;
typedef CGAL::Search_traits_2<K> Traits;
typedef CGAL::Kd_tree<Traits> Tree;

int main() {
    // Generate 2000 random points
    std::vector<Point_2> points;
    points.reserve(2000);
    Random_points_iterator rpg(10.0);
    for (int i = 0; i < 2000; ++i) {
        points.push_back(*rpg++);
    }

    // Construct KD tree
    Tree tree(points.begin(), points.end());

    // Define the query point A
    double query_x, query_y;
    std::cout << "Enter the x-coordinate of point A: ";
    std::cin >> query_x;
    std::cout << "Enter the y-coordinate of point A: ";
    std::cin >> query_y;
    Point_2 query_point(query_x, query_y);

    // Search for the nearest neighbor
    auto search_result = tree.search(query_point,1);

    // Get the nearest neighbor
    /*if (!search_result.empty()) {
        Point_2 nearest_neighbor = *search_result.begin();
        std::cout << "Nearest neighbor coordinates: (" << nearest_neighbor.x() << ", " << nearest_neighbor.y() << ")\n";
    }
    else {
        std::cout << "No nearest neighbor found.\n";
    }*/

    return 0;
}

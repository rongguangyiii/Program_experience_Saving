#include "flowSolve/include/globalData.h"
#include <iostream>
#include <filesystem>  


void GlobalData::createFolder(const std::string& basePath)
{
	std::filesystem::path folderPath = basePath;
	if (!std::filesystem::exists(folderPath))
		std::filesystem::create_directories(folderPath);
}

bool GlobalData::isPointInPolygon_moreInner(const CoordPtr& p, const std::vector<CoordPtr>& polygon) {
    const double EPSILON = 1e-10;
    bool inside = false;
    size_t j = polygon.size() - 1;

    for (size_t i = 0; i < polygon.size(); i++) {
        // Check if point is on vertex
        if ((std::fabs(p->x - polygon[i]->x) < EPSILON) && (std::fabs(p->y - polygon[i]->y) < EPSILON)) {
            return true; // Point is on a vertex
        }

        // Check if point is on edge
        double minX = std::min(polygon[i]->x, polygon[j]->x);
        double maxX = std::max(polygon[i]->x, polygon[j]->x);
        double minY = std::min(polygon[i]->y, polygon[j]->y);
        double maxY = std::max(polygon[i]->y, polygon[j]->y);

        if ((minX <= p->x && p->x <= maxX) && (minY <= p->y && p->y <= maxY)) {
            double dx = polygon[j]->x - polygon[i]->x;
            double dy = polygon[j]->y - polygon[i]->y;
            if (std::fabs(dy) < EPSILON) { // horizontal line
                if (std::fabs(p->y - polygon[i]->y) < EPSILON) {
                    return true; // Point is on horizontal edge
                }
            }
            else if (std::fabs(dx) < EPSILON) { // vertical line
                if (std::fabs(p->x - polygon[i]->x) < EPSILON) {
                    return true; // Point is on vertical edge
                }
            }
            else {
                double slope = dy / dx;
                double intercept = polygon[i]->y - slope * polygon[i]->x;
                if (std::fabs(p->y - (slope * p->x + intercept)) < EPSILON) {
                    return true; // Point is on non-vertical edge
                }
            }
        }

        // Check if point is inside using ray-casting algorithm
        if (((polygon[i]->y <= p->y) && (p->y < polygon[j]->y)) ||
            ((polygon[j]->y <= p->y) && (p->y < polygon[i]->y))) {
            double intersectionX = polygon[i]->x + (p->y - polygon[i]->y) * (polygon[j]->x - polygon[i]->x) / (polygon[j]->y - polygon[i]->y);
            if (intersectionX < p->x) {
                inside = !inside;
            }
        }

        j = i;
    }
    return inside;
}

// 计算质心
Coord GlobalData::calculateCentroid(const std::vector<CoordPtr>& points) {
    double sumX = 0, sumY = 0;
    for (const auto& point : points) {
        sumX += point->x;
        sumY += point->y;
    }
    size_t n = points.size();
    return Coord(PointId(0,0), sumX / n, sumY / n);
}

// 计算给定点相对于质心的角度（0 到 2π）
double GlobalData::calAngle_02pi(const Coord& centroid, const CoordPtr& point) {
    double angle = std::atan2(point->y - centroid.y, point->x - centroid.x);
    // 将负角度转换到 0 到 2π 范围内
    return (angle >= 0) ? angle : (2 * GlobalData::PI + angle);
}

// 计算给定点相对于质心的角度（-π 到 π）
double GlobalData::calAngle_pipi(const Coord& centroid, const CoordPtr& point) {
    return std::atan2(point->y - centroid.y, point->x - centroid.x);
}

// 按逆时针方向（从0到2π）排序
void GlobalData::sortPointsCounterclockwise(std::vector<CoordPtr>& points) {
    // 计算质心
    Coord centroid = calculateCentroid(points);

    // 使用质心来计算每个点的角度，并根据这些角度对点进行排序
    std::sort(points.begin(), points.end(), [&centroid](const CoordPtr& a, const CoordPtr& b) {
        double angleA = calAngle_02pi(centroid, a);
        double angleB = calAngle_02pi(centroid, b);
        return angleA < angleB;  // 按从小到大的顺序排序
        });
}

// 根据给定角度范围筛选点
std::vector<CoordPtr> GlobalData::filterPointsByAngleRange(const std::vector<CoordPtr>& points, double angleStart, double angleEnd) {
    Coord centroid = calculateCentroid(points);
    std::vector<CoordPtr> filteredPoints;
    const double epsilon = 1e-6;

    for (const auto& point : points) {
        double angle = calAngle_pipi(centroid, point);
        // 将角度调整到 0 到 2π 范围内
        if (angle < 0) {
            angle += 2 * PI;
        }

        // 将输入的角度范围也调整到 0 到 2π 范围内
        if (angleStart < 0) {
            angleStart += 2 * PI;
        }
        if (angleEnd < 0) {
            angleEnd += 2 * PI;
        }

        // 判断点的角度是否在指定范围内，使用 EPSILON 进行浮点数比较
        if ((angle > angleStart || std::fabs(angle - angleStart) < epsilon) &&
            (angle < angleEnd || std::fabs(angle - angleEnd) < epsilon)) {
            filteredPoints.push_back(point);
        }
    }

    return filteredPoints;
}

// 计算给定角度的差异，结果在 [0, π] 范围内
double GlobalData::angleDifference(double angle1, double angle2) {
    double diff = std::fabs(angle1 - angle2);
    return std::fmin(diff, 2 * PI - diff);  // 确保差值在 [0, π] 范围内
}

// 根据给定角度找到最近的点
CoordPtr& GlobalData::PointAtAngle(std::vector<CoordPtr>& points, double targetAngle) {
    Coord centroid = calculateCentroid(points);

    size_t closestIndex = 0;
    double minAngleDifference = 1e6;  // 初始设定为一个很大的值

    for (size_t i = 0; i < points.size(); ++i) {
        double angle = calAngle_pipi(centroid, points[i]);
        // 将角度调整到 0 到 2π 范围内
        if (angle < 0) {
            angle += 2 * PI;
        }

        // 计算当前点与目标角度的差异
        double currentDifference = angleDifference(angle, targetAngle);

        // 更新最小差异和最近点
        if (currentDifference < minAngleDifference) {
            minAngleDifference = currentDifference;
            closestIndex = i;
        }
    }

    return points[closestIndex];  // 返回对points中最近点的引用
}

/*-------------------------------------------------------------------
*  Function: isPointInPolygon_moreInner
*  Purpose: 判断点是否在多边形内部版本1：包括与边界重合或顶点重合的点都认为在内部
*  Arguments:
*    p - 目标判断点
*    polygon - 有序封闭多边形
*  Returns:
*    boll
-------------------------------------------------------------------*/
bool Tools::isPointInPolygon_moreInner(const Coord& p, const std::vector<Coord>& polygon) {
    const double EPSILON = 1e-10;
    bool inside = false;
    size_t j = polygon.size() - 1;

    for (size_t i = 0; i < polygon.size(); i++) {
        // Check if point is on vertex
        if ((std::fabs(p.x - polygon[i].x) < EPSILON) && (std::fabs(p.y - polygon[i].y) < EPSILON)) {
            return true; // Point is on a vertex
        }

        // Check if point is on edge
        double minX = std::min(polygon[i].x, polygon[j].x);
        double maxX = std::max(polygon[i].x, polygon[j].x);
        double minY = std::min(polygon[i].y, polygon[j].y);
        double maxY = std::max(polygon[i].y, polygon[j].y);

        if ((minX <= p.x && p.x <= maxX) && (minY <= p.y && p.y <= maxY)) {
            double dx = polygon[j].x - polygon[i].x;
            double dy = polygon[j].y - polygon[i].y;
            if (std::fabs(dy) < EPSILON) { // horizontal line
                if (std::fabs(p.y - polygon[i].y) < EPSILON) {
                    return true; // Point is on horizontal edge
                }
            }
            else if (std::fabs(dx) < EPSILON) { // vertical line
                if (std::fabs(p.x - polygon[i].x) < EPSILON) {
                    return true; // Point is on vertical edge
                }
            }
            else {
                double slope = dy / dx;
                double intercept = polygon[i].y - slope * polygon[i].x;
                if (std::fabs(p.y - (slope * p.x + intercept)) < EPSILON) {
                    return true; // Point is on non-vertical edge
                }
            }
        }

        // Check if point is inside using ray-casting algorithm
        if (((polygon[i].y <= p.y) && (p.y < polygon[j].y)) ||
            ((polygon[j].y <= p.y) && (p.y < polygon[i].y))) {
            double intersectionX = polygon[i].x + (p.y - polygon[i].y) * (polygon[j].x - polygon[i].x) / (polygon[j].y - polygon[i].y);
            if (intersectionX < p.x) {
                inside = !inside;
            }
        }

        j = i;
    }
    return inside;
}

/*-------------------------------------------------------------------
*  Function: isPointInPolygon_moreOuter
*  Purpose: 判断点是否在多边形内部版本2：包括与边界重合或顶点重合的点都认为在外部
*  Arguments:
*    p - 目标判断点
*    polygon - 有序封闭多边形
*  Returns:
*    boll
-------------------------------------------------------------------*/
bool Tools::isPointInPolygon_moreOuter(const Coord& p, const std::vector<Coord>& polygon) {
    const double EPSILON = 1e-10;
    bool inside = false;
    size_t j = polygon.size() - 1;

    for (size_t i = 0; i < polygon.size(); i++) {
        // Check if point is on vertex
        if ((std::fabs(p.x - polygon[i].x) < EPSILON) && (std::fabs(p.y - polygon[i].y) < EPSILON)) {
            return false; // Point is on a vertex
        }

        // Check if point is on edge
        double minX = std::min(polygon[i].x, polygon[j].x);
        double maxX = std::max(polygon[i].x, polygon[j].x);
        double minY = std::min(polygon[i].y, polygon[j].y);
        double maxY = std::max(polygon[i].y, polygon[j].y);

        if ((minX <= p.x && p.x <= maxX) && (minY <= p.y && p.y <= maxY)) {
            double dx = polygon[j].x - polygon[i].x;
            double dy = polygon[j].y - polygon[i].y;
            if (std::fabs(dy) < EPSILON) { // horizontal line
                if (std::fabs(p.y - polygon[i].y) < EPSILON) {
                    return false; // Point is on horizontal edge
                }
            }
            else if (std::fabs(dx) < EPSILON) { // vertical line
                if (std::fabs(p.x - polygon[i].x) < EPSILON) {
                    return false; // Point is on vertical edge
                }
            }
            else {
                double slope = dy / dx;
                double intercept = polygon[i].y - slope * polygon[i].x;
                if (std::fabs(p.y - (slope * p.x + intercept)) < EPSILON) {
                    return false; // Point is on non-vertical edge
                }
            }
        }

        // Check if point is inside using ray-casting algorithm
        if (((polygon[i].y <= p.y) && (p.y < polygon[j].y)) ||
            ((polygon[j].y <= p.y) && (p.y < polygon[i].y))) {
            double intersectionX = polygon[i].x + (p.y - polygon[i].y) * (polygon[j].x - polygon[i].x) / (polygon[j].y - polygon[i].y);
            if (intersectionX < p.x) {
                inside = !inside;
            }
        }

        j = i;
    }
    return inside;
}

/*-------------------------------------------------------------------
*  Function: areCollinear(Coordinate p1, Coordinate p2, Coordinate p3)
*  Purpose: 检查三角形的三个点是否共线
*  Arguments:
*   Coordinate p1,p2,p3 - 三角形的三个顶点
*  Returns:
*    bool - 是否共线的布尔值
-------------------------------------------------------------------*/
bool Tools::isCollinear(const Coord& p1, const Coord& p2, const Coord& p3)
{
    // 使用斜率法判断三点是否共线  
    // 如果斜率相同或任意两点重合（在容差范围内），则它们共线  
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    double dx1 = p2.x - p1.x;
    double dy1 = p2.y - p1.y;
    double dx2 = p3.x - p1.x;
    double dy2 = p3.y - p1.y;

    // 使用斜率来判断是否共线  
    if (std::abs(dx1 * dy2 - dx2 * dy1) < epsilon) {
        return true;
    }
    // 检查是否接近0（即检查是否有重合的点）  
    if (std::abs(dx1) < epsilon && std::abs(dx2) < epsilon) {
        // 如果x坐标都接近相同，则检查y坐标是否也接近相同  
        return std::abs(dy1) < epsilon && std::abs(dy2) < epsilon;
    }

    return false;
}

/*-------------------------------------------------------------------
*  Function: reorderPointsIndex()
*  Purpose: 根据点的坐标，按逆时针排序点的索引
*  Description: 工具函数，用于重新按逆时针排序点的索引，
*  Arguments:
*    PointLists - 点坐标和索引的pair Vector
*    modelCenter - 多边形的中心点
*  Returns:
*    null
-------------------------------------------------------------------*/
void Tools::reorderPointsIndex(std::vector<std::pair<Coord, size_t>>& PointLists, const Coord& modelCenter)
{
    std::size_t wallNodeNum = PointLists.size();
    std::vector<std::pair<double, size_t>> angleAndIndex;
    angleAndIndex.resize(wallNodeNum);
    for (size_t iNode = 0; iNode < wallNodeNum; ++iNode)
    {
        auto& perNode = PointLists.at(iNode).first;
        //Coord vec = perNode - modelCenter;
        Coord vec{ (perNode.x - modelCenter.x), (perNode.y - modelCenter.y), (perNode.z - modelCenter.z) };
        double angle = atan2(vec.y, vec.x);
        angleAndIndex[iNode] = std::make_pair(angle, iNode);
    }
    std::sort(angleAndIndex.begin(), angleAndIndex.end(),
        [](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b) {return a.first < b.first; });

    std::vector<size_t> reorderedIndices;
    std::vector<Coord> NodeVec_temp;
    std::vector<std::pair<Coord, size_t>> NodePairVec_temp2;
    NodeVec_temp.resize(wallNodeNum);
    reorderedIndices.resize(wallNodeNum);
    for (size_t i = 0; i < wallNodeNum; ++i)
    {
        NodeVec_temp[i] = PointLists[angleAndIndex[i].second].first;
        reorderedIndices[i] = PointLists[angleAndIndex[i].second].second;
        NodePairVec_temp2.push_back(std::make_pair(NodeVec_temp[i], reorderedIndices[i]));
    }
    PointLists = NodePairVec_temp2;//会修改原本点的顺序;
}

double Tools::Distance(const CoordPtr& point1, const CoordPtr& point2) {
    // 使用距离公式 d = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    double dx = point2->x - point1->x;
    double dy = point2->y - point1->y;
    return std::sqrt(dx * dx + dy * dy);
}

std::vector<CoordPtr> Tools::cood2CoodPtr(const std::vector<Coord>& points)
{
    std::vector<CoordPtr> pointsPtr(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        pointsPtr[i] = std::make_shared<Coord>(points[i]);
    }
    return pointsPtr;
}

std::vector<Coord> Tools::coodPtr2Cood(const std::vector<CoordPtr>& pointsPtr)
{
    std::vector<Coord> points(pointsPtr.size());
    for (size_t i = 0; i < pointsPtr.size(); ++i) {
        points[i] = *pointsPtr[i];
    }
    return points;
}
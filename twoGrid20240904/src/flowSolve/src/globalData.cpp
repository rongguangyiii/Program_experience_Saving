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

// ��������
Coord GlobalData::calculateCentroid(const std::vector<CoordPtr>& points) {
    double sumX = 0, sumY = 0;
    for (const auto& point : points) {
        sumX += point->x;
        sumY += point->y;
    }
    size_t n = points.size();
    return Coord(PointId(0,0), sumX / n, sumY / n);
}

// �����������������ĵĽǶȣ�0 �� 2�У�
double GlobalData::calAngle_02pi(const Coord& centroid, const CoordPtr& point) {
    double angle = std::atan2(point->y - centroid.y, point->x - centroid.x);
    // �����Ƕ�ת���� 0 �� 2�� ��Χ��
    return (angle >= 0) ? angle : (2 * GlobalData::PI + angle);
}

// �����������������ĵĽǶȣ�-�� �� �У�
double GlobalData::calAngle_pipi(const Coord& centroid, const CoordPtr& point) {
    return std::atan2(point->y - centroid.y, point->x - centroid.x);
}

// ����ʱ�뷽�򣨴�0��2�У�����
void GlobalData::sortPointsCounterclockwise(std::vector<CoordPtr>& points) {
    // ��������
    Coord centroid = calculateCentroid(points);

    // ʹ������������ÿ����ĽǶȣ���������Щ�ǶȶԵ��������
    std::sort(points.begin(), points.end(), [&centroid](const CoordPtr& a, const CoordPtr& b) {
        double angleA = calAngle_02pi(centroid, a);
        double angleB = calAngle_02pi(centroid, b);
        return angleA < angleB;  // ����С�����˳������
        });
}

// ���ݸ����Ƕȷ�Χɸѡ��
std::vector<CoordPtr> GlobalData::filterPointsByAngleRange(const std::vector<CoordPtr>& points, double angleStart, double angleEnd) {
    Coord centroid = calculateCentroid(points);
    std::vector<CoordPtr> filteredPoints;
    const double epsilon = 1e-6;

    for (const auto& point : points) {
        double angle = calAngle_pipi(centroid, point);
        // ���Ƕȵ����� 0 �� 2�� ��Χ��
        if (angle < 0) {
            angle += 2 * PI;
        }

        // ������ĽǶȷ�ΧҲ������ 0 �� 2�� ��Χ��
        if (angleStart < 0) {
            angleStart += 2 * PI;
        }
        if (angleEnd < 0) {
            angleEnd += 2 * PI;
        }

        // �жϵ�ĽǶ��Ƿ���ָ����Χ�ڣ�ʹ�� EPSILON ���и������Ƚ�
        if ((angle > angleStart || std::fabs(angle - angleStart) < epsilon) &&
            (angle < angleEnd || std::fabs(angle - angleEnd) < epsilon)) {
            filteredPoints.push_back(point);
        }
    }

    return filteredPoints;
}

// ��������ǶȵĲ��죬����� [0, ��] ��Χ��
double GlobalData::angleDifference(double angle1, double angle2) {
    double diff = std::fabs(angle1 - angle2);
    return std::fmin(diff, 2 * PI - diff);  // ȷ����ֵ�� [0, ��] ��Χ��
}

// ���ݸ����Ƕ��ҵ�����ĵ�
CoordPtr& GlobalData::PointAtAngle(std::vector<CoordPtr>& points, double targetAngle) {
    Coord centroid = calculateCentroid(points);

    size_t closestIndex = 0;
    double minAngleDifference = 1e6;  // ��ʼ�趨Ϊһ���ܴ��ֵ

    for (size_t i = 0; i < points.size(); ++i) {
        double angle = calAngle_pipi(centroid, points[i]);
        // ���Ƕȵ����� 0 �� 2�� ��Χ��
        if (angle < 0) {
            angle += 2 * PI;
        }

        // ���㵱ǰ����Ŀ��ǶȵĲ���
        double currentDifference = angleDifference(angle, targetAngle);

        // ������С����������
        if (currentDifference < minAngleDifference) {
            minAngleDifference = currentDifference;
            closestIndex = i;
        }
    }

    return points[closestIndex];  // ���ض�points������������
}

/*-------------------------------------------------------------------
*  Function: isPointInPolygon_moreInner
*  Purpose: �жϵ��Ƿ��ڶ�����ڲ��汾1��������߽��غϻ򶥵��غϵĵ㶼��Ϊ���ڲ�
*  Arguments:
*    p - Ŀ���жϵ�
*    polygon - �����ն����
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
*  Purpose: �жϵ��Ƿ��ڶ�����ڲ��汾2��������߽��غϻ򶥵��غϵĵ㶼��Ϊ���ⲿ
*  Arguments:
*    p - Ŀ���жϵ�
*    polygon - �����ն����
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
*  Purpose: ��������ε��������Ƿ���
*  Arguments:
*   Coordinate p1,p2,p3 - �����ε���������
*  Returns:
*    bool - �Ƿ��ߵĲ���ֵ
-------------------------------------------------------------------*/
bool Tools::isCollinear(const Coord& p1, const Coord& p2, const Coord& p3)
{
    // ʹ��б�ʷ��ж������Ƿ���  
    // ���б����ͬ�����������غϣ����ݲΧ�ڣ��������ǹ���  
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    double dx1 = p2.x - p1.x;
    double dy1 = p2.y - p1.y;
    double dx2 = p3.x - p1.x;
    double dy2 = p3.y - p1.y;

    // ʹ��б�����ж��Ƿ���  
    if (std::abs(dx1 * dy2 - dx2 * dy1) < epsilon) {
        return true;
    }
    // ����Ƿ�ӽ�0��������Ƿ����غϵĵ㣩  
    if (std::abs(dx1) < epsilon && std::abs(dx2) < epsilon) {
        // ���x���궼�ӽ���ͬ������y�����Ƿ�Ҳ�ӽ���ͬ  
        return std::abs(dy1) < epsilon && std::abs(dy2) < epsilon;
    }

    return false;
}

/*-------------------------------------------------------------------
*  Function: reorderPointsIndex()
*  Purpose: ���ݵ�����꣬����ʱ������������
*  Description: ���ߺ������������°���ʱ��������������
*  Arguments:
*    PointLists - �������������pair Vector
*    modelCenter - ����ε����ĵ�
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
    PointLists = NodePairVec_temp2;//���޸�ԭ�����˳��;
}

double Tools::Distance(const CoordPtr& point1, const CoordPtr& point2) {
    // ʹ�þ��빫ʽ d = sqrt((x2 - x1)^2 + (y2 - y1)^2)
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
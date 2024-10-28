#pragma once
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include "grid/include/Coord.h"

struct LSGradCoord {
    double x, y, value;
    LSGradCoord(double x_, double y_, double value_) : x(x_), y(y_), value(value_) {}
};

class LSGrad {
public:
    // 构造函数，传入邻居节点的坐标和当前节点的坐标
    LSGrad(const std::vector<LSGradCoord>& neighbors, const LSGradCoord& current)
        : neighbors_(neighbors), current_(current) {}

    // 计算梯度
    void computeGradient(double& gx, double& gy) {
        Eigen::Matrix2d A;
        Eigen::Vector2d b;
        A.setZero();
        b.setZero();

        // 遍历所有邻居节点
        for (const auto& neighbor : neighbors_) {
            double deltaX = neighbor.x - current_.x;
            double deltaY = neighbor.y - current_.y;
            double deltaVar = neighbor.value - current_.value;

            // 距离的权重（这里的 p 设置为 1）
            double distance = std::sqrt(deltaX * deltaX + deltaY * deltaY);
            double weight = 1.0 / std::pow(distance, 1);

            // 填充矩阵 A 和向量 b
            A(0, 0) += weight * deltaX * deltaX;
            A(0, 1) += weight * deltaX * deltaY;
            A(1, 0) += weight * deltaY * deltaX;
            A(1, 1) += weight * deltaY * deltaY;

            // 修正后的 b 向量
            b(0) += weight * deltaX * deltaVar;
            b(1) += weight * deltaY * deltaVar;
        }

        // 计算 A 的逆矩阵并求解梯度
        Eigen::Vector2d grad = A.inverse() * b;
        gx = grad(0);
        gy = grad(1);
    }

private:
    std::vector<LSGradCoord> neighbors_;
    LSGradCoord current_;
};

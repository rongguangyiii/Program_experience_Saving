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
    // ���캯���������ھӽڵ������͵�ǰ�ڵ������
    LSGrad(const std::vector<LSGradCoord>& neighbors, const LSGradCoord& current)
        : neighbors_(neighbors), current_(current) {}

    // �����ݶ�
    void computeGradient(double& gx, double& gy) {
        Eigen::Matrix2d A;
        Eigen::Vector2d b;
        A.setZero();
        b.setZero();

        // ���������ھӽڵ�
        for (const auto& neighbor : neighbors_) {
            double deltaX = neighbor.x - current_.x;
            double deltaY = neighbor.y - current_.y;
            double deltaVar = neighbor.value - current_.value;

            // �����Ȩ�أ������ p ����Ϊ 1��
            double distance = std::sqrt(deltaX * deltaX + deltaY * deltaY);
            double weight = 1.0 / std::pow(distance, 1);

            // ������ A ������ b
            A(0, 0) += weight * deltaX * deltaX;
            A(0, 1) += weight * deltaX * deltaY;
            A(1, 0) += weight * deltaY * deltaX;
            A(1, 1) += weight * deltaY * deltaY;

            // ������� b ����
            b(0) += weight * deltaX * deltaVar;
            b(1) += weight * deltaY * deltaVar;
        }

        // ���� A �����������ݶ�
        Eigen::Vector2d grad = A.inverse() * b;
        gx = grad(0);
        gy = grad(1);
    }

private:
    std::vector<LSGradCoord> neighbors_;
    LSGradCoord current_;
};

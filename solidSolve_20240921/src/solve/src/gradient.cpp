
//#include "solve/include/gradient.h"
//#include <iostream>
//
//int main() {
//    // �������ĵ�
//    LSGradCoord center = { 0.0, 0.0, 0.0 };
//
//    // ���������Ͷ�Ӧ��fֵ
//    std::vector<LSGradCoord> points = {
//        {1.0, 0.0, 1.0},
//        {0.0, 1.0, 1.0},
//        {-1.0, 0.0, -1.0},
//        {0.0, -1.0, -1.0}
//    };
//
//    LSGrad gc(points, center);
//    double gx=0, gy=0;
//    gc.computeGradient(gx, gy);
//
//    std::cout << "Gradient gx: " << gx << ", gy: " << gy << std::endl;
//
//    return 0;
//}

#include "RHcondition/include/RHcondition.h"
#include <iostream>

int main() {
    double rho1 = 1.225; // ����ǰ�ܶȣ���λ kg/m^3
    double u1 = 340; // ����ǰ�ٶ� u ��������λ m/s
    double v1 = 0; // ����ǰ�ٶ� v ��������λ m/s
    double p1 = 101325; // ����ǰѹ������λ Pa
    double us = 500; // �����ٶȣ���λ m/s

    RHcondition rh(rho1, u1, v1, p1, us);
    rh.calculatePostShockState();
    rh.printPostShockState();
    return 0;
}
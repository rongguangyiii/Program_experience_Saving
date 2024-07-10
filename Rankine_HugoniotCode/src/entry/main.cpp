
#include "RHcondition/include/RHcondition.h"
#include <iostream>

int main() {
    double rho1 = 1.225; // 激波前密度，单位 kg/m^3
    double u1 = 340; // 激波前速度 u 分量，单位 m/s
    double v1 = 0; // 激波前速度 v 分量，单位 m/s
    double p1 = 101325; // 激波前压力，单位 Pa
    double us = 500; // 激波速度，单位 m/s

    RHcondition rh(rho1, u1, v1, p1, us);
    rh.calculatePostShockState();
    rh.printPostShockState();
    return 0;
}
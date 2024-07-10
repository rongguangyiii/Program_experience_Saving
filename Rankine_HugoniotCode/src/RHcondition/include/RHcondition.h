
#ifndef RHCONDITION_H
#define RHCONDITION_H
#include <iostream>
#include <map>
#include <cmath>
#include <string>

class RHcondition {
private:
    std::map<std::string, double> upStream_;
    std::map<std::string, double> downStream_;
    double shockSpeed_;
    double gamma_;
    double R_;

    // 目标函数F和其导数
    double F(double P2);
    double dF(double P2);

public:
    RHcondition(double rho1_up = 0.0, double u1_up = 0.0, double v1_up = 0.0, double p1_up = 0.0, double us = 0.0, double gammain = 1.4, double Rin = 287.058);
    RHcondition(const RHcondition& cond);
    RHcondition& operator=(const RHcondition& cond);
    void calculatePostShockState();
    void printPostShockState();
};



#endif // !RHCONDITION_H
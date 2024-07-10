#include "RHcondition/include/RHcondition.h"

RHcondition::RHcondition(double rho1_up, double u1_up, double v1_up, double p1_up, double us, double gammain, double Rin)
    : shockSpeed_(us), gamma_(gammain), R_(Rin) {
    upStream_["rho"] = rho1_up;
    upStream_["u"] = u1_up;
    upStream_["v"] = v1_up;
    upStream_["p"] = p1_up;
}

RHcondition::RHcondition(const RHcondition& cond)
    : upStream_(cond.upStream_), downStream_(cond.downStream_), shockSpeed_(cond.shockSpeed_), gamma_(cond.gamma_), R_(cond.R_) {}

RHcondition& RHcondition::operator=(const RHcondition& cond) {
    if (this != &cond) {
        upStream_ = cond.upStream_;
        downStream_ = cond.downStream_;
        shockSpeed_ = cond.shockSpeed_;
        gamma_ = cond.gamma_;
        R_ = cond.R_;
    }
    return *this;
}

double RHcondition::F(double P2) {
    double rho1 = upStream_["rho"];
    double u1 = upStream_["u"];
    double P1 = upStream_["p"];
    double term1 = P2 / (gamma_ - 1);
    double term2 = 0.5 * ((rho1 * u1 * u1 * (P2 - P1)) / (P1 * (gamma_ + 1) + P2 * (gamma_ - 1)));
    double term3 = P1 / (gamma_ - 1) + 0.5 * rho1 * u1 * u1;
    return term1 + term2 - term3;
}

double RHcondition::dF(double P2) {
    double rho1 = upStream_["rho"];
    double u1 = upStream_["u"];
    double P1 = upStream_["p"];
    double term1 = 1 / (gamma_ - 1);
    double num = rho1 * u1 * u1 * (P1 * (gamma_ + 1) + P2 * (gamma_ - 1));
    double denom = P1 * (gamma_ + 1) + P2 * (gamma_ - 1);
    double term2 = 0.5 * (num / (denom * denom));
    return term1 + term2;
}

void RHcondition::calculatePostShockState() {
    double P2 = upStream_["p"] * 2; // 初始猜测值
    double tol = 1e-8;  // 收敛容限
    int maxIter = 1000; // 最大迭代次数
    int iter = 0;

    while (iter < maxIter) {
        double P2_new = P2 - F(P2) / dF(P2);
        if (std::fabs(P2_new - P2) < tol) {
            P2 = P2_new;
            break;
        }
        P2 = P2_new;
        iter++;
    }

    if (iter == maxIter) {
        std::cerr << "Newton-Raphson方法未收敛" << std::endl;
        return;
    }

    double rho1 = upStream_["rho"];
    double u1 = upStream_["u"];
    double P1 = upStream_["p"];
    downStream_["p"] = P2;
    downStream_["rho"] = rho1 * (u1 / (u1 - shockSpeed_ + shockSpeed_ * std::sqrt(1 + (gamma_ + 1) * (P2 / P1 - 1) / (2 * gamma_))));
    downStream_["u"] = shockSpeed_ * (1 - 1 / std::sqrt(1 + (gamma_ + 1) * (P2 / P1 - 1) / (2 * gamma_)));
    downStream_["v"] = upStream_["v"]; // 假设v分量不变
}

void RHcondition::printPostShockState() {
    std::cout << "激波后的状态:" << std::endl;
    std::cout << "密度 rho2: " << downStream_["rho"] << std::endl;
    std::cout << "速度 u2: " << downStream_["u"] << std::endl;
    std::cout << "速度 v2: " << downStream_["v"] << std::endl;
    std::cout << "压力 p2: " << downStream_["p"] << std::endl;
}
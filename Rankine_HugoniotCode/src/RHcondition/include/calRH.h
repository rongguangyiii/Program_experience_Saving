
#ifndef CALRH_H
#define CALRH_H
#include <iostream>
#include <Eigen/Dense>
#include <map>
#include <cmath>
#include <string>
#include <vector>

class CALRH {
private:
    std::map<std::string, double> downStream_;//����Ŀ�꣺���α߽�״̬
    std::map<std::string, double> UL_;//����״ֱ̬�����Ƶõ����α߽�״̬
    std::map<std::string, double> UR_;//����״̬
    Eigen::Vector3d normal_;
    Eigen::Vector3d shockSpeed_;
    double gamma_;

    // Ŀ�꺯��F���䵼��
    double F(double a, double u, double q, double Mx) const;
    double DF(double a, double u, double q, double Mx) const;

public:
    CALRH(double rho1_up = 0.0, double u1_up = 0.0, double v1_up = 0.0, double p1_up = 0.0, double gammain = 1.4);
    CALRH(const CALRH& cond);
    CALRH& operator=(const CALRH& cond);
    void calculateRH();
    double basicSonic(const double & m_rho, const double& m_p, const double& m_gamma)const;

};



#endif // !CALRH_H
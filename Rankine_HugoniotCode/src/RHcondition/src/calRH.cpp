#include "RHcondition/include/calRH.h"
#include <cmath>


CALRH::CALRH(double rho1_up, double u1_up, double v1_up, double p1_up, double gammain)
	:gamma_(gammain){
	UL_["rho"] = rho1_up;
	UL_["u"] = u1_up;
	UL_["v"] = v1_up;
	UL_["p"] = p1_up;
}

CALRH::CALRH(const CALRH& cond)
	: UL_(cond.UL_), downStream_(cond.downStream_), shockSpeed_(cond.shockSpeed_), gamma_(cond.gamma_) {}

CALRH& CALRH::operator=(const CALRH& cond) {
	if (this != &cond) {
		UL_ = cond.UL_;
		downStream_ = cond.downStream_;
		shockSpeed_ = cond.shockSpeed_;
		gamma_ = cond.gamma_;
	}
	return *this;
}

void CALRH::calculateRH()
{
	double c_l = basicSonic(UL_.at("rho"), UL_.at("p"), gamma_);
	double c_r = basicSonic(UR_.at("rho"), UR_.at("p"), gamma_);
	Eigen::Vector3d uLv(UL_.at("u"), UL_.at("v"), UL_.at("w"));
	Eigen::Vector3d uRv(UR_.at("u"), UR_.at("v"), UR_.at("w"));
	double vn1 = uLv.dot(normal_);//计算UL的法向速度大小,注意点乘变成了标量！！
	double vn2 = uRv.dot(normal_);
	Eigen::Vector3d vt1 = uLv - vn1 * normal_;

	double Ms = (UL_.at("rho") * vn1 - UR_.at("rho") * vn2) / (UL_.at("rho") * UR_.at("rho"));//激波速度
	double M0 = (vn1 - Ms) / c_l;
	double J3 = 2.0 / (gamma_ - 1.0) * c_r - uRv.dot(normal_);//左传黎曼不变量
	M0 = std::pow(M0, 2);//在这平方

	double tol = 1e-6;  // 收敛容限
	size_t maxIter = 1000; // 最大迭代次数
	size_t iter = 0;
	while (iter < maxIter)
	{
		double M_new = M0 - F(c_l, vn1, J3, M0) / DF(c_l, vn1, J3, M0);
		if (M_new < 0)
			M_new = 0.5 * M0;
		if (std::fabs(M_new - M0) < tol)
		{
			M0 = M_new;
			break;
		}
		M0 = M_new;
		iter++;
	}
	if (iter == maxIter) {
		std::cerr << "Newton-Raphson方法未收敛" << std::endl;
		return;
	}

	Eigen::Vector3d vs = (vn1 - pow(M0, 0.5) * c_l) * normal_;//激波速度
	double gm1 = gamma_ - 1.0;
	double gm2 = gamma_ + 1.0;
	double ratio_p = 2.0 * gamma_ * M0 - gm1;
	ratio_p = ratio_p / gm2;
	double ratio_r = gm2 * M0;
	ratio_r = ratio_r / (gm1 * M0 + 2.0);
	double ratio_m = M0 + 2.0 / gm1;
	ratio_m = ratio_m / (2.0 * gamma_ / gm1 * M0 - 1.0);

	shockSpeed_ = vs;
	downStream_["rho"] = ratio_r * UL_.at("rho");
	downStream_["p"] = ratio_p * UL_.at("p");
	double cc = basicSonic(downStream_.at("rho"), downStream_.at("p"), gamma_);
	double vr = sqrt(ratio_m) * cc;
	Eigen::Vector3d vt2 = vt1;
	Eigen::Vector3d ubnd = vt2 + vr * normal_ + vs;
	downStream_["u"] = ubnd.x();
	downStream_["v"] = ubnd.y();
	downStream_["w"] = ubnd.z();

}

double CALRH::F(double a, double u, double q, double Mx)const
{
	double r1 = gamma_ - 1;
	double r2 = gamma_ + 1;
	double r3 = 1 / r1;
	double r4 = 2.0 * gamma_ * Mx - r1;
	double r5 = r1 * Mx + 2.0;
	double fx = 0.0;
	if (r4 * r5 / Mx < 0)
		fx = (Mx - 1.0) / pow(Mx, 0.5);
	else
		fx = r3 * sqrt(r4 * r5 / Mx) + (Mx - 1.0) / pow(Mx, 0.5);

	return fx - 0.5 * r2 / a * (q + u);
}

double CALRH::DF(double a, double u, double q, double Mx)const
{
	double r1 = gamma_ - 1.0;
	double r2 = gamma_ + 1.0;
	double r3 = 1.0 / r1;
	double x2 = Mx - 1.0;
	double x3 = 0.5 * x2 / pow(Mx, 1.5);
	double x4 = 2.0 * gamma_ * Mx - r1;
	double x5 = r1 * Mx + 2.0;

	double df = 2.0 * gamma_ * (r1 * Mx + 2.0) / Mx + r1 * x4 / Mx - x4 * x5 / pow(Mx, 2.0);
	if (x4 * x5 / Mx < 0)
		df = -x3 + 1.0 / pow(Mx, 0.5);
	else
		df = df / (2.0 * r1 * sqrt(x4 * x5 / Mx)) - x3 + 1.0 / pow(Mx, 0.5);
	return df;
}

double CALRH::basicSonic(const double& m_rho, const double& m_p, const double& m_gamma) const
{
	return sqrt(m_gamma * m_p / m_rho);
}
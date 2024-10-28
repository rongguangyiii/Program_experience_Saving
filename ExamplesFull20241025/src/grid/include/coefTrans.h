#pragma once

class CoefTrans
{
public:
	CoefTrans() :ksi_t_(0.0), ksi_x_(0.0), ksi_y_(0.0),
		eta_t_(0.0), eta_x_(0.0), eta_y_(0.0), jacob_(0.0) {}
	CoefTrans(double ksx, double ksy, double etx, double eaty, double Jac, double kst=0.0, double ett=0.0)
		: ksi_x_(ksx), ksi_y_(ksy), eta_x_(etx), eta_y_(eaty), jacob_(Jac), ksi_t_(kst), eta_t_(ett) {}
	~CoefTrans() {};
	void setCoordTrans(double inksi_x, double inksi_y, double ineta_x, double ineat_y, double inJacobian, double inksi_t=0.0, double ineta_t=0.0);
	double ksi_x() const { return ksi_x_; }
	double ksi_y() const { return ksi_y_; }
	double eta_x() const { return eta_x_; }
	double eta_y() const { return eta_y_; }
	double jacob() const { return jacob_; }
	double ksi_t() const { return ksi_t_; }
	double eta_t() const { return eta_t_; }
	// Overload assignment operator
	CoefTrans& operator=(const CoefTrans& other);

private:
	double ksi_t_, ksi_x_, ksi_y_;
	double eta_t_, eta_x_, eta_y_;
	double jacob_;
};


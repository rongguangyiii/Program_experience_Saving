#include "gridGenerate/include/coefTrans.h"

void CoefTrans::setCoordTrans(double inksi_x, double inksi_y, double ineta_x, double ineta_y, double inJacobian, double inksi_t, double ineta_t)
{
	ksi_x_ = inksi_x;
	ksi_y_ = inksi_y;
	eta_x_ = ineta_x;
	eta_y_ = ineta_y;
	jacob_ = inJacobian;
	ksi_t_ = inksi_t;
	eta_t_ = ineta_t;
}

CoefTrans& CoefTrans::operator=(const CoefTrans& other) {
	if (this == &other) return *this; // Handle self-assignment
	ksi_x_ = other.ksi_x_;
	ksi_y_ = other.ksi_y_;
	eta_x_ = other.eta_x_;
	eta_y_ = other.eta_y_;
	jacob_ = other.jacob_;
	ksi_t_ = other.ksi_t_;
	eta_t_ = other.eta_t_;
	return *this;
}

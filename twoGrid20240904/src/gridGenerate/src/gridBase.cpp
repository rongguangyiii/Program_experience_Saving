#include "gridGenerate/include/gridBase.h"


void GridBase::setCoefTrans(const CoordPtr& p, const CoefTrans& coef)
{
	coefTrans_[p->id_] = coef;
}


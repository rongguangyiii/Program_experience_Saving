
#include "flowSolve/include/ctrlFlow.h"
#include "gridGenerate/include/kdTreee.h"
#include <exception>  // ���� std::exception ���ͷ�ļ�
#include "triangle/include/triBase.h"

void CtrlFlow::solve() {

	markBackgrid();
	//markBackWallBound();
	treatBound();
	//movepoints_ = gridVec_.at(1)->getPoints();//�������Զ��仯�𣿻ᣡ
	movepoints_.clear();
	for (const auto& it : gridVec_.at(1)->getPoints())
		movepoints_.push_back(*it);
	calCoefNoTime(gridVec_.at(0));
	calCoefTime(gridVec_.at(1));
	/* ------------------------------------------------------------------------------------------------------------------
	 ˵�����ڼ��㶯����ʱ�����ڵ�1��ʱ���ƽ�ʹ�õ�����͵�0����������һ���ġ������������漸�д���ļ�����ǵ�0��allcoefTrans_�Ķ���ϵ����
		  ͬʱҲ������Ϊ���ڼ����1����allcoefTrans_���ڳ���ʵ�ֵ�ʱ��Ϊ�˷��㣬��allcoefTrans_old_�����˵�0�������ݣ���1����������
		  allcoefTrans_��ʾ������allcoefTrans_old_������nʱ�̵����ݣ�n+1ʱ�̵�������allcoefTrans_��ʾ��
	 -------------------------------------------------------------------------------------------------------------------*/
	allcoefTrans_old_ = allcoefTrans_;
	//isfirst_ = false;

}

void CtrlFlow::move(double vel, double dt)
{
	//step 1 �ƶ�����ǰ�ȸ���һ��nʱ�̵�����Ͷ���ϵ��,���ڼ���ksi_t,eta_t
	//movepoints_ = gridVec_.at(1)->getPoints();//q:�������Զ��仯��a:��
	movepoints_.clear();
	for (const auto & it: gridVec_.at(1)->getPoints()){
		movepoints_.push_back(*it);
	}
	allcoefTrans_old_ = allcoefTrans_;
	genOldGridNeighbors();//������һ�����е���ھӹ�ϵ��
	//step 2 �ƶ�����
	moveGrid(vel, dt);
	// ��Ҫ�趨һ���о�(To Do...)��������ÿһ��ʱ���ƽ�����Ҫ���¸����ھӺͱ�ǽڵ㣬�߽硣
	if (isSwallowPoints()) {
		updateNeighbors();//ֻ�ڸ�����������ھӹ�ϵ���͸��Խڵ�ĳ�ʼtype_��
		markBackgrid();//�ƶ��󣬷����������������Ҫ���±�Ǳ�������ڵ����͡�
		treatBound();//��������������ھӹ�ϵ��
	}
	else if (isSpitPoints())
	{
		updateNeighbors();
		markBackgrid();
		treatBound();
		foundSpitPoints();
		calSpitPointsCoesfTrans();
	}
	//step 3 ���¼������ϵ��
	calCoefNoTime(gridVec_.at(0));
	calCoefTime(gridVec_.at(1));
}

void CtrlFlow::moveGrid(double vel, double dt)
{
	//vel�����������������ٶȣ�dt����ǰ����ʱ�䲽,��
	auto& mgrid = gridVec_.at(1);
	const auto& mPointsvec = mgrid->getPoints();
	size_t mx = mgrid->getXnum();
	size_t my = mgrid->getYnum();
	auto& m_neighbors = mgrid->getReneighbors();
	for (size_t i = 0; i < mx; i++)
	{
		double tempv = 3*(vel / mx) + (vel / mx) * i;//���Էֲ�
		for (size_t j = 0; j < my; j++)
		{
			auto& perP = mPointsvec[i * my + j];
			perP->x = perP->x + tempv * dt;
			//perP->y = perP->y + tempv * dt;
		}
	}
}

void CtrlFlow::markBackgrid()
{
	std::vector<CoordPtr> boundpoints = foundMoveGridBound();
	GlobalData::sortPointsCounterclockwise(boundpoints);
	auto& back = gridVec_.at(0)->getPoints();
	for (auto& it : back) {
		if (GlobalData::isPointInPolygon_moreInner(it, boundpoints))
		{
			it->type_ = PointType::b_notcal;
			it->tag_ = PointTag::shade;
		}
	}

}

void CtrlFlow::treatBound()
{
	auto& mgrid = gridVec_.at(1);
	const auto& mPointsvec = mgrid->getPoints();
	size_t mx = mgrid->getXnum();
	size_t my = mgrid->getYnum();
	auto& m_neighbors = mgrid->getReneighbors();

	auto& curgrid = gridVec_.at(0);
	size_t xpoints = curgrid->getXnum();
	size_t ypoints = curgrid->getYnum();
	auto& points = curgrid->getPoints();
	auto& neighbors = curgrid->getReneighbors();

	for (size_t i = 0; i < xpoints; i++) {
		for (size_t j = 0; j < ypoints; j++) {

			auto& curpoint = points[i * ypoints + j];
			if (curpoint->type_ == PointType::b_notcal)
				continue;

			std::vector<CoordPtr>& neigvec = neighbors[curpoint->id_];

			// ������ǰ�ھ�
			for (auto it = neigvec.begin(); it != neigvec.end(); ++it) {
				if ((*it)->type_ == PointType::b_notcal) {
					// �ҵ���Ҫ�滻���ھ�
					CoordPtr oldNeighbor = *it;

					// ѡ���µ��ھӽ����滻
					CoordPtr newNeighbor;
					if (curpoint->x < mPointsvec[j]->x) {
						newNeighbor = mPointsvec[j];
						curpoint->type_ = PointType::b_Lbound;
					}
					else if (curpoint->x > mPointsvec[j]->x) {
						newNeighbor = mPointsvec[(mx - 1) * my + j];
						curpoint->type_ = PointType::b_Rbound;
					}
					else {
						throw std::runtime_error("Error message in treat bound!");
					}

					// �滻�ھӣ������ھӷŵ���ɾ���ھӵ�λ��
					*it = newNeighbor;

					// ������һ��������ӵ���ھӵ�˫���ϵ
					m_neighbors[newNeighbor->id_].push_back(curpoint);//��һ�������иõ�û�б���ȷ����
				}
			}


		}
	}

}

std::vector<CoordPtr> CtrlFlow::foundMoveGridBound()
{
	auto& curgrid = gridVec_.at(1);
	size_t xpoints = curgrid->getXnum();
	size_t ypoints = curgrid->getYnum();
	auto& points = curgrid->getPoints();

	//���㱳������Ĳ�����Ŀ���Ǽӿ���Ʒ�Χ��ȷ�������ƶ����������������ϵı߽�㲻Ҫ̫�ӽ���
	size_t bx = gridVec_.at(0)->getXnum();
	size_t by = gridVec_.at(0)->getYnum();
	double xsta = gridVec_.at(0)->getXstart();
	double ysta = gridVec_.at(0)->getYstart();
	double xend = gridVec_.at(0)->getXend();
	double yend = gridVec_.at(0)->getYend();
	double xstep = (xend - xsta) / (bx-1);
	double ystep = (yend - ysta) / (by-1);
	Coord Xrangepoint(xstep / 10, 0);
	Coord Yrangepoint(0, ystep / 10);
	std::vector<CoordPtr> moveboud;
	for (size_t j = 0; j < ypoints; j++) {
		Coord curpoint = *points[0 * ypoints + j];
		moveboud.push_back(std::make_shared<Coord>(curpoint -= Xrangepoint));
	}
	for (size_t j = 0; j < ypoints; j++) {
		Coord curpoint = *points[(xpoints - 1) * ypoints + j];
		moveboud.push_back(std::make_shared<Coord>(curpoint += Xrangepoint));
	}

	//�����ظ����
	for (size_t i = 1; i < xpoints - 1; i++) {
		Coord curpoint = *points[i * ypoints + 0];
		moveboud.push_back(std::make_shared<Coord>(curpoint -= Yrangepoint));

	}
	for (size_t i = 1; i < xpoints - 1; i++) {
		Coord curpoint = *points[i * ypoints + (ypoints - 1)];
		moveboud.push_back(std::make_shared<Coord>(curpoint += Yrangepoint));
	}
	return moveboud;
}

void CtrlFlow::calCoefNoTime(const std::shared_ptr<GridBase>& curgrid)
{
	const auto& Pointsvec = curgrid->getPoints();
	size_t numx = curgrid->getXnum();
	size_t numy = curgrid->getYnum();
	auto& neighbors = curgrid->getReneighbors();

	double ksi_t = 0, ksi_x = 0, ksi_y = 0, eta_x = 0, eta_y = 0, eta_t = 0;
	double jacob = 0.0;
	double XK = 0, XI = 0, YK = 0, YI = 0, det = 0, XT = 0, YT = 0;
	double dt = 0.01;
	//ֻ�����ڲ���Ķ���ϵ�����߽��Ŀǰ��������,�ҽ������㲻��ȷ��Ҳû����һʱ�̵�����
	for (size_t i = 1; i < numx - 1; i++) {
		for (size_t j = 1; j < numy - 1; j++)
		{
			const auto& it = Pointsvec[i * numy + j];
			//ֻ����
			if (it->type_ == PointType::b_notcal)
				continue;
			Tools::checkKeyExists(neighbors, it->id_, "neighbors");
			auto& perneighbors = neighbors[it->id_];
			//Ӧ�����ھӵ����������ϵ��,������Ҫ��һ����������ȷ���Ƿ������ȷ���ؼ����ھӵ�λ���Ƿ���ȷ��

			CoordPtr& pR = GlobalData::PointAtAngle(perneighbors, 0);
			CoordPtr& pU = GlobalData::PointAtAngle(perneighbors, GlobalData::PI / 2);
			CoordPtr& pL = GlobalData::PointAtAngle(perneighbors, GlobalData::PI);
			CoordPtr& pD = GlobalData::PointAtAngle(perneighbors, GlobalData::PI * 3 / 2);

			XK = (pR->x - pL->x) / 2.0;
			YK = (pR->y - pL->y) / 2.0;
			XI = (pU->x - pD->x) / 2.0;
			YI = (pU->y - pD->y) / 2.0;

			det = XK * YI - XI * YK;  // ����Am������ʽ

			ksi_x = YI / det;
			ksi_y = -XI / det;
			eta_x = -YK / det;
			eta_y = XK / det;
			jacob = 1.0 / det;//�ſ˱�Ϊ����Am������ʽ�ĵ���
			allcoefTrans_[it->id_] = CoefTrans(ksi_x, ksi_y, eta_x, eta_y, jacob);
			curgrid->setCoefTrans(it, allcoefTrans_[it->id_]);
		}
	}

}

void CtrlFlow::calCoefTime(const std::shared_ptr<GridBase>& curgrid)
{
	const auto& Pointsvec = curgrid->getPoints();
	const auto& oldPointsvec = movepoints_;
	size_t numx = curgrid->getXnum();
	size_t numy = curgrid->getYnum();
	auto& neighbors = curgrid->getReneighbors();

	double ksi_t = 0, ksi_x = 0, ksi_y = 0, eta_x = 0, eta_y = 0, eta_t = 0;
	double jacob = 0.0;
	double XK = 0, XI = 0, YK = 0, YI = 0, JJ = 0, XT = 0, YT = 0;
	double dt = 0.01;
	//�ƶ�������Ҫ�������ұ߽�
	for (size_t i = 0; i < numx; i++) {
		for (size_t j = 1; j < numy - 1; j++)
		{
			//Ӧ�����ھӵ����������ϵ��
			const auto& it = Pointsvec[i * numy + j];
			const auto& oldit = oldPointsvec[i * numy + j];
			Tools::checkKeyExists(neighbors, it->id_, "neighbors");
			auto& perneighbors = neighbors[it->id_];
			CoordPtr& pR = GlobalData::PointAtAngle(perneighbors, 0);
			CoordPtr& pU = GlobalData::PointAtAngle(perneighbors, GlobalData::PI / 2);
			CoordPtr& pL = GlobalData::PointAtAngle(perneighbors, GlobalData::PI);
			CoordPtr& pD = GlobalData::PointAtAngle(perneighbors, GlobalData::PI * 3 / 2);

			XK = (pR->x - pL->x) / 2.0;
			YK = (pR->y - pL->y) / 2.0;
			XI = (pU->x - pD->x) / 2.0;
			YI = (pU->y - pD->y) / 2.0;

			XT = (it->x - oldit.x) / dt;
			YT = (it->y - oldit.y) / dt;

			JJ = XK * YI - XI * YK;  // ����Am������ʽ

			ksi_t = (XI * YT - XT * YI) / JJ;
			ksi_x = YI / JJ;
			ksi_y = -XI / JJ;
			eta_x = -YK / JJ;
			eta_y = XK / JJ;
			eta_t = (XT * YK - XK * YT) / JJ;
			jacob = 1.0 / JJ;//�ſ˱�Ϊ����Am������ʽ�ĵ���
			allcoefTrans_[it->id_] = CoefTrans(ksi_x, ksi_y, eta_x, eta_y, jacob, ksi_t, eta_t);
			curgrid->setCoefTrans(it, allcoefTrans_[it->id_]);
		}
	}

}

void CtrlFlow::outplot_B(std::ofstream& out, const std::string& zone_title)
{
	//step1 �Ȱ�Ҫ����ĵ�ɾѡ����
	auto& mgrid = gridVec_.at(1);
	const auto& mPointsvec = mgrid->getPoints();
	size_t mx = mgrid->getXnum();
	size_t my = mgrid->getYnum();

	auto& bgrid = gridVec_.at(0);
	size_t bx = bgrid->getXnum();
	size_t by = bgrid->getYnum();
	auto& bPointsVec = bgrid->getPoints();
	std::vector<CoordPtr> outPoints;

	//step 1 ȷ�����
	size_t LB = 0;
	size_t RB = 0;
	for (size_t i = 0; i < bx; i++) {
		if (bPointsVec[i * by + 0]->type_ == PointType::b_Lbound)
			LB = i;
		if (bPointsVec[i * by + 0]->type_ == PointType::b_Rbound)
			RB = i;
	}
	//step 2 ����
	for (size_t i = 0; i <= LB; i++) {
		for (size_t j = 0; j < by; j++)
		{
			const auto& it = bPointsVec[i * by + j];
			outPoints.emplace_back(it);
		}
	}
	//step 3 ���movegrid�е���߽磬
	for (size_t j = 0; j < my; j++)
	{
		const auto& it = mPointsvec[0 * my + j];
		outPoints.emplace_back(it);
	}
	//step 4 ���movegrid�е��ұ߽磬
	for (size_t j = 0; j < my; j++)
	{
		const auto& it = mPointsvec[(mx - 1) * my + j];
		outPoints.emplace_back(it);
	}
	//step 4 �Ұ�ߣ�
	for (size_t i = RB; i < bx; i++) {
		for (size_t j = 0; j < by; j++)
		{
			const auto& it = bPointsVec[i * by + j];
			outPoints.emplace_back(it);
		}
	}

	//step 5 �γɵ�Ԫ
	std::vector<std::vector<size_t>> outelements;
	for (size_t i = 0; i < LB + 1; i++)
	{
		for (size_t j = 0; j < by - 1; j++)
		{
			size_t p1 = i * by + j + 1;
			size_t p2 = p1 + 1;
			size_t p3 = p1 + by;
			size_t p4 = p3 + 1;
			outelements.push_back({ p1, p3, p4, p2 });
		}

	}

	size_t outRow = 0;
	size_t row_num = bx - 1 - RB;
	for (size_t i = 0; i < row_num + 1; i++)
	{
		outRow = LB + 2 + i;
		for (size_t j = 0; j < by - 1; j++)
		{
			size_t p1 = outRow * by + j + 1;
			size_t p2 = p1 + 1;
			size_t p3 = p1 + by;
			size_t p4 = p3 + 1;
			outelements.push_back({ p1, p3, p4, p2 });
		}

	}

	out << "TITLE = \"" << zone_title << "\"\n";
	out << "VARIABLES = \"X\", \"Y\"\n";
	out << "ZONE T=\"" << zone_title << "\", N=" << outPoints.size() << ", E=" << outelements.size() << ", F=FEPOINT, ET=QUADRILATERAL\n";

	for (const auto& point : outPoints) {
		out << std::fixed << std::setprecision(6) << point->x << " " << point->y << "\n";
	}

	for (const auto& elem : outelements) {
		out << elem[0] << " " << elem[1] << " " << elem[2] << " " << elem[3] << "\n";
	}

	std::cout << zone_title << " output complete!\n";
}

void CtrlFlow::outplot_M(std::ofstream& out, const std::string& zone_title) {
	gridVec_.at(1)->outputToTecplot(out, zone_title);
}

void CtrlFlow::markBackWallBound()
{
	//�ú���Ŀǰû��
	auto& mgrid = gridVec_.at(1);
	const auto& mPointsvec = mgrid->getPoints();
	size_t mx = mgrid->getXnum();
	size_t my = mgrid->getYnum();
	auto& m_neighbors = mgrid->getReneighbors();
	std::vector<Coord> mgpoints;
	for (const auto &perit:mPointsvec)
	{
		mgpoints.push_back(*perit);
	}
	PointKDTreeManager mgPointsKDtree(mgpoints,"2D");

	auto& curgrid = gridVec_.at(0);
	size_t xpoints = curgrid->getXnum();
	size_t ypoints = curgrid->getYnum();
	double sta = curgrid->getXstart();
	double end = curgrid->getXend();
	auto& points = curgrid->getPoints();
	auto& neighbors = curgrid->getReneighbors();
	double step = (end - sta) / (xpoints - 1);

	std::vector<CoordPtr> outerP;

	for (size_t i = 0; i < xpoints; i++) {
		for (size_t j = 0; j < ypoints; j++) {

			auto& curpoint = points[i * ypoints + j];
			if (curpoint->type_ == PointType::b_notcal)
				continue;
			const auto RinPoints = mgPointsKDtree.findPointsWithinRadius(*curpoint, step /** std::sqrt(2)*/);
			for (size_t rh :RinPoints)
			{
				if (mgpoints[rh].type_==PointType::s_filedbound&& curpoint->x < mPointsvec[j]->x)//������!
				{
					curpoint->type_ = PointType::b_Lbound;
					outerP.push_back(curpoint);
				}
				else if (mgpoints[rh].type_ == PointType::s_filedbound && curpoint->x > mPointsvec[j]->x)
				{
					curpoint->type_ = PointType::b_Rbound;
					outerP.push_back(curpoint);

				}
			}
			std::vector<CoordPtr> outerP;
		}
	}

	std::vector<CoordPtr> boundpoints = foundMoveGridBound();

	TriBase tribase(outerP, boundpoints,"ring_points");
	tribase.boostMesh();
}

void CtrlFlow::updateNeighbors()
{
	for (const auto& perit : gridVec_)
	{
		perit->genNeiborNode();
	}
}

bool CtrlFlow::isSwallowPoints()
{
	isSwallow_ = false;
	size_t bx = gridVec_.at(0)->getXnum();
	size_t by = gridVec_.at(0)->getYnum();
	double xsta = gridVec_.at(0)->getXstart();
	double ysta = gridVec_.at(0)->getYstart();
	double xend = gridVec_.at(0)->getXend();
	double yend = gridVec_.at(0)->getYend();
	double xstep = (xend - xsta) / (bx - 1);
	double ystep = (yend - ysta) / (by - 1);

	for (auto& mit : gridVec_.at(1)->getPoints()) {

		if (mit->type_ != PointType::s_filedbound)
			continue;
		auto& m_neighbors = gridVec_.at(1)->getReneighbors();
		if (m_neighbors.find(mit->id_) == m_neighbors.end())
			std::cout << "error!\n";
		auto& nneighbors = m_neighbors[mit->id_];
		for (auto& nit : nneighbors)
		{
			if (nit->id_.gridId != 0)//�ҵ����������
				continue;
			auto dis = Tools::Distance(nit, mit);
			if (dis < (xstep / 10))
			{
				isSwallow_ = true;
				return isSwallow_;
			}
		}
	}
	return isSwallow_;

}

bool CtrlFlow::isSpitPoints()
{
	isSpit_ = false;
	size_t bx = gridVec_.at(0)->getXnum();
	size_t by = gridVec_.at(0)->getYnum();
	double xsta = gridVec_.at(0)->getXstart();
	double ysta = gridVec_.at(0)->getYstart();
	double xend = gridVec_.at(0)->getXend();
	double yend = gridVec_.at(0)->getYend();
	double xstep = (xend - xsta) / (bx - 1);
	double ystep = (yend - ysta) / (by - 1);

	for (auto& mit : gridVec_.at(1)->getPoints()) {

		if (mit->type_ != PointType::s_filedbound)
			continue;
		auto& m_neighbors = gridVec_.at(1)->getReneighbors();
		if (m_neighbors.find(mit->id_) == m_neighbors.end())
			std::cout << "error!\n";
		auto& nneighbors = m_neighbors[mit->id_];
		for (auto& nit : nneighbors)
		{
			if (nit->id_.gridId != 0)//�ҵ����������
				continue;
			auto dis = Tools::Distance(nit, mit);
			if (dis >(xstep + xstep / 10))
			{
				isSpit_ = true;
				return isSpit_;
			}
		}
	}
	return isSpit_;
}

void CtrlFlow::foundSpitPoints()
{
	std::vector<CoordPtr> boundpoints = foundMoveGridBound();
	GlobalData::sortPointsCounterclockwise(boundpoints);
	auto& back = gridVec_.at(0)->getPoints();//�����ڱ��������ҡ�
	for (auto& it : back) {
		if (GlobalData::isPointInPolygon_moreInner(it, boundpoints))
		{
			//std::cout << "������ !\n";
		}
		else if(!GlobalData::isPointInPolygon_moreInner(it, boundpoints)&&(it->tag_ == PointTag::shade))
		{
			if (std::fabs(it->y - gridVec_.at(0)->getYstart()) < 1E-6 || 
				std::fabs(it->y - gridVec_.at(0)->getYend()) < 1E-6)
			{
				it->tag_ = PointTag::surface;//�߽�㲻�������
				//std::cout << "�߽�� !\n";
			}
			else if (it->y > gridVec_.at(0)->getYstart() && it->y < gridVec_.at(0)->getYend())//�߽�㲻�������
			{
				spitoutPoints_.push_back(it);
				it->tag_ = PointTag::spit;//�ҵ��Ժ��ٽ�������tag_
				//std::cout << "Ŀ��� !\n";

			}
			else
			{
				//std::cout << "ERROR����\n";
			}
		}
		else
		{
			//std::cout << "�����棬�ҵ���һ��û������!\n";
		}
	}

	//step2. �����ҵ���spit�㣬�ҵ���Ӧ��moveGrid�ϵ����ӵ㣬��Ϊ�³����ĵ���������Ǹ��µġ�
	const auto& BPointsvec = gridVec_.at(0)->getPoints();
	auto& B_neighbors = gridVec_.at(0)->getReneighbors();
	//std::vector<PointId> tempId;
	for (const auto& perid : spitoutPoints_)
	{
		const auto& perNeighbor = B_neighbors[perid->id_];
		for (const auto& it : perNeighbor)
		{
			//����spitoutPoint���ھӵ㣬�Ƿ���movegrid�ĵ㡣
			if (it->type_ == PointType::s_filedbound)
			{
				//tempId.push_back(it->id_);//�ҵ���id
				//�洢��spit��ֱ�����ӵģ��˶������Ӧ����һʱ�̵�λ�á�
				OldconnectPoints_.push_back(std::make_shared<Coord>(getoldCoord(it->id_)));

			}
		}
	
	}

}

void CtrlFlow::calSpitPointsCoesfTrans()
{
	//step2.Ϊ�˺�������spit���u_new����Ҫ�����ⲿ�ֵ������ϵ����
	double ksi_t = 0, ksi_x = 0, ksi_y = 0, eta_x = 0, eta_y = 0, eta_t = 0;
	double jacob = 0.0;
	double XK = 0, XI = 0, YK = 0, YI = 0, JJ = 0, XT = 0, YT = 0;
	double dt = 0.01;
	for (size_t i = 0; i < spitoutPoints_.size(); ++i)
	{
		//����һʱ��û�б��³����ģ���ӦmoveGrid����ھӹ�ϵ���������ϵ��
		const auto& newit = spitoutPoints_[i];
		const auto& oldit = OldconnectPoints_[i];//movegrid���ƶ�ǰ��λ��
		Tools::checkKeyExists(oldmoveGridNeighbors, oldit->id_, "oldmoveGridNeighbors");
		auto perneighbors = Tools::cood2CoodPtr(oldmoveGridNeighbors[oldit->id_]);
		CoordPtr& pR = GlobalData::PointAtAngle(perneighbors, 0);
		CoordPtr& pU = GlobalData::PointAtAngle(perneighbors, GlobalData::PI / 2);
		CoordPtr& pL = GlobalData::PointAtAngle(perneighbors, GlobalData::PI);
		CoordPtr& pD = GlobalData::PointAtAngle(perneighbors, GlobalData::PI * 3 / 2);

		XK = (pR->x - pL->x) / 2.0;
		YK = (pR->y - pL->y) / 2.0;
		XI = (pU->x - pD->x) / 2.0;
		YI = (pU->y - pD->y) / 2.0;
		XT = (newit->x - oldit->x) / dt;
		YT = (newit->y - oldit->y) / dt;

		JJ = XK * YI - XI * YK;  // ����Am������ʽ

		ksi_t = (XI * YT - XT * YI) / JJ;
		ksi_x = YI / JJ;
		ksi_y = -XI / JJ;
		eta_x = -YK / JJ;
		eta_y = XK / JJ;
		eta_t = (XT * YK - XK * YT) / JJ;
		jacob = 1.0 / JJ;//�ſ˱�Ϊ����Am������ʽ�ĵ���
		spitoutPointsCoefTrans_[oldit->id_] = CoefTrans(ksi_x, ksi_y, eta_x, eta_y, jacob, ksi_t, eta_t);
	}

}

Coord CtrlFlow::getoldCoord(PointId index)
{
	for (const auto& it : movepoints_)
	{
		if (it.id_ == index)
		{
			return it;
		}
	}

	// ���û���ҵ���Ӧ�ĵ㣬���׳��쳣
	throw std::runtime_error("Point with given index not found in movepoints_");
}

void CtrlFlow::genOldGridNeighbors()
{
	for (const auto& perit : gridVec_)
	{
		const auto& curNeighbors = perit->getReneighbors();
		for (const auto&it: curNeighbors)
		{
			oldmoveGridNeighbors[it.first] = Tools::coodPtr2Cood(it.second);
		}
	}

}
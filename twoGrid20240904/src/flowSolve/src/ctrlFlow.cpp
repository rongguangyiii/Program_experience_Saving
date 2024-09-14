
#include "flowSolve/include/ctrlFlow.h"
#include "gridGenerate/include/kdTreee.h"
#include <exception>  // 包含 std::exception 类的头文件
#include "triangle/include/triBase.h"

void CtrlFlow::solve() {

	markBackgrid();
	//markBackWallBound();
	treatBound();
	//movepoints_ = gridVec_.at(1)->getPoints();//这样会自动变化吗？会！
	movepoints_.clear();
	for (const auto& it : gridVec_.at(1)->getPoints())
		movepoints_.push_back(*it);
	calCoefNoTime(gridVec_.at(0));
	calCoefTime(gridVec_.at(1));
	/* ------------------------------------------------------------------------------------------------------------------
	 说明：在计算动网格时，由于第1步时间推进使用的网格和第0步的网格是一样的。所以这里上面几行代码的计算的是第0步allcoefTrans_的度量系数，
		  同时也可以认为是在计算第1步的allcoefTrans_。在程序实现的时候，为了方便，用allcoefTrans_old_代表了第0步的内容，第1步的内容用
		  allcoefTrans_表示。后续allcoefTrans_old_代表了n时刻的内容，n+1时刻的内容用allcoefTrans_表示。
	 -------------------------------------------------------------------------------------------------------------------*/
	allcoefTrans_old_ = allcoefTrans_;
	//isfirst_ = false;

}

void CtrlFlow::move(double vel, double dt)
{
	//step 1 移动网格前先复制一份n时刻的坐标和度量系数,用于计算ksi_t,eta_t
	//movepoints_ = gridVec_.at(1)->getPoints();//q:这样会自动变化吗？a:会
	movepoints_.clear();
	for (const auto & it: gridVec_.at(1)->getPoints()){
		movepoints_.push_back(*it);
	}
	allcoefTrans_old_ = allcoefTrans_;
	genOldGridNeighbors();//保存上一步所有点的邻居关系。
	//step 2 移动网格
	moveGrid(vel, dt);
	// 需要设定一个判据(To Do...)，并不是每一步时间推进都需要重新更新邻居和标记节点，边界。
	if (isSwallowPoints()) {
		updateNeighbors();//只在各自网格更新邻居关系，和各自节点的初始type_。
		markBackgrid();//移动后，发生有吞吐网格点需要重新标记背景网格节点类型。
		treatBound();//在相邻网格更新邻居关系。
	}
	else if (isSpitPoints())
	{
		updateNeighbors();
		markBackgrid();
		treatBound();
		foundSpitPoints();
		calSpitPointsCoesfTrans();
	}
	//step 3 重新计算度量系数
	calCoefNoTime(gridVec_.at(0));
	calCoefTime(gridVec_.at(1));
}

void CtrlFlow::moveGrid(double vel, double dt)
{
	//vel代表整个网格的最大速度，dt代表当前计算时间步,。
	auto& mgrid = gridVec_.at(1);
	const auto& mPointsvec = mgrid->getPoints();
	size_t mx = mgrid->getXnum();
	size_t my = mgrid->getYnum();
	auto& m_neighbors = mgrid->getReneighbors();
	for (size_t i = 0; i < mx; i++)
	{
		double tempv = 3*(vel / mx) + (vel / mx) * i;//线性分布
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

			// 遍历当前邻居
			for (auto it = neigvec.begin(); it != neigvec.end(); ++it) {
				if ((*it)->type_ == PointType::b_notcal) {
					// 找到需要替换的邻居
					CoordPtr oldNeighbor = *it;

					// 选择新的邻居进行替换
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

					// 替换邻居：将新邻居放到被删除邻居的位置
					*it = newNeighbor;

					// 更新另一层网格被添加点的邻居的双向关系
					m_neighbors[newNeighbor->id_].push_back(curpoint);//另一层网格中该点没有被正确排序
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

	//计算背景网格的步长，目的是加宽控制范围以确定网格移动后控制两个网格层上的边界点不要太接近。
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

	//避免重复添加
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
	//只计算内部点的度量系数，边界点目前都不计算,且交界点计算不正确，也没有上一时刻的数据
	for (size_t i = 1; i < numx - 1; i++) {
		for (size_t j = 1; j < numy - 1; j++)
		{
			const auto& it = Pointsvec[i * numy + j];
			//只计算
			if (it->type_ == PointType::b_notcal)
				continue;
			Tools::checkKeyExists(neighbors, it->id_, "neighbors");
			auto& perneighbors = neighbors[it->id_];
			//应该用邻居点来计算度量系数,这里需要用一个简单网格来确定是否计算正确，关键在邻居点位置是否正确。

			CoordPtr& pR = GlobalData::PointAtAngle(perneighbors, 0);
			CoordPtr& pU = GlobalData::PointAtAngle(perneighbors, GlobalData::PI / 2);
			CoordPtr& pL = GlobalData::PointAtAngle(perneighbors, GlobalData::PI);
			CoordPtr& pD = GlobalData::PointAtAngle(perneighbors, GlobalData::PI * 3 / 2);

			XK = (pR->x - pL->x) / 2.0;
			YK = (pR->y - pL->y) / 2.0;
			XI = (pU->x - pD->x) / 2.0;
			YI = (pU->y - pD->y) / 2.0;

			det = XK * YI - XI * YK;  // 矩阵Am的行列式

			ksi_x = YI / det;
			ksi_y = -XI / det;
			eta_x = -YK / det;
			eta_y = XK / det;
			jacob = 1.0 / det;//雅克比为矩阵Am的行列式的倒数
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
	//移动网格需要计算左右边界
	for (size_t i = 0; i < numx; i++) {
		for (size_t j = 1; j < numy - 1; j++)
		{
			//应该用邻居点来计算度量系数
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

			JJ = XK * YI - XI * YK;  // 矩阵Am的行列式

			ksi_t = (XI * YT - XT * YI) / JJ;
			ksi_x = YI / JJ;
			ksi_y = -XI / JJ;
			eta_x = -YK / JJ;
			eta_y = XK / JJ;
			eta_t = (XT * YK - XK * YT) / JJ;
			jacob = 1.0 / JJ;//雅克比为矩阵Am的行列式的倒数
			allcoefTrans_[it->id_] = CoefTrans(ksi_x, ksi_y, eta_x, eta_y, jacob, ksi_t, eta_t);
			curgrid->setCoefTrans(it, allcoefTrans_[it->id_]);
		}
	}

}

void CtrlFlow::outplot_B(std::ofstream& out, const std::string& zone_title)
{
	//step1 先把要输出的点删选出来
	auto& mgrid = gridVec_.at(1);
	const auto& mPointsvec = mgrid->getPoints();
	size_t mx = mgrid->getXnum();
	size_t my = mgrid->getYnum();

	auto& bgrid = gridVec_.at(0);
	size_t bx = bgrid->getXnum();
	size_t by = bgrid->getYnum();
	auto& bPointsVec = bgrid->getPoints();
	std::vector<CoordPtr> outPoints;

	//step 1 确定编号
	size_t LB = 0;
	size_t RB = 0;
	for (size_t i = 0; i < bx; i++) {
		if (bPointsVec[i * by + 0]->type_ == PointType::b_Lbound)
			LB = i;
		if (bPointsVec[i * by + 0]->type_ == PointType::b_Rbound)
			RB = i;
	}
	//step 2 左半边
	for (size_t i = 0; i <= LB; i++) {
		for (size_t j = 0; j < by; j++)
		{
			const auto& it = bPointsVec[i * by + j];
			outPoints.emplace_back(it);
		}
	}
	//step 3 添加movegrid中的左边界，
	for (size_t j = 0; j < my; j++)
	{
		const auto& it = mPointsvec[0 * my + j];
		outPoints.emplace_back(it);
	}
	//step 4 添加movegrid中的右边界，
	for (size_t j = 0; j < my; j++)
	{
		const auto& it = mPointsvec[(mx - 1) * my + j];
		outPoints.emplace_back(it);
	}
	//step 4 右半边，
	for (size_t i = RB; i < bx; i++) {
		for (size_t j = 0; j < by; j++)
		{
			const auto& it = bPointsVec[i * by + j];
			outPoints.emplace_back(it);
		}
	}

	//step 5 形成单元
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
	//该函数目前没用
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
				if (mgpoints[rh].type_==PointType::s_filedbound&& curpoint->x < mPointsvec[j]->x)//会误判!
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
			if (nit->id_.gridId != 0)//找到背景网格点
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
			if (nit->id_.gridId != 0)//找到背景网格点
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
	auto& back = gridVec_.at(0)->getPoints();//首先在背景点中找。
	for (auto& it : back) {
		if (GlobalData::isPointInPolygon_moreInner(it, boundpoints))
		{
			//std::cout << "在里面 !\n";
		}
		else if(!GlobalData::isPointInPolygon_moreInner(it, boundpoints)&&(it->tag_ == PointTag::shade))
		{
			if (std::fabs(it->y - gridVec_.at(0)->getYstart()) < 1E-6 || 
				std::fabs(it->y - gridVec_.at(0)->getYend()) < 1E-6)
			{
				it->tag_ = PointTag::surface;//边界点不参与计算
				//std::cout << "边界点 !\n";
			}
			else if (it->y > gridVec_.at(0)->getYstart() && it->y < gridVec_.at(0)->getYend())//边界点不参与计算
			{
				spitoutPoints_.push_back(it);
				it->tag_ = PointTag::spit;//找到以后再进行修正tag_
				//std::cout << "目标点 !\n";

			}
			else
			{
				//std::cout << "ERROR！！\n";
			}
		}
		else
		{
			//std::cout << "在外面，且点上一次没被覆盖!\n";
		}
	}

	//step2. 根据找到的spit点，找到对应在moveGrid上的连接点，因为吐出来的点就是用他们更新的。
	const auto& BPointsvec = gridVec_.at(0)->getPoints();
	auto& B_neighbors = gridVec_.at(0)->getReneighbors();
	//std::vector<PointId> tempId;
	for (const auto& perid : spitoutPoints_)
	{
		const auto& perNeighbor = B_neighbors[perid->id_];
		for (const auto& it : perNeighbor)
		{
			//查找spitoutPoint的邻居点，是否有movegrid的点。
			if (it->type_ == PointType::s_filedbound)
			{
				//tempId.push_back(it->id_);//找到点id
				//存储与spit点直接连接的，运动网格对应点上一时刻的位置。
				OldconnectPoints_.push_back(std::make_shared<Coord>(getoldCoord(it->id_)));

			}
		}
	
	}

}

void CtrlFlow::calSpitPointsCoesfTrans()
{
	//step2.为了后续计算spit点的u_new，需要计算这部分点的坐标系数。
	double ksi_t = 0, ksi_x = 0, ksi_y = 0, eta_x = 0, eta_y = 0, eta_t = 0;
	double jacob = 0.0;
	double XK = 0, XI = 0, YK = 0, YI = 0, JJ = 0, XT = 0, YT = 0;
	double dt = 0.01;
	for (size_t i = 0; i < spitoutPoints_.size(); ++i)
	{
		//用上一时刻没有被吐出来的，对应moveGrid点的邻居关系来计算度量系数
		const auto& newit = spitoutPoints_[i];
		const auto& oldit = OldconnectPoints_[i];//movegrid上移动前的位置
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

		JJ = XK * YI - XI * YK;  // 矩阵Am的行列式

		ksi_t = (XI * YT - XT * YI) / JJ;
		ksi_x = YI / JJ;
		ksi_y = -XI / JJ;
		eta_x = -YK / JJ;
		eta_y = XK / JJ;
		eta_t = (XT * YK - XK * YT) / JJ;
		jacob = 1.0 / JJ;//雅克比为矩阵Am的行列式的倒数
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

	// 如果没有找到对应的点，则抛出异常
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
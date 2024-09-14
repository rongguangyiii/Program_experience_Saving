
#include "gridGenerate/include/gridBoundS.h"


void GridBoundS::generatePoints() {
	double x_step = (x_end_ - x_start_) / (x_points_ - 1);
	double y_step = (y_end_ - y_start_) / (y_points_ - 1);
	size_t x_ctrlNum = ctrl_x_;
	size_t y_ctrlNum = ctrl_y_;
	// Control curve points generation
	for (size_t j = 0; j < y_ctrlNum; ++j) {
		size_t num = x_ctrlNum * y_points_ + j;
		points_[num] = std::make_shared<Coord>(PointId(level_, num), x_start_ + x_ctrlNum * x_step, y_start_ + j * y_step, 0);
	}
	double A = 0.02;
	double B = ctrl_y_ * y_step;
	double C = 2.0 * GlobalData::PI / ((y_end_ - y_start_) - 2 * B);
	for (size_t j = y_ctrlNum; j < y_points_ - y_ctrlNum; j++) {
		double xt = x_start_ + x_ctrlNum * x_step;
		double yy = y_start_ + j * y_step;
		//double xx = xt + A * std::sin(C * (yy - B));
		double xx = xt + 0.05 * std::sin(2.5 * GlobalData::PI * (yy - 0.1));
		size_t num = x_ctrlNum * y_points_ + j;
		points_[num] = std::make_shared<Coord>(PointId(level_, num), xx, yy, 0);
	}
	for (size_t j = y_points_ - y_ctrlNum; j < y_points_; ++j) {
		size_t num = x_ctrlNum * y_points_ + j;
		points_[num] = std::make_shared<Coord>(PointId(level_, num), x_start_ + x_ctrlNum * x_step, y_start_ + j * y_step, 0);
	}

	// Left of control line
	for (size_t i = 0; i < x_ctrlNum; i++) {
		for (size_t j = 0; j < y_points_; j++) {
			size_t temp = x_ctrlNum * y_points_ + j;
			size_t index = i * y_points_ + j;
			double xx = x_start_ + i * (points_[temp]->x - x_start_) / x_ctrlNum;
			double yy = y_start_ + j * y_step;
			points_[index] = std::make_shared<Coord>(PointId(level_, index), xx, yy, 0);
		}
	}

	// Right of control line
	size_t inum = 0;
	for (size_t i = x_ctrlNum + 1; i < x_points_; ++i) {
		inum++;
		for (size_t j = 0; j < y_points_; ++j) {
			size_t temp = x_ctrlNum * y_points_ + j;
			double xx = points_[temp]->x + inum * (x_end_ - points_[temp]->x) / ((x_points_ - x_ctrlNum) - 1);
			double yy = y_start_ + j * y_step;
			size_t index = i * y_points_ + j;
			points_[index] = std::make_shared<Coord>(PointId(level_, index), xx, yy, 0);
		}
	}
	//初始化自己的坐标，
	for (size_t i = 0; i < x_points_; i++)
	{
		for (size_t j = 0; j < y_points_; j++)
		{
			x[i][j] = points_[i * y_points_ + j]->x;
			y[i][j] = points_[i * y_points_ + j]->y;
		}
	}
}

void GridBoundS::generateElements() {
	for (size_t i = 0; i < x_points_ - 1; ++i) {
		for (size_t j = 0; j < y_points_ - 1; ++j) {
			size_t p1 = i * y_points_ + j + 1;
			size_t p2 = p1 + 1;
			size_t p3 = p1 + y_points_;
			size_t p4 = p3 + 1;
			elements_.push_back({ p1, p3, p4, p2 });
		}
	}
}

void GridBoundS::genNeiborNode()
{
	size_t i = 0, j = 0;
	size_t xNodeNum = x_points_, yNodeNum = y_points_;
	for (size_t iNode = 0; iNode < points_.size(); ++iNode)
	{
		auto& currentNode = points_.at(iNode);
		PointType currentNodeType = PointType::inner;
		std::vector<CoordPtr> neiborNode;

		i = iNode / yNodeNum;
		j = iNode % yNodeNum;
		std::vector<size_t> neiborNodeIndex = { iNode + yNodeNum, iNode + 1, iNode - yNodeNum, iNode - 1 };//设置邻居节点的索引，右，上，左，下
		if (j == 0)
		{
			neiborNodeIndex.erase(neiborNodeIndex.begin() + 3);
			currentNodeType = PointType::s_filedbound;
		}

		if (i == 0)
		{
			neiborNodeIndex.erase(neiborNodeIndex.begin() + 2);
			currentNodeType = PointType::s_filedbound;
		}

		if (j == yNodeNum - 1)
		{
			neiborNodeIndex.erase(neiborNodeIndex.begin() + 1);
			currentNodeType = PointType::s_filedbound;
		}

		if (i == xNodeNum - 1)
		{
			neiborNodeIndex.erase(neiborNodeIndex.begin());
			currentNodeType = PointType::s_filedbound;
		}

		//if ((i == xNodeNum - 1 || i == 0) && (j == 0 || j == yNodeNum - 1))
		//{
		//	currentNodeType = PointType::coner;
		//}
		for (const auto it : neiborNodeIndex)
			neiborNode.emplace_back(points_[it]);
		currentNode->type_ = currentNodeType;
		Reneighbors_[currentNode->id_] = neiborNode;

	}
}

//重载操作符=
std::shared_ptr<GridBase> GridBoundS::operator=(const std::shared_ptr<GridBase>& other) {
	points_ = other->getPoints();
	elements_ = other->getelements();
	Reneighbors_ = other->getReneighbors();
	coefTrans_ = other->getTransCoef();
	x_points_ = other->getXnum();
	y_points_ = other->getYnum();
	x_start_ = other->getXstart();
	x_end_ = other->getXend();
	y_start_ = other->getYstart();
	y_end_ = other->getYend();
	return std::make_shared<GridBoundS>(*this);
}

void GridBoundS::outputToTecplot(std::ofstream& out, const std::string& zone_title) const {
	out << "TITLE = \"" << zone_title << "\"\n";
	out << "VARIABLES = \"X\", \"Y\"\n";
	out << "ZONE T=\"" << zone_title << "\", N=" << points_.size() << ", E=" << elements_.size() << ", F=FEPOINT, ET=QUADRILATERAL\n";

	for (const auto& point : points_) {
		out << std::fixed << std::setprecision(6) << point->x << " " << point->y << "\n";
	}

	for (const auto& elem : elements_) {
		const auto& eleconn = elem.getConnectivty();
		if (eleconn.size() == 3)
		{
			out << eleconn[0] << " " << eleconn[1] << " " << eleconn[2] << " " << eleconn[2] << "\n";
		}
		else if (eleconn.size() == 4)
		{
			out << eleconn[0] << " " << eleconn[1] << " " << eleconn[2] << " " << eleconn[3] << "\n";
		}
		else
		{
			std::cout << "error in outplot!\n";
		}
	}

	std::cout << zone_title << " output complete!\n";
}

void GridBoundS::coordTrans() {
	double XK, XI, YK, YI, JJ;
	double ksi_t = 0, ksi_x = 0, ksi_y = 0, eta_x = 0, eta_y = 0, eta_t = 0;
	double jacob = 0.0;
	const size_t mx = x_points_;
	const size_t my = y_points_;
	for (size_t i = 0; i < mx; ++i) {
		for (size_t j = 0; j < my; ++j) {
			if (i == 0 && j != 0 && j != my - 1) {  // 1. 左边界 (i == 0) 非角点 (0 < j < my)
				XK = 0.5 * (4.0 * x[i + 1][j] - x[i + 2][j] - 3.0 * x[i][j]);
				YK = 0.5 * (4.0 * y[i + 1][j] - y[i + 2][j] - 3.0 * y[i][j]);
				XI = (x[i][j + 1] - x[i][j - 1]) / 2.0;
				YI = (y[i][j + 1] - y[i][j - 1]) / 2.0;
			}
			else if (i == mx - 1 && j != 0 && j != my - 1) {  // 2. 右边界 (i == mx) 非角点 (0 < j < my)
				XK = 0.5 * (-4.0 * x[i - 1][j] + x[i - 2][j] + 3.0 * x[i][j]);
				YK = 0.5 * (-4.0 * y[i - 1][j] + y[i - 2][j] + 3.0 * y[i][j]);
				XI = (x[i][j + 1] - x[i][j - 1]) / 2.0;
				YI = (y[i][j + 1] - y[i][j - 1]) / 2.0;
			}
			else if (i == 0 && j == 0) {  // 3. 左下角 (i == 0, j == 0)
				XK = 0.5 * (4.0 * x[i + 1][j] - x[i + 2][j] - 3.0 * x[i][j]);
				YK = 0.5 * (4.0 * y[i + 1][j] - y[i + 2][j] - 3.0 * y[i][j]);
				XI = 0.5 * (4.0 * x[i][j + 1] - x[i][j + 2] - 3.0 * x[i][j]);
				YI = 0.5 * (4.0 * y[i][j + 1] - y[i][j + 2] - 3.0 * y[i][j]);
			}
			else if (i == 0 && j == my - 1) {  // 4. 左上角 (i == 0, j == my)
				XK = 0.5 * (4.0 * x[i + 1][j] - x[i + 2][j] - 3.0 * x[i][j]);
				YK = 0.5 * (4.0 * y[i + 1][j] - y[i + 2][j] - 3.0 * y[i][j]);
				XI = 0.5 * (-4.0 * x[i][j - 1] + x[i][j - 2] + 3.0 * x[i][j]);
				YI = 0.5 * (-4.0 * y[i][j - 1] + y[i][j - 2] + 3.0 * y[i][j]);
			}
			else if (i == mx - 1 && j == my - 1) {  // 5. 右上角 (i == mx, j == my)
				XK = 0.5 * (-4.0 * x[i - 1][j] + x[i - 2][j] + 3.0 * x[i][j]);
				YK = 0.5 * (-4.0 * y[i - 1][j] + y[i - 2][j] + 3.0 * y[i][j]);
				XI = 0.5 * (-4.0 * x[i][j - 1] + x[i][j - 2] + 3.0 * x[i][j]);
				YI = 0.5 * (-4.0 * y[i][j - 1] + y[i][j - 2] + 3.0 * y[i][j]);
			}
			else if (i == mx - 1 && j == 0) {  // 6. 右下角 (i == mx, j == 0)
				XK = 0.5 * (-4.0 * x[i - 1][j] + x[i - 2][j] + 3.0 * x[i][j]);
				YK = 0.5 * (-4.0 * y[i - 1][j] + y[i - 2][j] + 3.0 * y[i][j]);
				XI = 0.5 * (4.0 * x[i][j + 1] - x[i][j + 2] - 3.0 * x[i][j]);
				YI = 0.5 * (4.0 * y[i][j + 1] - y[i][j + 2] - 3.0 * y[i][j]);
			}
			else if (j == 0) {  // 7. 下边界 (j == 0) 非角点 (0 < i < mx)
				XI = 0.5 * (4.0 * x[i][j + 1] - x[i][j + 2] - 3.0 * x[i][j]);
				YI = 0.5 * (4.0 * y[i][j + 1] - y[i][j + 2] - 3.0 * y[i][j]);
				XK = (x[i + 1][j] - x[i - 1][j]) / 2.0;
				YK = (y[i + 1][j] - y[i - 1][j]) / 2.0;
			}
			else if (j == my - 1) {  // 8. 上边界 (j == my) 非角点 (0 < i < mx)
				XI = 0.5 * (-4.0 * x[i][j - 1] + x[i][j - 2] + 3.0 * x[i][j]);
				YI = 0.5 * (-4.0 * y[i][j - 1] + y[i][j - 2] + 3.0 * y[i][j]);
				XK = (x[i + 1][j] - x[i - 1][j]) / 2.0;
				YK = (y[i + 1][j] - y[i - 1][j]) / 2.0;
			}
			else {  // 9. 内部点
				XK = (x[i + 1][j] - x[i - 1][j]) / 2.0;
				YK = (y[i + 1][j] - y[i - 1][j]) / 2.0;
				XI = (x[i][j + 1] - x[i][j - 1]) / 2.0;
				YI = (y[i][j + 1] - y[i][j - 1]) / 2.0;
			}

			JJ = XK * YI - XI * YK;  // 雅可比行列式

			ksi_x = YI / JJ;
			ksi_y = -XI / JJ;
			eta_x = -YK / JJ;
			eta_y = XK / JJ;
			jacob = 1.0 / JJ;
			coefTrans_[points_[i * y_points_ + j]->id_] = CoefTrans(ksi_x, ksi_y, eta_x, eta_y, jacob);
		}
	}
}

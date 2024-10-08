#include "grid/include/gridSin.h"
#include "grid/include/globalData.h"

void GridBoundS::generatePoints() {
	double x_step = (x_end_ - x_start_) / (x_points_ - 1);
	double y_step = (y_end_ - y_start_) / (y_points_ - 1);
	size_t x_ctrlNum = (x_points_ - 1) / 2;
	size_t y_ctrlNum = ctrl_y_;
	// Control curve points generation
	for (size_t j = 0; j < y_ctrlNum; ++j) {
		size_t num = x_ctrlNum * y_points_ + j;
		points_[num] = std::make_shared<Coord>(PointId(level_, num), x_start_ + x_ctrlNum * x_step, y_start_ + j * y_step, 0);
	}
	double A = 2.5 * x_step;
	double B = ctrl_y_ * y_step;
	double C = 2.0 * GlobalData::PI / ((y_end_ - y_start_) - 2 * B);
	for (size_t j = y_ctrlNum; j < y_points_ - y_ctrlNum; j++) {
		double xt = x_start_ + x_ctrlNum * x_step;
		double yy = y_start_ + j * y_step;
		double xx = xt + A * std::sin(C * (yy - B));
		//double xx = xt + 0.05 * std::sin(2.5 * GlobalData::PI * (yy - 0.1));
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

void GridBoundS::generatePoints_1() {
	double x_step = (x_end_ - x_start_) / (x_points_ - 1);
	double y_step = (y_end_ - y_start_) / (y_points_ - 1);
	size_t x_midNum = (x_points_ - 1) / 2;//中间曲线坐在点编号
	size_t y_ctrlNum = ctrl_y_;
	// Control curve points generation
	for (size_t j = 0; j < y_ctrlNum; ++j) {
		size_t num = x_midNum * y_points_ + j;
		points_[num] = std::make_shared<Coord>(PointId(level_, num), x_start_ + x_midNum * x_step, y_start_ + j * y_step, 0);
	}
	double A = 2.5 * x_step;
	double B = ctrl_y_ * y_step;
	double C = 2.0 * GlobalData::PI / ((y_end_ - y_start_) - 2 * B);
	for (size_t j = y_ctrlNum; j < y_points_ - y_ctrlNum; j++) {
		double xt = x_start_ + x_midNum * x_step;
		double yy = y_start_ + j * y_step;
		double xx = xt + A * std::sin(C * (yy - B));
		//double xx = xt + 0.05 * std::sin(2.5 * GlobalData::PI * (yy - 0.1));
		size_t num = x_midNum * y_points_ + j;
		points_[num] = std::make_shared<Coord>(PointId(level_, num), xx, yy, 0);
	}
	for (size_t j = y_points_ - y_ctrlNum; j < y_points_; ++j) {
		size_t num = x_midNum * y_points_ + j;
		points_[num] = std::make_shared<Coord>(PointId(level_, num), x_start_ + x_midNum * x_step, y_start_ + j * y_step, 0);
	}

	// Left of control line
	for (size_t i = 0; i < ctrl_x_; i++) {
		for (size_t j = 0; j < y_points_; j++) {
			size_t index = i * y_points_ + j;
			double xx = x_start_ + i * x_step;
			double yy = y_start_ + j * y_step;
			points_[index] = std::make_shared<Coord>(PointId(level_, index), xx, yy, 0);
		}
	}
	size_t it = 0;
	double leftstart = x_start_ + ctrl_x_ * x_step;
	for (size_t i = ctrl_x_; i < x_midNum; i++) {
		for (size_t j = 0; j < y_points_; j++) {
			size_t temp = x_midNum * y_points_ + j;
			size_t index = i * y_points_ + j;
			double xx = leftstart + it * (points_[temp]->x - leftstart) / (x_midNum - ctrl_x_);
			double yy = y_start_ + j * y_step;
			points_[index] = std::make_shared<Coord>(PointId(level_, index), xx, yy, 0);
		}
		it++;
	}

	// Right of control line
	size_t inum = 0;
	double rightend = x_end_ - ctrl_x_ * x_step;
	for (size_t i = x_midNum + 1; i < x_points_ - ctrl_x_; ++i) {
		inum++;
		for (size_t j = 0; j < y_points_; ++j) {
			size_t temp = x_midNum * y_points_ + j;
			double xx = points_[temp]->x + inum * (rightend - points_[temp]->x) / ((x_points_ - x_midNum - ctrl_x_) - 1);
			double yy = y_start_ + j * y_step;
			size_t index = i * y_points_ + j;
			points_[index] = std::make_shared<Coord>(PointId(level_, index), xx, yy, 0);
		}
	}
	for (size_t i = x_points_ - ctrl_x_; i < x_points_; i++) {
		for (size_t j = 0; j < y_points_; j++) {
			size_t index = i * y_points_ + j;
			double xx = x_start_ + i * x_step;
			double yy = y_start_ + j * y_step;
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

void GridBoundS::generatePoints_2() {
	double x_step = (x_end_ - x_start_) / (x_points_ - 1);
	double y_step = (y_end_ - y_start_) / (y_points_ - 1);
	size_t y_midNum = (y_points_ - 1) / 2;//中间曲线坐在点编号
	size_t y_ctrlNum = ctrl_y_;
	size_t x_ctrlNum = ctrl_x_;
	// Control curve points generation
	for (size_t i = 0; i < x_ctrlNum; ++i) {
		size_t num = i * y_points_ + y_midNum;
		points_[num] = std::make_shared<Coord>(PointId(level_, num), x_start_ + i * x_step, y_start_ + y_midNum * y_step, 0);
	}
	double A = 2.5 * y_step;
	double B = ctrl_x_ * x_step;
	double C = 2.0 * GlobalData::PI / ((x_end_ - x_start_) - 2 * B);
	for (size_t i = x_ctrlNum; i < x_points_ - x_ctrlNum; i++) {
		double yt = y_start_ + y_midNum * x_step;
		double xx = x_start_ + i * x_step;
		double yy = yt + A * std::sin(C * (xx - B));
		size_t num = i * y_points_ + y_midNum;
		points_[num] = std::make_shared<Coord>(PointId(level_, num), xx, yy, 0);
	}
	for (size_t i = x_points_ - x_ctrlNum; i < y_points_; ++i) {
		size_t num = i * y_points_ + y_midNum;
		points_[num] = std::make_shared<Coord>(PointId(level_, num), x_start_ + i * x_step, y_start_ + y_midNum * y_step, 0);
	}

	for (size_t i = 0; i < x_points_; i++)
	{
		// Down of control line
		for (size_t j = 0; j < ctrl_y_; j++) {
			size_t index = i * y_points_ + j;
			double xx = x_start_ + i * x_step;
			double yy = y_start_ + j * y_step;
			points_[index] = std::make_shared<Coord>(PointId(level_, index), xx, yy, 0);
		}
		size_t it = 0;
		double downstart = y_start_ + ctrl_y_ * y_step;
		for (size_t j = ctrl_y_; j < y_midNum; j++) {
			size_t temp = i * y_points_ + y_midNum;
			size_t index = i * y_points_ + j;
			double yy = downstart + it * (points_[temp]->y - downstart) / (y_midNum - ctrl_y_);
			double xx = x_start_ + i * x_step;
			points_[index] = std::make_shared<Coord>(PointId(level_, index), xx, yy, 0);
			it++;
		}

		// Up of control line
		size_t inum = 0;
		double upend = y_end_ - ctrl_y_ * y_step;
		for (size_t j = y_midNum + 1; j < y_points_ - ctrl_y_; j++) {
			inum++;
			size_t temp = i * y_points_ + y_midNum;
			double yy = points_[temp]->y + inum * (upend - points_[temp]->y) / ((y_points_ - y_midNum - ctrl_y_) - 1);
			double xx = x_start_ + i * x_step;
			size_t index = i * y_points_ + j;
			points_[index] = std::make_shared<Coord>(PointId(level_, index), xx, yy, 0);
		}
		for (size_t j = y_points_ - ctrl_y_; j < y_points_; j++) {
			size_t index = i * y_points_ + j;
			double xx = x_start_ + i * x_step;
			double yy = y_start_ + j * y_step;
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
	//coefTrans_ = other->getTransCoef();
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

void GridBoundS::outputToTecplot_VCFEM(std::ofstream& out, const std::string& zone_title) const {
	out << "TITLE = \"" << zone_title << "\"\n";
	out << "VARIABLES = \"X\", \"Y\"\n";
	out << "ZONE T=\"" << zone_title << "\", N=" << points_.size() << ", E=" << elements_.size() << ", F=FEPOINT, ET=QUADRILATERAL\n";

	//1. CONNECTIVITY
	out << "CONNECTIVITY \n";
	for (size_t i = 0; i < elements_.size(); i++)
	{
		const auto& eleconn = elements_[i].getConnectivty();
		const size_t num = eleconn.size();
		if (eleconn.size() == 3) {
			std::cout << "error!!!\n";
		}
		else if (eleconn.size() == 4) {
			out << i + 1 << " " << num << " " << eleconn[0] << " " << eleconn[1] << " " << eleconn[2] << " " << eleconn[3] << " " << -1 << " " << 1 << "\n";
		}
		else {
			std::cout << "error in outplot!\n";
		}
	}
	//2. COORDINATES
	out << "COORDINATES \n";
	size_t nn = 0;
	for (const auto& point : points_) {
		nn++;
		out << std::fixed << std::setprecision(6) << nn << " " << point->x << " " << point->y << "\n";
	}
	//3. INCLUSIONS
	out << "INCLUSIONS \n";
	for (size_t i = 0; i < elements_.size(); i++) {

		out << i + 1 << " " << 0 << "\n";
	}

	out << "INFORMATION AND VALUES OF PRESCRIBED INCREMENT DISPLACEMENT \n";

	//左边x方向固定u=0-------------------------------------------------------------
	out << 1 << ",    " << 1 << ",    " << y_points_ << "\n";
	out << 0 << "\n";
	for (size_t i = 0; i < y_points_; i++) {
		out << i + 1 << ",    ";
	}
	out << std::endl;
	//左边y方向固定一个点v=0-------------------------------------------------------------
	out << 2 << ",    " << 2 << ",    " << 1 << "\n";
	out << 0 << "\n";
	out << 1 << "\n";
	//左边x方向均匀拉伸u=0.01-------------------------------------------------------------
	out << 3 << ",    " << 1 << ",    " << y_points_ << "\n";
	out << 0.01 << "\n";
	size_t sta = (x_points_ - 1) * y_points_;
	for (size_t i = 0; i < y_points_; i++) {
		out << sta + i + 1 << ",    ";
	}

	std::cout << zone_title << " output complete!\n";
}

void GridBoundS::coefTrans() {
	coeftransVec_.clear();
	double ksi_x = 0, ksi_y = 0, eta_x = 0, eta_y = 0, jacob = 0;
	double XK = 0, XI = 0, YK = 0, YI = 0, JJ = 0;
	const size_t mx = x_points_;
	const size_t my = x_points_;
	size_t index = 0;
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
			index = i * y_points_ + j;
			JJ = XK * YI - XI * YK;  // 雅可比行列式
			ksi_x = YI / JJ;
			ksi_y = -XI / JJ;
			eta_x = -YK / JJ;
			eta_y = XK / JJ;
			jacob = 1.0 / JJ;
			coeftransVec_.emplace_back(CoefTrans(ksi_x, ksi_y, eta_x, eta_y, jacob));
		}
	}
}
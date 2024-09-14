#pragma once
#include <string>
#include <unordered_map>
#include <set>
#include "gridGenerate/include/coord.h"
#include "gridGenerate/include/gridBase.h"
#include "flowSolve/include/globalData.h"


class GridBoundS : public GridBase {
public:
    GridBoundS(double x_start, double x_end, double y_start, double y_end, size_t x_points, size_t y_points, size_t ctrlx, size_t ctrly)
        : x_start_(x_start), x_end_(x_end), y_start_(y_start), y_end_(y_end), level_(1),
        x_points_(x_points), y_points_(y_points), ctrl_x_(ctrlx), ctrl_y_(ctrly), x(x_points_, std::vector<double>(y_points_, 0.0)),
        y(x_points_, std::vector<double>(y_points_, 0.0))
    {
        points_.resize(x_points_ * y_points_, std::make_shared<Coord>(0.0,0.0));
        generatePoints();
        generateElements();
        genNeiborNode();
    }
    ~GridBoundS() override = default;
    void generatePoints() override;
    void generateElements() override;
    void genNeiborNode()override;
    void coordTrans()override;

    //ÖØÔØ²Ù×÷·û=
    std::shared_ptr<GridBase> operator=(const std::shared_ptr<GridBase>& other) override;

    void outputToTecplot(std::ofstream& out, const std::string& zone_title) const override;

    size_t getXnum() const { return x_points_; }
    size_t getYnum() const { return y_points_; }
    double getXstart() const { return x_start_; }
    double getXend() const { return x_end_; }
    double getYstart() const { return y_start_; }
    double getYend() const { return y_end_; }
private:
    double x_start_, x_end_, y_start_, y_end_;
    size_t x_points_, y_points_, ctrl_x_, ctrl_y_,level_;
    std::vector<std::vector<double>> x;
    std::vector<std::vector<double>> y;
};

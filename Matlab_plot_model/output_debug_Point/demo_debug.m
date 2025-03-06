clear;clc;
colors = [ 0 0 0;       % 黑色1  
    1 0 0;       % 红色2  
    0 0 1;       % 蓝色3 
    0 1 0;       % 绿色4
    1 1 0;       % 黄色5
    1 0.5 0.5;   % 粉色（浅红色）6  
    0.5 0 0.5;   % 紫色7
    0 1 1;       % 青色8
    0.5 0.8 1;   % 天蓝色9 
    0.5 0.25 0;  % 棕色10
    0.5 0.5 0.5;]; % 灰色11
markerStyles = {'.';% 1 
    'o';%圆圈2
    '*';%星号3
    '+';%加号4
    '^';%上三角5
    'p';%五边形6
    'hexagon';%六边形7
    'pentagram';%五角星8
    'square';%方块8
    'diamond';};%菱形9  
fig_handle = figure;
isGenGrid=true;%是否生成网格点，默认生成
isflag=false;%是否标记网格点顺序，默认不标记
%----------------------------------------DEBUG:
datafilename_debug=["Data_Backgrid_grid.dat";"Data_CurvePC_x3_grid.dat";
                    "Data_Curve_x3_grid.dat";"Data_Curve3_Hole_grid.dat";
                     "Data_Curve2_divideorder_grid.dat";
                    "Data_SolidNode_check_grid.dat"; "Data_OuterPC002_grid.dat"; 
                    "Data_SolidNode002_grid.dat"; "Data_HoleNode002_grid.dat"; ];

fig_handle = ViewGrid(fig_handle, datafilename_debug(1), colors(1,:), markerStyles{1}, isGenGrid, isflag);%背景网格点
fig_handle = ViewGrid(fig_handle, datafilename_debug(2), colors(4,:), markerStyles{1}, isGenGrid, isflag);%模型点云
fig_handle = ViewGrid(fig_handle, datafilename_debug(3), colors(2,:), markerStyles{3}, isGenGrid, isflag);
fig_handle = ViewGrid(fig_handle, datafilename_debug(4), colors(3,:), markerStyles{6}, isGenGrid, isflag);
fig_handle = ViewGrid(fig_handle, datafilename_debug(5), colors(9,:), markerStyles{2}, ~isGenGrid, ~isflag);
%附面层最外层多边形点云：Data_OuterPC002_grid
fig_handle = ViewGrid(fig_handle, datafilename_debug(6), colors(4,:), markerStyles{1}, ~isGenGrid, isflag);
%第二次标记的固体点：Data_SolidNode002_grid_0
fig_handle = ViewGrid(fig_handle, datafilename_debug(7), colors(3,:), markerStyles{3}, ~isGenGrid, isflag);
%背景网格洞边界：Data_HoleNode002_grid_0
fig_handle = ViewGrid(fig_handle, datafilename_debug(8), colors(6,:), markerStyles{8}, ~isGenGrid, isflag);



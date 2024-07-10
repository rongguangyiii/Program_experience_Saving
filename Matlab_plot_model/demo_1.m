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

datafilename=["Data_Backgrid_grid_0.dat"; "Data_modelCloud_grid_.dat"; "Data_BoundNode_grid_1.dat";...
              "Data_BoundLayerNode_grid_1.dat"; "Data_SolidNode002_grid_0.dat"; 'Data_HoleNode001_grid_0.dat';...
              %---------------------DEBUG
             'Data_HoleNode002_grid_0.dat';'Data_HoleNode_grid_1.dat'; 'Data_HoleNode_grid_74.dat';];
fig_handle = figure;
isGenGrid=true;%是否生成网格点，默认生成
isflag=false;%是否标记网格点顺序，默认不标记
fig_handle = ViewGrid(fig_handle, datafilename(1), colors(1,:), markerStyles{1}, isGenGrid, isflag);%背景网格点
fig_handle = viewMode(fig_handle, datafilename(2), colors(9,:), markerStyles{1}, isGenGrid);%模型点云
fig_handle = ViewGrid(fig_handle, datafilename(3), colors(6,:), markerStyles{2}, ~isGenGrid, isflag);%物面边界点
fig_handle = ViewGrid(fig_handle, datafilename(4), colors(10,:), markerStyles{6}, ~isGenGrid, isflag);%边界层节点
fig_handle = ViewGrid(fig_handle, datafilename(5), colors(7,:), markerStyles{3}, ~isGenGrid, isflag);%固体节点
fig_handle = ViewGrid(fig_handle, datafilename(6), colors(3,:), markerStyles{8}, ~isGenGrid, isflag);%第一次背景网格洞边界
fig_handle = ViewGrid(fig_handle, datafilename(7), colors(4,:), markerStyles{8}, ~isGenGrid, ~isflag);%第二次背景网格洞边界
fig_handle = ViewGrid(fig_handle, datafilename(8), colors(6,:), markerStyles{2}, ~isGenGrid, isflag);%非结构网格洞边界


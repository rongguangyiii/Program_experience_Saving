% 定义圆心坐标和半径  
center_x = 0.5;  
center_y = 0.0;  
radius = 0.06;  
num_points = 61; % 定义点的数量  
  
% 计算圆上100个点的坐标  
theta = linspace(0, 2*pi, num_points); % 生成100个从0到2*pi的等差角  
x = center_x + radius * cos(theta); % 根据圆的标准方程计算x坐标  
y = center_y + radius * sin(theta); % 根据圆的标准方程计算y坐标  
  
% 将坐标写入coor.dat文件  
fileID = fopen('coor.dat', 'w');  
if fileID == -1  
    error('无法打开文件coor.dat进行写入');  
end  
  
% 写入前缀信息  
fprintf(fileID, 'Poly\n'); % 写入第一行：Poly  
fprintf(fileID, '%d\n', num_points); % 修正：写入第二行：点数和具体的数量  
  
% 遍历每个点，将坐标写入文件  
for i = 1:num_points  
    fprintf(fileID, '%.12f %.12f\n', x(i), y(i)); % 写入一个点的x和y坐标，然后换行  
end  
  
fclose(fileID); % 关闭文件  
  
% 绘制圆  
figure; % 创建一个新的图形窗口  
plot(x, y, 'b.'); % 使用蓝色点绘制圆  
axis equal; % 设置x轴和y轴具有相同的刻度，使得圆看起来是圆的  
grid on; % 显示网格  
title('Circle with Center (0.5, 0.0) and Radius 0.2'); % 设置图形标题  
xlabel('X Coordinate'); % 设置x轴标签  
ylabel('Y Coordinate'); % 设置y轴标签
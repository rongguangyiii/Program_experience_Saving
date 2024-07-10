clear;clc;
% 定义y的范围
pointnum = 101;
y = linspace(-0.5, 0.5, pointnum); % 使用linspace创建从-3到3的1000个点

% 计算每个函数对应的x值
x0 = 1.0 * y .^ 2 - 1.2; % x=1y^2-1.2
x1 = 2.0 * y .^ 2 - 0.9; % x=2y^2-0.9
x2 = 4.0 * y .^ 2 - 0.6; % x=4y^2-0.6
x3 = 6.0 * y .^ 2 - 0.3; % x=6y^2-0.3

% 绘制图形  
figure; % 创建一个新的图形窗口  
plot(x0, y, 'c.', 'LineWidth', 2); 
hold on; % 保持当前图形，以便在同一张图上绘制其他函数  
plot(x1, y, 'g.', 'LineWidth', 2); 
plot(x2, y, 'r.', 'LineWidth', 2);  
plot(x3, y, 'b-', 'LineWidth', 2);   

% 添加图例  
legend('x0=1y^2-1.2', 'x1=2y^2-0.9', 'x2=4y^2-0.6','x3=6y^2-0.3', 'Location', 'NorthWest');  

% 显示网格（可选）  
grid on;  
  
% 释放hold状态（可选，但建议在使用hold on后总是使用hold off）  
hold off;

% 定义输出文件名
outputFileName = 'mode.input';
% 打开文件准备写入
fid = fopen(outputFileName, 'w');
% 写入数据
fprintf(fid, 'Curve : x0=1y^2-1.2\n');
fprintf(fid, '%d\n', pointnum);
for i = 1:length(y)
    fprintf(fid, '%f %f\n', x0(i), y(i));
end

fprintf(fid, 'Curve : x1=2y^2-0.9\n');
fprintf(fid, '%d\n', pointnum);
for i = 1:length(y)
    fprintf(fid, '%f %f\n', x1(i), y(i));
end

fprintf(fid, 'Curve : x2=4y^2-0.6\n');
fprintf(fid, '%d\n', pointnum);
for i = 1:length(y)
    fprintf(fid, '%f %f\n', x2(i), y(i));
end

fprintf(fid, 'Curve : x3=6y^2-0.3\n');
fprintf(fid, '%d\n', pointnum);
for i = 1:length(y)
    fprintf(fid, '%f %f\n', x3(i), y(i));
end
% 关闭文件
fclose(fid);

% 输出成功提示
disp('Data written to file successfully.');


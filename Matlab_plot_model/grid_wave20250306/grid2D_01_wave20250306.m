% 二维网格的波浪动画
clear all;
close all;

% 初始化参数
n = 36;              % 每个方向的点数
x = linspace(0, 4, n);  % x坐标范围0到4
y = linspace(0, 4, n);  % y坐标范围0到4
[X, Y] = meshgrid(x, y); % 创建初始均匀网格

% 初始化高度Z
Z = zeros(size(X));

% 创建图形窗口
figure;
h = surf(X, Y, Z, 'EdgeColor', 'b', 'FaceColor', 'none');
axis([0 4 0 4 -2 2]); % 设置坐标轴范围
grid on;
title('二维网格波浪动画');
xlabel('X轴');
ylabel('Y轴');
zlabel('Z轴');
view(-37.5, 30); % 设置视角

% 动画参数
t = 0;              % 时间参数
dt = 0.1;          % 时间步长
amplitude = 1;     % 波浪幅度
freq_x = 1;        % x方向波浪频率
freq_y = 1;        % y方向波浪频率

% 动画循环
while true
    % 计算新的Z坐标，x和y方向都添加波浪
    Z = amplitude * (sin(freq_x * X + t) + sin(freq_y * Y + t));
    
    % 更新图形
    set(h, 'ZData', Z);
    drawnow;
    
    % 更新时间
    t = t + dt;
    
    % 暂停一小段时间，控制动画速度
    pause(0.05);
end
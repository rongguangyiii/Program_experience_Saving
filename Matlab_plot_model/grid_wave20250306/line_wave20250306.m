% 创建36个点的直线并实现波浪动画
clear all;
close all;

% 初始化参数
n = 36;              % 点的数量
x = 1:n;            % x坐标 (1到36)
y = zeros(1,n);     % 初始y坐标 (直线)

% 创建图形窗口
figure;
h = plot(x, y, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
axis([0 n+1 -2 2]); % 设置坐标轴范围
grid on;
title('波浪动画');
xlabel('X轴');
ylabel('Y轴');

% 动画参数
t = 0;              % 时间参数
dt = 0.1;          % 时间步长
amplitude = 1;     % 波浪幅度
frequency = 0.2;   % 波浪频率

% 动画循环
while t<40
    % 计算新的y坐标，使用正弦函数创建波浪效果
    y = amplitude * sin(frequency * x + t);
    
    % 更新图形
    set(h, 'YData', y);
    drawnow;
    
    % 更新时间
    t = t + dt;
    
    % 暂停一小段时间，控制动画速度
    pause(0.05);
end
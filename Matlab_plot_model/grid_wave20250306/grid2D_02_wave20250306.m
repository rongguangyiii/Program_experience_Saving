% MATLAB 代码：生成带有波浪运动的网格并保存为视频
clear all;
close all;

% 定义网格范围和点数
x = linspace(0, 10, 36); % x方向21个点，从0到10
y = linspace(0, 10, 36); % y方向21个点，从0到10
[X, Y] = meshgrid(x, y); % 创建初始均匀正交网格

% 添加波浪运动参数
amplitude = 0.5; % 波浪幅度
frequency = 1;   % 波浪频率

% 创建视频文件
video = VideoWriter('wavy_grid.mp4', 'MPEG-4'); % 创建视频对象，保存为mp4格式
video.FrameRate = 20; % 设置帧率（每秒20帧）
open(video); % 打开视频文件


fig = figure('Position', [0, 0, 1920, 1080]); % 窗口大小设为1280x720像素

% 动画循环
for t = 0:0.1:10
    % 对x和y方向添加波浪扰动
    X_wave = X + amplitude * sin(frequency * Y + t); % x方向波浪
    Y_wave = Y + amplitude * sin(frequency * X + t); % y方向波浪
    
    % 绘制图形
    figure(1);
    clf; % 清除当前图形
    plot(X_wave, Y_wave, 'b-', X_wave', Y_wave', 'b-'); % 绘制波浪网格
    axis equal; % 保持x和y轴比例一致
    axis([-1 11 -1 11]); % 设置坐标轴范围
    grid on;
    title(['波浪网格 (t = ' num2str(t, '%.1f') ')']);
    xlabel('X轴');
    ylabel('Y轴');
    
    % 获取当前帧并写入视频
    frame = getframe(gcf); % 捕获当前图形作为帧
    writeVideo(video, frame); % 写入视频文件
    
    drawnow; % 更新显示
end

% 关闭视频文件
close(video);
close all;

disp('视频已保存为 wavy_grid.mp4');
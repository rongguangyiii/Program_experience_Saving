
function fig_handle = ViewGrid(fig_handle,datafilename,colors,marktag,isGenGrid,flag)% 文件名
if(~isGenGrid) 
    return;
end
% filename = 'D:\ZARAN3\bin\Debug\';
filename = 'D:\SF_dev\Ver_ShockfitMoveline\out\Debug\tempData\'+datafilename;
fileID = fopen(filename, 'r');
if fileID == -1
    error('无法打开文件 %s', filename);
end

% 读取前三行文件说明
fgetl(fileID); % 第一行
fgetl(fileID); % 第二行
numPoints=fscanf(fileID, '%d',1);
points = zeros(numPoints, 2);
% 读取点的坐标
for i = 1:numPoints
    points(i,:) = fscanf(fileID, '%f %f',2);
end
figure(fig_handle);
plot(points(:,1), points(:,2),marktag,'Color', colors);

if(flag)
    for i = 1:numPoints
        text(points(i,1),points(i,2), sprintf('%d', i-1), ...
            'HorizontalAlignment', 'right', ... % 水平对齐方式，根据需要调整
            'VerticalAlignment', 'bottom', ... % 垂直对齐方式，根据需要调整
            'FontSize', 8); % 字体大小，可以根据需要调整
    end
end
title('Model');
xlabel('X');
ylabel('Y');
grid on; % 打开网格
hold on;
% 关闭文件
fclose(fileID);
end


function fig_handle = viewMode(fig_handle,datafilename,colors,marktag,isGenGrid)
if(~isGenGrid) 
    return;
end
 filename = 'D:\ZARAN3\model\'+datafilename;  
%filename = 'D:\ZARAN3\bin\Debug\tempData\'+datafilename;   
fileID = fopen(filename, 'r');  
if fileID == -1  
    error('无法打开文件 %s', filename);  
end  
  
numModels = fscanf(fileID, '%d');  
modelDescription = fgetl(fileID); 
numPoints = fscanf(fileID, '%d', 1);  
fprintf('模型数: %d\n', numModels);  
fprintf('模型说明: %s\n', modelDescription); 
fprintf('点数: %d\n', numPoints);    

points = zeros(numPoints, 2);  
   
for i = 1:numPoints  
    points(i,:) = fscanf(fileID, '%f %f',2);  
end  
fclose(fileID);  
  
% 绘制点  
figure(fig_handle); % 创建一个新的图形窗口  
% plot(points(:,1), points(:,2),'Color', colors,'Marker',marktag); % 标记点和画线连起来
% plot(points(:,1), points(:,2),'o','MarkerFaceColor', colors);%MarkerFaceColor表示标记的线框颜色。
plot(points(:,1), points(:,2),marktag,'Color', colors);%Color表示标记本身的颜色
title('Model');  
xlabel('X');  
ylabel('Y');  
grid on; % 打开网格
hold on;
end


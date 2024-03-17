clear;clc;


AA = xlsread('E:\研究生二年级工作安排\动态器件\Dynamic-memristor-based-reservoir-computing-v1.0.0\Stuka255-Dynamic-memristor-based-reservoir-computing-9193344\频率湿度平滑\第一次尝试.xlsx','用来画图')% 数据表\sheet\行列
% B = xlsread('E:\研究生二年级工作安排\动态器件\Dynamic-memristor-based-reservoir-computing-v1.0.0\Stuka255-Dynamic-memristor-based-reservoir-computing-9193344\dt-湿度.xlsx','sheet8')% 数据表\sheet\行列
% C = xlsread('E:\研究生二年级工作安排\动态器件\Dynamic-memristor-based-reservoir-computing-v1.0.0\Stuka255-Dynamic-memristor-based-reservoir-computing-9193344\dt-湿度.xlsx','sheet9')% 数据表\sheet\行列
A=AA(1:6,1:51);
B=AA(8:13,1:51);
C=AA(15:20,1:51);
s = 0.4; % 柱子宽度
% map = TheColor('dream',1);
figure(4)
yvalues = {'0.01','0.1','1','5','10','50','100','500','1000'};
% xvalues = {'10000','1000','200','100','20','10','5','2','1','0.1'};
% yvalues = {'0.3','0.5','0.7','0.9','1.1','1.3','1.5','1.7','1.9'};
xvalues = {'10000','1000','200','100','20','10','5','2','1','0.1'};
% yvalues = {'0.1','1','5','10','20','30','40','50','100'};
h = heatmap(xvalues,yvalues,A);
set(gcf,'color','white');
h.Title = 'good simu(NRMSE<0.02)';
% h.YLabel = 'humidity';
h.XLabel = 'sample rate(point/sec)';
h.YLabel = 'current amplify';
colormap summer
h.FontName = 'Arial';
h.FontSize = 12;
% c = colorbar;
% c.Label.String = 'Machine learning effects(78 experiments)';


figure(5)
yvalues = {'0.01','0.1','1','5','10','50','100','500','1000'};
% xvalues = {'10000','1000','200','100','20','10','5','2','1','0.1'};
% yvalues = {'0.3','0.5','0.7','0.9','1.1','1.3','1.5','1.7','1.9'};
xvalues = {'10000','1000','200','100','20','10','5','2','1','0.1'};
% yvalues = {'0.1','1','5','10','20','30','40','50','100'};
h = heatmap(xvalues,yvalues,B);
set(gcf,'color','white');
h.Title = 'good simu with fluctuation(NRMSE<0.1)';
% h.XLabel = 'Sample rate(point/sec)';
h.XLabel = 'sample rate(point/sec)';
h.YLabel = 'current amplify';
% h.YLabel = 'humidity';
colormap summer
% c = colorbar;
% c.Label.String = 'Machine learning effects(78 experiments)';
h.FontName = 'Arial';
h.FontSize = 12;


figure(6)
yvalues = {'0.01','0.1','1','5','10','50','100','500','1000'};
% xvalues = {'10000','1000','200','100','20','10','5','2','1','0.1'};
% yvalues = {'0.3','0.5','0.7','0.9','1.1','1.3','1.5','1.7','1.9'};
xvalues = {'10000','1000','200','100','20','10','5','2','1','0.1'};
h = heatmap(xvalues,yvalues,C);
set(gcf,'color','white');
h.Title = 'recognisable(NRMSE<0.3)';
% h.XLabel = 'Sample rate(point/sec)';
h.XLabel = 'sample rate(point/sec)';
h.YLabel = 'current amplify';
% h.YLabel = 'humidity';
colormap summer
% c = colorbar;
% c.Label.String = 'Machine learning effects(78 experiments)';
h.FontName = 'Arial';
h.FontSize = 12;

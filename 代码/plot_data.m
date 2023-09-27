clear;clc;close all;
%% 复杂度比较
% 输入序列的长度 N
t=1:7;
N=2.^t;

% DFT 算法的复杂度 C
C_DFT = N.^2;
% FFT 算法的复杂度 C，这里假设为 O(N log N)
C_FFT = N .* log2(N);
figure1=figure;
% 绘制折线图
plot(N, C_DFT, '-o', 'LineWidth',2,'Color','#6699CC');  % 自定义颜色为淡蓝色
hold on;
plot(N, C_FFT, '-o', 'LineWidth',2,'Color','#CC3333');    % 自定义颜色为深蓝色
hold off;

% 添加标签和标题
xlabel('输入序列的长度 (N)');
ylabel('复杂度 (C)');
title('DFT vs FFT 算法的复杂度');

% 添加图例
legend('DFT', 'FFT');
saveas(figure1, 'complexity_plot.png');

%% 问题一复杂度
filename = 'C_results.txt'; % 文件名
fileID = fopen(filename, 'r'); % 打开文件以读取模式

% 读取文件中的数据
data = textscan(fileID, 't=%d, RMSE=%.6f, C=%f');
t_values = data{1}; % 提取t的值
RMSE_values=data{2};% 提取RMSE
C_values = data{3}; % 提取C的值
fclose(fileID);

figure2=figure;
plot(t_values,C_values,'o-','LineWidth',2)
xlabel('t');
ylabel('C');
title('q=16,复杂度C=ql随t的变化');
xticks(t_values);
grid on;
% 在图上显示每个数据点的数值
for i = 1:length(t_values)
    text(t_values(i), C_values(i), sprintf('%.2f', C_values(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
saveas(figure2, '复杂度C=ql随t的变化.png');

figure3=figure;
plot(t_values,RMSE_values,'o-','LineWidth',2)
xlabel('t');
ylabel('RMSE');
title('最小误差RMSE随t的变化');
xticks(t_values);
grid on
% 在图上显示每个数据点的数值
for i = 1:length(t_values)
    text(t_values(i), RMSE_values(i), sprintf('%f', RMSE_values(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
saveas(figure3, '最小误差RMSE随t的变化.png');

%% 问题二 t-iter-bestFitness

iter=1:2000;
data=xlsread('results.xlsx',1,'B12:F2011');
figure4=figure;
subplot(3,2,1);
plot(iter,data(:,1),'Color','#1450A3','LineWidth',1.8);legend('t=1');
hold on
subplot(3,2,2);
plot(iter,data(:,2),'Color','#337CCF','LineWidth',1.8);legend('t=2');
subplot(3,2,3);
plot(iter,data(:,3),'Color','#FFC436','LineWidth',1.8);legend('t=3');
subplot(3,2,4);
plot(iter,data(:,4),'Color','#4682A9','LineWidth',1.8);legend('t=4');
subplot(3,2,5);
plot(iter,data(:,5),'Color','#749BC2','LineWidth',1.8);legend('t=5');
sgtitle('适应度曲线');

% set(figure4, 'Position', [100, 100, 800, 1000]);
saveas(figure4, '适应度曲线.png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 基于二维阵的高频地波雷达电离层回波方向估计
% % Author:RadarSun(ZhicongSun)
% % Data:from 2020.03.05 to 2020.xx.xx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%全局变量初始化
clear; close all;
derad = pi/180;
num_of_echoes = 1;			%信源个数

signal.freq = 4.7*10^6;		%发射（接收）信号频率4.7Mhz
signal.theta = 0;			%发射（接收）信号初始相位
signal.lamda = (3*10^8)/signal.freq;%发射（接收）信号波长

array.num = 8;			%天线阵元总个数
array.x_num = 4;		%X方向阵元个数
array.y_num = 4;		%Y方向阵元个数
array.spacing = signal.lamda/2;%阵元间距

%%阵元位置信息 8x8共16个阵元
axis_range.x = [array.spacing  2*array.spacing  3*array.spacing  4*array.spacing ...
                5*array.spacing 6*array.spacing 7*array.spacing 8*array.spacing 0 0 0 0 0 0 0 0];
axis_range.y = [0 0 0 0 0 0 0 0 0  array.spacing  2*array.spacing  3*array.spacing ...
                4*array.spacing 5*array.spacing 6*array.spacing 7*array.spacing];
axis_range.z = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

figure(1);
scatter3(axis_range.x, axis_range.y, axis_range.z);
axis([0, 9*array.spacing, 0, 9*array.spacing, 0, 9*array.spacing]);
title("L阵接收阵列示意图");
xlabel('X');
ylabel('Y：海岸线');
zlabel('Z');

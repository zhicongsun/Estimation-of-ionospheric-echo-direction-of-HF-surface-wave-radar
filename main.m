%{
    * 基于二维阵的高频地波雷达电离层回波方向估计
    * Description: my graduation projec
    * Author: ZhicongSun from HITWH
    * Email:hitsunzhicong@163.com
    * Github address: https://github.com/RadarSun
    * Data: from 2020.03.05 to 2020.03.12
%}

%%
%%全局变量初始化
close all; 
clear global;
clear all; 

global g_para g_signal g_array g_axis_range;
global g_echos;

%静态参数
g_para.rad = pi/180;        %常量
g_signal.freq = 4.7*10^6;     %发射（接收）信号频率4.7Mhz
g_signal.lamda = (3*10^8)/g_signal.freq;%发射（接收）信号波长

%默认的可变参数
testmode = 'dbf';
g_array.num = 16;                %天线阵元总个数
g_array.x_num = 8;		        %X方向阵元个数
g_array.y_num = 8;		        %Y方向阵元个数
g_array.span = g_signal.lamda/2; %阵元间距
g_array.x_pos = g_array.span : g_array.span : (g_array.x_num)*g_array.span;
g_array.y_pos = 0 : g_array.span : (g_array.y_num-1)*g_array.span;

g_echos.theta.num = 45;
g_echos.phi.num = 45;
g_echos.theta.rad = g_echos.theta.num*g_para.rad;
g_echos.phi.rad = g_echos.phi.num*g_para.rad;
g_echos.num = 1;        %回波数
g_echos.snr = 10;
g_echos.snapshot = 1000; %节拍数
g_echos.t = [0:99]/1000;
g_echos.signal=sqrt(10^(g_echos.snr/10))*exp(j*2*pi*g_signal.freq*g_echos.t);  %构造有用信号 

if((strcmp(testmode,'array_num')==0)&(strcmp(testmode,'array_span')==0))
    plotarray();
end
switch testmode
    case 'array_num'
        g_array.num = 8;                %天线阵元总个数
        g_array.x_num = 4;		        %X方向阵元个数
        g_array.y_num = 4;		        %Y方向阵元个数
        g_array.span = g_signal.lamda/2; %阵元间距
        g_array.x_pos = g_array.span : g_array.span : (g_array.x_num)*g_array.span;
        g_array.y_pos = 0 : g_array.span : (g_array.y_num-1)*g_array.span;
        beamforming('normal');
        g_array.num = 16;                %天线阵元总个数
        g_array.x_num = 8;		        %X方向阵元个数
        g_array.y_num = 8;		        %Y方向阵元个数
        g_array.span = g_signal.lamda/2; %阵元间距
        g_array.x_pos = g_array.span : g_array.span : (g_array.x_num)*g_array.span;
        g_array.y_pos = 0 : g_array.span : (g_array.y_num-1)*g_array.span;
        beamforming('normal');
        g_array.num = 32;                %天线阵元总个数
        g_array.x_num = 16;		        %X方向阵元个数
        g_array.y_num = 16;		        %Y方向阵元个数
        g_array.span = g_signal.lamda/2; %阵元间距
        g_array.x_pos = g_array.span : g_array.span : (g_array.x_num)*g_array.span;
        g_array.y_pos = 0 : g_array.span : (g_array.y_num-1)*g_array.span;
        beamforming('normal');
    case 'array_span'
        g_array.span = g_signal.lamda/4; %阵元间距
        g_array.x_pos = g_array.span : g_array.span : (g_array.x_num)*g_array.span;
        g_array.y_pos = 0 : g_array.span : (g_array.y_num-1)*g_array.span;
        beamforming('normal');
        g_array.span = g_signal.lamda/2; %阵元间距
        g_array.x_pos = g_array.span : g_array.span : (g_array.x_num)*g_array.span;
        g_array.y_pos = 0 : g_array.span : (g_array.y_num-1)*g_array.span;
        beamforming('normal');
        g_array.span = g_signal.lamda/2*3; %阵元间距
        g_array.x_pos = g_array.span : g_array.span : (g_array.x_num)*g_array.span;
        g_array.y_pos = 0 : g_array.span : (g_array.y_num-1)*g_array.span;
        beamforming('normal');
    case 'snapshot_num'
    case 'doa'
        g_echos.num = 3;            %回波数 
        g_echos.theta.num = [15 45 75] ;
        g_echos.phi.num = [70 46 10];
        g_echos.theta.rad = g_echos.theta.num*g_para.rad;
        g_echos.phi.rad = g_echos.phi.num*g_para.rad;
        g_echos.signal = rand(g_echos.num,g_echos.snapshot);
        doa();
        g_echos.num = 5;            %回波数 
        g_echos.theta.num = [15 35 45 60 75] ;
        g_echos.phi.num = [70 60 45 35 15];
        g_echos.theta.rad = g_echos.theta.num*g_para.rad;
        g_echos.phi.rad = g_echos.phi.num*g_para.rad;
        g_echos.signal = rand(g_echos.num,g_echos.snapshot);
        doa();
%         g_echos.signal = [exp(j*2*pi*g_signal.freq*g_echos.t);...
%                           exp(j* (2*pi*(g_signal.freq+1)*g_echos.t+30*g_para.rad) );...
%                           exp(j* (2*pi*(g_signal.freq+2)*g_echos.t+60*g_para.rad) )];
    case 'dbf'
        dbf_mode = 'normal_and_capon';
        switch dbf_mode
            case 'normal_and_capon'
                %% 对比二维普通波束形成和二维capon
                %可分辨
                g_echos.theta.num = [34,60];
                g_echos.phi.num = [60,35];
                g_echos.theta.rad = g_echos.theta.num*g_para.rad;
                g_echos.phi.rad = g_echos.phi.num*g_para.rad;
                g_echos.num = length(g_echos.theta.num);
                g_echos.snr = 10;
                g_echos.snapshot = 100; %节拍数
                g_echos.t = [0:99]/1000;
                g_echos.signal=rand(g_echos.num,g_echos.snapshot);
                beamforming('normal');
                beamforming('capon');
                %不可分辨
                g_echos.theta.num = [34,50];
                g_echos.phi.num = [40,34];
                g_echos.theta.rad = g_echos.theta.num*g_para.rad;
                g_echos.phi.rad = g_echos.phi.num*g_para.rad;
                g_echos.num = length(g_echos.theta.num);
                g_echos.snr = 10;
                g_echos.snapshot = 100; %节拍数
                g_echos.t = [0:99]/1000;
                g_echos.signal=rand(g_echos.num,g_echos.snapshot);
                beamforming('normal');
                beamforming('capon');

            case 'capon' 
                %% 测试二维capon
                g_echos.theta.num = [34,60];
                g_echos.phi.num = [60,35];
                g_echos.theta.rad = g_echos.theta.num*g_para.rad;
                g_echos.phi.rad = g_echos.phi.num*g_para.rad;
                g_echos.num = length(g_echos.theta.num);
                g_echos.snr = 10;
                g_echos.snapshot = 100; %节拍数
                g_echos.t = [0:99]/1000;
                g_echos.signal=rand(g_echos.num,g_echos.snapshot);
                beamforming('capon');
        
            case 'normal'
                %% 测普通数字波束形成
                g_echos.theta.num = 45;
                g_echos.phi.num = 50;
                g_echos.theta.rad = g_echos.theta.num*g_para.rad;
                g_echos.phi.rad = g_echos.phi.num*g_para.rad;
                g_echos.num = length(g_echos.theta.num);
                g_echos.snr = 10;
                g_echos.snapshot = 100; %节拍数
                g_echos.t = [0:99]/1000;
                g_echos.signal=rand(g_echos.num,g_echos.snapshot);
                beamforming('normal');

            case 'normal_theta_and_direction'
                %% 测指向方向对普通波束形成波束宽度的影响
                g_echos.theta.num = 40;
                g_echos.phi.num = 60;
                g_echos.theta.rad = g_echos.theta.num*g_para.rad;
                g_echos.phi.rad = g_echos.phi.num*g_para.rad;
                g_echos.num = length(g_echos.theta.num);
                g_echos.snr = 10;
                g_echos.snapshot = 100; %节拍数
                g_echos.t = [0:99]/1000;
                g_echos.signal=rand(g_echos.num,g_echos.snapshot);
                beamforming('normal');
                g_echos.theta.num = 10;
                g_echos.phi.num = 60;
                g_echos.theta.rad = g_echos.theta.num*g_para.rad;
                g_echos.phi.rad = g_echos.phi.num*g_para.rad;
                g_echos.num = length(g_echos.theta.num);
                g_echos.snr = 10;
                g_echos.snapshot = 100; %节拍数
                g_echos.t = [0:99]/1000;
                g_echos.signal=rand(g_echos.num,g_echos.snapshot);
                beamforming('normal');
                g_echos.theta.num = 60;
                g_echos.phi.num = 10;
                g_echos.theta.rad = g_echos.theta.num*g_para.rad;
                g_echos.phi.rad = g_echos.phi.num*g_para.rad;
                g_echos.num = length(g_echos.theta.num);
                g_echos.snr = 10;
                g_echos.snapshot = 100; %节拍数
                g_echos.t = [0:99]/1000;
                g_echos.signal=rand(g_echos.num,g_echos.snapshot);
                beamforming('normal');

            case 'normal_resolution'
                %% 测普通波束形成分辨率问题
                g_echos.theta.num = [40,55];
                g_echos.phi.num = [55,40];
                g_echos.theta.rad = g_echos.theta.num*g_para.rad;
                g_echos.phi.rad = g_echos.phi.num*g_para.rad;
                g_echos.num = length(g_echos.theta.num);
                g_echos.snr = 10;
                g_echos.snapshot = 100; %节拍数
                g_echos.t = [0:99]/1000;
                g_echos.signal=rand(g_echos.num,g_echos.snapshot);
                beamforming('normal');
                g_echos.theta.num = [40,60];
                g_echos.phi.num = [40,40];
                g_echos.theta.rad = g_echos.theta.num*g_para.rad;
                g_echos.phi.rad = g_echos.phi.num*g_para.rad;
                g_echos.num = length(g_echos.theta.num);
                g_echos.snr = 10;
                g_echos.snapshot = 100; %节拍数
                g_echos.t = [0:99]/1000;
                g_echos.signal=rand(g_echos.num,g_echos.snapshot);
                beamforming('normal');
            otherwise
                disp('Please input dbf mode!');
        end
    otherwise
        disp('Please input global mode!');
end

function  plotarray()
%{
    Function description:
            画8x8阵列图，参数都是固定的
    Syntax：
    Log description：
            2020.03.17  建立函数
%}  
    global   g_array g_axis_range;
    g_axis_range.x = [g_array.span  2*g_array.span  3*g_array.span  4*g_array.span ...    %阵元位置信息 8x8共16个阵元
        5*g_array.span 6*g_array.span 7*g_array.span 8*g_array.span 0 0 0 0 0 0 0 0];
    g_axis_range.y = [0 0 0 0 0 0 0 0 0  g_array.span  2*g_array.span  3*g_array.span ...
        4*g_array.span 5*g_array.span 6*g_array.span 7*g_array.span];
    g_axis_range.z = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    figure('Name','L阵示意图','NumberTitle','off','Color','white','Position',[200 200 400 400]);
    scatter3(g_axis_range.x, g_axis_range.y, g_axis_range.z);
    axis([0, 9*g_array.span, 0, 9*g_array.span, 0, 9*g_array.span]);
    title("L阵接收阵列示意图");
    xlabel('X/m');
    ylabel('Y/m：海岸线');
    zlabel('Z');
end

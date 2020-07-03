clc;
clear all;
close all;
imag=sqrt(-1);
element_num=8;%阵元数为8
d_lambda=0.5;%间距为半波长
theta=-90:0.5:90;%扫描范围
theta0=0;%来波方向
theta1=50;%干扰方向
L=1024;%采样单元数
for i=1:L
    amp0=10*randn(1);
    amp1=50*randn(1);
    ampn=0.5;
  s(:,i)=amp0*exp(imag*2*pi*d_lambda*sin(theta0*pi/180)*[0:element_num-1]');
  j(:,i)=amp1*exp(imag*2*pi*d_lambda*sin(theta1*pi/180)*[0:element_num-1]');
  n(:,i)=ampn*exp(imag*2*pi*randn(1)*[0:element_num-1]');
end
Rx=1/L*(s+j+n)*(s+j+n)';%接收信号自相关矩阵
Rnj=1/L*(j+n)*(j+n)';%%干拢+噪声的自相关矩阵
e=exp(imag*2*pi*d_lambda*sin(theta0*pi/180)*[0:element_num-1]');
Wopt_Rx=inv(Rx)*e/(e'*inv(Rx)*e);%采用接收信号的权矢量
Wopt_Rnj=inv(Rnj)*e/(e'*inv(Rnj)*e);%采用干拢+噪声信号的权矢量
for j=1:length(theta)
    a=exp(imag*2*pi*d_lambda*sin(theta(j)*pi/180)*[0:element_num-1]');
    f1(j)=Wopt_Rx'*a;
    f2(j)=Wopt_Rnj'*a;
end
F1=20*log10(abs(f1)/max(max(abs(f1))));
F2=20*log10(abs(f2/max(max(abs(f2)))));
figure;
plot(theta,F1,'b',theta,F2,'r');
grid on;
hold on;
plot(theta0,-50:0,'.');
plot(theta1,-50:0,'.');
xlabel('theta/°');
ylabel('F(1,2)/dB');
title('不同方法估计协方差矩阵的Capon波束形成');
axis([-90 90 -60 0]);

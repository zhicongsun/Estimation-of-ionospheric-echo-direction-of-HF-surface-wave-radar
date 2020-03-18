function capon()
    clc;
    clear all;
    close all;
    imag=sqrt(-1);
    element_num=8;%阵元数
    d_lambda=0.5;%阵元间距与波长的关系
    theta=-90:0.5:90; %搜索范围
    theta0=0; %三个信号源的来波方向
    theta1=20;
    theta2=60;
    L=1000;%采样单元数
    for i=1:L
        amp0=10*randn(1);
        amp1=200*randn(1);
        amp2=200*randn(1);
        ampn=3;   x(:,i)=amp0*exp(imag*2*pi*d_lambda*sin(theta0*pi/180)*[0:element_num-1]')+amp1*exp(imag*2*pi*d_lambda*sin(theta1*pi/180)*[0:element_num-1]')+amp2*exp(imag*2*pi*d_lambda*sin(theta2*pi/180)*[0:element_num-1]')+ampn*(randn(element_num,1)+imag*randn(element_num,1));
    end
    Rx=1/L*x*x';
    R=inv(Rx);
    steer=exp(imag*2*pi*d_lambda*sin(theta0*pi/180)*[0:element_num-1]');
    w=R*steer/(steer'*R*steer);%最优权矢量
    for j=1:length(theta)
        a=exp(imag*2*pi*d_lambda*sin(theta(j)*pi/180)*[0:element_num-1]');
        f(j)=w'*a;
        p(j)=1/(a'*R*a);
    end
    F=20*log10(abs(f)/(max(max(abs(f)))));
    subplot(1,2,1)
    plot(theta,F);
    grid on;
    hold on;
    plot(theta0,-50:0,'.');
    plot(theta1,-50:0,'.');
    plot(theta2,-50:0,'.');
    xlabel('theta/°');
    ylabel('F/dB');
    title('Capon beamforming 方向图');
    axis([-90 90 -50 0]);
    P=20*log10(abs(p)/(max(max(abs(p)))));
    subplot(1,2,2)
    plot(theta,P);
    grid on;
    hold on;
    xlabel('theta/°');
    ylabel('P/dB');
    title('Capon beamforming 功率谱');
    axis([-90 90 -90 0]);
end
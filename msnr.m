function msnr()
    clc;
    clear all;
    close all;
    imag=sqrt(-1);
    element_num=8;%阵元数为8
    d_lambda=0.5;%间距为半波长
    theta=-90:0.5:90;%扫描范围
    theta0=0;%来波方位
    theta1=20;%干扰方向
    L=512;%采样点数
    for i=1:L
        amp0=10*randn(1);
        amp1=200*randn(1);
        ampn=1;    s(:,i)=amp0*exp(imag*2*pi*0.5*sin(theta0*pi/180)*[0:element_num-1]');   j(:,i)=amp1*exp(imag*2*pi*0.5*sin(theta1*pi/180)*[0:element_num-1]');    n(:,i)=ampn*exp(randn(element_num,1)+imag*randn(element_num,1));
    end
    Rs=1/L*s*s';%信号自相关矩阵
    Rnj=1/L*(j*j'+n*n'); %干扰+噪声的自相关矩阵
    [V,D]=eig(Rs,Rnj); %（Rs,Rnj）的广义特征值和特征向量
    [D,I]=sort(diag(D)); %特征向量排序
    Wopt=V(:,I(8));%最优权矢量
    for j=1:length(theta)
    a=exp(imag*2*pi*d_lambda*sin(theta(j)*pi/180)*[0:element_num-1]');
        f(j)=Wopt'*a;
        p(j)=a'*Rs*a+a'*Rnj*a;
    end
    F=20*log10(abs(f)/max(max(abs(f))));
    P=20*log10(abs(p)/max(max(abs(p))));
    subplot(1,2,1)
    plot(theta,F);
    grid on;
    hold on;
    plot(theta0,-80:0,'.');
    plot(theta1,-80:0,'.');
    xlabel('theta/0');
    ylabel('F in dB');
    title('max-SNR 方向图');
    axis([-90 90 -80 0]);
    hold on;
    subplot(1,2,2);
    plot(theta,P,'r');
    grid on;
    xlabel('theta/0');
    ylabel('功率 in dB');
    title('max-SNR 功率谱');
    grid on;
    axis([-90 90 -80 0]);
end
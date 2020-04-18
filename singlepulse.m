%function [abs_f,abs_p]=singlepulse(theta0,element_num)
    %{
        Function description:
                一维普通波束形成
        Syntax：
                Input:
                        theta0：目标角度,单位为rad,可多目标,如: theta0 = [30,20]/180*pi
                        element_num：阵元个数
                        d_lamda: 阵元间距与lamda的比例，单位1表示波长，0.5表示半波长
                Output:
                        abs_f：方向图的幅值
                        abs_p：功率谱幅值
        Log description：
                2020.03.17  建立函数
                2020.03.20  为了对比不同算法而将波形显示剔除，加入输入输出参数
    %}  
    clear all;
    clear global;
    close all;
    %if nargin == 0
        theta0 = 30;
        element_num = 8;
    %end
    imag=sqrt(-1);
    d_lamda = 0.5;
    thetaB = 50.8/(element_num*d_lamda);
    thetaK = 0.5*thetaB;
    thetaddd = 2*asind((1/element_num/0.5));

    thetaddd = 2*asind( (1/element_num/0.5) + sin((theta0+thetaK)/180*pi))
    % snapshot = 200;
    % snr = 10;
    
    w1=exp(imag*2*pi*d_lamda*[0:element_num-1]'*sin((theta0+thetaK)/180*pi));
    w2=exp(imag*2*pi*d_lamda*[0:element_num-1]'*sin((theta0-thetaK)/180*pi));
    % if length(theta0) ~= 1
    %     wtemp = w(:,1);
    % else
    %     wtemp = w;
    % end
    % S = rand(length(theta0),snapshot);%KxP
    % X = w*S;%MxP
    % X1 = awgn(X,snr,'measured');
    % Rxx = X1*X1'/snapshot; 

    theta=linspace(-90,90,181);
    j=1:length(theta);
    a=exp(imag*2*pi*d_lamda*[0:element_num-1]'*sin(theta(j)/180*pi));
    f1(j)=w1'*a;f1 =abs(f1);
    f2(j)=w2'*a;f2 =abs(f2);
    sumf = abs(f1)+abs(f2);
    delf = abs(f1)-abs(f2);
    K = delf./sumf;

    figure();
    plot(theta,f1);hold on;
    plot(theta,f2);
    % plot(theta(-10+91:10+91),f1(-10+91:10+91));hold on;
    % plot(theta(-10+91:10+91),f2(-10+91:10+91));
    legend('f1','f2');
    figure();
    plot(theta,sumf);hold on;
    plot(theta,delf);
    % plot(theta(-10+91:10+91),sumf(-10+91:10+91));hold on;
    % plot(theta(-10+91:10+91),delf(-10+91:10+91));
    legend('sum','del');
    figure();
    plot(theta,K);


%end    
function set_rectangle(f,theta)

end   
    
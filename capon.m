function [abs_f,abs_p]=capon(theta0,element_num,d_lamda)
    %{
    Function description:
            一维Capon算法
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
    if nargin == 0
        theta0 = 30;
    end
    imag=sqrt(-1)
    snapshot = 200;
    snr = 10;
    theta=linspace(-90,90,200);
    a=exp(imag*2*pi*d_lamda*[0:element_num-1]'*sin(theta0/180*pi));%MxK
    if length(theta0) ~= 1
        steer = a(:,1);
    else
        steer = a;
    end
    S = rand(length(theta0),snapshot);%KxP
    X = a*S;%MxP
    X1 = awgn(X,snr,'measured');
    Rxx = X1*X1'/snapshot; 
    R = inv(Rxx);
    C = [a(:,1) a(:,2)];
    F = [1 0]';
    A=R*C*(inv(C'*R*C))*F;
    for  j=1:length(theta)
        w=exp(imag*2*pi*d_lamda*[0:element_num-1]'*sin(theta(j)/180*pi));%Mx1
        f(j) = A'*w;
        p(j) = 1/(w'*R*w);
    end
    abs_f=abs(f);
    abs_f_max=max(abs_f);
    abs_f = 10*log10(abs_f/abs_f_max);
    abs_p=abs(p);
    abs_p_max=max(abs_p);
    abs_p = 10*log10(abs_p/abs_p_max);
end    
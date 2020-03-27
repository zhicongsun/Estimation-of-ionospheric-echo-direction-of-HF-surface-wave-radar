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
    ampn=3;   
    x(:,i)=amp0*exp(imag*2*pi*d_lambda*sin(theta0*pi/180)*[0:element_num-1]')...
    +amp1*exp(imag*2*pi*d_lambda*sin(theta1*pi/180)*[0:element_num-1]')...
    +amp2*exp(imag*2*pi*d_lambda*sin(theta2*pi/180)*[0:element_num-1]')...
    +ampn*(randn(element_num,1)+imag*randn(element_num,1));
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



function doa(input)
    %{
        Function description:
                2D-DOA估计算法选择器
        Syntax：
                input:doa算法类型
        Log description：
                2020.03.10 添加doa选择器，加入music2d
    %}
        if nargin<1
            input = "music2d";
        end
        switch input
            case 'music2d'
                music2d();
            otherwise 
                error('Error! Invalid input! ');
        end
    end
    
    
    function music2d()
    %{
        Function description:
                2D-MUSIC DOA算法
        Syntax：
        Log description：
                2020.03.10 完成music2d函数，存在非理想谱峰的问题
    %}
        global g_signal;
        global g_array;
        global g_echos;
        global g_para;
        
        %行是阵元个数，列是回波个数
        Ax = exp( -j*2*pi/g_signal.lamda*(cos(g_echos.theta.rad).*cos(g_echos.phi.rad)).'*g_array.x_pos);
        %k =g_array.x_pos.'/g_signal.lamda
        Ay = exp( -j*2*pi/g_signal.lamda*(cos(g_echos.theta.rad).*sin(g_echos.phi.rad)).'*g_array.y_pos);
        A = [Ax ; Ay];%行总阵元个数
        X = A*g_echos.signal;
        X1 = awgn(X,g_echos.snr,'measured');
        %X1= X;%测试无噪声
        Rxx = X1*X1'/g_echos.snapshot;
        
        [eigenvector,eigenvalue] = eig(Rxx);
        [EVA,I] = sort(diag(eigenvalue).');
        eigenvector = fliplr(eigenvector(:,I));
        Un = eigenvector(:,g_echos.num+1:end);
        for ang1 = 1:90
            for ang2 = 1:90
                thet(ang1) = ang1-1;%0~89度
                phim1 = thet(ang1)*g_para.rad; 
                f(ang2) = ang2-1;
                phim2 = f(ang2)*g_para.rad;
                ax = exp(-j*2*pi/g_signal.lamda*g_array.x_pos.*cos(phim1)*cos(phim2));
                ay = exp(-j*2*pi/g_signal.lamda*g_array.y_pos.*cos(phim1)*sin(phim2));
                a = [ax;ay];
                SP(ang1,ang2) = 1/(a'*Un*Un'*a);
            end
        end
        SP=abs(SP);
        %[linemax,linemax_id] = max(SP);
        %[rowmax,rowmax_id] = max(linemax);
        %rowmax_id
        %SP(rowmax_id,linemax_id);
        SPmax=max(max(SP));
        SP=SP/SPmax; 
        figure('Color','white');
        h = mesh(thet,f,SP);
        set(h,'Linewidth',2)
        xlabel('俯仰(degree)');
        ylabel('方位(degree)');
        zlabel('magnitude(dB)');
    
    end

    function normaldbf()
        %{
            Function description:
                    普通DBF算法，用于2D-DOA估计
            Syntax：
            Log description：
                    2020.03.11  初步实现DBF，单个目标波信合理；问题：多目标出现额外波峰，三维图XY颠倒
                    2020.03.12  解决mesh存在的XY坐标颠倒问题，用sum处理多目标导向矢量，解决额外波峰
        %}
            global g_signal;
            global g_array;
            global g_echos;
            global g_para;
            tic;
            
            theta.rad=linspace(-pi/2,pi/2,181);
            theta.num = theta.rad/pi*180;
            phi.rad=linspace(-pi/2,pi/2,181);
            phi.num = phi.rad/pi*180;
        
            %普通数字波束形成
            wx=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*sin(g_echos.theta.rad).*cos(g_echos.phi.rad));
            wx = sum(wx,2);%对行求和,即对各个信源的导向向量累加,作用类似wx = wx(:,1) + wx(:,2);
            wy=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*sin(g_echos.theta.rad).*sin(g_echos.phi.rad));
            wy = sum(wy,2);
            W = kron(wx,wy);%相当于将wx*wy矩阵的每一列移位到第一列末尾，形成一维向量
            X = W*g_echos.signal;
            X1 = awgn(X,g_echos.snr,'measured');
            Rxx = X1*X1'/g_echos.snapshot;   
            for k = 1:length(theta.rad)
                for  i = 1:length(phi.rad)
                    ax=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*sin(theta.rad(k)).*cos(phi.rad(i)));
                    ay=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*sin(theta.rad(k)).*sin(phi.rad(i)));
                    A = kron(ax,ay);
                    F(k,i) = W'*A;
                    P(k,i) = A'*Rxx*A;
                end        
            end
            toc;
            disp(['普通波束形成算法用时：',num2str(toc),'s']);
        
            %画方向图
            abs_F=abs(F);
            %abs_F_norm=abs_F/max(max(abs_F));
            %abs_F_dB=20*log10(abs_F);
            %abs_F_norm_dB=20*log10(abs_F_norm);
            absf_theta = sum(abs_F,2);
            absf_phi = sum(abs_F,1);
            figure('Name','方向图','NumberTitle','off','Color','white','Position',[600 50 550 750]);
            subplot(311);
            f = meshc(phi.num,theta.num,abs_F);% mesh(X,Y,Z)X,Y分别对应Z的列、行，mesh太坑了，要先XY置换下
            title(['二维方向图,阵元数:' num2str(g_array.num) ',阵元间距:' num2str(g_array.span/g_signal.lamda) 'lamda,指向[theta,phi]=[' num2str(g_echos.theta.num) ',' num2str(g_echos.phi.num) ']' ]);
            set(f,'Linewidth',2);
            xlabel('方位角phi(degree)');
            ylabel('俯仰角theta(degree)');
            zlabel('abs of F');
            subplot(312);
            plot(theta.num,absf_theta);
            title('俯仰角 的方向图');
            xlabel('theta(degree)');
            ylabel('abs of Ptheta');
            subplot(313);
            plot(phi.num,absf_phi);
            title('方位角 的方向图');
            xlabel('phi(degree)');
            ylabel('abs of Pphi');
            
            %画功率谱
            abs_P = abs(P);
            absp_theta = sum(abs_P,2);%2是对列相加，得到列向量；1则是得到行向量
            absp_phi = sum(abs_P,1);
            %谱峰搜索
            [temp_ptheta,theta_estimation] = max(absp_theta);
            [temp_pphi,phi_estimation] = max(absp_phi);
            theta_estimation = theta_estimation-91;
            phi_estimation = phi_estimation-91;
            disp(['估计角度为[',num2str(theta_estimation),',',num2str(phi_estimation),']']);
            
            figure('Name','功率谱','NumberTitle','off','Color','white');
            subplot(311);
            p = meshc(phi.num(91:181),theta.num(91:181),abs_P(91:181,91:181));
            title(['二维功率谱,阵元间距:' num2str(g_array.span/g_signal.lamda) 'lamda, 实际方向为[theta,phi]=[' num2str(g_echos.theta.num) ',' num2str(g_echos.phi.num) ']' ...
                '估计方向为[' num2str(theta_estimation) ',' num2str(phi_estimation) ']']);
            set(p,'Linewidth',2);
            xlabel('方位角phi(degree)');
            ylabel('俯仰角theta(degree)');
            zlabel('abs of P');
            subplot(312);
            plot(theta.num,absp_theta);
            title('俯仰角 的功率谱');
            xlabel('theta(degree)');
            ylabel('abs of Ptheta');
            subplot(313);
            plot(phi.num,absp_phi);
            title('方位角 的功率谱');
            xlabel('phi(degree)');
            ylabel('abs of Pphi');
        
        end 
        
        

        function [abs_p,peak_ang,snr,rmse]=music(theta0,element_num,d_lamda,multipath_mode)
            %{
                    Function description:
                            一维普通波束形成
                    Syntax：
                            Input:
                                    theta0：目标角度,单位为rad,可多目标,如: theta0 = [30,20]/180*pi
                                    element_num：阵元个数
                                    d_lamda: 阵元间距与lamda的比例，单位1表示波长，0.5表示半波长
                            Output:
                                    abs_p：功率谱幅值
                    Log description：
                            2020.03.25  建立函数
            %}  
            derad = pi/180;        
            radeg = 180/pi;
            twpi = 2*pi;
            d=0:d_lamda:(element_num-1)*d_lamda;     
            iwave = 3;              
            snr = 10;               
            n = 500;                 
            A=exp(-j*twpi*d.'*sin(theta0));
            if strcmp(multipath_mode,'multi_path')==1
                S0 = randn(iwave-1,n);
                S = [S0(1,:);S0];
            else
                S=randn(iwave,n);
            end
            X=A*S;
            snr0=0:3:100;
            for isnr=1:20
                X1=awgn(X,snr0(isnr),'measured');
                Rxx=X1*X1'/n;
                InvS=inv(Rxx); 
                [EV,D]=eig(Rxx); 
                EVA=diag(D)';
                [EVA,I]=sort(EVA);
                EVA=fliplr(EVA);
                EV=fliplr(EV(:,I));
                % MUSIC
                for iang = 1:361
                        angle(iang)=(iang-181)/2;
                        phim=derad*angle(iang);
                        a=exp(-j*twpi*d*sin(phim)).';
                        L=iwave;    
                        En=EV(:,L+1:element_num);
                        P(iang)=(a'*a)/(a'*En*En'*a);
                end
                abs_p=abs(P);
                abs_p_max=max(abs_p);
                abs_p=10*log10(abs_p/abs_p_max);
                peak_ang = [];
                derivative = diff(abs_p);
                for iang = 2:360
                    if( (abs_p(iang)>abs_p(iang-1)) && (abs_p(iang)>abs_p(iang+1)) && (derivative(iang-1)>0.5) )
                        peak_ang = [peak_ang iang];
                    end
                end
                peak_ang = (peak_ang-181)/2;
                size(peak_ang)
                rmse(isnr) = sqrt( sum((theta0/pi*180-peak_ang).^2)/iwave);
            end
                snr = snr0(1:20);
            end
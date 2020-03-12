% function test()
%     global g_signal;
%     global g_array;
%     global g_echos;
%     global g_para;
% 
%     element_num=32;%阵元数为8
% 
%     theta=linspace(-pi/2,pi/2,200);
%     phi = linspace(-pi/2,pi/2,200);
%     theta0=80/180*pi;%来波方向
%     phi0 = 70/180*pi;%来波方向
% 
%     wx=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*cos(theta0)*cos(phi0)) 
%     for  i=1:length(phi)
%         ax=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*cos(theta0)*cos(phi(i)));
%        px(i)=wx'*ax;
%     end
% 
%     wy=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*cos(theta0)*sin(phi0));
%     for  k=1:length(theta)
%         ay=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*cos(theta(k))*sin(phi0));
%        py(k)=wy'*ay;
%     end
% 
%     figure();
%     patternmag=abs(px);
%     patternmagnorm=patternmag/max(max(patternmag));
%     patterndB=20*log10(patternmag);
%     patterndBnorm=20*log10(patternmagnorm);
%     subplot(2,1,1);
%     plot(phi*180/pi,patternmag);
%     grid on;
%     xlabel('theta/radian')
%     ylabel('amplitude/dB')
%     title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(phi0*180/pi) '度']);
%     hold on;
%     subplot(2,1,2);
%     plot(phi,patterndBnorm,'r');
%     grid on;
%     xlabel('theta/radian')
%     ylabel('amplitude/dB')
%     title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(phi0*180/pi) '度']);
%     axis([-1.5 1.5 -50 0]);
% 
%     figure();
%     patternmag=abs(py);
%     patternmagnorm=patternmag/max(max(patternmag));
%     patterndB=20*log10(patternmag);
%     patterndBnorm=20*log10(patternmagnorm);
%     subplot(2,1,1);
%     plot(theta*180/pi,patternmag);
%     grid on;
%     xlabel('theta/radian')
%     ylabel('amplitude/dB')
%     title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(theta0*180/pi) '度']);
%     hold on;
%     subplot(2,1,2);
%     plot(theta,patterndBnorm,'r');
%     grid on;
%     xlabel('theta/radian')
%     ylabel('amplitude/dB')
%     title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(theta0*180/pi) '度']);
%     axis([-1.5 1.5 -50 0]);
% end

% function test()
%     global g_signal;
%     global g_array;
%     global g_echos;
%     global g_para;
% 
%     element_num=8;%阵元数为8
% 
%     theta=linspace(-pi/2,pi/2,200);
%     phi = linspace(-pi/2,pi/2,200);
%     theta0=[20 80]/180*pi;%来波方向
%     phi0 = 47/180*pi;%来波方向
% 
%     wx=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*cos(theta0)*cos(phi0)) 
%     wx = wx(:,1) + wx(:,2);
%     size(wx)
%     for  i=1:length(phi)
%         ax=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*cos(theta0)*cos(phi(i)));
%         ax = ax(:,1)+ax(:,2);
%        px(i)=wx'*ax;
%     end
% 
%     wy=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*cos(theta0)*sin(phi0));
%     wy = wy(:,1) + wy(:,2);
%     %wy = wy(:,2);
%     for  k=1:length(theta)
%         ay=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*cos(theta(k))*sin(phi0));
%         %ay = ay(:,1).*ay(:,2);
%        py(k)=wy'*ay;
%     end
% 
%     figure();
%     patternmag=abs(px);
%     patternmagnorm=patternmag/max(max(patternmag));
%     patterndB=20*log10(patternmag);
%     patterndBnorm=20*log10(patternmagnorm);
%     subplot(2,1,1);
%     plot(theta*180/pi,patternmag);
%     grid on;
%     xlabel('theta/radian')
%     ylabel('amplitude/dB')
%     title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(phi0*180/pi) '度']);
%     hold on;
%     subplot(2,1,2);
%     plot(theta,patterndBnorm,'r');
%     grid on;
%     xlabel('theta/radian')
%     ylabel('amplitude/dB')
%     title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(phi0*180/pi) '度']);
%     axis([-1.5 1.5 -50 0]);
% 
%     figure();
%     patternmag=abs(py);
%     patternmagnorm=patternmag/max(max(patternmag));
%     patterndB=20*log10(patternmag);
%     patterndBnorm=20*log10(patternmagnorm);
%     subplot(2,1,1);
%     plot(theta*180/pi,patternmag);
%     grid on;
%     xlabel('theta/radian')
%     ylabel('amplitude/dB')
%     title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(theta0*180/pi) '度']);
%     hold on;
%     subplot(2,1,2);
%     plot(theta,patterndBnorm,'r');
%     grid on;
%     xlabel('theta/radian')
%     ylabel('amplitude/dB')
%     title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(theta0*180/pi) '度']);
%     axis([-1.5 1.5 -50 0]);
% end

% clc;
% clear all;
% close all;
% element_num=8;%阵元数为8
% d_lamda=1/2;%阵元间距d与波长lamda的关系
% 
% theta=linspace(-pi/2,pi/2,200);
% phi = linspace(-pi/2,pi/2,200);
% theta0=80/180*pi;%来波方向
% phi0 = 47/180*pi;%来波方向
% 
% wx=exp(-j*2*pi*d_lamda*[0:element_num-1]'*cos(theta0)*cos(phi0)) 
% wx = wx(:,1)+wx(:,2);
% size(wx)
% for  i=1:length(phi)
%     ax=exp(-j*2*pi*d_lamda*[0:element_num-1]'*cos(theta0)*cos(phi(i)));
%     ax = ax(:,1)+ax(:,2);
%    px(i)=wx'*ax;
% end
% 
% wy=exp(-j*2*pi*d_lamda*[0:element_num-1]'*cos(theta0)*sin(phi0));
% wy = wy(:,1) + wy(:,2);
% %wy = wy(:,2);
% for  k=1:length(theta)
%     ay=exp(-j*2*pi*d_lamda*[0:element_num-1]'*cos(theta(k))*sin(phi0));
%     %ay = ay(:,1).*ay(:,2);
%    py(k)=wy'*ay;
% end
% 
% figure();
% patternmag=abs(px);
% patternmagnorm=patternmag/max(max(patternmag));
% patterndB=20*log10(patternmag);
% patterndBnorm=20*log10(patternmagnorm);
% subplot(2,1,1);
% plot(theta*180/pi,patternmag);
% grid on;
% xlabel('theta/radian')
% ylabel('amplitude/dB')
% title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(phi0*180/pi) '度']);
% hold on;
% subplot(2,1,2);
% plot(theta,patterndBnorm,'r');
% grid on;
% xlabel('theta/radian')
% ylabel('amplitude/dB')
% title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(phi0*180/pi) '度']);
% axis([-1.5 1.5 -50 0]);
% 
% figure();
% patternmag=abs(py);
% patternmagnorm=patternmag/max(max(patternmag));
% patterndB=20*log10(patternmag);
% patterndBnorm=20*log10(patternmagnorm);
% subplot(2,1,1);
% plot(theta*180/pi,patternmag);
% grid on;
% xlabel('theta/radian')
% ylabel('amplitude/dB')
% title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(theta0*180/pi) '度']);
% hold on;
% subplot(2,1,2);
% plot(theta,patterndBnorm,'r');
% grid on;
% xlabel('theta/radian')
% ylabel('amplitude/dB')
% title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(theta0*180/pi) '度']);
% axis([-1.5 1.5 -50 0]);
% 


clc;
clear all;
close all;
element_num=32;%阵元数为8
d_lamda=1/2;%阵元间距d与波长lamda的关系

theta=linspace(-pi/2,pi/2,200);
phi = linspace(-pi/2,pi/2,200);
theta0=80/180*pi;%来波方向
phi0 = 70/180*pi;%来波方向

wx=exp(-j*2*pi*d_lamda*[0:element_num-1]'*cos(theta0)*cos(phi0)) ;
for  i=1:length(phi)
    ax=exp(-j*2*pi*d_lamda*cos(theta0)*cos(phi(i))*[0:element_num-1]');
   px(i)=wx'*ax;
end

wy=exp(-j*2*pi*d_lamda*cos(theta0)*sin(phi0)*[0:element_num-1]');
for  k=1:length(theta)
    ay=exp(-j*2*pi*d_lamda*cos(theta(k))*sin(phi0)*[0:element_num-1]');
   py(k)=wy'*ay;
end

figure();
patternmag=abs(px);
patternmagnorm=patternmag/max(max(patternmag));
patterndB=20*log10(patternmag);
patterndBnorm=20*log10(patternmagnorm);
subplot(2,1,1);
plot(phi*180/pi,patternmag);
grid on;
xlabel('theta/radian')
ylabel('amplitude/dB')
title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(phi0*180/pi) '度']);
hold on;
subplot(2,1,2);
plot(phi,patterndBnorm,'r');
grid on;
xlabel('theta/radian')
ylabel('amplitude/dB')
title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(phi0*180/pi) '度']);
axis([-1.5 1.5 -50 0]);

figure();
patternmag=abs(py);
patternmagnorm=patternmag/max(max(patternmag));
patterndB=20*log10(patternmag);
patterndBnorm=20*log10(patternmagnorm);
subplot(2,1,1);
plot(theta*180/pi,patternmag);
grid on;
xlabel('theta/radian')
ylabel('amplitude/dB')
title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(theta0*180/pi) '度']);
hold on;
subplot(2,1,2);
plot(theta,patterndBnorm,'r');
grid on;
xlabel('theta/radian')
ylabel('amplitude/dB')
title([num2str(element_num) '阵元均匀线阵方向图','来波方向为' num2str(theta0*180/pi) '度']);
axis([-1.5 1.5 -50 0]);




% function test()
%     global g_signal;
%     global g_array;
%     global g_echos;
%     global g_para;
%     
%     theta = linspace(-pi/2,pi/2,201);
%     phi = linspace(-pi/2,pi/2,201);
%     %行是阵元个数，列是回波个数 8x3(假设有3个回波)
%     Wx = exp( -j*2*pi/g_signal.lamda*g_array.x_pos.'*(cos(g_echos.theta.rad).*cos(g_echos.phi.rad)));
%     Wy = exp( -j*2*pi/g_signal.lamda*g_array.y_pos.'*(cos(g_echos.theta.rad).*sin(g_echos.phi.rad)));
% 
%     for i = 1:length(theta)
%         for k = 1:length(phi)
%             ax = exp( -j*2*pi/g_signal.lamda*g_array.x_pos.'*(cos(theta(i)).*cos(phi(k))) );
%             ay = exp( -j*2*pi/g_signal.lamda*g_array.y_pos.'*(cos(theta(i)).*sin(phi(k))) );
%             Px(i,k) = Wx'*ax;% (3x8)*(8x1)=(3x1)
%             Py(i,k) = Wy'*ay;
%         end
%     end
% 
%     Px(101,:);
%     Py(:,10)';
%     size(Px(101,:));
%     size(Py(:,101));
% 
%     patternmag_x=abs(Px(101,:));
%     patternmagnorm_x=patternmag_x/max(max(patternmag_x));
%     patterndB_x=20*log10(patternmag_x);
%     patterndBnorm_x=20*log10(patternmagnorm_x);
%     figure();
%     subplot(2,1,1);
%     plot(phi*180/pi,patternmag_x);
%     grid on;
%     xlabel('pi/度')
%     ylabel('amplitude/dB')
%     title(['方位角方向图','来波方向为' num2str(g_echos.phi.num) '度']);
%     hold on;
%     subplot(2,1,2);
%     plot(phi,patterndBnorm_x,'r');
%     grid on;
%     xlabel('pi/rad')
%     ylabel('amplitude/dB')
%     title(['方位角方向图','来波方向为' num2str(g_echos.phi.num) '度']);
%     axis([-1.5 1.5 -50 0]);
%     hold on;
%     
%     patternmag_y=abs(Py(:,10));
%     patternmagnorm_y=patternmag_y/max(max(patternmag_y));
%     patterndB_y=20*log10(patternmag_y);
%     patterndBnorm_y=20*log10(patternmagnorm_y);
%     figure();
%     subplot(2,1,1);
%     plot(theta*180/pi,patternmag_y);
%     grid on;
%     xlabel('theta/度')
%     ylabel('amplitude/dB')
%     title(['俯仰方向图','来波方向为' num2str(g_echos.theta.num) '度']);
%     hold on;
%     subplot(2,1,2);
%     plot(theta,patterndBnorm_y,'r');
%     grid on;
%     xlabel('theta/rad')
%     ylabel('amplitude/dB')
%     title(['俯仰角方向图','来波方向为' num2str(g_echos.theta.num) '度']);
%     axis([-1.5 1.5 -50 0]);
%     hold on;
% end


%{
clear all
close all
clc

twpi = 2*pi;
rad = pi/180;
deg = 180/pi;

kelm = 8;
snr  = 10;
iwave = 3;
theta = [10 30 50];
fe = [15 25 35];
n = 100;
%dd = 3*10^8/4.7*10^6;
dd =5;
d = 0:dd:(kelm-1)*dd;
d1 = dd:dd:(kelm-1)*dd;
% d = dd:dd:(kelm)*dd;
% d1 = 0:dd:(kelm-1)*dd;
Ax = exp(-j*twpi*d.'*(cos(theta*rad).*cos(fe*rad)));
Ay = exp(-j*twpi*d1.'*(cos(theta*rad).*sin(fe*rad)));
A = [Ax;Ay];

freq = 4.7*10^6;
t = (0:99)/1000;
S = [sin(2*pi*freq*t) ; sin(2*pi*freq*(t+10)) ;sin(2*pi*freq*(t+20)) ];
%S = randn(iwave,n);

X = A*S;
X1 = awgn(X,snr,'measured');
Rxx = X1*X1'/n;
[EV,D] = eig(Rxx);
[EVA,I] = sort(diag(D).');
EV = fliplr(EV(:,I));
Un = EV(:,iwave+1:end);
for ang1 = 1:90
    for ang2 = 1:90
        thet(ang1) = ang1-1;
        phim1 = thet(ang1)*rad;
        f(ang2) = ang2-1;
        phim2 = f(ang2)*rad;
        a1 = exp(-j*twpi*d.'*cos(phim1)*cos(phim2));
        a2 = exp(-j*twpi*d1.'*cos(phim1)*sin(phim2));
        a = [a1;a2];
        SP(ang1,ang2) = 1/(a'*Un*Un'*a);
    end
end
SP=abs(SP);
SPmax=max(max(SP));
SP=SP/SPmax; 
h = mesh(thet,f,SP);
set(h,'Linewidth',2)
xlabel('elevation(degree)')
ylabel('azimuth(degree)')
zlabel('magnitude(dB)')
%}

%{    
    function output = test()
    %myFun - Description
    %
    % Syntax: output = myFun(input)
    %
    % Long description
    global g_array;
    output = g_array.x_num;
    end
   
%}

%{
    g_echos(2).theta.num = 30;
    g_echos(2).theta.rad = 30*g_para.rad;
    g_echos(2).phi.num = 30;
    g_echos(2).phi.rad = 30*g_para.rad;
    g_echos(3).theta.num = 60;
    g_echos(3).theta.rad = 60*g_para.rad;
    g_echos(3).phi.num = 60;
    g_echos(3).phi.rad = 60*g_para.rad;
%}
% 二维角度解析

function test2()
    global g_signal;
    global g_array;
    global g_echos;
    global g_para;
    
    theta=linspace(-pi/2,pi/2,201);
    phi=linspace(-pi/2,pi/2,201);
    theta0=60/180*pi;%来波方向
    phi0 = 30/180*pi;

    wx=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*cos(theta0).*cos(phi0));
    wy=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*cos(theta0).*sin(phi0));
    W = kron(wx,wy);
    size(W)
    for  k=1:length(theta)
        for  i=1:length(phi)
            ax=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*cos(theta(k)).*cos(phi(i)));
            ay=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*cos(theta(k)).*sin(phi(i)));
            A = kron(ax,ay);
            P(k,i) = W'*A;
        end        
    end
    size(A)
    patternmag=abs(P);
    patternmagnorm=patternmag/max(max(patternmag));
    patterndB=20*log10(patternmag);
    patterndBnorm=20*log10(patternmagnorm);
    figure();
     h = mesh(theta(101:200)/pi*180,phi(101:200)/pi*180,patternmag(101:200,101:200));
     set(h,'Linewidth',2)
    %plot3(theta/pi*180,phi/pi*180,patternmag)
    xlabel('俯仰(degree)');
    ylabel('方位(degree)');
    zlabel('magnitude(dB)');


%     subplot(211);
%     plot(phi*180/pi,patternmag);
%     grid on;
%     xlabel('phi/degree')
%     ylabel('amplitude/dB')
%     title(['phi方向图','来波方向为' num2str(phi0*180/pi) '度']);
%     subplot(212);
%     plot(phi,patterndBnorm,'r');
%     grid on;
%     xlabel('phi/radian')
%     ylabel('amplitude/dB')
%     title(['phi方向图','来波方向为' num2str(phi0*180/pi) '度']);
%     axis([-1.5 1.5 -50 0]);

end 

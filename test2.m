% 单个角度解析
function test2()
    global g_signal;
    global g_array;
    global g_echos;
    global g_para;
    
    theta=linspace(-pi/2,pi/2,200);
    phi=linspace(-pi/2,pi/2,200);
    theta0=15/180*pi;%来波方向
    phi0 = 25/180*pi;

    wx=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*cos(theta0).*cos(phi0));
    for  i=1:length(phi)
        ax=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*cos(theta0).*cos(phi(i)));
        px(i)=wx'*ax;
    end

    wy=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*cos(theta0).*sin(phi0));
    for  k=1:length(theta)
        ay=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*cos(theta(k)).*sin(phi0));
        py(k)=wy'*ay;
    end

    patternmag_x=abs(px);
    patternmagnorm_x=patternmag_x/max(max(patternmag_x));
    patterndB_x=20*log10(patternmag_x);
    patterndBnorm_x=20*log10(patternmagnorm_x);
    figure();
    subplot(211);
    plot(phi*180/pi,patternmag_x);
    grid on;
    xlabel('phi/degree')
    ylabel('amplitude/dB')
    title(['phi方向图','来波方向为' num2str(phi0*180/pi) '度']);
    subplot(212);
    plot(phi,patterndBnorm_x,'r');
    grid on;
    xlabel('phi/radian')
    ylabel('amplitude/dB')
    title(['phi方向图','来波方向为' num2str(phi0*180/pi) '度']);
    axis([-1.5 1.5 -50 0]);

    patternmag_y=abs(py);
    patternmagnorm_y=patternmag_y/max(max(patternmag_y));
    patterndB_y=20*log10(patternmag_y);
    patterndBnorm_y=20*log10(patternmagnorm_y);
    figure();
    subplot(211);
    plot(theta*180/pi,patternmag_y);
    grid on;
    xlabel('theta/degree')
    ylabel('amplitude/dB')
    title(['theta方向图','来波方向为' num2str(theta0*180/pi) '度']);
    subplot(212);
    plot(theta,patterndBnorm_y,'r');
    grid on;
    xlabel('theta/radian')
    ylabel('amplitude/dB')
    title(['theta方向图','来波方向为' num2str(theta0*180/pi) '度']);
    axis([-1.5 1.5 -50 0]);
end 

function beamforming(input)
 %{
    Function description:
            数字波束形成算法选择器
    Syntax：
    Long description：
            2020.03.11 加入normaldbf
%}
    if nargin<1
        input = "normal";
    end
    switch input
        case 'normal'
            normaldbf();
        otherwise 
            error('Error! Invalid input! ');
    end
end

function normaldbf()
%{
    Function description:
            普通DBF算法，用于2D-DOA估计
    Syntax：
    Long description：
            2020.03.11 初步实现DBF，单个目标波信合理，多目标出现额外波峰
%}
    global g_signal;
    global g_array;
    global g_echos;
    global g_para;
    
    theta.rad=linspace(-pi/2,pi/2,201);
    theta.num = theta.rad/pi*180;
    phi.rad=linspace(-pi/2,pi/2,201);
    phi.num = phi.rad/pi*180;

    wx=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*cos(g_echos.theta.rad).*cos(g_echos.phi.rad));
    %wx = wx(:,1) + wx(:,2);
    wy=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*cos(g_echos.theta.rad).*sin(g_echos.phi.rad));
    %wy = wy(:,1) + wy(:,2);
    W = kron(wx,wy);
    for k = 1:length(theta.rad)
        for  i = 1:length(phi.rad)
            ax=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*cos(theta.rad(k)).*cos(phi.rad(i)));
            ay=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*cos(theta.rad(k)).*sin(phi.rad(i)));
            A = kron(ax,ay);
            P(k,i) = W'*A;
        end        
    end
    patternmag=abs(P);
    patternmagnorm=patternmag/max(max(patternmag));
    patterndB=20*log10(patternmag);
    patterndBnorm=20*log10(patternmagnorm);
    figure();
    %subplot(211);
    h = mesh(theta.num(101:201),phi.num(101:201),patternmag(101:201,101:201));
    set(h,'Linewidth',2);
    xlabel('俯仰(degree)');
    ylabel('方位(degree)');
    zlabel('magnitude(dB)');
%     subplot(212);
%     mesh(theta.num(101:201),phi.num(101:201),patterndBnorm(101:201,101:201));
%     set(h,'Linewidth',2);
%     xlabel('俯仰(degree)');
%     ylabel('方位(degree)');
%     zlabel('magnitude(dB)');


%     subplot(211);
%     plot(phi.rad*180/pi,patternmag);
%     grid on;
%     xlabel('phi.rad/degree')
%     ylabel('amplitude/dB')
%     title(['phi方向图','来波方向为' num2str(g_echos.phi.rad.rad*180/pi) '度']);
%     subplot(212);
%     plot(phi.rad,patterndBnorm,'r');
%     grid on;
%     xlabel('phi.rad/radian')
%     ylabel('amplitude/dB')
%     title(['phi方向图','来波方向为' num2str(g_echos.phi.rad.rad*180/pi) '度']);
%     axis([-1.5 1.5 -50 0]);

end 

function beamforming(input)
 %{
    Function description:
            数字波束形成算法选择器
    Syntax：
    Log description：
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

    wx=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*cos(g_echos.theta.rad).*cos(g_echos.phi.rad));
    wx = sum(wx,2);%对行求和,即对各个信源的导向向量累加,作用类似wx = wx(:,1) + wx(:,2);
    wy=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*cos(g_echos.theta.rad).*sin(g_echos.phi.rad));
    wy = sum(wy,2);
    W = kron(wx,wy);%相当于将wx*wy矩阵的每一列移位到第一列末尾，形成一维向量
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
    
%     absp_theta = abs(P(:,(g_echos.phi.num+91)));
%     absp_phi = abs(P((g_echos.theta.num+91),:));
    absp_theta = sum(patternmag,1);
    absp_phi = sum(patternmag,2);
    
    figure('Name','方向图','NumberTitle','off','Color','white','Position',[600 50 550 750]);
    subplot(311);
    %patternmag = patternmag';% mesh太坑了，要先XY置换下
    %h = meshc(phi.num,theta.num,patternmag);% mesh(X,Y,Z)X,Y分别对应Z的列、行
    h = meshc(phi.num(91:181),theta.num(91:181),patternmag(91:181,91:181));% mesh(X,Y,Z)X,Y分别对应Z的列、行
    title(['二维方向图 来波方向为theta=[' num2str(g_echos.theta.num) ']' ',phi=[' num2str(g_echos.phi.num) ']' ]);
    set(h,'Linewidth',2);
    xlabel('方位角phi(degree)');
    ylabel('俯仰角theta(degree)');
    zlabel('abs of P');
    subplot(312);
    plot(theta.num(91:181),absp_theta(91:181));
    %plot(theta.num,absp_theta);
    title('俯仰角 的方向图');
    xlabel('theta(degree)');
    ylabel('abs of Ptheta');
    subplot(313);
    plot(phi.num(91:181),absp_phi(91:181));   
    %plot(phi.num,absp_phi);
    title('方位角 的方向图');
    xlabel('phi(degree)');
    ylabel('abs of Pphi');
    toc;
    disp(['普通波束形成算法用时：',num2str(toc),'s']);
end 

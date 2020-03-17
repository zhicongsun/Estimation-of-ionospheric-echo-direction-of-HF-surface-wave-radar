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

    %普通数字波束形成
    wx=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*sin(g_echos.theta.rad).*cos(g_echos.phi.rad));
    wx = sum(wx,2);%对行求和,即对各个信源的导向向量累加,作用类似wx = wx(:,1) + wx(:,2);
    wy=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*sin(g_echos.theta.rad).*sin(g_echos.phi.rad));
    wy = sum(wy,2);
    W = kron(wx,wy);%相当于将wx*wy矩阵的每一列移位到第一列末尾，形成一维向量
    X = g_echos.signal*W;
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
    title(['二维方向图 来波方向为theta=[' num2str(g_echos.theta.num) ']' ',phi=[' num2str(g_echos.phi.num) ']' ]);
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
    figure('Name','功率谱','NumberTitle','off','Color','white');
    subplot(311);
    p = meshc(phi.num(91:181),theta.num(91:181),abs_P(91:181,91:181));
    title(['二维功率谱 来波方向为theta=[' num2str(g_echos.theta.num) ']' ',phi=[' num2str(g_echos.phi.num) ']' ]);
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

    %谱峰搜索
    [temp_ptheta,theta_estimation] = max(absp_theta);
    [temp_pphi,phi_estimation] = max(absp_phi);
    theta_estimation = theta_estimation-91;
    phi_estimation = phi_estimation-91;
    disp(['估计角度为[',num2str(theta_estimation),',',num2str(phi_estimation),']']);

end 

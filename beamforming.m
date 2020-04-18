function beamforming(input)
 %{
    Function description:
            数字波束形成算法选择器
    Syntax：
    Log description：
            2020.03.11  加入normaldbf
            2020.03.19  加入capondbf
%}
    if nargin<1
        input = "normal";
    end
    switch input
        case 'normal'
            normaldbf();
        case 'capon'
            capondbf();
        otherwise 
            error('Error! Invalid input! ');
    end
end

function capondbf()
%{
    Function description:
            普通DBF算法，用于2D-DOA估计
    Syntax：
    Log description：
            2020.03.19  实现算法，相比普通波束形成提高分辨力
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

    %Capon数字波束形成
    ax=exp(j*2*pi/g_signal.lamda*g_array.x_pos.'*(sin(g_echos.theta.rad).*cos(g_echos.phi.rad)));
    ay=exp(j*2*pi/g_signal.lamda*g_array.y_pos.'*(sin(g_echos.theta.rad).*sin(g_echos.phi.rad)));
    A = khatri_rao(ay,ax);%M^2xK
    X = A*g_echos.signal;%M^2xK*KxP=M^2xP 
    X1 = awgn(X,g_echos.snr,'measured');
    Rxx = X1*X1'/g_echos.snapshot;
    R = inv(Rxx);
    if g_echos.num ~= 1
        steerA = A(:,2);%选第一个方向为目标方向，其他为干扰
    else
        steerA = A;
    end
    oneA = R*steerA/(steerA'*R*steerA);
    for k = 1:length(theta.rad)
        for  i = 1:length(phi.rad)
            wx=exp(j*2*pi/g_signal.lamda*g_array.x_pos.'*sin(theta.rad(k)).*cos(phi.rad(i)));
            wy=exp(j*2*pi/g_signal.lamda*g_array.y_pos.'*sin(theta.rad(k)).*sin(phi.rad(i)));
            W = khatri_rao(wy,wx);
            F(k,i) = W'*oneA;%1xM^2 * M^2x1 
            P(k,i) = 1/(W'*R*W);%1xM^2 * M^2xM^2 *M^2x1
        end        
    end
    toc;
    disp(['Capon波束形成算法用时：',num2str(toc),'s']);

    %画方向图
    abs_F=abs(F);
    absf_theta = sum(abs_F,2);
    absf_phi = sum(abs_F,1);
    abs_F_max = max(max(abs_F));
    abs_F = 10*log10(abs_F/abs_F_max);
    absf_theta_max = max(absf_theta);
    absf_theta = 10*log10(absf_theta/absf_theta_max);
    absf_phi_max = max(absf_phi);
    absf_phi = 10*log10(absf_phi/absf_phi_max);
    figure('Name','方向图','NumberTitle','off','Color','white','Position',[600 50 550 750]);
    subplot(311);
    f = meshc(theta.num,phi.num,abs_F);% mesh(X,Y,Z)X,Y分别对应Z的列、行，mesh太坑了，要先XY置换下
    title(['二维方向图,阵元数:' num2str(g_array.num) ',阵元间距:' num2str(g_array.span/g_signal.lamda) 'lamda,指向[theta,phi]=[' num2str(g_echos.theta.num) ',' num2str(g_echos.phi.num) ']' ]);
    set(f,'Linewidth',2);
    xlabel('方位角phi/degree');
    ylabel('俯仰角theta/degree');
    zlabel('abs of F');
    subplot(312);
    plot(theta.num,absf_theta);
    title('俯仰角 的方向图');
    xlabel('theta/degree');
    ylabel('abs of Ptheta');
    subplot(313);
    plot(phi.num,absf_phi);
    title('方位角 的方向图');
    xlabel('phi/degree');
    ylabel('abs of Pphi');
    
    %画功率谱
    abs_P = abs(P);
    absp_theta = sum(abs_P,2);%2是对列相加，得到列向量；1则是得到行向量
    absp_phi = sum(abs_P,1);
    abs_P_max = max(max(abs_P));
    abs_P = 10*log10(abs_P/abs_P_max);
    absp_theta_max = max(absp_theta);
    absp_theta = 10*log10(absp_theta/absp_theta_max);
    absp_phi_max = max(absp_phi);
    absp_phi = 10*log10(absp_phi/absp_phi_max);
    %谱峰搜索
    [temp_ptheta,theta_estimation] = max(absp_theta);
    [temp_pphi,phi_estimation] = max(absp_phi);
    theta_estimation = theta_estimation-91;
    phi_estimation = phi_estimation-91;
    disp(['估计角度为[',num2str(theta_estimation),',',num2str(phi_estimation),']']);
    
    figure('Name','功率谱','NumberTitle','off','Color','white');
    subplot(311);
    p = meshc(theta.num,phi.num,abs_P);
    if g_echos.num ~=1
        title(['二维功率谱,阵元间距:' num2str(g_array.span/g_signal.lamda) 'lamda, 实际方向为[theta,phi]=[' num2str(g_echos.theta.num) ',' num2str(g_echos.phi.num) ']'] );
    else
        title(['二维功率谱,阵元间距:' num2str(g_array.span/g_signal.lamda) 'lamda, 实际方向为[theta,phi]=[' num2str(g_echos.theta.num) ',' num2str(g_echos.phi.num) ']' ...
        '估计方向为[' num2str(theta_estimation) ',' num2str(phi_estimation) ']']);
    end
    set(p,'Linewidth',2);
    xlabel('方位角phi/degree');
    ylabel('俯仰角theta/degree');
    zlabel('abs of P');
    subplot(312);
    plot(theta.num,absp_theta);
    title('俯仰角 的功率谱');
    xlabel('theta/degree');
    ylabel('abs of Ptheta');
    subplot(313);
    plot(phi.num,absp_phi);
    title('方位角 的功率谱');
    xlabel('phi(degree)');
    ylabel('abs of Pphi');

   
    figure('Color','white');
    imagesc(phi.num,theta.num,abs_P,'CDataMapping','scaled');
    title('二维Capon算法波束扫描功率谱');
    xlabel('方位角phi/degree');
    ylabel('俯仰角theta/degree');    
    colorbar;

end 
    

function normaldbf()
%{
    Function description:
            二维Capon算法，用于2D-DOA估计
    Syntax：
    Log description：
            2020.03.11  初步实现DBF，单个目标波信合理；问题：多目标出现额外波峰，三维图XY颠倒
            2020.03.12  解决mesh存在的XY坐标颠倒问题，用sum处理多目标导向矢量，解决额外波峰
            2020.03.19  Ax和Ay用张量积处理，修正了数学模型
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
    ax=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*(sin(g_echos.theta.rad).*cos(g_echos.phi.rad)));
    ay=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*(sin(g_echos.theta.rad).*sin(g_echos.phi.rad)));
    A = khatri_rao(ay,ax);%M^2xK
    X = A*g_echos.signal;%M^2xK*KxP=M^2xP 
    X1 = awgn(X,g_echos.snr,'measured');
    Rxx = X1*X1'/g_echos.snapshot;
    if g_echos.num ~= 1
        oneA = A(:,1);%选第一个方向画方向图
    else
        oneA = A;
    end
    for k = 1:length(theta.rad)
        for  i = 1:length(phi.rad)
            wx=exp(-j*2*pi/g_signal.lamda*g_array.x_pos.'*sin(theta.rad(k)).*cos(phi.rad(i)));
            wy=exp(-j*2*pi/g_signal.lamda*g_array.y_pos.'*sin(theta.rad(k)).*sin(phi.rad(i)));
            W = khatri_rao(wy,wx);
            F(k,i) = W'*oneA;%1xM^2 * M^2x1 
            P(k,i) = W'*Rxx*W;%1xM^2 * M^2xM^2 *M^2x1
        end        
    end
    toc;
    disp(['普通波束形成算法用时：',num2str(toc),'s']);

    %画方向图
    abs_F=abs(F);
    absf_theta = sum(abs_F,2);
    absf_phi = sum(abs_F,1);
    abs_F_max = max(max(abs_F));
    abs_F = 10*log10(abs_F/abs_F_max);
    absf_theta_max = max(absf_theta);
    absf_theta = 10*log10(absf_theta/absf_theta_max);
    absf_phi_max = max(absf_phi);
    absf_phi = 10*log10(absf_phi/absf_phi_max);

    figure('Name','方向图','NumberTitle','off','Color','white','Position',[600 50 550 750]);
    subplot(311);
    f = meshc(theta.num,phi.num,abs_F);% mesh(X,Y,Z)X,Y分别对应Z的列、行，mesh太坑了，要先XY置换下
    if  g_echos.num ~=1
        title(['二维方向图,阵元数:' num2str(g_array.num) ',阵元间距:' num2str(g_array.span/g_signal.lamda) 'lamda,指向[theta,phi]=[' num2str(g_echos.theta.num(1)) ',' num2str(g_echos.phi.num(1)) ']' ]);
    else
        title(['二维方向图,阵元数:' num2str(g_array.num) ',阵元间距:' num2str(g_array.span/g_signal.lamda) 'lamda,指向[theta,phi]=[' num2str(g_echos.theta.num) ',' num2str(g_echos.phi.num) ']' ]);
    end
    set(f,'Linewidth',2);
    xlabel('方位角phi/degree');
    ylabel('俯仰角theta/degree');
    zlabel('abs of F');
    subplot(312);
    plot(theta.num,absf_theta);
    title('俯仰角 的方向图');
    xlabel('theta/degree');
    ylabel('abs of Ftheta');
    subplot(313);
    plot(phi.num,absf_phi);
    title('方位角 的方向图');
    xlabel('phi/degree/');
    ylabel('abs of Fphi');
    
    %画功率谱
    abs_P = abs(P);
    absp_theta = sum(abs_P,2);%2是对列相加，得到列向量；1则是得到行向量
    absp_phi = sum(abs_P,1);
    abs_P_max = max(max(abs_P));
    abs_P = 10*log10(abs_P/abs_P_max);
    absp_theta_max = max(absp_theta);
    absp_theta = 10*log10(absp_theta/absp_theta_max);
    absp_phi_max = max(absp_phi);
    absp_phi = 10*log10(absp_phi/absp_phi_max);

    %谱峰搜索
    [temp_ptheta,theta_estimation] = max(absp_theta);
    [temp_pphi,phi_estimation] = max(absp_phi);
    theta_estimation = theta_estimation-91;
    phi_estimation = phi_estimation-91;
    disp(['估计角度为[',num2str(theta_estimation),',',num2str(phi_estimation),']']);
    
    figure('Name','功率谱','NumberTitle','off','Color','white','Position',[600 50 550 750]);
    subplot(311);
    p = meshc(theta.num,phi.num,abs_P);
    if g_echos.num ~=1
        title(['二维功率谱,阵元间距:' num2str(g_array.span/g_signal.lamda) 'lamda, 实际方向为[theta,phi]=[' num2str(g_echos.theta.num) ',' num2str(g_echos.phi.num) ']'] );
    else
        title(['二维功率谱,阵元间距:' num2str(g_array.span/g_signal.lamda) 'lamda, 实际方向为[theta,phi]=[' num2str(g_echos.theta.num) ',' num2str(g_echos.phi.num) ']' ...
        '估计方向为[' num2str(theta_estimation) ',' num2str(phi_estimation) ']']);
    end
    set(p,'Linewidth',2);
    xlabel('方位角phi/degree');
    ylabel('俯仰角theta/degree');
    zlabel('abs of P');
    subplot(312);
    plot(theta.num,absp_theta);
    title('俯仰角 的功率谱');
    xlabel('theta/degree');
    ylabel('abs of Ptheta');
    subplot(313);
    plot(phi.num,absp_phi);
    title('方位角 的功率谱');
    xlabel('phi/degree');
    ylabel('abs of Pphi');
    
    figure('Color','white');
    imagesc(phi.num,theta.num,abs_P,'CDataMapping','scaled');
    title('二维延迟相加法波束扫描功率谱');
    xlabel('方位角phi/degree');
    ylabel('俯仰角theta/degree');    
    colorbar;

end 

function mm=khatri_rao(A,B)
    mm=[];
    n=size(A,1);
    for im=1:n
         mm=[mm;B*diag(A(im,:))];
    end
end

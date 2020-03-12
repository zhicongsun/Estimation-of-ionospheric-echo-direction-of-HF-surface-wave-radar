function doa(input)
%{
    Function description:
            2D-DOA估计算法选择器
    Syntax：
            input:doa算法类型
    Long description：
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
    Long description：
            2020.03.10 完成music2d函数，存在非理想谱峰的问题
%}
    global g_signal;
    global g_array;
    global g_echos;
    global g_para;
    
    %行是阵元个数，列是回波个数
    Ax = exp( -j*2*pi/g_signal.lamda*g_array.x_pos.'*(cos(g_echos.theta.rad).*cos(g_echos.phi.rad)));
    Ay = exp( -j*2*pi/g_signal.lamda*g_array.y_pos.'*(cos(g_echos.theta.rad).*sin(g_echos.phi.rad)));
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
            ax = exp(-j*2*pi*g_array.x_pos.'*cos(phim1)*cos(phim2));
            ay = exp(-j*2*pi*g_array.y_pos.'*cos(phim1)*sin(phim2));
            a = [ax;ay];
            SP(ang1,ang2) = 1/(a'*Un*Un'*a);
        end
    end
    SP=abs(SP);
    [linemax,linemax_id] = max(SP);
    [rowmax,rowmax_id] = max(linemax);
    rowmax_id
    SP(rowmax_id,linemax_id);
    SPmax=max(max(SP));
    SP=SP/SPmax; 
    figure(2);
    h = mesh(thet,f,SP);
    set(h,'Linewidth',2)
    xlabel('俯仰(degree)');
    ylabel('方位(degree)');
    zlabel('magnitude(dB)');
    %figure(3);
    %subplot(1,2,1);
    %plot(thet,SP);
    %xlabel('方位/degree');
    %subplot(1,2,2);
    %plot(f,SP);
    %xlabel('俯仰/degree');

end
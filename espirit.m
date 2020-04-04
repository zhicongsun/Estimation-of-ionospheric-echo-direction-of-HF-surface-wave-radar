function [rmse]=espirit(theta0,element_num)
    %*******************************************************
    % 共轭ESPRIT算法
    %
    % Inputs
    %    d_lamda       sensor separation in wavelength
    %    Rxx(K,K)      array output covariance matrix
    %    echos_num     estimated number of sources ==L=iwave
    %   
    % Output
    %    estimate  estimated angles in degrees
    %              estimated powers
    %*******************************************************
    source_number=length(theta0);%信元数
    sub_sensor_number=element_num-1;%子阵元数
    theta0_sort = sort(theta0);
    snapshot_number=1024; 
    d_lamda = 0.5;
    A=exp(j*d_lamda*2*pi*(0:element_num-1).'*sin(theta0/180*pi));
    %估计信源个数
    snr=10;
    s=sqrt(10.^(snr/10))*randn(source_number,snapshot_number);%仿真信号
    x=A*s+(1/sqrt(2))*(randn(element_num,snapshot_number)+1j*randn(element_num,snapshot_number));
    Rxx = x*x'/snapshot_number;
    [~,value]=eig(Rxx);
    value = diag(value);
    [value_sort,~] = sort(value,'descend');
    for i = 1:(size(value)-2)
        gama(i) = value_sort(i)/value_sort(i+1);
    end
    [~,esti_source_num] = max(gama);
    disp(['信噪比为10dB下估计信源数目：' num2str(esti_source_num)]);

    snr0=-10:1:10;
    rmse = zeros(1,20);
    store_doa = zeros(20,source_number);
    for isnr=1:20
        s=sqrt(10.^(snr0(isnr)/10))*randn(source_number,snapshot_number);%仿真信号
        x=A*s+(1/sqrt(2))*(randn(element_num,snapshot_number)+1j*randn(element_num,snapshot_number));
        x1=x(1:sub_sensor_number,:);
        x2 = x(2:sub_sensor_number+1,:);

        %对两个子阵的模型进行合并
        X=[x1;x2];
        R=X*X'/snapshot_number;

        %对R进行奇异值分解
        [U,S,V]=svd(R);
        R=R-S(2*sub_sensor_number,2*sub_sensor_number)*eye(2*sub_sensor_number);
        [U,S,V]=svd(R);
        Us=U(:,1:source_number);
        Us1=Us(1:sub_sensor_number,:);
        Us2=Us((sub_sensor_number+1):2*sub_sensor_number,:);
        %形成矩阵Us12
        Us12=[Us1,Us2];
        %对“Us12'*Us12”进行特征分解，得到矩阵E
        [F,Sa,Va]=svd(Us12'*Us12);
        %将E分解为四个小矩阵
        %F11=F(1:source_number,1:source_number);
        F12=F(1:source_number,(1+source_number):(2*source_number));
        %F21=F((1+source_number):(2*source_number),1:source_number);
        F22=F((1+source_number):(2*source_number),(1+source_number):(2*source_number));
        %按照公式得到旋转不变矩阵M
        E=-(F12*(inv(F22)));
        %对得到的旋转不变矩阵进行特征分解
        [V,d_lamda]=eig(E);
        d_lamda=(diag(d_lamda)).';
        doa=asin(angle(d_lamda)/pi)*180/pi;
        doa=sort(doa);
        rmse(isnr) = sqrt( sum(((theta0_sort-doa).^2))/source_number );
        i = 1:source_number;
        store_doa(isnr,i) = doa(i);
    end 
    figure('Color','white');
    plot(snr0(1:20),store_doa(1:20,1:source_number).','o-');
    grid on;
    xlabel('SNR/dB');
    ylabel('DOA 估计/度');
    title('ESPRIT 算法在不同信噪比下的DOA估计');
end

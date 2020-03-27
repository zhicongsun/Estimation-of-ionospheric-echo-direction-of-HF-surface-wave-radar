function [snr,rmse]=espirit(theta0,element_num,d_lamda)
    twpi = 2*pi;
    d=0:d_lamda:(element_num-1)*d_lamda;     
    iwave = length(theta0);    
    n = 200;
    A=exp(-j*twpi*d.'*sin(theta0));
    S=randn(iwave,n);
    snr0=0:1:30;
    for isnr=1:20
        X0=A*S;
        X=awgn(X0,snr0(isnr),'measured');
        Rxx=X*X'/n;
        [EV,D]=eig(Rxx);
        EVA=diag(D)';
        [EVA,I]=sort(EVA);
        EVA=fliplr(EVA);
        EV=fliplr(EV(:,I)); 
        estimates=(tls_esprit(d_lamda,Rxx,iwave));
        theta0_sort = sort(theta0)/pi*180;
        doaes(isnr,:)=sort(estimates(1,:));
        rmse(isnr) = sqrt( sum((theta0_sort-doaes(isnr,:)).^2)/iwave);
    end
    snr = snr0(1:20);
end

function estimate =  tls_esprit(d_lamda,Rxx,echos_num)
    %*******************************************************
    % This function calculates TLS-ESPRIT estimator for
    % uniform linear array.
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
    twpi = 2.0*pi;
    derad = pi / 180.0;
    radeg = 180.0 / pi;
    
    % eigen decomposition of Rxx
    [K,KK] = size(Rxx);
    [V,D]=eig(Rxx);
    EVA = real(diag(D)');
    [EVA,I] = sort(EVA);
    EVA=fliplr(EVA);
    EV=fliplr(V(:,I));
    
    %  composition of E_{xy} and E_{xy}^H E_{xy} = E_xys
    Exy = [EV(1:K-1,1:echos_num)...
            EV(2:K,1:echos_num)];
    E_xys = Exy'*Exy;
    
    % eigen decomposition of E_xys
    [V,D]=eig(E_xys);
    EVA_xys = real(diag(D)');
    [EVA_xys,I] = sort(EVA_xys);
    EVA_xys=fliplr(EVA_xys);
    EV_xys=fliplr(V(:,I));
    
    % decomposition of eigenvectors
    Gx = EV_xys(1:echos_num,echos_num+1:echos_num*2);
    Gy = EV_xys(echos_num+1:echos_num*2,echos_num+1:echos_num*2);
    
    % calculation of  Psi = - Gx [Gy]^{-1}
    Psi = - Gx/Gy;
    
    % eigen decomposition of Psi
    [V,D]=eig(Psi);
    EGS = diag(D).';
    [EGS,I] = sort(EGS);
    EGS=fliplr(EGS);
    EVS=fliplr(V(:,I));
     
    % DOA estimates
    ephi = atan2(imag(EGS), real(EGS));
    ange = - asin( ephi / twpi / d_lamda ) * radeg;
    estimate(1,:)=ange;
    
    % power estimates
    T = inv(EVS);
    powe = T*diag(EVA(1:echos_num) - EVA(K))*T';
    powe = abs(diag(powe).')/K;
    estimate(2,:)=powe;
end
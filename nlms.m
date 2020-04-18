function [beam,en] = nlms()

    %生成回波信号
    element_num=16;                           % the number of antenna
    theta0=[0 20];                            % DOA
    f = 4.7e6;
    f_int = 1e2;
    K=length(theta0);                         % the number of sources 
    d_lamda=0.5;                              % antenna spacing
    snapshot=500;                             % samples 
    Meann=0;  varn=1;                         % mean of noise,variance of noise                               
    SNR=20;                                   % signal-to-noise ratio 
    INR=20;                                   % interference-to-noise ratio 
    rvar1=sqrt(varn) * 10^(SNR/20);           % power of signal 
    rvar2=sqrt(varn) * 10^(INR/20);           % power of interference 
    S=[rvar1*exp(j*2*pi*(f*0.001*[0:snapshot-1]));...
        rvar2*exp(j*2*pi*(f_int*0.001*[0:snapshot-1]+rand))]; % generate the source signals 
    A=exp(j*2*pi*d_lamda*[0:element_num-1].'*sin(theta0*pi/180));   % the direction matrix  
    N=sqrt(varn/2)*(randn(element_num,snapshot)+j*randn(element_num,snapshot));       % the noise 
    X=A*S+N;                                        % the received data 

    % LMS algorithm
    de =S(1, :); 
    bata = 1;afa = 0.0001;
    w = zeros(element_num, 1); 
    e = zeros(1,element_num);
    for k = 1:snapshot 
        e(k)  = de(k) - w'*X(:,k);        
        w = w + bata * (X(:,k)*conj(e(k)))./(afa+X(:,k)'*X(:,k)); 
    end 

    % beamforming using the LMS method 
    theta=linspace(-90,90,200);
    beam=zeros(1,length(theta)); 
    for i = 1 : length(theta) 
    a=exp(j*2*pi*d_lamda*[0:element_num-1].'*sin(-pi/2 + pi*(i-1)/200)); 
    beam(i)=20*log10(abs(w'*a)); 
    end 
    
    % % plotting 
    % figure('Color','white');
    % angle=-90:180/200:(90-180/200); 
    % subplot(2,1,1);
    % plot(angle,beam); grid on
    % xlabel('方向角/degree'); 
    % ylabel('幅度响应/dB'); 
    % title('NLMS');
    for k = 1:snapshot 
        en(k)=(abs(e(k))).^2; 
    end 
    % subplot(2,1,2);
    % semilogy(en); hold on
    % xlabel('迭代次数'); 
    % ylabel('MSE');  

end
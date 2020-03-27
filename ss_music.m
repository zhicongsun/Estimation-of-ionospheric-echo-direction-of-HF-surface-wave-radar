function SP = ss_music(theta0,array_num,d_lamda,multipath_mode)
derad = pi/180;         % deg -> rad
twpi = 2*pi;
subarray_num = array_num-1;               
d=0:d_lamda:(array_num-1)*d_lamda;     
iwave = 3;              
n = 200;
A=exp(-j*twpi*d.'*sin(theta0));
if strcmp(multipath_mode,'multi_path')==1
    S0 = randn(iwave-1,n);
    S = [S0(1,:);S0];
else
    S=randn(iwave,n);
end
X0=A*S;
X=awgn(X0,10,'measured');
Rxxm=X*X'/n;
issp = 1;  

% spatial smoothing music
if issp == 1
  Rxx = ssp(Rxxm,subarray_num);
elseif issp == 2
  Rxx = mssp(Rxxm,subarray_num);
else
  Rxx = Rxxm;
  subarray_num = array_num;
end
  
% Rxx
[EV,D]=eig(Rxx);
EVA=diag(D)'; 
[EVA,I]=sort(EVA);
EVA=fliplr(EVA);
EV=fliplr(EV(:,I));

for iang = 1:361
        angle(iang)=(iang-181)/2;
        phim=derad*angle(iang);
        a=exp(-j*twpi*d(1:subarray_num)*sin(phim)).';
        L=iwave;     
        En=EV(:,L+1:subarray_num);
        SP(iang)=(a'*a)/(a'*En*En'*a);
end
   
SP=abs(SP);
SPmax=max(SP);
SP=10*log10(SP/SPmax);

end

function crs = mssp(cr, K)
    % modified spatial smoothing
    [M,MM]=size(cr);
    N=M-K+1;
    J = fliplr(eye(M));
    crfb = (cr + J*cr.'*J)/2;
    crs = zeros(K,K);
    for  in =1:N
      crs = crs + crfb(in:in+K-1,in:in+K-1);
    end
    crs = crs / N;
end

function crs = ssp(cr, K)
    % spatial smoothing 
    [M,MM]=size(cr);
    N=M-K+1;
    crs = zeros(K,K);
    for  in =1:N
            crs = crs + cr(in:in+K-1,in:in+K-1);
    end
    crs = crs / N;
end
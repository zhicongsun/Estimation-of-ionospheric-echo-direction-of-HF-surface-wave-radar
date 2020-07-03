load FrequencyData
mesh(v,r,squeeze(db(Frequenydata_result(2,:,:))))
% signal(1)=squeeze(db(Frequenydata_result(1,43,229)));
% signal(2)=squeeze(db(Frequenydata_result(2,43,229)));
% signal(3)=squeeze(db(Frequenydata_result(3,43,229)));
% signal(4)=squeeze(db(Frequenydata_result(6,43,229)));
% signal(5)=squeeze(db(Frequenydata_result(7,43,229)));
% signal(6)=squeeze(db(Frequenydata_result(8,43,229)));
signal(1)=squeeze((Frequenydata_result(1,43,229)));
signal(2)=squeeze((Frequenydata_result(2,43,229)));
signal(3)=squeeze((Frequenydata_result(3,43,229)));
signal(4)=squeeze((Frequenydata_result(6,43,229)));
signal(5)=squeeze((Frequenydata_result(7,43,229)));
signal(6)=squeeze((Frequenydata_result(8,43,229)));

X=[signal(1); signal(2); signal(3)];
Y=[signal(4); signal(5); signal(6)];
Rxx=X*X';
Ryy=Y*Y';

XY=khatri_rao(Y,X);
Rrr=XY*XY';
R=inv(Rrr);

theta.rad=linspace(-pi/2,pi/2,181);
theta.num = theta.rad/pi*180;
phi.rad=linspace(-pi/2,pi/2,181);
phi.num = phi.rad/pi*180;

freq = 4.7*10^6;     %发射（接收）信号频率4.7Mhz
lamda = (3*10^8)/freq;%发射（接收）信号波长

num = 6;                %天线阵元总个数
x_num = 3;		        %X方向阵元个数
y_num = 3;		        %Y方向阵元个数
span = lamda/2; %阵元间距
x_pos = span : span : x_num*span;
y_pos = 0 : span : (y_num-1)*span;

for k = 1:length(theta.rad)
    for  i = 1:length(phi.rad)
        wx=exp(j*2*pi/lamda*x_pos.'*sin(theta.rad(k)).*cos(phi.rad(i)));
        wy=exp(j*2*pi/lamda*y_pos.'*sin(theta.rad(k)).*sin(phi.rad(i)));
        W = khatri_rao(wy,wx);
        P(k,i) = 1/(W'*R*W);
    end        
end
abs_P = abs(P);
abs_P_max = max(max(abs_P));
abs_P = 10*log10(abs_P/abs_P_max);

%%画二维功率谱
figure('Color','white');
imagesc(phi.num,theta.num,abs_P);
hold on;
title('二维Capon算法波束扫描功率谱');
xlabel('方位角phi/degree');
ylabel('俯仰角theta/degree');    
colormap(jet);
colorbar;

function mm=khatri_rao(A,B)  %张量积
    mm=[];
    n=size(A,1);%得到A的行数
    for im=1:n
         mm=[mm;B*diag(A(im,:))];%A的第n行的每个元素和B的对应行相乘
    end%在本课题中A和B都是列向量，该函数生成的是有M^2个元素的列向量
end
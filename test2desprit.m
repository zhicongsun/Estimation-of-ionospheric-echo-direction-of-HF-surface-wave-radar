%From Improvement of 2-D DOA Estimation Based on L-shaped Array 
%使用相关性来匹配算法2
close all
clear all
clc
N=600;  %N表示迭代次数
Mx=13;%水平方向天线阵元个数
My=12;%垂直方向天线阵元个数
f=2e9;
lamda=3e8/f;
dx=lamda/2;
dy=lamda/2;
j=sqrt(-1);
theta=(-90:0.3:89.7);
fai=(0:0.3:179.7);
k=2*pi/lamda;
doa=2;
snr=10;
theta_target1=10;  %垂直方向
theta_target2=50;
theta_target3=50;
fai_target1=30; %水平方向
fai_target2=80;
fai_target3=30;
degrad=pi/180;
a_s1x=zeros(Mx,1);
a_s2x=zeros(Mx,1);
a_s3x=zeros(Mx,1);
a_s1y=zeros(My,1);
a_s2y=zeros(My,1);
a_s3y=zeros(My,1);
for mx=0:Mx-1
    a_s1x(mx+1)=exp(-j*k*mx*dx*cos(fai_target1*degrad)*sin(theta_target1*degrad));
    a_s2x(mx+1)=exp(-j*k*mx*dx*cos(fai_target2*degrad)*sin(theta_target2*degrad));
    a_s3x(mx+1)=exp(-j*k*mx*dx*cos(fai_target3*degrad)*sin(theta_target3*degrad));  
end
for my=0:My-1
    a_s1y(my+1)=exp(-j*k*my*dy*sin(fai_target1*degrad)*sin(theta_target1*degrad));
    a_s2y(my+1)=exp(-j*k*my*dy*sin(fai_target2*degrad)*sin(theta_target2*degrad));
    a_s3y(my+1)=exp(-j*k*my*dy*sin(fai_target3*degrad)*sin(theta_target3*degrad)); 
end
%% 构造均匀线阵ESPRIT的信号模型,a_s3x,a_s3y
A_x=[a_s1x,a_s2x];
A_y=[a_s1y,a_s2y];
% thetax=[cos(fai_target1*degrad)*sin(theta_target1*degrad),cos(fai_target2*degrad)*sin(theta_target2*degrad),cos(fai_target3*degrad)*sin(theta_target3*degrad)];
% thetay=[sin(fai_target1*degrad)*sin(theta_target1*degrad),sin(fai_target2*degrad)*sin(theta_target2*degrad),sin(fai_target3*degrad)*sin(theta_target3*degrad)];
% mux=k*dx*thetax;
% Phimx = diag(exp(j*mux));
% Aphimx = A_x*Phimx;
% AAx=[A_x;Aphimx];
% thetay=[sin(fai_target1*degrad)*sin(theta_target1*degrad),sin(fai_target2*degrad)*sin(theta_target2*degrad),sin(fai_target3*degrad)*sin(theta_target3*degrad)];
% muy=k*dy*thetay
% Phimy = diag(exp(j*muy));
% Aphimy = A_y*Phimy;
% AAy=[A_y;Aphimy];
% cnt=0;
% for nn=1:100
% cnt=cnt+1
% Rzx=R_zx-min(diag(Delt_x))*eye(2*Mx);


S=rand(doa,N);
X=A_x*S;
Y=A_y*S;
% Xy=[X;Y];
% XY=awgn(Xy,snr);
% X=XY(1:Mx,:);
% Y=Xy(Mx+1:Mx+My,:);
X=awgn(X,snr,'measured');
Y=awgn(Y,snr,'measured');
Y1=zeros(My,N);
X1=X(1:Mx-1,:);
X2=X(2:Mx,:);
Y1(1,:)=X(1,:);
Y1(2:My,:)=Y(1:My-1,:);
Y2=Y(1:My,:);
c1=X1*X1'/N;
c2=X1*X2'/N;
c3=X1*Y1'/N;
c4=X1*Y2'/N;
[Ec1,Deltc1]=svd(c1);
J=zeros(My,My);
L=zeros(My,My);
L(1,1)=1;
for i=1:My-1
    J(i+1,i)=1;
end
Ldelt=diag(Deltc1).';
Minsigma=sum(Ldelt(doa+1:My))/(My-doa);
C1=c1-Minsigma*eye(My);
C2=c2-Minsigma*J;
C3=c3-Minsigma*L;
C4=c4;
C=[C1;C2;C3;C4];
%% 对整体的特征矩阵进行特征分解
[Ec,Deltc] = svd(C);
E0=Ec(1:My,1:doa);
E1=Ec(My+1:2*My,1:doa);
E2=Ec(2*My+1:3*My,1:doa);
E3=Ec(3*My+1:4*My,1:doa);
Fx=inv(E0'*E0)*E0'*E1;
Fy=inv(E2'*E2)*E2'*E3;
%% 分别对x和y的特征向量进行特征分解
[Exf,Lamdax] =eig(Fx);
Lx=diag(Lamdax)';
Mux= -imag(log(Lx));
Doathetax = Mux/(k*dx);
DoaBelta=asind(Doathetax);
[Eyf,Lamday] = eig(Fy);
Ly=diag(Lamday)';
Muy=-imag(log(Ly));
Doathetay = Muy/(k*dy);
Doaelfa=asind(Doathetay);
%% 用相关性来匹配
G=Eyf'*Exf;
Gxy=abs(G);
for i=1:doa
    [Y(i),I(i)]=max(G(i,:));
end
for i=1:doa
    doatheta(i)=asind(sqrt(Doathetax(I(i))*Doathetax(I(i))+Doathetay(i)*Doathetay(i)))
    doafai(i)=atand(Doathetay(i)/Doathetax(I(i)))
end

% for j=1:nn
%     scatter(doatheta(:,j),doafai(:,j),'b.');
%     hold on
% end
% xlabel('方位角'); 
% ylabel('俯仰角');
% axis([0 180 0 90]);
% grid on;



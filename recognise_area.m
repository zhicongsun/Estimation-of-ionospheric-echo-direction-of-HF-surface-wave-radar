function recognise_area()
    I=imread('abs_p.jpg');
    I=imcrop(I,[round(103.5),round(50.5),605,532]);
    I_gray=rgb2gray(I);
    BG = imopen(I_gray,strel('disk',15));
    I2 = imsubtract(I_gray,BG);
    %求出阈值
    level = graythresh(I2)
    %与阈值进行比较，即进行黑白化
    bw = im2bw(I2,level);
    figure;
    imshow(bw)
    bw=imresize(bw,[2*181,2*181]);

    %% 1.2车牌定位
    [y,x,z]=size(bw);		%将I5的行列像素点总和分别放入变量y,x中
    I6=double(bw);
    Y1=zeros(y,1);			%定义一个y行1列的矩阵
    for i=1:y
        for j=1:x
            if(I6(i,j,1) ==1)
                Y1(i,1)=Y1(i,1)+1;
            end
        end
    end
    [temp MaxY]=max(Y1);	%计算行列像素值总和。并将最大的列所在的位置存放到MaxY中
    PY1=MaxY;				%先将最大的像素总和位置赋值给PY1
    while((Y1(PY1,1) >=50) && (PY1>1))
        PY1=PY1-1;
    end		%while循环得到和最大像素点连通区域的上边界
    PY2=MaxY;	%先将最大的像素总和位置赋值给PY2
    while((Y1(PY2,1) >=50) && (PY2<y))
        PY2=PY2+1;
    end			%while循环得到的和最大像素点连通的区域的下边界

    %求的车牌的列起始位置和终止位置
    X1=zeros(1,x);	%定义一个1行x列的矩阵
    for j=1:x
        for i=PY1:PY2
            if(I6(i,j,1)==1)
                X1(1,j)=X1(1,j)+1;
            end
        end
    end  %通过循环，求得到像素值总和最大 的位置
    PX1=1;	%先将PX1赋值为1
    while((X1(1,PX1)<3) && (PX1<x))
        PX1=PX1+1;
    end		%用while循环得到左边界
    PX2=x;	%将PX2赋值为最大值x
    while((X1(1,PX2)<3) && (PX2>PX1))
        PX2=PX2-1;
    end  %用while循环得到右边界
    PX1=PX1-1;
    PX2=PX2+1;
    PY1=PY1+15/915*(PY2-PY1);
    PX1=PX1+15/640*(PX2-PX1);
    PX2=PX2-15/940*(PX2-PX1);

    %根据经验值对车牌图像的位置进行微调
    I=imresize(I,[2*181,2*181]);
    PX1=round(PX1);
    PX2=round(PX2);
    PY1=round(PY1);
    PY2=round(PY2);

    dw=drawRect(I,[PX1,PY1],[PX2-PX1,PY2-PY1],3);
    figure('Color','white');
    imshow(dw);
    set(get(gca, 'XLabel'), 'String', '方位角：-90° ~ 90°');
    set(get(gca, 'YLabel'), 'String', '俯仰角：-90° ~ 90°');
    set(get(gca, 'Title'), 'String', '功率谱');
    px1=PX1/2-91;
    px2=PX2/2-91;
    py1=PY1/2-91;
    py2=PY2/2-91;
    text(2,8,['方位角：' num2str(px1) '~' num2str(px2)],'Color','white','FontSize',14);
    text(2,8+14,['俯仰角：' num2str(py1) '~' num2str(py2)],'Color','white','FontSize',14);
end

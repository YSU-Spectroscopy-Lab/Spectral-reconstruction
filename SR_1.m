clear;clc;close all;warning off
%% 平均最后n个光谱数据文件
file_list1 = dir('F:\课题\数据\武\SO2标准谱\2023.3.22SO2标准谱\bd\*.txt');
file_list2 = dir('F:\课题\数据\武\SO2标准谱\2023.3.22SO2标准谱\50\*.txt');
file_list3 = dir('F:\课题\数据\武\SO2标准谱\2023.3.22SO2标准谱\300\*.txt');
u0=203.5;v0=230;%测量波段
u1=202;v1=232;%拟合波段
%标注浓度，测量浓度
S_C=50;M_C=300;xishu=S_C/M_C;
len1 = size(file_list1, 1);
len2 = size(file_list2, 1);
len3 = size(file_list3, 1);
Q1=[];Q2=[];Q3=[];
%背底谱采用最后10个文件
for i1=(len1-9):1:len1
    file_list1(i1).name;
    %     fileID=fopen(file_list(i).name);(路径不全,使用strcat连接文本)
    fileID1=fopen(strcat('F:\课题\数据\武\SO2标准谱\2023.3.22SO2标准谱\bd\',file_list1(i1).name));
    A1=textscan(fileID1,'%f%f');
    %     figure
    %     plot(A1{1,1},A1{1,2});
    %     hold on;循环打开n个光谱图
    Q1=[Q1;A1];
    fclose all;
end
%平均了最后10个文件的背底谱数据
beidibochang=Q1{1,1};
a1=zeros(length(beidibochang),1);
for i2=1:length(Q1)
    a1=a1+Q1{i2,2};
end
beidiguangqiang=a1/length(Q1);
beidishuju=[beidibochang beidiguangqiang];
%差分谱1采用最后10个
for i3=(len2-9):1:len2
    file_list2(i3).name;
    fileID2=fopen(strcat('F:\课题\数据\武\SO2标准谱\2023.3.22SO2标准谱\50\',file_list2(i3).name));
    A2=textscan(fileID2,'%f%f');
    %     figure
    %     plot(A2{1,1},A2{1,2});
    %     hold on;
    Q2=[Q2;A2];
    fclose all;
end
celiangbochang1=Q2{1,1};
a2=zeros(length(celiangbochang1),1);
for i4=1:length(Q2)
    a2=a2+Q2{i4,2};
end
%平均了最后10个文件的测量谱1数据
celiangguangqiang1=a2/length(Q2);
celiangshuju1=[celiangbochang1 celiangguangqiang1];
%差分谱2采用最后10个
for i5=(len3-9):1:len3
    file_list3(i5).name;
    fileID3=fopen(strcat('F:\课题\数据\武\SO2标准谱\2023.3.22SO2标准谱\300\',file_list3(i5).name));
    A3=textscan(fileID3,'%f%f');
    %     figure
    %     plot(A3{1,1},A3{1,2});
    %     hold on;
    Q3=[Q3;A3];
    fclose all;
end
celiangbochang2=Q3{1,1};
a3=zeros(length(celiangbochang2),1);
for i6=1:length(Q3)
    a3=a3+Q3{i6,2};
end
%平均了最后10个文件的测量谱2数据
celiangguangqiang2=a3/length(Q3);
celiangshuju2=[celiangbochang2 celiangguangqiang2];
fclose all;
%% 测量波长与拟合波长选取
array=beidishuju(:,1);
value1=u0;value2=v0;value3=u1;value4=v1;
u00 = findClosestNum(array, value1);
v00 = findClosestNum(array, value2);
u11 = findClosestNum(array, value3);
v11 = findClosestNum(array, value4);
wave1=[u00,v00];
wave2=[u11,v11];
fprintf('吸收波长: ');
fprintf('%f ', wave1);
fprintf('\n');
fprintf('拟合波长: ');
fprintf('%f ', wave2);
fprintf('\n');  % 可选，添加换行符
%% 差分谱1
xishoupu2=celiangguangqiang1./beidiguangqiang;
u01 = find(celiangbochang1==u00);
v01 = find(celiangbochang1==v00);
yongdexishoupu1=xishoupu2(u01:v01,1);
yongdebochang1=celiangbochang1(u01:v01,1);%用的吸收谱波长选取

u12 = find(celiangbochang1==u11);
v12 = find(celiangbochang1==v11);
nihebochang1=celiangbochang1(u12:v12,1);
nihexishoupu1=xishoupu2(u12:v12,1);
manbianxishou=polyfit(nihebochang1,nihexishoupu1,4);
manbianxishou1=polyval(manbianxishou,nihebochang1);
manbianpu1=[nihebochang1,manbianxishou1];
u22 = find(nihebochang1==u00);
v22 = find(nihebochang1==v00);
yongdemanbian1=manbianpu1(u22:v22,:);%用的慢变拟合波段选取

kuaibianxishoupu2=yongdexishoupu1./yongdemanbian1(:,2);
chafenpu1=log(kuaibianxishoupu2);%取对数
%% 差分谱2
xishoupu2=celiangguangqiang2./beidiguangqiang;
yongdexishoupu2=xishoupu2(u01:v01,1);
yongdebochang2=celiangbochang1(u01:v01,1);%用的吸收谱波长选取

nihebochang2=celiangbochang1(u12:v12,1);
nihexishoupu2=xishoupu2(u12:v12,1);
manbianxishou=polyfit(nihebochang1,nihexishoupu2,4);
manbianxishou2=polyval(manbianxishou,nihebochang1);
manbianpu2=[nihebochang2,manbianxishou2];%用的慢变拟合波段选取

yongdemanbian2=manbianpu2(u22:v22,:);
kuaibianxishoupu2=yongdexishoupu2./yongdemanbian2(:,2);
chafenpu2=log(kuaibianxishoupu2);%取对数
%% 图1两个差分谱
y1=chafenpu1;
y2=chafenpu2;
% plot(yongdebochang1,chafenpu1,'-ro')
% hold on
% plot(yongdebochang2,chafenpu2,'-bo')
% title('两组差分谱');
% xlabel('波长λ(nm)');
% ylabel('差分吸收光谱(a.u)');
% legend('标准浓度差分谱','待测浓度差分谱')
%% 图2 插值
bochang=u0:0.02:v0;
pp1=csape(yongdebochang1,y1,'second');
pp2=csape(yongdebochang1,y2,'second');
y_0=fnval(pp1,bochang);
y_1=fnval(pp2,bochang);
zz0=find(isnan(y_0));
zz1=find(isnan(y_1));
y_0(zz0)=[];
y_1(zz1)=[];
bochang(zz0)=[];
lunwen00=[bochang',y_0'];
lunwen01=[bochang',y_1'];
% figure
% plot(bochang',y_0,'bo')
% hold on
% plot(bochang',y_1,'ro')
y_0=y_0';
y_1=y_1';
%% 差分谱选取排序
m=length(bochang);
a=1:m;
Y=[y_0 y_1 a'];
Y0=[y_0 y_1];
Y1=sortrows(Y,1);%标准谱从小到大排序，测量谱跟着变动
Y01=sortrows(Y0,1);
%% 图3 标准谱线性，测量谱变动
X=[];
X(1)=0;
for p1=2:length(Y1)
    X(p1)=100*(Y1(p1,1)-Y1(p1-1,1))+X(p1-1);
    p1=p1+1;
end
XX=X';
lunwen02=[X' Y01];
% figure
% plot(X,Y01,'-o');
%% 图4 线性标准谱，拟合测量谱
%拟合变动后的测量谱
cc=Y1(:,2)';
yn=polyfit(X,cc,8);
ynn=polyval(yn,X);
% %n倍的标准谱，与测量谱相比较
% dd=6*Y1(:,1);
% dd1=polyfit(X,dd',4);
% ynn1=polyval(dd1,X);
%Y1排序
%Y2标准谱，测量谱，
%Y3标准谱，拟合后的测量谱，
Y2=[Y1(:,1) cc'];
Y3=[Y1(:,1) ynn'];
% figure
% plot(X,Y3);
%% 图5 标准谱横坐标，测量谱做纵坐标
% figure
% plot(Y1(:,1),cc','r');
% hold on
% plot(Y1(:,1),ynn','b');%拟合完
%% 图6 逆重构
Y31=[Y1(:,1) ynn' Y1(:,3)];
Y32=sortrows(Y31,3);
Y33=[Y32(:,1) Y32(:,2)];
bochang=bochang';
% figure
% plot(bochang,y_1,'-ro');hold on
% plot(bochang,Y33(:,1),'-ko');hold on
% plot(bochang,Y33(:,2),'-bo');
% title('三组差分谱');
% xlabel('波长λ(nm)');
% ylabel('差分吸收光谱(a.u)');
% legend('待测浓度差分谱','标准浓度差分谱','逆重构差分谱')
a_1=1:length(cc);
U=[a_1' Y1(:,1) ynn'];
%% 定义滑动窗参数
MAX=[];
for N=80:20:200
% N=50;
Nn=20;
nn=(length(U)-N)/Nn;
roll=floor(nn);
%% 滑动开始
I=[];
I1=[];
for p=0:roll
    j=0;
     for i = 1:length(U)       
            if U(i,1)>=p*Nn && U(i,1)<N+p*Nn
                j=j+1;
                E(j)=i;         
            end   
     end
     c=[];
     d=[];
     c1=[];
     for o=1:j  
     newy0=U(E(o),2);newy1=U(E(o),3);
     c=[c,newy0];
     d=[d,newy1];
     c1=polyfit(c,d,1); 
     end 
       k=sum(c.*d)/sum(c.^2);
       k1=c1(1);
%      k=abs(k_1);
     I=[I k];
     I1=[I1 k1];
     K=I';
     K1=I1';
     K2=[K K1];
end
%% 图7 滑动结果分布
figure
% plot(c,d,'ro')
% plot(I*50,'-bo');%图6
% hold on
III=I1*S_C;
III=III';
plot(I1*S_C,'-ro');grid on
format compact
% x11=(sum(I(1,:))/length(I))
% y11=50*x11
% z1=(y1-200)/2.0

[maxVal, ~] = max(III); % 获取最大值
[minVal, ~] = min(III); % 获取最大值
chazhi=maxVal-minVal;
% thresh = 0.17; % 设置阈值
% % 找出比最大值小一定值的点
% pointsBelowThresh = III(III > maxVal-thresh);
% % 输出结果
% disp(pointsBelowThresh)
% Ti=length(pointsBelowThresh);
% idx=find(III>=min(pointsBelowThresh))
% xin_a=N+idx(1)*Nn
% xin_b=idx(length(idx))*Nn+N
% xin_ynn=cc(1,xin_a:xin_b);
% xin_X= X(1,xin_a:xin_b);
% % xin_U=U(xin_a:xin_b,2)
% xin_newy0=U(xin_a:xin_b,2);xin_newy1=U(xin_a:xin_b,3);
% xin_c=[];
% xin_d=[];
% xin_c=[xin_c,xin_newy0];
% xin_d=[xin_d,xin_newy1];
% xin_c1=polyfit(xin_c,xin_d,1)
% jieguo=xin_c1(1)*50
% SD=std(III);
% r=SD/(sqrt(length(III)))
% disp(maxVal)
maxVal;
fprintf('浓度: ');
fprintf('%f ', maxVal);
MAX=[MAX maxVal];
error=abs((maxVal-M_C)/M_C)*100;
fprintf('误差: ');
fprintf('%f ', error);
fprintf('\n');
if abs(MAX(1)-maxVal)>0.1
    break
end
end
fclose all;
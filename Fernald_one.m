%%本程序反演米散射单个廓线信号：输出原始信号，距离修正信号、消光系数、退偏比
close all;
clear;
clc;
t1=clock;
%% 读取数据
path = '原始数据mat文件\';
strpath = '20211029';
load([path strpath '.mat'])       %原始数据
%% 可调参数
i = 2;%选择某条廓线进行反演
ZS = 12;          %%展示高度，单位km
delt_z = 0.015;               %距离分辨率
HT = round(ZS/delt_z);
HT_EX = round(3/delt_z);
Threshold = 0.005;                   %信噪比阈值
Gain_ratio = 2;                  %增益比
Noise = 200;   %背景噪声采样点数量
MovingPoints = 11;%滑动平均数
%% 参数定义
[m,n] = size(PS);             %读数据长度
P0 = PS(1:m/2,i);           %0-5000为平行
S0 = PS(m/2+1:m,i);        %5001-10000为垂直
Height = delt_z*linspace(1,m/2,m/2);
Height = Height';
lambda = 532;%发射激光的波长
Beta_m = 1.54*10^(-3)*(532/lambda)^4*exp(-Height/7);%大气分子后向散射系数
% % beta_m = 1.54*10^(-3)*(532/lambda)^4*exp(-R/7); % 大气分子后向散射系数。(单位:km时，beta_m的单位也是km^-1)
% % beta_a = (2.47*10^(-3)*exp(-R/2)+5.13*10^(-6)*exp(-(R-20).^2/36))*(532/lambda);
Sa = 80;%气溶胶消光散射比
Sm = 8*pi/3;%大气分子消光散射比
Asr = 1.01;%气溶胶散射比
Alpha_m = Beta_m*Sm;
[mm,nn] = size(P0);
%% 扣除背景噪声
P1 = P0-mean(P0(mm-Noise:mm));  %扣除平行通道背景噪声
NOP = mean(P0(mm-Noise:mm));       %平行通道底噪
NOP = abs(NOP);
S1 = S0-mean(S0(mm-Noise:mm));  %扣除垂直通道背景噪声
NOS = mean(S0(mm-Noise:mm));       %垂直通道底噪
NOS = abs(NOS);
%% 滑动平均与距离修正
P2 = smooth(P1,MovingPoints);
SNR_P = P2/NOP;                          %平行通道信噪比
wuxiaoP = find(SNR_P<Threshold);         %索引低于信噪比的位置
% P3 = P2.*Height.^2./overlap;               %平行通道距离平方修正几何重叠因子订正
P3=P2;

S2 = smooth(S1,MovingPoints);
SNR_S = S2./NOS;                           %垂直通道信噪比
wuxiaoS = find(SNR_S<Threshold);         %索引低于信噪比的位置
% S3 = S2.*Height.^2./overlap;               %垂直通道距离平方修正几何重叠因子订正
S3=S2;

DR = Gain_ratio.*S3./P3;                %计算退偏比
DR_Denoising = DR;
DR_Denoising(wuxiaoS) = NaN;            %退偏比低信噪比置空
%% 确定参考高度,选取 （距离平方/分子后散） 的最小值
Minrange = 4;
Maxrange = 9;
PB_m = P2./Beta_m;              %PB_m = P2/Beta_m
Min_height = Minrange/delt_z;
Max_height = Maxrange/delt_z;
A_PB_m = PB_m(Min_height:Max_height);
MAX = max(max(A_PB_m));
A_PB_m(A_PB_m<=0) = MAX;
[x y] = find(A_PB_m==min(min(A_PB_m)));
Zc = floor(x+Min_height-1);
%% 计算后向散射fernald方法
Beta_a_Backward = 0*ones(Zc,1);
Beta_ac = (Asr-1)*Beta_m(Zc);   %标定高度的气溶胶后向散射系数
Beta_a_Backward(Zc) = Beta_ac;
Alpha_ac = Beta_ac*Sa;          %标定高度的气溶胶消光系数
for i = Zc:-1:2
    i;
    X(i-1) = P2(i-1);           %原始数据已经做了距离平方纠正
    A_B = (Sa-Sm)*(Beta_m(i)+Beta_m(i-1))*delt_z;%计算A1
    X(i) = P2(i);
    Beta_a_Backward(i-1) = -Beta_m(i-1)+X(i-1)*exp(A_B)/(X(i)/(Beta_a_Backward(i)+Beta_m(i))+Sa*(X(i)+X(i-1)*exp(A_B))*delt_z); 
end
Alpha_a_Backward = Beta_a_Backward*Sa;       %参考高度以下的消光
%% 计算前向散射
Beta_a_forward = 0*ones(n,1);
Beta_a_forward(Zc) = Beta_ac;
for j = Zc:1:n-1
    j;
    X(j+1) = P2(j+1);           %原始数据已经做了距离平方纠正
    A_F = (Sa-Sm)*(Beta_m(j)+Beta_m(j+1))*delt_z;%计算A1
    X(j) = P2(j);
    Beta_a_forward(j+1) = -Beta_m(j+1)+X(j+1)*exp(-A_F)/(X(j)/(Beta_a_forward(j)+Beta_m(j))-Sa*(X(j)+X(j+1)*exp(-A_F))*delt_z); 
end
Beta_a = [Beta_a_Backward;Beta_a_forward(Zc+1:n)];        %Beta_a_forward为参考高度以上的消光
Alpha_a = Beta_a*Sa; 
Alpha_a_Denoising = Alpha_a;
Alpha_a_Denoising(wuxiaoP) = NaN;            %退偏比低信噪比置空
%% 画图
ZI = 1000;
figure(1)
plot(Height(1:ZI),P0(1:ZI),Height(1:ZI),S0(1:ZI));
xlabel('Range(km)'),ylabel('signal amplitude'),title('Original signal');
grid on

figure(2)
semilogy(Height(1:ZI),P2(1:ZI));
xlabel('Range(km)');ylabel('signal amplitude');title('Square-Range-Correction');
grid on
hold on
semilogy(Height(1:ZI),S2(1:ZI));
legend('Parallel channel','vertical channel')

figure(3)
semilogx(Alpha_a_Backward(10:Zc),Height(10:Zc),'k','linewidth',2);
ylabel('Range(km)'),xlabel('km^-1'),title('Aerosol extinction coefficient');
grid on
hold on
semilogx(Alpha_m(10:Zc),Height(10:Zc),'g','linewidth',2);
legend('Aerosol','Molecules')

figure(4)
plot(DR_Denoising(5:Zc),Height(5:Zc),'linewidth',2)
ylabel('Range(km)'),xlabel('ratio'),title('Depolarization ratio');
grid on

%%����������ɢ�䵥�������źţ����ԭʼ�źţ����������źš�����ϵ������ƫ��
close all;
clear;
clc;
t1=clock;
%% ��ȡ����
path = 'ԭʼ����mat�ļ�\';
strpath = '20211029';
load([path strpath '.mat'])       %ԭʼ����
%% �ɵ�����
i = 2;%ѡ��ĳ�����߽��з���
ZS = 12;          %%չʾ�߶ȣ���λkm
delt_z = 0.015;               %����ֱ���
HT = round(ZS/delt_z);
HT_EX = round(3/delt_z);
Threshold = 0.005;                   %�������ֵ
Gain_ratio = 2;                  %�����
Noise = 200;   %������������������
MovingPoints = 11;%����ƽ����
%% ��������
[m,n] = size(PS);             %�����ݳ���
P0 = PS(1:m/2,i);           %0-5000Ϊƽ��
S0 = PS(m/2+1:m,i);        %5001-10000Ϊ��ֱ
Height = delt_z*linspace(1,m/2,m/2);
Height = Height';
lambda = 532;%���伤��Ĳ���
Beta_m = 1.54*10^(-3)*(532/lambda)^4*exp(-Height/7);%�������Ӻ���ɢ��ϵ��
% % beta_m = 1.54*10^(-3)*(532/lambda)^4*exp(-R/7); % �������Ӻ���ɢ��ϵ����(��λ:kmʱ��beta_m�ĵ�λҲ��km^-1)
% % beta_a = (2.47*10^(-3)*exp(-R/2)+5.13*10^(-6)*exp(-(R-20).^2/36))*(532/lambda);
Sa = 80;%���ܽ�����ɢ���
Sm = 8*pi/3;%������������ɢ���
Asr = 1.01;%���ܽ�ɢ���
Alpha_m = Beta_m*Sm;
[mm,nn] = size(P0);
%% �۳���������
P1 = P0-mean(P0(mm-Noise:mm));  %�۳�ƽ��ͨ����������
NOP = mean(P0(mm-Noise:mm));       %ƽ��ͨ������
NOP = abs(NOP);
S1 = S0-mean(S0(mm-Noise:mm));  %�۳���ֱͨ����������
NOS = mean(S0(mm-Noise:mm));       %��ֱͨ������
NOS = abs(NOS);
%% ����ƽ�����������
P2 = smooth(P1,MovingPoints);
SNR_P = P2/NOP;                          %ƽ��ͨ�������
wuxiaoP = find(SNR_P<Threshold);         %������������ȵ�λ��
% P3 = P2.*Height.^2./overlap;               %ƽ��ͨ������ƽ�����������ص����Ӷ���
P3=P2;

S2 = smooth(S1,MovingPoints);
SNR_S = S2./NOS;                           %��ֱͨ�������
wuxiaoS = find(SNR_S<Threshold);         %������������ȵ�λ��
% S3 = S2.*Height.^2./overlap;               %��ֱͨ������ƽ�����������ص����Ӷ���
S3=S2;

DR = Gain_ratio.*S3./P3;                %������ƫ��
DR_Denoising = DR;
DR_Denoising(wuxiaoS) = NaN;            %��ƫ�ȵ�������ÿ�
%% ȷ���ο��߶�,ѡȡ ������ƽ��/���Ӻ�ɢ�� ����Сֵ
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
%% �������ɢ��fernald����
Beta_a_Backward = 0*ones(Zc,1);
Beta_ac = (Asr-1)*Beta_m(Zc);   %�궨�߶ȵ����ܽ�����ɢ��ϵ��
Beta_a_Backward(Zc) = Beta_ac;
Alpha_ac = Beta_ac*Sa;          %�궨�߶ȵ����ܽ�����ϵ��
for i = Zc:-1:2
    i;
    X(i-1) = P2(i-1);           %ԭʼ�����Ѿ����˾���ƽ������
    A_B = (Sa-Sm)*(Beta_m(i)+Beta_m(i-1))*delt_z;%����A1
    X(i) = P2(i);
    Beta_a_Backward(i-1) = -Beta_m(i-1)+X(i-1)*exp(A_B)/(X(i)/(Beta_a_Backward(i)+Beta_m(i))+Sa*(X(i)+X(i-1)*exp(A_B))*delt_z); 
end
Alpha_a_Backward = Beta_a_Backward*Sa;       %�ο��߶����µ�����
%% ����ǰ��ɢ��
Beta_a_forward = 0*ones(n,1);
Beta_a_forward(Zc) = Beta_ac;
for j = Zc:1:n-1
    j;
    X(j+1) = P2(j+1);           %ԭʼ�����Ѿ����˾���ƽ������
    A_F = (Sa-Sm)*(Beta_m(j)+Beta_m(j+1))*delt_z;%����A1
    X(j) = P2(j);
    Beta_a_forward(j+1) = -Beta_m(j+1)+X(j+1)*exp(-A_F)/(X(j)/(Beta_a_forward(j)+Beta_m(j))-Sa*(X(j)+X(j+1)*exp(-A_F))*delt_z); 
end
Beta_a = [Beta_a_Backward;Beta_a_forward(Zc+1:n)];        %Beta_a_forwardΪ�ο��߶����ϵ�����
Alpha_a = Beta_a*Sa; 
Alpha_a_Denoising = Alpha_a;
Alpha_a_Denoising(wuxiaoP) = NaN;            %��ƫ�ȵ�������ÿ�
%% ��ͼ
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

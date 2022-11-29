clear;
close all;
clc;

s= [1	1 0 0 0 0 0 0 0 0 0;
   1	1.4142136 1 0 0 0 0 0 0 0 0;
   1	2	2 1 0 0 0 0 0 0 0;
   1	2.6131259	3.4142136	2.6131259 1 0 0 0 0 0 0;
   1	3.236068	5.236068	5.236068	3.236068 1 0 0 0 0 0;
   1	3.8637033	7.4641016	9.1416202	7.4641016	3.8637033 1 0 0 0 0;
   1	4.4939592	10.0978347	14.5917939	14.5917939	10.0978347	4.4939592 1 0 0 0;
   1	5.1258309	13.1370712	21.846151	25.6883559	21.846151	13.1370712	5.1258309 1 0 0;
   1	5.7587705	16.5817187	31.1634375	41.9863857	41.9863857	31.1634375	16.5817187	5.7587705 1 0;
   1	6.3924532	20.4317291	42.8020611	64.8823963	74.2334292	64.8823963	42.8020611	20.4317291	6.3924532 1;];


PASS_F = 300;%通带截止频率300hz
STOP_F = 200;%阻带截止频率
rp = 1;%通带衰减为1dB
rs = 20;%阻带衰减为20dB
fs = 1000;%采样频率为1000HZ
T = 1 / fs;



W_PASS = 2 * pi * PASS_F * T;
W_STOP = 2 * pi * STOP_F * T;%计算通带阻带截止频率的数字频率
omega_stop = cot(W_STOP / 2);
omega_pass = cot(W_PASS / 2);%对通带阻带截止频率进行预畸

N = ceil(log10((10^(rs / 10) - 1) / (10^(rp / 10) - 1)) / (2 * log10(omega_stop / omega_pass)));
omega_c = omega_stop / (10^(rs / 10)  - 1)^(1/(2*N));%计算截至频率（3db频率）
a = s(N,1:(N+1));%通过查表法获得系统函数的归一化系数
b = omega_c^N;%计算分母系数
sb=zeros(1,N+1);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%计算分母的系数%%%%%
for i = 1:(N+1)
    sb(i) = a(i) * omega_c^(N+1-i);%计算拉氏变换下的传递函数的分母系数
end

lenth = N+1;
c=[];
for i=1:lenth
    z_1 = triangle(i-1,'ADD');%计算1+z的幂系数
    z_2 = triangle(lenth-i,'SUB');%计算1-z的幂系数
    mult = multiplication(z_1,length(z_1),z_2,length(z_2));
    c=[c;mult];%存储每次迭代产生的系数6*6的矩阵（每一项展开后各z的幂对应的系数）
end

den=zeros(1,lenth);
for i=1:lenth%计算H(z)分母的系数
    for j=1:lenth
    den(j)=den(j)+(sb(i)*c(i,lenth+1-j));
    end
end
dens=den;
for i=1:lenth
    dens(i)=dens(i)/den(1);%获得分母的归一化系数dens，系数按照幂数从高到低排列
end
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%计算分子的系数%%%%%
num=zeros(1,lenth);

z_3=triangle(N,'SUB');%利用杨辉三角计算分子(1-z^-1)^N的系数
for i=1:lenth
    num(i)=b*(z_3(lenth+1-i));
end
nums=num;
for i=1:lenth
    nums(i)=nums(i)/den(1);%获得分母的归一化系数dens，系数按照幂数从高到低排列
end
%%%%%%%%%%画出图像%%%%%%%%
disp('获得的巴特沃斯滤波器的分子系数为');
disp(nums);
disp('获得的巴特沃斯滤波器的分母系数为');
disp(dens);
[h,f] = freqz(nums,dens,1024,fs);
figure(1);
plot(f,20*log10(abs(h)),200,-20,'*',300,0,'bo');
axis([0,750,-200,20])
grid on;
title(['通带截止频率为',num2str(PASS_F),'hz','阻带截至频率为',num2str(STOP_F),'hz'])
xlabel('频率/Hz')
ylabel('振幅/dB')
strings={'巴特沃斯滤波器的分子系数为';num2str(nums);'获得的巴特沃斯滤波器的分母系数为';num2str(dens)};
annotation('textbox',[0.2,0.2,0.3,0.1],'LineStyle','-','LineWidth',2,'String',strings)
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%将各多项式对应系数相乘，获得对应的各幂数的系数%%%
function result = multiplication(z1,n,z2,nextn)
    temp = [z1,zeros(1,nextn-1)];
    result = zeros(1,n+nextn-1);

    for j=1:nextn
        for i=j:(n+nextn-1)
            result(i)=result(i)+temp(i-j+1)*z2(j);
        end
    end
    return
end
%%%%%计算杨辉三角的系数%%%
function result = triangle(n,fun)
    result=zeros(1,n+1);%建立数据的数组
    result(1)=1;
    if n==0
        return
    elseif n==1
        if(strcmp(fun,'ADD'))
            result(2)=1;
        else
            result(1)=-1;%如果是减号，则第二项系数是-1
            result(2)=1;
        end
        return
    end

    len = n + 1;%数组长度
    temp = zeros(1,len);%定义中间变量
    temp(1) = 1;
    temp(2) = 1;

    for i=3:(n+1)
        result(i) = 1;
        for j=2:i
            result(j) = temp(j-1)+temp(j);
        end
        if i==(n+1)%最后一次不需要为中间变量赋值
           if(strcmp(fun,'SUB'))
              for k=1:len
                  result(k)=result(k)*((-1)^(len-k));
              end
           end
           return
        end
        for j=2:i
            temp(j)=result(j);
        end
    end
end

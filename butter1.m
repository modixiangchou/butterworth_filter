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


PASS_F = 300;%ͨ����ֹƵ��300hz
STOP_F = 200;%�����ֹƵ��
rp = 1;%ͨ��˥��Ϊ1dB
rs = 20;%���˥��Ϊ20dB
fs = 1000;%����Ƶ��Ϊ1000HZ
T = 1 / fs;



W_PASS = 2 * pi * PASS_F * T;
W_STOP = 2 * pi * STOP_F * T;%����ͨ�������ֹƵ�ʵ�����Ƶ��
omega_stop = cot(W_STOP / 2);
omega_pass = cot(W_PASS / 2);%��ͨ�������ֹƵ�ʽ���Ԥ��

N = ceil(log10((10^(rs / 10) - 1) / (10^(rp / 10) - 1)) / (2 * log10(omega_stop / omega_pass)));
omega_c = omega_stop / (10^(rs / 10)  - 1)^(1/(2*N));%�������Ƶ�ʣ�3dbƵ�ʣ�
a = s(N,1:(N+1));%ͨ��������ϵͳ�����Ĺ�һ��ϵ��
b = omega_c^N;%�����ĸϵ��
sb=zeros(1,N+1);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%�����ĸ��ϵ��%%%%%
for i = 1:(N+1)
    sb(i) = a(i) * omega_c^(N+1-i);%�������ϱ任�µĴ��ݺ����ķ�ĸϵ��
end

lenth = N+1;
c=[];
for i=1:lenth
    z_1 = triangle(i-1,'ADD');%����1+z����ϵ��
    z_2 = triangle(lenth-i,'SUB');%����1-z����ϵ��
    mult = multiplication(z_1,length(z_1),z_2,length(z_2));
    c=[c;mult];%�洢ÿ�ε���������ϵ��6*6�ľ���ÿһ��չ�����z���ݶ�Ӧ��ϵ����
end

den=zeros(1,lenth);
for i=1:lenth%����H(z)��ĸ��ϵ��
    for j=1:lenth
    den(j)=den(j)+(sb(i)*c(i,lenth+1-j));
    end
end
dens=den;
for i=1:lenth
    dens(i)=dens(i)/den(1);%��÷�ĸ�Ĺ�һ��ϵ��dens��ϵ�����������Ӹߵ�������
end
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%������ӵ�ϵ��%%%%%
num=zeros(1,lenth);

z_3=triangle(N,'SUB');%����������Ǽ������(1-z^-1)^N��ϵ��
for i=1:lenth
    num(i)=b*(z_3(lenth+1-i));
end
nums=num;
for i=1:lenth
    nums(i)=nums(i)/den(1);%��÷�ĸ�Ĺ�һ��ϵ��dens��ϵ�����������Ӹߵ�������
end
%%%%%%%%%%����ͼ��%%%%%%%%
disp('��õİ�����˹�˲����ķ���ϵ��Ϊ');
disp(nums);
disp('��õİ�����˹�˲����ķ�ĸϵ��Ϊ');
disp(dens);
[h,f] = freqz(nums,dens,1024,fs);
figure(1);
plot(f,20*log10(abs(h)),200,-20,'*',300,0,'bo');
axis([0,750,-200,20])
grid on;
title(['ͨ����ֹƵ��Ϊ',num2str(PASS_F),'hz','�������Ƶ��Ϊ',num2str(STOP_F),'hz'])
xlabel('Ƶ��/Hz')
ylabel('���/dB')
strings={'������˹�˲����ķ���ϵ��Ϊ';num2str(nums);'��õİ�����˹�˲����ķ�ĸϵ��Ϊ';num2str(dens)};
annotation('textbox',[0.2,0.2,0.3,0.1],'LineStyle','-','LineWidth',2,'String',strings)
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%��������ʽ��Ӧϵ����ˣ���ö�Ӧ�ĸ�������ϵ��%%%
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
%%%%%����������ǵ�ϵ��%%%
function result = triangle(n,fun)
    result=zeros(1,n+1);%�������ݵ�����
    result(1)=1;
    if n==0
        return
    elseif n==1
        if(strcmp(fun,'ADD'))
            result(2)=1;
        else
            result(1)=-1;%����Ǽ��ţ���ڶ���ϵ����-1
            result(2)=1;
        end
        return
    end

    len = n + 1;%���鳤��
    temp = zeros(1,len);%�����м����
    temp(1) = 1;
    temp(2) = 1;

    for i=3:(n+1)
        result(i) = 1;
        for j=2:i
            result(j) = temp(j-1)+temp(j);
        end
        if i==(n+1)%���һ�β���ҪΪ�м������ֵ
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

% ��ʼ��
clear;close all;clc;

% ���������ļ�
Rrs=load('data\Rrs.txt');      % ң�з�����
aw=load('data\aw.txt');        % ��ˮ������ϵ��
bw=load('data\bw.txt');        % ����ɶ��
bands=load('data\bands.txt');  % ���β�����ֵ


[n_sample,n_band] = size(Rrs);  % ���ݵ���������

% ��������õ��Ĳ�����ֵ
% ���û�иò��Σ������������δ��棬��ֱ���޸ĺ����ֵ
B412 = 412;
B443 = 443;
B490 = 488;
B560 = 555;
B645 = 645;
B665 = 667;
S0=0.015;

% ָ��ÿ�����ε��к�
i_412=find(bands==B412);
i_443=find(bands==B443);
i_490=find(bands==B490);
i_560=find(bands==B560);
i_645=find(bands==B645);
i_665=find(bands==B665);

rrs = Rrs ./ (0.52+1.7.*Rrs); %step 0 ����rrs

g0 = 0.089; g1 = 0.1245;
u = (-g0 + (g0^2 + 4*g1.*rrs).^0.5 ) / (2 * g1); %step 1 ����u


B0 = zeros(n_sample,1);     % ���� ��0 ��ֵ B560 �� B665
bbp_B0 = zeros(n_sample,1); % ������������ bbp(��0)
aph_QAA=ones(n_sample,1);
bbp_QAA=ones(n_sample,1);
a_QAA=ones(n_sample,1);
adg_QAA=ones(n_sample,1);
slope_adg=ones(1,n_band);

% �ж�ÿ������ Rrs(665) < 0.0015 λ�õ�ֵ
for n = 1:n_sample
    if Rrs(n,i_665) < 0.0015
        B0(n) = B645;
        i_B0 = i_645;
        p1 = rrs(n,i_443) + rrs(n,i_490);
        p2 = rrs(n,i_645) + 5 * rrs(n,i_665) ./rrs(n,i_490) .* rrs(n,i_665);
        x = log10(p1./p2);  % ��Ϊlog10
        h0 = -1.146;
        h1 = -1.366;
        h2 = -0.469;
        a_645 = aw(i_645) + 10^(h0+h1*x+h2*x.^2); % step 2
        a_B0 = a_645;
        
    else
        B0(n) = B645;
        i_B0 = i_645;
        a_645=aw(i_645) + 0.39*(Rrs(n,i_645)/(Rrs(n,i_443)+Rrs(n,i_490)))^1.14; % step 2
        a_B0 = a_645;
    end
    
    bbp_B0(n) = u(n,i_B0) * a_B0 / (1-u(n,i_B0)) - bw(i_B0); % step 3 ������� bbp(��0)
end

e=exp(-0.9*rrs(:,i_443)./rrs(:,i_645));
g=2.0*(1-1.2*e); % step 4

bbp=bbp_B0.*(B0./bands).^g; % step 5 ����bbp

a = (1-u).*(bw+bbp)./u; % step 6 ����a

Zeta=0.74+0.2/(0.8+rrs(n,i_443)/rrs(n,i_645)); % step 7
    
S=S0+0.002/(0.6+rrs(n,i_443)/rrs(n,i_645)); % step 8
    
Xi=exp(S*(443-411)); % step 8

adg443=zeros(n_sample,1);
adg443(:,1)=((a(:,i_412)-Zeta*a(:,i_443))-(aw(:,i_412)-Zeta*aw(:,i_443)))/(Xi-Zeta); % step 9


adg=adg443.*exp(-S*(bands-443)); % step 8

%for n=n_sample:1
%   adg_B0(n)=adg443(n,1)*exp(-S*(B0(n)-443)); % step 8
%end

aph=a-adg-aw; % step 10
    
%aph_QAA(i,:)=aph;
%bbp_QAA(i,:)=bbp;
%a_QAA(i,:)=a;
%adg_QAA(i,:)=adg;
%slope_adg(n)=S;



%figure,mesh(bbp);
%figure,mesh(a);
figure,mesh(adg),title('adg');

figure,mesh(aph),title('aph')



Save2file('data\out_bbp.txt',bbp);
Save2file('data\out_a.txt',a);
Save2file('data\out_aph.txt',aph);
Save2file('data\out_adg.txt',adg);



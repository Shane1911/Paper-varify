clc
close all
clear 
% �������������������������ķָ��ߡ������������������� %
% Tep=30;
Tep_b=7.335;
Tep_r=6.24;
r_b=3.65;
R_i=29.65;
alpha_b=3.2e-6;
alpha_r=12.5e-6;
expansion_b=alpha_b*r_b.*Tep_b;
sprintf('������������Ϊ: %0.2d mm',expansion_b)
expansion_r=alpha_r*R_i.*Tep_r;
sprintf('��Ȧ��������Ϊ: %0.2d mm',expansion_r)
y=[];
for i=1:length(Tep_b)
    y(i,:)=[expansion_b(i) expansion_r(i)];
end
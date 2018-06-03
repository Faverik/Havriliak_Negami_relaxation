clear all;
close all;
clc;
name0='F';
doc_num=60;
strok=4;
freq_num=58;
% ��������� ���������� ��� ����� � ���������� ���������, �� ������������
% �������� �������� � ��������� �������� ������ �� ���
%name0='F'; %��������� ����� �������� �����
% ������ ��������� �������� ����� �� ���������� � ��������� �����. ������:
% TEMP23.txt
%doc_num=60;

All=cell(1, (doc_num+1));
for q=1:(doc_num+1) %59 - ������������ ����� ���������+1
   name = sprintf('%s%d.TXT',name0, q-1);
   BigM = dlmread(name, '\t', strok, 0); % 4 - � ����� ������ ���������
   All{q}=BigM; %������� �������, � ������� ������ ��������� �� ����� ������� ��� ���� ��������� �������
end
   

xx=1;
p=3;  %������� �������� ��������
for t=1:1:61  %������ �������� �� ������ ������� ������ p-��� �������
    TempAll=All{t};
    for p=7         %����� �������� ������ 5,6 � 9, � �������, ���� ������ �������� ��� ������
       for y=1:1:freq_num  % �� ���������� ������
           TempEPS(y,xx)=TempAll(y, p);
       end
       xx=xx+1;
    end
end   
All_data=TempEPS;

dlmwrite('matricaEPS2_00.xls', TempEPS, '\t');
    

clc;
close all;
clear all;
% x=[220, 250, 280, 290, 300, 310,320, 330, 340, 380];
% x=x';
% for i=1:10
% % ryetry(i)=10-(266000000/((x(i)-253)*8.314));
% logt(i)=12-(2660*150/(((x(i))-150)*8.314));
% end
% %ryetry=ryetry';
% logt=logt';
% plot(x, logt);

ftype = fittype('logta-(Ea/(x*8.314))'); % вид функции (полином всех фиттингов)
Data1=[123
128
133
138
143
148
153
158
163
168
173
178
183
188
193
198
203
208
213
218
223
];
Data3=[-0.310179881
-0.027355191
0.429780755
0.947311028
1.467067876
1.957852268
2.371058927
2.743383883
3.111825784
3.422992267
3.429910449
3.758157168
4.045848633
4.334094512
4.595353248
4.84556851
5.031885801
5.246750315
5.427495389
6.029389856
6.236577223
];

[f, l ,k]=fit(Data1,Data3,ftype, 'StartPoint', [5000, 10], 'TolFun' , 0.0000000000000000000000000000000001);

%plot(f, Data1,Data3);

coeffvals1 = coeffvalues(f);
Ea =coeffvals1(1);

logta=coeffvals1(2);

for i=1:size(Data1)
    tau(i)=logta-(Ea/((Data1(i)*8.314)));
end
tau=tau';
Data1=1000*(Data1').^-1;

figure('Color','w');
hold on;

set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','Arilal Black');  %формат подписей осей
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arilal Black'); 


plot(Data1, Data3,'ok', Data1, tau);


f
title('ishodnij');
%print('-dpng','-r720','ishodnijArren');
text_all=sprintf('%s', Ea);
%text(3, 5, text_all, 'EdgeColor','k', 'LineWidth', 2);








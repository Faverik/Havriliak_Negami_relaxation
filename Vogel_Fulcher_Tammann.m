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

ftype = fittype('logta-(D*T0/((x-T0)))'); % вид функции (полином всех фиттингов)
Data1=[273
278
283
288
293
298
303
308
313
318
323
328
333
338
];
Data3=[0.3393675634
0.929771828
1.365604166
2.108090449
2.64721501
3.051113528
3.479967135
3.732261809
4.122235813
4.346345993
4.620020407
4.770041739
4.979209127
5.125708349
];

[f, l ,k]=fit(Data1,Data3,ftype, 'StartPoint', [5, 200, 10], 'TolFun' , 0.0000000000000000000000000000000001);

%plot(f, Data1,Data3);

coeffvals1 = coeffvalues(f);
D =coeffvals1(1);
T0 =coeffvals1(2);
logta=coeffvals1(3);
Data11=[
273
278
283
288
293
298
303
308
313
318
323
328
333
338
];

for i=1:size(Data11)
    tau(i)=logta-(D*T0/((Data11(i)-T0)));
end
tau=tau';
Data11=1000*(Data11').^-1;
Data1=1000*(Data1').^-1;
figure('Color','w');
hold on;

set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','Arilal Black');  %формат подписей осей
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arilal Black'); 


plot(Data1, Data3,'ok', Data11, tau);
Tg=232
Ea=(D*T0*8.314)/(1-(T0/Tg))


f
title('ishodnij');
%print('-dpng','-r720','ishodnijVFT');
text_all=sprintf('%s', Ea);
%text(3, 5, text_all, 'EdgeColor','k', 'LineWidth', 2);








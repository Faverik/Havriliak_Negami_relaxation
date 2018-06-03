% To use the program, place it in the data folder of the Novocontrol with
% txt files. 
% To get the correct result, you must enter variables and change the values
% of coeff_data which is fitting bounds. For more info print "help fit" in
% command window. If necessary, you can contact questions by e-mail:
% nefaeron@gmail.com
% This programm is a part of my PhD, so I made it very quickly, 
% if you have any difficulties with using it, write to me. 
% Soon I will try to create an improved version and, perhaps, 
% add a graphical interface. Possibly, the next version will be 
% made with a Python 3.0 and will become more accessible.

% Cleaning
clc;
close all;
clear all;
% Variables
peak=3; % peaks number (2 or 3)
name0='F'; %text part of filename (F25.txt of T12.txt)
doc_num=60; % max number of document
strok=4; %row number, with text
sample_name='alumin'; % name of sample
num_doc_1=38; %number of first document
num_doc_2=39; %number of last document

All_data=ReadEPSall(name0, doc_num, strok);
eps_end=size(All_data, 2);

% In case of 3 peaks
if peak==3
ftype = fittype('abs(imag(deps1/(complex(1, ((1*10^tay1)*x*2*pi)^alfa1))^betta1))+abs(imag(deps2/(complex(1, ((1*10^tay2)*x*2*pi)^alfa2))^betta2))+abs(imag(deps3/(complex(1, ((1*10^tay3)*x*2*pi)^alfa3))^betta3))'); % âèä ôóíêöèè (ïîëèíîì âñåõ ôèòòèíãîâ)
xx=2;
coeff_table(1:12,1)=[0.177212341815809;0.454861730781002;0.392354999514640;0.969043302439048;0.820831659025781;0.700000226327451;0.336345910800798;0.00689361759145735;0.0101876524032858;1.5;-1;-4];
for plot_num=(num_doc_1*2-1):2:(num_doc_2*2-1)
    Data1=All_data(1:end, (eps_end-(plot_num+1)+1));
    Data3=All_data(1:end, (eps_end-(plot_num)+1));
    for i=1:(size(Data1)-1)
        Ddata1(i,1)=Data1(i+1)-Data1(i);
        Ddata3(i,1)=Data3(i+1)-Data3(i);
        diffur(i,1)=-(Ddata3(i)/Ddata1(i)-1);
    end

Data111=Data1;
Data111(1)=[];
figure;
loglog(Data111, diffur);
% ncoeffs = numcoeffs(ftype) % åñëè íàäî óçíàòü ÷èñëî êîýôôèöèåíòîâ â ô-öè
% coeffs = coeffnames(ftype) % åñëè íàäî óçíàòü èìåíà êîýôôèöèåíòîâ â ô-öè
    
[f, l ,k]=fit(Data1,Data3,ftype, 'StartPoint',...
    [coeff_table(1, xx-1), coeff_table(2, xx-1), coeff_table(3, xx-1),...
    coeff_table(4, xx-1), coeff_table(5, xx-1), coeff_table(6, xx-1),...
    coeff_table(7, xx-1), coeff_table(8, xx-1), coeff_table(9, xx-1),...
    coeff_table(10, xx-1), coeff_table(11, xx-1), coeff_table(12, xx-1)],...
    'Lower', [0.05, 0.25, 0.05, 0.7, 0.7, 0.7,0.000006, 0.000006, 0.000006, -2, -6, -15],...
    'Upper',[1, 1, 1, 1, 1, 1,5, 5, 5, coeff_table(10, xx-1), coeff_table(11, xx-1), coeff_table(12, xx-1)], 'TolFun' , 0.0000000000000000000000001, 'MaxIter', 2000);

coeffvals1 = coeffvalues(f);
alfa1 =coeffvals1(1);
alfa2 =coeffvals1(2);
alfa3 =coeffvals1(3);
betta1 =coeffvals1(4);
betta2 =coeffvals1(5);
betta3 =coeffvals1(6);
deps1 =coeffvals1(7);
deps2 =coeffvals1(8);
deps3 =coeffvals1(9);
tay1 =coeffvals1(10);
tay2 =coeffvals1(11);
tay3 =coeffvals1(12);
       
Fit_plot=zeros(length(Data1), 1);
for i=1:length(Data1)
Fit_plot(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1))...
    +abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2))...
    +abs(imag(deps3/(complex(1, ((1*10^tay3)*Data1(i)*2*pi)^alfa3))^betta3));
end


for ut=1:numcoeffs(ftype)
coeff_table(ut, xx)=coeffvals1(ut);
end
%coeff_table(numcoeffs(ftype)+1, xx)=doc_num-num_doc_1-(xx-2)+1;
figure('Color','w');
    
set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','Arilal Black');  %ôîðìàò ïîäïèñåé îñåé
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arilal Black'); 
semilogx(Data1, Data3, 'ko', Data1, Fit_plot);
xlim([10^-1,10^6]);
xlabel('frequency');
ylabel('\epsilon''''', 'fontsize',20);
hold on;

Fit_plot1=zeros(1, length(Data1));
Fit_plot2=zeros(1, length(Data1));
Fit_plot3=zeros(1, length(Data1));

for i=1:length(Data1)
Fit_plot1(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1));
end

for i=1:length(Data1)
Fit_plot2(i)=abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2));
end

for i=1:length(Data1)
Fit_plot3(i)=abs(imag(deps3/(complex(1, ((1*10^tay3)*Data1(i)*2*pi)^alfa3))^betta3));
end

plot(Data1, Fit_plot1,'--', Data1, Fit_plot2,'--', Data1, Fit_plot3,'--');
hold off;
name_graf=sprintf('%s_graph%d',sample_name , xx-1 );
%print('-dpng','-r720',name_graf);
xx=xx+1;

end
coeff_table(:, 1)=[];

coeff_table2(1:12,1)=coeff_table(1:12, xx-2);
yy=2;
for plot_num=1:2:2*(num_doc_2-num_doc_1+1)
    eps_middle=eps_end-((num_doc_2*2)-1);
        
    Data1=All_data(1:end, (eps_middle+plot_num-1));
    Data3=All_data(1:end, (eps_middle+plot_num));
    


% ncoeffs = numcoeffs(ftype) % åñëè íàäî óçíàòü ÷èñëî êîýôôèöèåíòîâ â ô-öè
% coeffs = coeffnames(ftype) % åñëè íàäî óçíàòü èìåíà êîýôôèöèåíòîâ â ô-öè

    
[f, l ,k]=fit(Data1,Data3,ftype, 'StartPoint',...
    [coeff_table2(1, yy-1), coeff_table2(2, yy-1), coeff_table2(3, yy-1),...
    coeff_table2(4, yy-1), coeff_table2(5, yy-1), coeff_table2(6, yy-1),...
    coeff_table2(7, yy-1), coeff_table2(8, yy-1), coeff_table2(9, yy-1),...
    coeff_table2(10, yy-1), coeff_table2(11, yy-1), coeff_table2(12, yy-1)],...
    'Lower', [0.10, 0.10, 0.10, 0.7, 0.7,0.7,0.00000006,0.00000006, 0.00000006, coeff_table2(10, yy-1), coeff_table2(11, yy-1), coeff_table2(12, yy-1)],...
    'Upper',[1, 1, 1, 1,1, 1,3, 3, 3, 4, -0.5, -4], 'TolFun' , 0.000000000000000000000000000000001, 'MaxIter', 2000);

coeffvals1 = coeffvalues(f);
alfa1 =coeffvals1(1);
alfa2 =coeffvals1(2);
alfa3 =coeffvals1(3);
betta1 =coeffvals1(4);
betta2 =coeffvals1(5);
betta3 =coeffvals1(6);
deps1 =coeffvals1(7);
deps2 =coeffvals1(8);
deps3 =coeffvals1(9);
tay1 =coeffvals1(10);
tay2 =coeffvals1(11);
tay3 =coeffvals1(12);

       
Fit_plot=zeros(length(Data1), 1);
for i=1:length(Data1)
Fit_plot(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1))+abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2))+abs(imag(deps3/(complex(1, ((1*10^tay3)*Data1(i)*2*pi)^alfa3))^betta3));
end


for ut=1:numcoeffs(ftype)
coeff_table2(ut, yy)=coeffvals1(ut);
end
%coeff_table(numcoeffs(ftype)+1, xx)=doc_num-num_doc_1-(xx-2)+1;
figure('Color','w');
    
set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','Arilal Black');  %ôîðìàò ïîäïèñåé îñåé
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arilal Black'); 
semilogx(Data1, Data3, 'ko', Data1, Fit_plot);
xlim([10^-1,10^6]);
xlabel('frequency');
ylabel('\epsilon''''', 'fontsize',20);
hold on;

for i=1:length(Data1)
Fit_plot1(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1));
end

for i=1:length(Data1)
Fit_plot2(i)=abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2));
end

for i=1:length(Data1)
Fit_plot3(i)=abs(imag(deps3/(complex(1, ((1*10^tay3)*Data1(i)*2*pi)^alfa3))^betta3));
end

plot(Data1, Fit_plot1,'--', Data1, Fit_plot2,'--', Data1, Fit_plot3,'--');
hold off;
name_graf=sprintf('%s_graph%d',sample_name , yy-1 );
%print('-dpng','-r720',name_graf);
yy=yy+1;

end
coeff_table2(:, 1)=[];
Temp_coeff_table=coeff_table2;
close all;
xx=1;
uu=2;
average_coeff=zeros(numcoeffs(ftype), size(coeff_table, 2));

for i=1:size(coeff_table, 2)
    for j=1:numcoeffs(ftype)
average_coeff(j, i)=(coeff_table2(j, (num_doc_2-num_doc_1-i+2))+coeff_table(j, i))/2;
    end
end
coeff_table3(1:12, 1)=average_coeff(1:12, 1);
for plot_num=(num_doc_1*2-1):2:(num_doc_2*2-1)
    
    %Data1=All_data(1:end, (eps_end-(plot_num+2)+1));
    Data1=All_data(1:end, (eps_end-(plot_num+1)+1));
    Data3=All_data(1:end, (eps_end-(plot_num)+1));
% ncoeffs = numcoeffs(ftype) % åñëè íàäî óçíàòü ÷èñëî êîýôôèöèåíòîâ â ô-öè
% coeffs = coeffnames(ftype) % åñëè íàäî óçíàòü èìåíà êîýôôèöèåíòîâ â ô-öè
[f, l ,k]=fit(Data1,Data3,ftype, 'StartPoint',...
    [average_coeff(1, xx), average_coeff(2, xx), average_coeff(3, xx),...
    average_coeff(4, xx), average_coeff(5, xx), average_coeff(6, xx),...
    average_coeff(7, xx), average_coeff(8, xx), average_coeff(9, xx), average_coeff(10, xx), average_coeff(11, xx), average_coeff(12, xx)],...
    'Lower', [0.10, 0.10,0.10, 0.7,0.7, 0.7, 0.00000006,0.00000006, 0.00000006, -1, -5, -15],...
    'Upper',[1, 1, 1, 1,1, 1,3, 3, 3, coeff_table3(10, uu-1), coeff_table3(11, uu-1), coeff_table3(12, uu-1)], 'TolFun' , 0.000000000000000000000000000000001, 'MaxIter', 2000);

coeffvals1 = coeffvalues(f);
alfa1 =coeffvals1(1);
alfa2 =coeffvals1(2);
alfa3 =coeffvals1(3);
betta1 =coeffvals1(4);
betta2 =coeffvals1(5);
betta3 =coeffvals1(6);
deps1 =coeffvals1(7);
deps2 =coeffvals1(8);
deps3 =coeffvals1(9);
tay1 =coeffvals1(10);
tay2 =coeffvals1(11);
tay3 =coeffvals1(12);
       
Fit_plot=zeros(length(Data1), 1);
for i=1:length(Data1)
Fit_plot(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1))+abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2))+abs(imag(deps3/(complex(1, ((1*10^tay3)*Data1(i)*2*pi)^alfa3))^betta3));
end
for ut=1:numcoeffs(ftype)
coeff_table3(ut, uu)=coeffvals1(ut);
end
coeff_table(numcoeffs(ftype)+1, xx)=doc_num-num_doc_1-(xx-2)+1;
figure('Color','w');
    
set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','Arilal Black');  %ôîðìàò ïîäïèñåé îñåé
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arilal Black'); 
semilogx(Data1, Data3, 'ko', Data1, Fit_plot);
xlim([10^-1,10^6]);
xlabel('frequency');
ylabel('\epsilon''''', 'fontsize',20);
hold on;

for i=1:length(Data1)
Fit_plot1(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1));
end

for i=1:length(Data1)
Fit_plot2(i)=abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2));
end

for i=1:length(Data1)
Fit_plot3(i)=abs(imag(deps3/(complex(1, ((1*10^tay3)*Data1(i)*2*pi)^alfa3))^betta3));
end

plot(Data1, Fit_plot1,'--', Data1, Fit_plot2,'--', Data1, Fit_plot3,'--');
hold off;
name_graf=sprintf('%s_graph%d',sample_name , yy-1 );
%print('-dpng','-r720',name_graf);
xx=xx+1;
uu=uu+1;

end
end
%In case of 2 peaks
if peak==2
    ftype = fittype('abs(imag(deps1/(complex(1, ((1*10^tay1)*x*2*pi)^alfa1))^betta1))+abs(imag(deps2/(complex(1, ((1*10^tay2)*x*2*pi)^alfa2))^betta2))'); % âèä ôóíêöèè (ïîëèíîì âñåõ ôèòòèíãîâ)

xx=2;
%coeff_table=zeros(numcoeffs(ftype), eps_end/2);
coeff_table(1:8,1)=[0.200000000000022;0.999999968570415;0.998537426987678;0.700000034599178;0.0383380820780626;1.47273633387733e-03;-1.08051851123043;-5.9906476390894];

for plot_num=(num_doc_1*2-1):2:(num_doc_2*2-1)
    
    %Data1=All_data(1:end, (eps_end-(plot_num+2)+1));
    Data1=All_data(1:end, (eps_end-(plot_num+1)+1));
    Data3=All_data(1:end, (eps_end-(plot_num)+1));
    


% ncoeffs = numcoeffs(ftype) % åñëè íàäî óçíàòü ÷èñëî êîýôôèöèåíòîâ â ô-öè
% coeffs = coeffnames(ftype) % åñëè íàäî óçíàòü èìåíà êîýôôèöèåíòîâ â ô-öè

    
[f, l ,k]=fit(Data1,Data3,ftype, 'StartPoint',...
    [coeff_table(1, xx-1), coeff_table(2, xx-1), coeff_table(3, xx-1),...
    coeff_table(4, xx-1), coeff_table(5, xx-1), coeff_table(6, xx-1),...
    coeff_table(7, xx-1), coeff_table(8, xx-1)],...
    'Lower', [0.10, 0.10, 0.7, 0.7,0.00000006, 0.00000006, -5 , -10],...
    'Upper',[1, 1, 1, 1,3, 3, coeff_table(7, xx-1), coeff_table(8, xx-1)], 'TolFun' , 0.000000000000000000000000000000001, 'MaxIter', 2000);

coeffvals1 = coeffvalues(f);
alfa1 =coeffvals1(1);
alfa2 =coeffvals1(2);

betta1 =coeffvals1(3);
betta2 =coeffvals1(4);

deps1 =coeffvals1(5);
deps2 =coeffvals1(6);

tay1 =coeffvals1(7);
tay2 =coeffvals1(8);

       
Fit_plot=zeros(length(Data1), 1);
for i=1:length(Data1)
Fit_plot(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1))+abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2));
end


for ut=1:numcoeffs(ftype)
coeff_table(ut, xx)=coeffvals1(ut);
end
%coeff_table(numcoeffs(ftype)+1, xx)=doc_num-num_doc_1-(xx-2)+1;
figure('Color','w');
    
set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','Arilal Black');  %ôîðìàò ïîäïèñåé îñåé
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arilal Black'); 
semilogx(Data1, Data3, 'ko', Data1, Fit_plot);
xlim([10^-1,10^6]);
xlabel('frequency');
ylabel('\epsilon''''', 'fontsize',20);
hold on;

for i=1:length(Data1)
Fit_plot1(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1));
end
for i=1:length(Data1)
Fit_plot2(i)=abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2));
end
plot(Data1, Fit_plot1,'--', Data1, Fit_plot2,'--');
hold off;
name_graf=sprintf('%s_graph%d',sample_name , xx-1 );
%print('-dpng','-r720',name_graf);
xx=xx+1;
end
coeff_table(:, 1)=[];
coeff_table2(1:8,1)=coeff_table(1:8, xx-2);
yy=2;
for plot_num=1:2:2*(num_doc_2-num_doc_1+1)
    eps_middle=eps_end-((num_doc_2*2)-1);
        
    Data1=All_data(1:end, (eps_middle+plot_num-1));
    Data3=All_data(1:end, (eps_middle+plot_num));
    


% ncoeffs = numcoeffs(ftype) % åñëè íàäî óçíàòü ÷èñëî êîýôôèöèåíòîâ â ô-öè
% coeffs = coeffnames(ftype) % åñëè íàäî óçíàòü èìåíà êîýôôèöèåíòîâ â ô-öè

    
[f, l ,k]=fit(Data1,Data3,ftype, 'StartPoint',...
    [coeff_table2(1, yy-1), coeff_table2(2, yy-1), coeff_table2(3, yy-1),...
    coeff_table2(4, yy-1), coeff_table2(5, yy-1), coeff_table2(6, yy-1),...
    coeff_table2(7, yy-1), coeff_table2(8, yy-1)],...
    'Lower', [0.10, 0.10, 0.7, 0.7,0.00000006, 0.00000006, coeff_table2(7, yy-1), coeff_table2(8, yy-1)],...
    'Upper',[1, 1, 1, 1,3, 3, 1.5, 0.8], 'TolFun' , 0.000000000000000000000000000000001, 'MaxIter', 2000);

coeffvals1 = coeffvalues(f);
alfa1 =coeffvals1(1);
alfa2 =coeffvals1(2);

betta1 =coeffvals1(3);
betta2 =coeffvals1(4);

deps1 =coeffvals1(5);
deps2 =coeffvals1(6);

tay1 =coeffvals1(7);
tay2 =coeffvals1(8);

       
Fit_plot=zeros(length(Data1), 1);
for i=1:length(Data1)
Fit_plot(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1))...
    +abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2));
end


for ut=1:numcoeffs(ftype)
coeff_table2(ut, yy)=coeffvals1(ut);
end
%coeff_table(numcoeffs(ftype)+1, xx)=doc_num-num_doc_1-(xx-2)+1;
figure('Color','w');
    
set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','Arilal Black');  %ôîðìàò ïîäïèñåé îñåé
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arilal Black'); 
semilogx(Data1, Data3, 'ko', Data1, Fit_plot);
xlim([10^-1,10^6]);
xlabel('frequency');
ylabel('\epsilon''''', 'fontsize',20);
hold on;

for i=1:length(Data1)
Fit_plot1(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1));
end

for i=1:length(Data1)
Fit_plot2(i)=abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2));
end
plot(Data1, Fit_plot1,'--', Data1, Fit_plot2,'--');
hold off;
name_graf=sprintf('%s_graph%d',sample_name , yy-1 );
%print('-dpng','-r720',name_graf);
yy=yy+1;

end
coeff_table2(:, 1)=[];
Temp_coeff_table=coeff_table2;
close all;
xx=1;
for i=1:size(coeff_table, 2)
    for j=1:numcoeffs(ftype)
average_coeff(j, i)=(coeff_table2(j, (num_doc_2-num_doc_1-i+2))+coeff_table(j, i))/2;
    end
end
for plot_num=(num_doc_1*2-1):2:(num_doc_2*2-1)
    
    %Data1=All_data(1:end, (eps_end-(plot_num+2)+1));
    Data1=All_data(1:end, (eps_end-(plot_num+1)+1));
    Data3=All_data(1:end, (eps_end-(plot_num)+1));
    


% ncoeffs = numcoeffs(ftype) % åñëè íàäî óçíàòü ÷èñëî êîýôôèöèåíòîâ â ô-öè
% coeffs = coeffnames(ftype) % åñëè íàäî óçíàòü èìåíà êîýôôèöèåíòîâ â ô-öè

    
[f, l ,k]=fit(Data1,Data3,ftype, 'StartPoint',...
    [average_coeff(1, xx), average_coeff(2, xx), average_coeff(3, xx),...
    average_coeff(4, xx), average_coeff(5, xx), average_coeff(6, xx),...
    average_coeff(7, xx), average_coeff(8, xx)],...
    'Lower', [0.10, 0.10, 0.7, 0.7,0.00000006, 0.00000006, -5, -25],...
    'Upper',[1, 1, 1, 1, 3, 3, 1.5, 0.5], 'TolFun' , 0.000000000000000000000000000000001, 'MaxIter', 2000);

coeffvals1 = coeffvalues(f);
alfa1 =coeffvals1(1);
alfa2 =coeffvals1(2);

betta1 =coeffvals1(3);
betta2 =coeffvals1(4);

deps1 =coeffvals1(5);
deps2 =coeffvals1(6);

tay1 =coeffvals1(7);
tay2 =coeffvals1(8);

       
Fit_plot=zeros(length(Data1), 1);
for i=1:length(Data1)
Fit_plot(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1))...
    +abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2));
end


for ut=1:numcoeffs(ftype)
coeff_table3(ut, xx)=coeffvals1(ut);
end
coeff_table(numcoeffs(ftype)+1, xx)=doc_num-num_doc_1-(xx-2)+1;
figure('Color','w');
    
set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','Arilal Black');  %ôîðìàò ïîäïèñåé îñåé
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arilal Black'); 
semilogx(Data1, Data3, 'ko', Data1, Fit_plot);
xlim([10^-1,10^6]);
xlabel('frequency');
ylabel('\epsilon''''', 'fontsize',20);
hold on;

for i=1:length(Data1)
Fit_plot1(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1));
end

for i=1:length(Data1)
Fit_plot2(i)=abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2));
end

plot(Data1, Fit_plot1,'--', Data1, Fit_plot2,'--');
hold off;
name_graf=sprintf('%s_graph%d',sample_name , xx-1 );
%print('-dpng','-r720',name_graf);
xx=xx+1;

end
end
% Incase of 2 peaks and conduction
if peak==4
ftype = fittype('abs(imag(deps1/(complex(1, ((1*10^tay1)*x*2*pi)^alfa1))^betta1))+abs(imag(deps2/(complex(1, ((1*10^tay2)*x*2*pi)^alfa2))^betta2))+imag(complex(1, (vsigma/(8.85e-12*x))))'); % âèä ôóíêöèè (ïîëèíîì âñåõ ôèòòèíãîâ)
xx=2;
coeff_table(1:9,1)=[0.454861730781002;0.392354999514640;0.820831659025781;0.700000226327451;0.00719361759145735;0.0121876524032858;-1;-9; 8e-15];
for plot_num=(num_doc_1*2-1):2:(num_doc_2*2-1)
    
    %Data1=All_data(1:end, (eps_end-(plot_num+2)+1));
    Data1=All_data(1:end, (eps_end-(plot_num+1)+1));
    Data3=All_data(1:end, (eps_end-(plot_num)+1));
    %Data_examples(:, xx-1)=Data3(:,1);


% ncoeffs = numcoeffs(ftype) % åñëè íàäî óçíàòü ÷èñëî êîýôôèöèåíòîâ â ô-öè
% coeffs = coeffnames(ftype) % åñëè íàäî óçíàòü èìåíà êîýôôèöèåíòîâ â ô-öè

    
[f, l ,k]=fit(Data1,Data3,ftype, 'StartPoint',...
    [coeff_table(1, xx-1), coeff_table(2, xx-1), coeff_table(3, xx-1),...
    coeff_table(4, xx-1), coeff_table(5, xx-1), coeff_table(6, xx-1),...
    coeff_table(7, xx-1), coeff_table(8, xx-1), coeff_table(9, xx-1)],...
    'Lower', [0.25, 0.05, 0.7, 0.7,0.000006, 0.000006, -6, -15, 1e-24],...
    'Upper',[1, 1, 1, 1,5, 5, coeff_table(7, xx-1), coeff_table(8, xx-1), 1e-5], 'TolFun' , 0.00000000000000000000000000000000000000000000001, 'MaxIter', 2000);

coeffvals1 = coeffvalues(f);
alfa1 =coeffvals1(1);
alfa2 =coeffvals1(2);
betta1 =coeffvals1(3);
betta2 =coeffvals1(4);
deps1 =coeffvals1(5);
deps2 =coeffvals1(6);
tay1 =coeffvals1(7);
tay2 =coeffvals1(8);
sigma=coeffvals1(9);

Fit_plot=zeros(length(Data1), 1);
for i=1:length(Data1)
Fit_plot(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1))...
    +abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2))...
    +imag(complex(1, (sigma/(8.85e-12*Data1(i)))));
end


for ut=1:numcoeffs(ftype)
coeff_table(ut, xx)=coeffvals1(ut);
end
%coeff_table(numcoeffs(ftype)+1, xx)=doc_num-num_doc_1-(xx-2)+1;
figure('Color','w');
    
set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','Arilal Black');  %ôîðìàò ïîäïèñåé îñåé
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arilal Black'); 
loglog(Data1, Data3, 'ko', Data1, Fit_plot);
xlim([10^-1,10^6]);
xlabel('frequency');
ylabel('\epsilon''''', 'fontsize',20);
hold on;

Fit_plot1=zeros(1, length(Data1));
Fit_plot2=zeros(1, length(Data1));
Fit_plot3=zeros(1, length(Data1));

for i=1:length(Data1)
Fit_plot1(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1));
end

for i=1:length(Data1)
Fit_plot2(i)=abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2));
end

for i=1:length(Data1)
Fit_plot3(i)=imag(complex(1, (sigma/(8.85e-12*Data1(i)))));
end

plot(Data1, Fit_plot1,'--', Data1, Fit_plot2,'--', Data1, Fit_plot3,'--');
hold off;
name_graf=sprintf('%s_graph%d',sample_name , xx-1 );
%print('-dpng','-r720',name_graf);
xx=xx+1;

end
coeff_table(:, 1)=[];
coeff_table2(1:9,1)=coeff_table(1:9, xx-2);
yy=2;
for plot_num=1:2:2*(num_doc_2-num_doc_1+1)
    eps_middle=eps_end-((num_doc_2*2)-1);
        
    Data1=All_data(1:end, (eps_middle+plot_num-1));
    Data3=All_data(1:end, (eps_middle+plot_num));
    


% ncoeffs = numcoeffs(ftype) % åñëè íàäî óçíàòü ÷èñëî êîýôôèöèåíòîâ â ô-öè
% coeffs = coeffnames(ftype) % åñëè íàäî óçíàòü èìåíà êîýôôèöèåíòîâ â ô-öè

    
[f, l ,k]=fit(Data1,Data3,ftype, 'StartPoint',...
    [coeff_table2(1, yy-1), coeff_table2(2, yy-1), coeff_table2(3, yy-1),...
    coeff_table2(4, yy-1), coeff_table2(5, yy-1), coeff_table2(6, yy-1),...
    coeff_table2(7, yy-1), coeff_table2(8, yy-1), coeff_table2(9, yy-1)],...
    'Lower', [0.10, 0.10, 0.7, 0.7,0.00000006, 0.00000006, coeff_table2(7, yy-1), coeff_table2(8, yy-1), 1*10^-8],...
    'Upper',[1, 1, 1, 1,3, 3, 1.5, 0.8, 1*10^-4], 'TolFun' , 0.000000000000000000000000000000001, 'MaxIter', 2000);

coeffvals1 = coeffvalues(f);
alfa1 =coeffvals1(1);
alfa2 =coeffvals1(2);

betta1 =coeffvals1(3);
betta2 =coeffvals1(4);

deps1 =coeffvals1(5);
deps2 =coeffvals1(6);

tay1 =coeffvals1(7);
tay2 =coeffvals1(8);

sigma =coeffvals1(9);

       
Fit_plot=zeros(length(Data1), 1);
for i=1:length(Data1)
Fit_plot(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1))...
    +abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2))...
    +imag(complex(1, (sigma/(8.85e-12*Data1(i)))));
end


for ut=1:numcoeffs(ftype)
coeff_table2(ut, yy)=coeffvals1(ut);
end
%coeff_table(numcoeffs(ftype)+1, xx)=doc_num-num_doc_1-(xx-2)+1;
figure('Color','w');
    
set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','Arilal Black');  %ôîðìàò ïîäïèñåé îñåé
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arilal Black'); 
semilogx(Data1, Data3, 'ko', Data1, Fit_plot);
xlim([10^-1,10^6]);
xlabel('frequency');
ylabel('\epsilon''''', 'fontsize',20);
hold on;


for i=1:length(Data1)
Fit_plot1(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1));
end

for i=1:length(Data1)
Fit_plot2(i)=abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2));
end

for i=1:length(Data1)
Fit_plot3(i)=imag(complex(1, (sigma/(8.85e-12*Data1(i)))));
end

plot(Data1, Fit_plot1,'--', Data1, Fit_plot2,'--', Data1, Fit_plot3,'--');
hold off;
name_graf=sprintf('%s_graph%d',sample_name , yy-1 );
%print('-dpng','-r720',name_graf);
yy=yy+1;

end
coeff_table2(:, 1)=[];
Temp_coeff_table=coeff_table2;
close all;
xx=1;
for i=1:size(coeff_table, 2)
    for j=1:numcoeffs(ftype)
average_coeff(j, i)=(coeff_table2(j, (num_doc_2-num_doc_1-i+2))+coeff_table(j, i))/2;
    end
end

for plot_num=(num_doc_1*2-1):2:(num_doc_2*2-1)
    
    %Data1=All_data(1:end, (eps_end-(plot_num+2)+1));
    Data1=All_data(1:end, (eps_end-(plot_num+1)+1));
    Data3=All_data(1:end, (eps_end-(plot_num)+1));
    


% ncoeffs = numcoeffs(ftype) % åñëè íàäî óçíàòü ÷èñëî êîýôôèöèåíòîâ â ô-öè
% coeffs = coeffnames(ftype) % åñëè íàäî óçíàòü èìåíà êîýôôèöèåíòîâ â ô-öè

    
[f, l ,k]=fit(Data1,Data3,ftype, 'StartPoint',...
    [average_coeff(1, xx), average_coeff(2, xx), average_coeff(3, xx),...
    average_coeff(4, xx), average_coeff(5, xx), average_coeff(6, xx),...
    average_coeff(7, xx), average_coeff(8, xx), average_coeff(9, xx)],...
    'Lower', [0.10, 0.10, 0.7, 0.7,0.00000006, 0.00000006, -5, -25, 1*10^-7],...
    'Upper',[1, 1, 1, 1, 3, 3, 1.5, 0.5, 1*10^-4], 'TolFun' , 0.000000000000000000000000000000001, 'MaxIter', 2000);

coeffvals1 = coeffvalues(f);
alfa1 =coeffvals1(1);
alfa2 =coeffvals1(2);

betta1 =coeffvals1(3);
betta2 =coeffvals1(4);

deps1 =coeffvals1(5);
deps2 =coeffvals1(6);

tay1 =coeffvals1(7);
tay2 =coeffvals1(8);
sigma=coeffvals1(9);

       
Fit_plot=zeros(length(Data1), 1);
for i=1:length(Data1)
Fit_plot(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1))...
    +abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2))...
    +imag(complex(1, (sigma/(8.85e-12*Data1(i)))));
end


for ut=1:numcoeffs(ftype)
coeff_table3(ut, xx)=coeffvals1(ut);
end
coeff_table(numcoeffs(ftype)+1, xx)=doc_num-num_doc_1-(xx-2)+1;
figure('Color','w');
    
set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','Arilal Black');  %ôîðìàò ïîäïèñåé îñåé
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arilal Black'); 
semilogx(Data1, Data3, 'ko', Data1, Fit_plot);
xlim([10^-1,10^6]);
xlabel('frequency');
ylabel('\epsilon''''', 'fontsize',20);
hold on;

for i=1:length(Data1)
Fit_plot1(i)=abs(imag(deps1/(complex(1, ((1*10^tay1)*Data1(i)*2*pi)^alfa1))^betta1));
end

for i=1:length(Data1)
Fit_plot2(i)=abs(imag(deps2/(complex(1, ((1*10^tay2)*Data1(i)*2*pi)^alfa2))^betta2));
end

for i=1:length(Data1)
Fit_plot3(i)=imag(complex(1, (sigma/(8.85e-12*Data1(i)))));
end

plot(Data1, Fit_plot1,'--', Data1, Fit_plot2,'--', Data1, Fit_plot3,'--');
hold off;
name_graf=sprintf('%s_graph%d',sample_name , xx-1 );
%print('-dpng','-r720',name_graf);
xx=xx+1;

end
end



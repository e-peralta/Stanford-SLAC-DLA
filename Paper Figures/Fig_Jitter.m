clear all
close all

%cols(0,:)=[1 1 1];  %use 'w'
cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple
%cols(5,:)=[0 0 0]; %use 'k'

folder='~/Dropbox/Research/Data/';

width=3.5;
height=2;
close all
hFig = figure(94);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')


filename='jitter1.csv';
data=importdata([folder,filename]);

factor=5.1878;  %um/keV

t1 = data(1:end-19,1); %
v1 = data(20:end,2)*factor; %
m1= mean(v1);
std1=std(v1)*ones(size(t1));
v1=v1-m1;

filename='jitter3.csv';
data=importdata([folder,filename]);

t2 = t1(end)+data(:,1); %
v2 = data(:,2)*factor*.8; %
m2= mean(v2);
std2=std(v2)*ones(size(t2));
v2=v2-m2;

filename='jitter4.csv';
data=importdata([folder,filename]);

t3 = [t2; t2(end)+data(:,1)]; %
v3 = data(:,2)*factor*.8; %
m3= mean(v3);
v3=v3-m3;
v3= [v2; v3];
std3=std(v3)*ones(size(t3));

std1(1)
std3(1)


h1a=plot(t1,v1,'o','Color',cols(4,:),'Linewidth',1,...
    'MarkerSize',4,'MarkerFaceColor',cols(2,:));
hold on
h1b=plot(t1,std1,':',t1,-std1,':','Color',cols(4,:),'Linewidth',1);
hold on

h2a=plot(t3,v3,'o','Color',cols(3,:),'Linewidth',1,...
    'MarkerSize',4,'MarkerFaceColor',cols(1,:));
hold on
h2b=plot(t3,std3,':',t3,-std3,':','Color',cols(3,:),'Linewidth',1);
hold on

axis tight
hx=xlabel('Shot number');
hy=ylabel('Main peak position [keV]');
ylim([-150 150])

tightfig;
hL=legend([h1a,h2a],{'Before circuit','After circuit'});
set(hL,'Location','northeast','Box','Off')

set(gca,'FontSize', 10);
set([hx,hy,hL],'FontSize', 10);
    

export_fig Jitter.eps -cmyk -r300 -painters%-painters
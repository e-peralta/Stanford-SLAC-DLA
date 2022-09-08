% figure
% plot(OTR(:,1),OTR(:,2))
% %%
% figure
% plot(timing2(:,1),timing2(:,2))
%%
%save timing IR OTR both

%%
filepath='~/Dropbox/Research/Data/';
filename='timing';
load([filepath,filename]);

col1=[230 97 1]/255;
col2=[253 184 99]/255;
col3=[178 171 210]/255;
col4=[94 60 153]/255;

cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple

% width=8;
% height=6;

width=3;
height=2;
close all

figure (1)
hFig = gcf;
set(hFig,'ActivePositionProperty','position')
%set(hFig,'Units','centimeters')
set(hFig,'Units','inches')

set(hFig, 'Position', [0 0 width height])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')


t0=both(:,1);
v0=both(:,2);

%v0=smooth(smooth(smooth(both(:,2))));

plot(t0,v0,'Color',cols(3,:),'LineWidth',1)
ax=gca;
hx=xlabel('Arrival time,{\it t} [ps]');
hy=ylabel('Amplitude [mV]');
axis tight

set([hx,hy], 'FontSize' , 10 );
set(ax, 'FontSize', 10  );

export_fig TemporalOverlap.eps -cmyk -r300 %-painters

%%


%close all

figure (2)
hFig = gcf;
width=3.5;
height=2;

set(hFig,'ActivePositionProperty','position')
%set(hFig,'Units','centimeters')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

t1=IR(:,1);
v1=smooth(smooth(smooth(IR(:,2))))*2-3;

[~, p1]=max(diff(v1));

t2=OTR(:,1);
v2=(smooth(smooth(smooth(OTR(:,2))))*2-9.9);
v2flat=(v2(end)-v2(1))/length(v2)*(1:length(v2));
v2=(v2-v2flat')*10;

[~, p2]=max(diff(v2));


% line([t1(p1+1) t1(p1+1)],[0 30],'LineWidth',2,'Color',col2,'LineStyle',':');
% hold on
line([t2(p2+1) t2(p2+1)],[2 17],'LineWidth',2,'Color',cols(1,:),'LineStyle',':');
hold on

h1=plot(t1,v1,'Color',cols(3,:),'LineWidth',1)
hold on
h2=plot(t2,v2,'Color',cols(2,:),'LineWidth',2)

hold on
h3=plot(t1(2:end),diff(v1),'--','Color',cols(1,:),'LineWidth',2)
hold on
h4=plot(t2(2:end),diff(v2),'--','Color',cols(4,:),'LineWidth',1)
axis tight
xlim([200 2500])
ax=gca;
hx=xlabel('Arrival time,{\it t} [ps]');
hy=ylabel('Amplitude [mV]');

set([hx,hy], 'FontSize' , 10 );
set(ax, 'FontSize', 10  ,'Box','On');

hL=legend([h1,h2,h3,h4],{'IR','OTR X 10','{\it d/dt} IR','{\it d/dt} OTR X 10'});
set(hL, 'FontSize', 10 ,'Box','Off' );

export_fig TemporalOverlap2.eps -cmyk -r300 %-painters

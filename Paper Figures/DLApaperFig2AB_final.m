
clear all
close all

%folder='V:\ARDB\E163\Data\130409\';
folder='~/Dropbox/Research/Data/';
filename='run2216_ScrnAvgs.mat';
eval(['load ',folder,filename]);

[b,a]=butter(8,.1,'low');

%Figure Size
xwidth=5.5;  %8.9 for normal, 5.5 for colorbar
ywidth=2.75;

%
% custom jet colormap
N=64;
col1=jet(N);
col1b=jet(N/2);
dblue3=.35:.03:col1(1,3);
col0=[zeros(length(dblue3),2),dblue3'];
%dblue3=.3:.05:col1(1,3);
dblue3=.3:.0625:col1(1,3);
col0b=[zeros(length(dblue3),2),dblue3'];
col2=jet(2*N);
col3=jet(4*N);
col4=jet(5*N);
col5=jet(6*N);
col6=jet(8*N);
dred=col6(end,1):-.01:.3;
col7=[dred',zeros(length(dred),2)];
col=[col0b;col1b(1:round(N/2*2/5),:);col2(2*round(N*2/5)+1:2*round(N*3/5),:);col1b(round(N/2*3/5)+1:N/2*3/4,:);col6(8*N*3/4+1:8*N,:);col7];

%col=[col0b;col1(1:round(N*3/4),:);col2(2*round(N*3/4)+1:2*round(N*4/5),:);col6(8*N*4/5+1:8*N,:);col7];

ROI=[400 1024 250 1024];
%ROI=[200 1024 250 1024];
%ROI=[700 1024 250 1024];
xD=ROI(1):ROI(2);
yD=ROI(3):ROI(4);
yDp=yD*17.65/1000;

contL1=[.045];
contL3=[.5];

scrnAvgOn=scrnAvgOn(yD(1):yD(end),xD(1):xD(end));
scrnAvgOff=scrnAvgOff(yD(1):yD(end),xD(1):xD(end));
xD=(xD-858*ones(size(xD)))*1.2;    %run2216


scrnAvgOff=imrotate(scrnAvgOff,-.5,'crop');

hFig = figure(97);

%set(gcf,'PaperPositionMode','auto')
%set(hFig,'ActivePositionProperty','position')
set(hFig,'ActivePositionProperty','outerposition')
set(hFig,'Units','centimeters')
%set(hFig, 'Position', [0 0 xwidth ywidth])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

hold on

%imagesc(xD,yD,scrnAvgOff)
spectraOff=mean(scrnAvgOff(350:475,:));

Y1=spectraOff;
c1=FitSpectrum5c(Y1,18,260);
[yfit1, ymain1, ysignal1, ysig1a, ysig1b]=DataFitPlotter(c1,18);
ysub1=ymain1(1:length(Y1))+ysig1a(1:length(Y1));
Y1b=Y1-ysub1;

iden=ones(size(scrnAvgOff,1),1);
imag=iden*ysub1;

%scrnAvgOff2=scrnAvgOff-imag;  % the part that subtracts the background
%scrnAvgOff=scrnAvgOff2;

scrnAvgOff=medfilt2(scrnAvgOff,[7,7]);
imagesc(xD,yDp,scrnAvgOff)
scrnAvgOff=medfilt2(scrnAvgOff,[25,25]);
%scrnAvgOff=filter(b,a,scrnAvgOff);

hold on
hC=contour(xD,yDp,scrnAvgOff,contL1,'LineColor','w','LineWidth',1);
% hold on


axis xy

%ylim([475 875]*17.65/1000)
%xlim([-120 120])  %For looking only at transmitted

colormap(col)

caxis([0 1])
colorbar('location','northoutside')

%set(gca, 'YTick', []);
%set(gca,'XGrid','on','GridLineStyle', '-.','XColor','w','LineWidth',2.5);
%set(gca,'XGrid','on','GridLineStyle', '-.','XColor','k','LineWidth',2.5);
%set(gcf, 'Color', 'k');
set(gca,'YColor','k');
set(gca,'XColor','k');
set(gca,'YTick',[9 12 15])

hy=ylabel('Position (mm)');
%hx=xlabel('\Delta E (keV)');
%enhance_plot(0,0,-1);
%title ('Averaged Spectrum Screen - Laser off')

 set( gca                       , ...
    'FontName'   , 'Arial' );
set(hy, ...
    'FontName'   , 'Arial');
set(gca             , ...
    'FontSize'   , 6           );
set(hy  , ...
    'FontSize'   , 7          );

set(gca, ...
  'TickDir'     , 'out'     , ...
  'Box'         , 'on'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );

export_fig ScrOff4C.tif -cmyk -r300 %-painters%-painters

%
scrnAvgOn=imrotate(scrnAvgOn,-.5,'crop');

hFig = figure(98);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','centimeters')
%set(hFig, 'Position', [0 0 xwidth ywidth])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

hold on
%figure (98)
%imagesc(xD,yD,scrnAvgOn)
spectraOn=mean(scrnAvgOn(350:475,:));


Y2=spectraOn;
c2=FitSpectrum5c(Y2,18,260);
[yfit2, ymain2, ysignal2, ysig2a, ysig2b]=DataFitPlotter(c2,18);
ysub2=ymain2(1:length(Y2))+ysig1a(1:length(Y2));

Y2b=Y2-ysub2;

iden=ones(size(scrnAvgOn,1),1);
imag=iden*ysub2;

%scrnAvgOn2=scrnAvgOn-imag;      % the part that subtracts the background
%scrnAvgOn=scrnAvgOn2;

scrnAvgOn=medfilt2(scrnAvgOn,[7,7]);
imagesc(xD,yDp,scrnAvgOn)

scrnAvgOn=medfilt2(scrnAvgOn,[25,25]);
%scrnAvgOn=filter(b,a,scrnAvgOn);



hold on
hC=contour(xD,yDp,scrnAvgOn,contL1,'LineColor','w','LineWidth',1);
hold on

%hl=clabel(hC);
%set(hl,'Color','w','LabelSpacing',250)
axis xy

%ylim([yD(50) yD(end)])
%ylim([375 975])

%ylim([475 875]*17.65/1000)
%xlim([-120 120])  %For looking only at transmitted

%axis equal
%colormap(hot)
%     colormap(jet)
%     caxis([.05 .35])
colormap(col)
caxis([0 1])
%colorbar
set(gca,'YColor','k');
set(gca,'XColor','k');
set(gca,'YTick',[9 12 15])

hy=ylabel('Position (mm)');

 set( gca                       , ...
    'FontName'   , 'Arial' );
set(hy, ...
    'FontName'   , 'Arial');
set(gca             , ...
    'FontSize'   , 6           );
set(hy  , ...
    'FontSize'   , 7          );

set(gca, ...
  'TickDir'     , 'out'     , ...
  'Box'         , 'on'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );

export_fig ScrnOn4.tif -cmyk -r300 %-painters%-painters

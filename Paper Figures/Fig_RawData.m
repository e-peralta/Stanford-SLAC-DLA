
clear all
close all

%folder='V:\ARDB\E163\Data\130409\';
folder='~/Dropbox/Research/Data/';
filename='run2216_ScrnAvgs.mat';
eval(['load ',folder,filename]);

[b,a]=butter(8,.1,'low');

%Figure Size
xwidth=4.5;  %8.9 for normal, 5.5 for colorbar [cm]
ywidth=1.5;

% Custom Morgenstemning colormap for grating transmission
N=48;%32;
col1=colormap(morgenstemning(N,'invert',1));
col2=colormap(morgenstemning(2*N,'invert',1));
col4=colormap(morgenstemning(4*N,'invert',1));
col8=colormap(morgenstemning(8*N,'invert',1));
col16=colormap(morgenstemning(16*N,'invert',1));
%col=[col1(1:round(N*2/5),:);col4(round(4*N*2/5)+1:round(4*N*3/5),:);col1(round(N*3/5)+1:N*3/4,:);col16(16*N*3/4+1:16*N,:)]; %SAME
%col=[col1(1:round(N*2/5),:);col2(round(2*N*2/5)+1:round(2*N*3/5),:);col1(round(N*3/5)+1:N*3/4,:);col16(16*N*3/4+1:16*N,:)];
col=[col2(1:round(2*N*3/5),:);col1(round(N*3/5)+1:N*3/4,:);col16(16*N*3/4+2:16*N,:)];
colsize=size(col)


%%

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


ywidth=2.5; %for colorbar figure

%set(gcf,'PaperPositionMode','auto')
%set(hFig,'ActivePositionProperty','position')
set(hFig,'ActivePositionProperty','outerposition')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 xwidth ywidth])

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

% the part that subtracts the background (comment for raw data)
scrnAvgOff2=scrnAvgOff-imag;  
scrnAvgOff=scrnAvgOff2;

scrnAvgOff=medfilt2(scrnAvgOff,[7,7]);
imagesc(xD,yDp,scrnAvgOff)
scrnAvgOff=medfilt2(scrnAvgOff,[25,25]);
%scrnAvgOff=filter(b,a,scrnAvgOff);

% ADD CONTOUR (remove for raw data)
hold on
hC=contour(xD,yDp,scrnAvgOff,contL1,'LineColor','k','LineWidth',1);


axis xy
%  axis tight (use for raw data)
ylim([475 875]*17.65/1000)
xlim([-120 120])  %For looking only at transmitted

colormap(col)

caxis([0.02 1])
caxis([0.02 .5])
colorbar
%colorbar('location','northoutside')

%set(gca, 'YTick', []);
%set(gca,'XGrid','on','GridLineStyle', '-.','XColor','w','LineWidth',2.5);
%set(gca,'XGrid','on','GridLineStyle', '-.','XColor','k','LineWidth',2.5);
%set(gcf, 'Color', 'k');
set(gca,'YColor','k');
set(gca,'XColor','k');
set(gca,'YTick',[9 12 15])
%set(gca,'YTick',6:3:18)

hy=ylabel('Position [mm]');
%hx=xlabel('\Delta E (keV)');
%enhance_plot(0,0,-1);
%title ('Averaged Spectrum Screen - Laser off')

ax=gca;

set(gca,'Color','none','Box','on','FontSize',10);
set(hy, 'FontSize', 10);

%  set( gca                       , ...
%     'FontName'   , 'Arial' );
% set(hy, ...
%     'FontName'   , 'Arial');
% set(gca             , ...
%     'FontSize'   , 6           );
% set(hy  , ...
%     'FontSize'   , 7          );
% 
% set(gca, ...
%   'TickDir'     , 'out'     , ...
%   'Box'         , 'on'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'LineWidth'   , 1         );

export_fig ScrOffC.eps -cmyk -r300 -painters%-painters

%% Laser On image
scrnAvgOn=imrotate(scrnAvgOn,-.5,'crop');

hFig = figure(98);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 xwidth ywidth])

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

scrnAvgOn2=scrnAvgOn-imag;      % the part that subtracts the background
scrnAvgOn=scrnAvgOn2;

scrnAvgOn=medfilt2(scrnAvgOn,[7,7]);
imagesc(xD,yDp,scrnAvgOn)

scrnAvgOn=medfilt2(scrnAvgOn,[25,25]);
%scrnAvgOn=filter(b,a,scrnAvgOn);


% ADD CONTOUR
hold on
hC=contour(xD,yDp,scrnAvgOn,contL1,'LineColor','k','LineWidth',1);


%hl=clabel(hC);
%set(hl,'Color','w','LabelSpacing',250)
axis xy
axis tight

%ylim([yD(50) yD(end)])
%ylim([375 975])

ylim([475 875]*17.65/1000)
xlim([-120 120])  %For looking only at transmitted

%axis equal
%colormap(hot)
%     colormap(jet)
%     caxis([.05 .35])
colormap(col)
%caxis([0.02 1])
caxis([0.02 .5])

%colorbar
set(gca,'YColor','k');
set(gca,'XColor','k');
set(gca,'YTick',[9 12 15])
%set(gca,'YTick',6:3:18)
hy=ylabel('Position [mm]');
ax=gca;

set(gca,'Color','none','Box','on','FontSize',10);
set(hy, 'FontSize', 10);

% set(gca, ...
%   'TickDir'     , 'out'     , ...
%   'Box'         , 'on'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'LineWidth'   , 1         );

export_fig ScrnOn.eps -cmyk -r300 -painters%-painters

 %% Spectrum Projections
%cols(0,:)=[1 1 1];  %use 'w'
cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple
%cols(5,:)=[0 0 0]; %use 'k'

% 
% hFig = figure(99);
% 
% %set(gcf,'PaperPositionMode','auto')
% set(hFig,'ActivePositionProperty','position')
% set(hFig,'Units','inches')
% set(hFig, 'Position', [0 0 xwidth ywidth])
% 
% set(gcf, 'Color', 'w');
% set(hFig,'Units','points')
% 
% hold on
% 
% h1=plot(xD,Y1,'LineWidth',1,'Color',cols(1,:));
% hold on
% h2=plot(xD,Y2,'LineWidth',1,'Color',cols(4,:));
% axis xy
% axis tight
% 
% %ylim([yD(50) yD(end)])
% %ylim([375 975])
% 
% %ylim([475 875]*17.65/1000)
% %xlim([-120 120])  %For looking only at transmitted
% 
% %axis equal
% %colormap(hot)
% %     colormap(jet),'LineWidth',1,'Color','k');
% %     caxis([.05 .35])
% colormap(col)
% caxis
% caxis([0.02 1])
% %colorbar
% % set(gca,'YColor','k');
% % set(gca,'XColor','k');
% %set(gca,'YTick',[9 12 15])
% %set(gca,'YTick',6:3:18)
% %hy=ylabel('Position [mm]');
% hy=ylabel('Charge density [arb. unit]');
% ax=gca;
% 
% set(gca,'Color','none','FontSize',10);
% %set(hy, 'FontSize', 10);
% hL=legend([h1,h2],{'Laser off','Laser on'});
% set(hL,'Location','NorthEast','Box','off')
% set([hy,hL], 'FontSize', 10);
% 
% 
% % set(gca, ...
% %   'TickDir'     , 'out'     , ...
% %   'Box'         , 'on'     , ...
% %   'TickLength'  , [.02 .02] , ...
% %   'LineWidth'   , 1         );
% 
% export_fig RawData.eps -cmyk -r300 -painters%-painters
% 

%%

%%
%clear all
close all
%load raw_model_sim_data
%load V:\ARDB\E163\Data\2013\130513\raw_model_sim_data
load ~/Dropbox/Research/Data/raw_model_sim_data 

%
cout=FitSpectrum5c(raw_on,18,340,50);
[yfit, ymain, ysignal, ysig1, ysig2]=DataFitPlotter(cout,18);

% plot(1:1024,raw_off,'x',1:1024,raw_on,'x',...
%    1:1024,fit_full_off,1:1024,fit_main_off,...
%    1:1024,fit_main_on,1:1024,fit_plateau_off,...
%    1:1024,fit_signal_off,1:1024,yfit,1:1024,ysignal)
% xlim([600 1000])

%
%figure(2)


xwidth=4.5;  %8.9 for normal, 5.5 for colorbar [cm]
ywidth=2.3;

hFig = figure(95);


%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 xwidth ywidth])

 set(gcf, 'Color', 'w');
set(hFig,'Units','points')

hold on


%h1=plot(x0,raw_off-fit_main_off-fit_plateau_off,'bx','Color',[.65 .65 1],'LineWidth',1)
h1a=plot(x0,raw_off-fit_main_off-fit_plateau_off,'LineWidth',1,'Color',cols(1,:),...
    'Marker','x','MarkerSize',8,'LineStyle','none');

% 
% [nh, xh] = hist(EnergyIn,[-400:5:400]*1e3);
% %%plot(xh*1e-3,nh/max(nh)*.2175,'bs','LineWidth',3)
% 
%h3=plot(x0,raw_on-fit_main_on-fit_plateau_off,'rx','Color',[1 .65 .65],'LineWidth',1)
h2a=plot(x0,raw_on-fit_main_on-fit_plateau_off,'LineWidth',1,'Color',cols(2,:),...
    'Marker','x','MarkerSize',8,'LineStyle','none');

h1b=plot(x0,fit_signal_off,'LineWidth',1.5,'Color',cols(3,:))
% 
h2b=plot(x0,model_signal_on,'LineWidth',1.5,'Color',cols(4,:))
%     
[nh, xh] = hist(EnergyOut,[-400:8:400]*1e3);
h3=plot(xh*1e-3,nh/max(nh)*.1,'k','LineWidth',1,'Marker','.','MarkerSize',10,'LineStyle','none');
% 
% set(h5                         , ...
%   'Marker'          , 'o'         , ...
%   'MarkerSize'      , 2           , ...
%   'MarkerEdgeColor' , 'k'      , ...
%   'MarkerFaceColor' , 'k' );


%h6=plot(x0(820:end),ysignal(820:end)-fit_main_on(820:end)-.8*fit_plateau_off(820:end),'LineWidth',1,'Color',[.75 0 0])
%h6=plot(x0,ysignal-fit_main_on-.8*fit_plateau_off,'LineWidth',1,'Color',[.75 0 0])
%h6=plot(x0,ysignal-fit_main_on-fit_plateau_off,'LineWidth',1,'Color',[.75 0 0])


hx=xlabel('Energy deviation, {\it\DeltaE} [keV]')
hy=ylabel('Charge density [arb. unit]')
%title('Adding Simulation Data')

axis([-120 120 -.005 .225])

box on
%enhance_plot;


%xwidth=600;
%ywidth=300;

%hL=legend('Laser Off','Spec. Fit','Laser On','Model','Simulation','Location','Northeast')
%hL=legend('Laser Off','Location','Northeast')
%hL=legend('Laser Off','Spec. Fit Off','Laser On','Spec. Fit On','Location','Northwest')


hL=legend([h1a,h1b,h2a,h2b,h3],{'Laser off','Spectrum fit','Laser on','Model','Simulation'});
set(hL,'Location','NorthEast','Box','off')

set([hx,hy,hL], 'FontSize', 10);
set(gca,'Color','none','Box','on','FontSize',10);

export_fig Modulation.eps -painters -cmyk -r300


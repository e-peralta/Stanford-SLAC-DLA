%%
clear all
close all
[b,a]=butter(8,.1,'low'); %%%change the second parameter to change the degree of LP-FIltering. ([0,1], 0 for filtering everything, 1 for no filtering at all)
%c=[14.69 .1 .7324 .2796];
%fShift0=c(1)*exp(-(lpStrength-c(2))^c(3)/c(4));

% MAKE SURE COLORMAP DOESNT EXCEED 256 COLORS OR FILE SIZE BLOWS UP 
%% Improved Morgenstemning colormap
N=102;
col1=colormap(morgenstemning(N,'invert',1));
col3=colormap(morgenstemning(3*N,'invert',1));
col0=[col1(1:round(N/8),:);col3(round(3*N/8)+1:round(3*N*7/8),:);col1(round(7*N/8)+1:end,:)];
col0size=size(col0)
%% Custom Morgenstemning colormap for grating transmission
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
%col=[col1(1:round(N*2/5),:);col2(round(4*N*2/5)+1:round(4*N*3/5),:);col1(r
%ound(N*3/5)+1:N*3/4,:);col3(8*N*3/4+1:8*N,:)];

%cols(0,:)=[1 1 1];  %use 'w'
cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple
%cols(5,:)=[0 0 0]; %use 'k'


%ROI=[200 1024 250 1024];
ROI=[1 1024 1 1024];
xD=ROI(1):ROI(2);
yD=ROI(3):ROI(4);
yDp=yD*17.65/1000;
xD=(xD-858*ones(size(xD)))*1.2;    %run2216

%% Load Vacuum spectra
%folder='V:\ARDB\eperalta\';
%folder='~/Documents/MATLAB/';
folder='~/Dropbox/Research/Data/';
%filename='vacuum_grating_pictures.mat';
filename='spectrum2D_121220.mat';
eval(['load ',folder,filename]);

% Shift & straighten image
vac=circshift(vacuum, [0 268]);
vac=imrotate(vac,-1.5,'crop'); %-1.35

%screen=imag{i}(y(1):y(end),x(1):x(end));
screen=vac;
spectra0=mean(screen);

%re-scale the image for a nice plot

spectraF=filter(b,a,spectra0);
% This initial fit is used for various normalizations below
%[cout, ~, ~,~, ybkgd, ~, ~, ~]=FitSpectrum5c(spectraF(400:end),1);
[cout, ~, ~,~, ybkgd, ~, ~, ~]=FitSpectrum(1,spectraF(400:end));

peakAmp0=cout(1);

screenD=screen-ybkgd(1)*ones(size(screen));
screenD(screenD<0)=0;
specVac=mean(screenD);

%remove x-rays
screenD(screenD>2*peakAmp0)=0;

scrnVac=screenD;%(yD(1):yD(end),xD(1):xD(end));
Navg=5;
scrnVac=medfilt2(scrnVac,[Navg,Navg]);%/1.616;

%% Load Grating Spectrum
%folder='V:\ARDB\E163\Data\130409\';
%folder='~/Documents/MATLAB/';
folder='~/Dropbox/Research/Data/';
filename='run2216_ScrnAvgs.mat';
eval(['load ',folder,filename]);

% ROI=[200 1024 250 1024];
% xD=ROI(1):ROI(2);
% yD=ROI(3):ROI(4);
% yDp=yD*17.65/1000;


%scrnAvgOff=scrnAvgOff(yD(1):yD(end),xD(1):xD(end));
%xD=(xD-858*ones(size(xD)))*1.2;    %run2216

scrnAvgOff=imrotate(scrnAvgOff,-.5,'crop');

spectraOff=mean(scrnAvgOff(350:475,:));
Y1=spectraOff;

screen=scrnAvgOff*1.185;    %to match the area of the projected spectrum

spectra0=mean(screen);

%re-scale the image for a nice plot
spectraF=filter(b,a,spectra0);
% This initial fit is used for various normalizations below
%[cout, ~, ~,~, ybkgd, ~, ~, ~]=FitSpectrum5c(spectraF(400:end),1);
[cout, ~, ~,~, ybkgd, ~, ~, ~]=FitSpectrum(1,spectraF(400:end));

peakAmp0=cout(1);

screenD=screen-ybkgd(1)*ones(size(screen));
screenD(screenD<0)=0;
specGla=mean(screenD);
%remove x-rays
screenD(screenD>2*peakAmp0)=0;

scrnGra=screenD;%(yD(1):yD(end),xD(1):xD(end));
Navg=5;
scrnGra=medfilt2(scrnGra,[Navg,Navg]);%/1.616;


imagesc(xD,yDp,scrnGra)
%%
width=4;
height=4;
close all
hFig = figure(94);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

%%
% subplot(3,1,1)
% 
% imagesc(xD,yDp,scrnVac)
% 
% set(gca,'YTick',[1 5 9 13 17])
% axis xy
% colormap(col)
% %caxis
% %caxis([0 1450])
% %caxis  %did this to figure out how to scale data so it would be [0 1]
% xlim([-750 100])
% hy1=ylabel('Position [mm]');
% set(gca,'XTickLabel',[],'XMinorTick','on');
% ax1=gca;

%%
%subplot(3,1,2)
subplot(2,1,1)

imagesc(xD,yDp,scrnGra)

set(gca,'YTick',[1 5 9 13 17])
axis xy
colormap(col)
%caxis([5 250])
%%caxis([0 .95])
%caxis  %did this to figure out how to scale data so it would be [0 1]
xlim([-750 200])
hy1=ylabel('Position [mm]');
set(gca,'XTickLabel',[],'XMinorTick','on');
ax1=gca;

hold on

line([xD(1) xD(end)], [575*17.65/1000 575*17.65/1000],...
    'LineStyle',':','LineWidth',1,'Color','k');
line([xD(1) xD(end)], [775*17.65/1000 775*17.65/1000],...
    'LineStyle',':','LineWidth',1,'Color','k');
%% Load Joshs Particle Sim Results

filename='AllTheTeeth_Corrected.tsv';
eval(['load ',folder,filename]);

Energy=AllTheTeeth_Corrected(:,1);
y1c=AllTheTeeth_Corrected(:,2);

E=(Energy-824*ones(size(Energy)))*1.2;


%% Load individual data set (not averaged like screen shot)

filename='run2216_f18-4_1.mat';
eval(['load ',folder,filename]);
list=[2, 12, 13, 20, 21, 28, 29, 44, 45, 52, 53, 60, 61, 76, 77, 84, 85, 92, 93, 100, 101, 108, 109, 116, 117, 124, 125, 140, 141, 148, 149, 156, 157, 172, 173, 180, 181, 188, 189, 204, 205, 212, 213, 220, 221, 228, 229, 236, 237, 244, 245, 268, 269, 276, 277, 284, 285, 300, 301];
%list2=[2,3];

ROI=Data.ROI;
Peak1pos=Data.Peak1pos;
Peak1amp=Data.Peak1amp;
spectra=Data.spectra;
regen=Data.regen;
Cout1=Data.Cout1;

subplot(2,1,2)

x=(ROI(1):ROI(2));
%x=200:1024;

for i=1:1%length(list)
   ind=list(i);
   xVar=x-Peak1pos(ind)*ones(1,length(x));
   xVar=xVar*1.2-323*ones(1,length(x));   %run2216
   [yfit, ymain, ysignal]=DataFitPlotter(Cout1(ind,:),18);
   if regen(ind)==0

      %hOff=plot(xVar,spectra(ind,x)./Peak1amp(ind),'b+','LineWidth',.8,'MarkerSize',3);
      %hOff=plot(xVar,spectra(ind,x)./Peak1amp(ind),'c+','LineWidth',.8,'MarkerSize',3);
      %h2=plot(xVar,spectra(ind,x)./Peak1amp(ind),'kx','LineWidth',.8,'MarkerSize',3.5);
      h2=plot(xVar,spectra(ind,x)./Peak1amp(ind),'LineWidth',3,'Color',cols(1,:));
      
      hold on
      %hOff_fit=plot(xF,yfit,'b','LineWidth',2);
   end
end

h1=plot(E,y1c./max(y1c)*1.15,'.','Color',cols(2,:),'MarkerSize',5)

hold on
c=[14.69 .1 .7324 .2796];
fShift0=c(1)*exp(-(0-c(2))^c(3)/c(4));
xVar=xVar-fShift0*ones(1,length(x));

%plot(xVar,yfit,'Color',[1 .7 0],'LineWidth',.75);%,'LineStyle','--')
%plot(xVar,yfit,'Color','r','LineWidth',.75);%,'LineStyle','--')
%h3=plot(xVar,yfit,'Color',[1 .7 0],'LineWidth',.75);%,'LineStyle','--')
h3=plot(xVar,yfit,'Color',cols(4,:),'LineWidth',1);%,'LineStyle','--')

hx=xlabel('Energy deviation, {\it\DeltaE} [keV]');
hy2=ylabel('Charge density [a.u.]');
xlim([-750 200])
ylim([0 1])
%hL=legend('Simulation','Data','Spec. Fit')
hL=legend([h2,h1,h3],{'Data','Simulation','Spec. Fit'});

ax2=gca;
%
set(gca,'XMinorTick','on')
%tightfig;
set(hL,'Location','NorthWest','Box','Off')
set([ax1,ax2],'FontSize', 10);
set([hx,hy1,hy2,hL],'FontSize', 10);
%%
tightfig;

% 
%     %line([c1(2)-c1(3)*scale1 c1(2)+c1(4)*scale2],...
%     %    [c1(5)+c1(1)/2 c1(5)+c1(1)/2], 'LineWidth',2,'Color','r');
%     fprintf('\nmode %g: FWHM= %g += %g, Chi^2=%g\n',i,FWHM,FWHMstd,ChiSq)
% end
% axis tight
% xlim([-750 100])
% hx3=xlabel('Energy deviation, {\it\DeltaE} [keV]');
% hy3=ylabel('Charge density [a.u.]');
% ax3=gca;
% %
% set(gca,'XMinorTick','on')
% %tightfig;
% hL=legend([h1,h1f,h2,h2f],{'Vacuum','Gauss fit','Fused Silica','Lorentz fit'});
% set(hL,'Location','NorthWest','Box','Off')
% set([ax1,ax2,ax3],'FontSize', 10);
% set([hx3,hy1,hy2,hy3,hL],'FontSize', 10);
% 
% % set(gca, ...
% %     'TickDir'     , 'in'     , ...
% %     'Box'         , 'on'     , ...
% %     'TickLength'  , [.02 .02] , ...
% %     'LineWidth'   , 1         );
% 
% %export_fig ScrVac.tif -cmyk -r300 %-painters%-painters
% %export_fig ScrGrat.tif -cmyk -r300 %-painters%-painters
% 
% %export_fig ScrVac2.tif -rgb -r300 %-painters%-painters
% %export_fig ScrGrat2.tif -rgb -r300 %-painters%-painters
% %
% %export_fig ScrTrans2.eps -cmyk -r300 -painters%-painters
% %%
%%
export_fig ScrTrans.eps -cmyk -r300 -painters%-painters

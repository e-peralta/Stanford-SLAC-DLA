%%
[b,a]=butter(8,.1,'low'); %%%change the second parameter to change the degree of LP-FIltering. ([0,1], 0 for filtering everything, 1 for no filtering at all)
c=[14.69 .1 .7324 .2796];
%fShift0=c(1)*exp(-(lpStrength-c(2))^c(3)/c(4));


ROI=[200 1024 250 1024];
xD=ROI(1):ROI(2);
yD=ROI(3):ROI(4);
yDp=yD*17.65/1000;

% Load Vacuum spectra
%load V:\ARDB\eperalta\vacuum_grating_pictures.mat

%load ~/Documents/MATLAB/vacuum_grating_pictures.mat
folder='~/Documents/MATLAB/';
% Straighthen image
filename='vacuum_grating_pictures.mat';
eval(['load ',folder,filename]);

figure (1)
subplot(2,2,1)
imagesc(vacuum)
subplot(2,2,2)
vac=imrotate(vacuum,-1.5,'crop'); %-1.35
imagesc(vac)
subplot(2,2,3)
vac=circshift(vac, [0 268]);
imagesc(vac)


%screen=imag{i}(y(1):y(end),x(1):x(end));
screen=vac;
spectra0=mean(screen);

%re-scale the image for a nice plot

spectraF=filter(b,a,spectra0);
% This initial fit is used for various normalizations below
[cout, ~, ~, ybkgd, ~, ~, ~,~]=FitSpectrum5c(spectraF,3,0,0,0);
peakAmp0=cout(1);
peakPos0=cout(2);

screenD=screen-ybkgd(1)*ones(size(screen));

screenD(screenD>2*peakAmp0)=0;   %remove x-rays
screenD(screenD<0)=0;
screenD=screenD./peakAmp0;
subplot(2,2,4)
imagesc(screenD)

%%

scrnVac=screenD(yD(1):yD(end),xD(1):xD(end));

% Load Grating Spectrum
%folder='V:\ARDB\E163\Data\130409\';
filename='run2216_ScrnAvgs.mat';
eval(['load ',folder,filename]);

scrnAvgOff=scrnAvgOff(yD(1):yD(end),xD(1):xD(end));
xD=(xD-858*ones(size(xD)))*1.2;    %run2216

scrnAvgOff=imrotate(scrnAvgOff,-.5,'crop');


spectraOff=mean(scrnAvgOff(350:475,:));

Y1=spectraOff;
%%
figure
c1=FitSpectrum5c(Y1,18,260,0,1);
%%

% % custom jet colormap
% N=64;
% col1=jet(N);
% col1b=jet(N/2);
% dblue3=.35:.03:col1(1,3);
% col0=[zeros(length(dblue3),2),dblue3'];
% %dblue3=.3:.05:col1(1,3);
% dblue3=.3:.0625:col1(1,3);
% col0b=[zeros(length(dblue3),2),dblue3'];
% col2=jet(2*N);
% col3=jet(4*N);
% col4=jet(5*N);
% col5=jet(6*N);
% col6=jet(8*N);
% dred=col6(end,1):-.01:.3;
% col7=[dred',zeros(length(dred),2)];
% 
% col=[col0b;col1b(1:round(N/2*2/5),:);col2(2*round(N*2/5)+1:2*round(N*3/5),:);col1b(round(N/2*3/5)+1:N/2*3/4,:);col6(8*N*3/4+1:8*N,:);col7];

N=32;
col1=colormap(morgenstemning(N,'invert',1));
col2=colormap(morgenstemning(2*N,'invert',1));
col4=colormap(morgenstemning(4*N,'invert',1));
col8=colormap(morgenstemning(8*N,'invert',1));
col16=colormap(morgenstemning(16*N,'invert',1));
%col=[col1(1:round(N*2/5),:);col4(round(4*N*2/5)+1:round(4*N*3/5),:);col1(round(N*3/5)+1:N*3/4,:);col16(16*N*3/4+1:16*N,:)]; %SAME
%col=[col1(1:round(N*2/5),:);col2(round(2*N*2/5)+1:round(2*N*3/5),:);col1(round(N*3/5)+1:N*3/4,:);col16(16*N*3/4+1:16*N,:)];
col=[col2(1:round(2*N*3/5),:);col1(round(N*3/5)+1:N*3/4,:);col16(16*N*3/4+1:16*N,:)];

%col=[col1(1:round(N*2/5),:);col2(round(4*N*2/5)+1:round(4*N*3/5),:);col1(round(N*3/5)+1:N*3/4,:);col3(8*N*3/4+1:8*N,:)];



%% Figure
%xwidth=8.9;  %8.9 for normal, 5.5 for colorbar
%ywidth=3.75;  %2.75 used on Fig 2
width=4;
height=2;

Navg=5;

close all
hFig = figure(94);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

hold on

 % Vacuum Spectrum
scrnVacp=medfilt2(scrnVac,[Navg,Navg])/1.616;
imagesc(xD,yDp,scrnVacp)

% % Grating Spectrum
% scrnAvgOffp=medfilt2(scrnAvgOff,[Navg,Navg])/1.051;
% imagesc(xD,yDp,scrnAvgOffp)

set(gca,'YTick',[5 9 13 17])
axis xy
%colorbar
colormap(col)

%colormap(morgenstemning(256,'invert',1))
%caxis([.03 1])
%caxis  %did this to figure out how to scale data so it would be [0 1]
xlim([-750 200])

hy=ylabel('Position [mm]');
%hx=xlabel('\Delta E (keV)');

% set( gca                       , ...
%     'FontName'   , 'Arial' );
% set(hy, ...
%     'FontName'   , 'Arial');
set(gca             , ...
    'FontSize'   , 10           );
set(hy  , ...
    'FontSize'   , 10          );

% set(gca, ...
%     'TickDir'     , 'in'     , ...
%     'Box'         , 'on'     , ...
%     'TickLength'  , [.02 .02] , ...
%     'LineWidth'   , 1         );

%export_fig ScrVac.tif -cmyk -r300 %-painters%-painters
%export_fig ScrGrat.tif -cmyk -r300 %-painters%-painters

%export_fig ScrVac2.tif -rgb -r300 %-painters%-painters
%export_fig ScrGrat2.tif -rgb -r300 %-painters%-painters
%%
%export_fig ScrVac.eps -cmyk -r300 -painters%-painters
%export_fig ScrGrat.eps -cmyk -r300 -painters%-painters

%%
clear all
close all
[b,a]=butter(8,.1,'low'); %%%change the second parameter to change the degree of LP-FIltering. ([0,1], 0 for filtering everything, 1 for no filtering at all)
%c=[14.69 .1 .7324 .2796];
%fShift0=c(1)*exp(-(lpStrength-c(2))^c(3)/c(4));


%ROI=[200 1024 250 1024];
ROI=[200 1024 250 1024];
xD=ROI(1):ROI(2);
yD=ROI(3):ROI(4);
yDp=yD*17.65/1000;

% Load Vacuum spectra
%load V:\ARDB\eperalta\vacuum_grating_pictures.mat

% folder='~/Documents/MATLAB/';
% filename='vacuum_grating_pictures.mat';
% eval(['load ',folder,filename]);

%folder='~/Documents/MATLAB/';
filename='vacuum_grating_pictures.mat';
eval(['load ',filename]);


% Straighthen image
 figure (1)
% subplot(2,2,1)
% imagesc(vacuum)
%subplot(2,2,2)
vac=imrotate(vacuum,-1.5,'crop'); %-1.35
%imagesc(vac)
%subplot(2,2,3)
vac=circshift(vac, [0 268]);
imagesc(vac)

%screen=imag{i}(y(1):y(end),x(1):x(end));
screen=vac;
spectra0=mean(screen);

%re-scale the image for a nice plot

spectraF=filter(b,a,spectra0);
% This initial fit is used for various normalizations below
[cout, ~, ~,~, ybkgd, ~, ~, ~]=FitSpectrum5c(spectraF(400:end),1);
peakAmp0=cout(1);

screenD=screen-ybkgd(1)*ones(size(screen));

%remove x-rays
screenD(screenD>2*peakAmp0)=0;
screenD(screenD<0)=0;
screenD=screenD./peakAmp0;
%subplot(2,2,4)
figure
imagesc(screenD)

%%

scrnVac=screenD(yD(1):yD(end),xD(1):xD(end));

% Load Grating Spectrum
%%folder='V:\ARDB\E163\Data\130409\';
%filename='run2216_ScrnAvgs.mat';
%eval(['load ',folder,filename]);

%scrnAvgOff=scrnAvgOff(yD(1):yD(end),xD(1):xD(end));
%scrnAvgOff=imrotate(scrnAvgOff,-.5,'crop');
%spectraOff=mean(scrnAvgOff(350:475,:));
%Y1=spectraOff;

Y1=mean(scrnVac);
%Y1=spectra0;
%%
%figure
%c1=FitSpectrum5c(Y1,18,260,0,1);
x=1:length(Y1);
figure
plot(x,Y1,'k')
hold on
for i=1:4
    [c1, c1std, yfit, ChiSq]=FitSpectrum5c(Y1,i);
    plot(x,yfit,'b')
    hold on
    switch i
        case 1
            scale1=sqrt(2*log(2));
            scale2=sqrt(2*log(2));
        case 2
            scale1=1/2;
            scale2=1/2;
        case 3
            scale1=asech(sqrt(1/2));
            scale2=asech(sqrt(1/2));
        case 4
            scale1=1/2;
            scale2=asech(sqrt(1/2));
    end
    FWHM=c1(3)*scale1+c1(4)*scale2;
    if c1std(3)*scale1>c1std(4)*scale2
        FWHMstd=c1std(3)*scale1;
    else
        FWHMstd=c1std(4)*scale2;
    end
    
    % line([x(round(c1(2)-c1(3)*sqrt(2*log(2)))) x(round(c1(2)+c1(4)*sqrt(2*log(2))))],...
    %         [c1(5)+c1(1)/2 c1(5)+c1(1)/2], 'LineWidth',2,'Color','r','LineStyle',':');
    line([c1(2)-c1(3)*scale1 c1(2)+c1(4)*scale2],...
        [c1(5)+c1(1)/2 c1(5)+c1(1)/2], 'LineWidth',2,'Color','r');
    fprintf('\nmode %g: FWHM= %g += %g, Chi^2=%g\n',i,FWHM,FWHMstd,ChiSq)
end

%%
% Custom Morgenstemning colormap
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

xD=(xD-858*ones(size(xD)))*1.2;    %run2216
%% Figure
width=4;
height=2;

Navg=5;

%close all
hFig = figure(94);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

hold on

% % Vacuum Spectrum
scrnVacp=medfilt2(scrnVac,[Navg,Navg])/1.616;
imagesc(xD,yDp,scrnVacp)

% Grating Spectrum
%scrnAvgOffp=medfilt2(scrnAvgOff,[Navg,Navg])/1.051;
%imagesc(xD,yDp,scrnAvgOffp)

set(gca,'YTick',[5 9 13 17])
axis xy
%colorbar
colormap(col)
size(col)
%colormap(morgenstemning(256,'invert',1))
%caxis([0 1])
%caxis  %did this to figure out how to scale data so it would be [0 1]
xlim([-750 200])

hy=ylabel('Position (mm)');
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
%export_fig ScrnGlas.eps -cmyk -r300 -painters%-painters
%export_fig ScrGrat.eps -cmyk -r300 %-painters%-painters

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
%col=[col1(1:round(N*2/5),:);col2(round(4*N*2/5)+1:round(4*N*3/5),:);col1(round(N*3/5)+1:N*3/4,:);col3(8*N*3/4+1:8*N,:)];

col1(1,:) = [1 1 1];                 %white
col2 = [252 229 1]./255;        %yellow
col3=[253 160 80]/255;          %light orange
col4=[220 87 1]/255;            %dark orange
col5 = [20 50 95]./255;         %navy blue
col6=[0 0 0];                   %black

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

% Load Vacuum spectra
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
[cout, ~, ~,~, ybkgd, ~, ~, ~]=FitSpectrum5c(spectraF(400:end),1);
peakAmp0=cout(1);

screenD=screen-ybkgd(1)*ones(size(screen));
screenD(screenD<0)=0;
specVac=mean(screenD);

%remove x-rays
screenD(screenD>2*peakAmp0)=0;

scrnVac=screenD;%(yD(1):yD(end),xD(1):xD(end));
Navg=5;
scrnVac=medfilt2(scrnVac,[Navg,Navg]);%/1.616;

%% Load Glass Spectrum
filename='spectrum_run2064_event2.dat';
glass=load([folder,filename]);

glass=circshift(glass, [0 175]);
glass=imrotate(glass,-1,'crop');
screen=glass*1.185;    %to match the area of the projected spectrum

spectra0=mean(screen);

%re-scale the image for a nice plot
spectraF=filter(b,a,spectra0);
% This initial fit is used for various normalizations below
[cout, ~, ~,~, ybkgd, ~, ~, ~]=FitSpectrum5c(spectraF(400:end),1);
peakAmp0=cout(1);

screenD=screen-ybkgd(1)*ones(size(screen));
screenD(screenD<0)=0;
specGla=mean(screenD);
%remove x-rays
screenD(screenD>2*peakAmp0)=0;

scrnGla=screenD;%(yD(1):yD(end),xD(1):xD(end));
Navg=5;
scrnGla=medfilt2(scrnGla,[Navg,Navg]);%/1.616;

%%
width=4;
height=6;
close all
hFig = figure(94);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

% scrnVacS=zeros(size(scrnVac,1)/4,size(scrnVac,2));
% scrnGlaS=zeros(size(scrnGla,1)/4,size(scrnGla,2));
% for i=1:size(scrnVac,1)/4
%     scrnVacS(i,:)=sum(scrnVac(4*i-3:4*i,:))/4;
%     scrnGlaS(i,:)=sum(scrnGla(4*i-3:4*i,:))/4;
% end
%%
subplot(3,1,1)
%%
imagesc(xD,yDp,scrnVac)

set(gca,'YTick',[1 5 9 13 17])
axis xy
colormap(col0)
%caxis
caxis([0 1450])
%caxis  %did this to figure out how to scale data so it would be [0 1]
xlim([-750 100])
hy1=ylabel('Position [mm]');
set(gca,'XTickLabel',[],'XMinorTick','on');
ax1=gca;

%%
subplot(3,1,2)
%
imagesc(xD,yDp,scrnGla)

set(gca,'YTick',[1 5 9 13 17])
axis xy
colormap(col0)
caxis([5 250])
%%caxis([0 .95])
%caxis  %did this to figure out how to scale data so it would be [0 1]
xlim([-750 100])
hy2=ylabel('Position [mm]');
set(gca,'XTickLabel',[],'XMinorTick','on');
ax2=gca;
%%
subplot(3,1,3)
%
Y1=mean(scrnVac);
sum(Y1)

% plot(xD,Y1/1000,'r')
% hold on

h1=plot(xD,specVac/1000,'LineWidth',2,'Color',cols(1,:));%col2)
hold on
for i=1:1 %extend to 4 to test the other fitting algorithms (1 was best)
    [c1, c1std, yfit, ChiSq]=FitSpectrum5c(Y1,i);
    h1f=plot(xD,yfit/1000,':','LineWidth',1,'Color',cols(3,:));%col4)
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
    c1(3)*scale1
    c1(4)*scale2
    FWHM=c1(3)*scale1+c1(4)*scale2;
    if c1std(3)*scale1>c1std(4)*scale2
        FWHMstd=c1std(3)*scale1;
    else
        FWHMstd=c1std(4)*scale2;
    end
    
    %line([c1(2)-c1(3)*scale1 c1(2)+c1(4)*scale2],...
    %    [c1(5)+c1(1)/2 c1(5)+c1(1)/2], 'LineWidth',2,'Color','r');
    fprintf('\nmode %g: FWHM= %g += %g, Chi^2=%g\n',i,FWHM,FWHMstd,ChiSq)
end


h2=plot(xD,specGla/1000,'LineWidth',2,'Color',cols(2,:));%col3)
hold on

Y1=mean(scrnGla);
sum(Y1)

% plot(xD,Y1/1000,'k')
% hold on
for i=2:2 %extend to 4 to test the other fitting algorithms (2 was best)
    %[c1, c1std, yfit, ChiSq]=FitSpectrum5c(Y1,i);
    [c1, c1std, yfit, ChiSq]=FitSpectrum5c(specGla,i);
    h2f=plot(xD,yfit/1000,':','LineWidth',1,'Color',cols(4,:));%col5)
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
    c1(3)*scale1
    c1(4)*scale2
    FWHM=c1(3)*scale1+c1(4)*scale2;
    if c1std(3)*scale1>c1std(4)*scale2
        FWHMstd=c1std(3)*scale1;
    else
        FWHMstd=c1std(4)*scale2;
    end
    p2=xD(round(c1(2)))
    %line([c1(2)-c1(3)*scale1 c1(2)+c1(4)*scale2],...
    %    [c1(5)+c1(1)/2 c1(5)+c1(1)/2], 'LineWidth',2,'Color','r');
    fprintf('\nmode %g: FWHM= %g += %g, Chi^2=%g\n',i,FWHM,FWHMstd,ChiSq)
end
axis tight
xlim([-750 100])
hx3=xlabel('Energy deviation, {\it\DeltaE} [keV]');
hy3=ylabel('Charge density [a.u.]');
ax3=gca;
%
set(gca,'XMinorTick','on')
tightfig;
hL=legend([h1,h1f,h2,h2f],{'Vacuum','Gauss fit','Fused Silica','Lorentz fit'});
set(hL,'Location','NorthWest','Box','Off')
set([ax1,ax2,ax3],'FontSize', 10);
set([hx3,hy1,hy2,hy3,hL],'FontSize', 10);

% set(gca, ...
%     'TickDir'     , 'in'     , ...
%     'Box'         , 'on'     , ...
%     'TickLength'  , [.02 .02] , ...
%     'LineWidth'   , 1         );

%export_fig ScrVac.tif -cmyk -r300 %-painters%-painters
%export_fig ScrGrat.tif -cmyk -r300 %-painters%-painters

%export_fig ScrVac2.tif -rgb -r300 %-painters%-painters
%export_fig ScrGrat2.tif -rgb -r300 %-painters%-painters
%
%export_fig ScrTrans2.eps -cmyk -r300 -painters%-painters
%%
%export_fig ScrBulk.eps -cmyk -r300 -painters%-painters

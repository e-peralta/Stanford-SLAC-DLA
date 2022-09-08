%% Making the colorbars
clear all
close all
%dC1=.0417:.0417:1;
dC1=linspace(0,1,25);
%dC2=.0625:.0625:1;
dC2=linspace(0,1,17);
Chot=[[dC1',zeros(length(dC1),2)];...
    [ones(length(dC1),1),dC1',zeros(length(dC1),1)];...
    [ones(length(dC2),2),dC2']];
Chot=flipud(Chot);

Ccold=[[zeros(length(dC1),2),dC1'];...
    [zeros(length(dC1),1),dC1',ones(length(dC1),1)];...
    [dC2',ones(length(dC2),2)]];
%Ccold=flipud(Ccold);
Cboth=[Ccold;Chot];

%% All in units of microns
%sigma values are actually FWHM
FWHM2sigG=2*sqrt(2*log(2));
dist=1000;
Amp=.5;
sigEy=8/2;
sigEz=129/2;
zE=125;
yE=200;
sigLy=372;  %FWHM
sigLz=300/2;
zL=625;
yL=700;
x0=1:.5:dist; 

%for n=0:9
n=7.7; %enters at n=1,exits at n=9

vt=100*n;

%florentz=@(c,x) c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2);
%fgauss=@(c,x) c(1)*exp(-4*log(2)*(x-c(2)).^2/c(3)^2);
fgauss=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2/2);  %takes c3 as sigma
fsecant=@(c,x) c(1)*sech(2*asech(1/sqrt(2))*(x-c(2))./c(3)).^2; %takes c3 as FWHM

fEz=fgauss([Amp,zE+vt,sigEz],x0);
fEy=fgauss([Amp,yE,sigEy],x0)';
fLz=fgauss([Amp,zL,sigLz],x0);
fLy=fsecant([Amp,yL-vt,sigLy],x0)';

imag2=-fEy*fEz*2;  %changed factor to 3.5 at perfect overlap to see e-beam
                     % 2.5 one at +- one step from it

imag1=fLy*fLz*2;
imag=imag2+imag1;

% figure (1)
% imagesc(x0,x0,imag2)
% colormap(Chot)
% axis xy
% caxis([0 1])
% 
% figure (2)
% imagesc(x0,x0,imag1)
% colormap(Ccold)
% axis xy
% caxis([-1 0])


hFig=figure(round(3+n));
xwidth=17;%24;
ywidth=3.5;%5;

%set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','centimeters')
set(hFig, 'Position', [0 0 xwidth ywidth])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

set(gca,'FontName','Arial');
set(gca, 'FontSize', 12 );
set(gca, ...
  'TickDir'     , 'out'     , ...
  'Box'         , 'on'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 2         );
hy=ylabel('Vertical transverse position, y [\mum]');
hx=xlabel('Longitudinal position, z [\mum]');
set(hy, 'FontSize',12);
set(hx, 'FontSize',12);
set(gca,'YColor','k');
set(gca,'XColor','k');


%figure (3)
hold on
%hImg=imagesc(x0,x0,imag)
hImg=imagesc(x0-350*ones(size(x0)),x0,imag)
colormap(Cboth)
axis xy
%axis equal
caxis([-.55 .55])

hold on
%hC=contour(x0,x0,imag1,[.25],'LineColor','w','LineWidth',.5,'LineStyle',':');
hC=contour(x0-350*ones(size(x0)),x0,imag1,[.25],'LineColor','w','LineWidth',.5,'LineStyle',':');

hold on
text(-180,230, ['t = ',num2str((n-.74)/3), ' ps']);

hold on
%line([350 900],[201 201],'LineWidth',.2,'Color','k');
%line([350 900],[199 199],'LineWidth',.2,'Color','k');
line([0 550],[201 201],'LineWidth',.2,'Color','k');
line([0 550],[199 199],'LineWidth',.2,'Color','k');

%hold on
%hC=contour(x0,x0,imag2,[-.25],'LineColor','b','LineWidth',.2,'LineStyle',':');

if ywidth==3.5
    %ylim([150 250])
    ylim([160 240])
elseif ywidth==5
    ylim([100 300])
elseif ywidth==20
    ylim([0 800]) 
end
%xlim([0 1000])
xlim([-200 500])
set(gca,'YTick',[160 200 240])

%pause(.5)

%end

%figure (3)

%
%export_fig Overlap770_2.eps -painters -rgb

%%
%Cross-correlation 
FWHM2sigG=1/(2*sqrt(2*log(2)));

fgauss=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2/2);                  %takes c3 as sigma
fsecant=@(c,x) c(1)*sech(2*asech(1/sqrt(2))*(x-c(2))./c(3)).^2; %takes c3 as FWHM
x=1:.1:100;
%y=fgauss([1,50,10*FWHM2sigG],x);
y1=fgauss([1,50,25*FWHM2sigG],x);
y2=fsecant([1,50,25],x);

for i=1:length(x)
yavg1(i)=quad(@(y)fgauss([1,50,25*FWHM2sigG],y),0,x(i))/x(i);
yavg2(i)=quad(@(y)fsecant([1,50,25],y),0,x(i))/x(i);

end
figure
plot(x,y1,x,yavg1,x,y2,x,yavg2)
axis tight

%%
%clear all
close all
%cols(0,:)=[1 1 1];  %use 'w'
cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple
cols(5,:)=[0 0 0]; %use 'k'

width=3.25;
height=2.5;
hFig = figure(97);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
%set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')

clear L tau avgGz avgGt
L=10:1:1000;

tau=[1.25 .625 2.5 1.25 1.25];%1.25;%1.24;
sigZ=[310 310 310 155 620];%310;%250;%310;

for j=1:1%length(tau)
for i=1:length(L)
fgauss=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2/2);                  %takes c3 as sigma
fsecant=@(c,x) c(1)*sech(2*asech(1/sqrt(2))*(x-c(2))./c(3)).^2; %takes c3 as FWHM

%avgGz(i)=2*quad(@(y)fgauss([1,0,sigZ],y),0,L(i)/2)/L(i);
%avgGt(i)=2*quad(@(y)fsecant([1,0,tau],y),0,L(i)*um/c/ps/2)/(L(i)*um/c/ps);
%avgG(i)=avgGz*avgGt;

Pulse=@(z,t) exp(-z.^2/sigZ(j)^2/2).*(sech(2*asech(1/sqrt(2))*t./tau(j)).^2);
avgG(i)=4*quad2d(Pulse,0,L(i)/2,0,L(i)*um/c/ps/2)/L(i)/(L(i)*um/c/ps);

end
eFf=avgG.*L/L(end);
%subplot(1,2,1)
plot(L/1000,avgG,'Color',cols(6-j,:),'Linewidth',1);
%plot(L/1000,eFf,'Color',cols(6-j,:),'Linewidth',1);

hold on

[maxV ind]=max(eFf);
L(ind)
end
ax1=gca;
set(gca,'Color','none');

% ylim([0 180])
% set(ax1,'YTick',20:20:180)
% set(ax1,'XTickLabel',{'0.03','','','','','','','0.1','','0.3','','0.5','','','','','1.0','','2','3','4','5'})
% hy=ylabel('Energy gain, \Delta{\itE} [keV]');


%ylim([0 1.12])
set(ax1,'YTick',.1:.1:1)
hy=ylabel('Energy form factor, f_{F,E} / 1 mm');

% xlim([.03 5])
set(ax1,'XTick',.1:.1:1)
% set(ax1,'XTickLabel',{'0.03','','','','','','','0.1','','0.3','','0.5','','','','','1.0','','2','3','4','5'})
hx=xlabel('Interaction length, L [mm]');

ax2=copyobj(ax1,hFig);
delete(get(ax2,'Children'));
set(ax2,'Color','None','Box','off','Ygrid','on','Xgrid','on');
set(ax2,'Xcolor',[.9 .9 .9],'Ycolor',[.9 .9 .9]);%,'XTickLabel',[],'YTickLabel',[]);
%get(ax3)    %Unsure how to remove labels from this third axes.. will have
%to do by hand
uistack(ax2,'bottom');

linkaxes([ax1,ax2],'xy')

%hL=legend('f_F(\tau_p,\sigma_{lz})','f_F(.5\tau_p,\sigma_{lz})','f_F(2\tau_p,\sigma_{lz})','f_F(\tau_p,.5\sigma_{lz})','f_F(\tau_p,2\sigma_{lz})');
%set(hL,'Location','southwest','Box','Off')

set([ax1,ax2,hx,hy],'FontSize', 10);

%export_fig FormFactor2.eps -cmyk -r300 -painters%-painters

%%
N=625;
%close all
%hold on
theta=.00001:.000001:.1*pi/180;
 G=(1+sin(theta))/(2+sin(theta))*sin(2*pi*N*sin(theta))./(pi*N*sin(theta));
 plot(theta*1e3,G)
 

%%
clear all
close all
%load raw_model_sim_data
%load V:\ARDB\E163\Data\2013\130513\raw_model_sim_data
load ~/Dropbox/Research/Data/raw_model_sim_data 

%%
cout=FitSpectrum5c(raw_on,18,340,50);
[yfit, ymain, ysignal, ysig1, ysig2]=DataFitPlotter(cout,18);

plot(1:1024,raw_off,'x',1:1024,raw_on,'x',...
   1:1024,fit_full_off,1:1024,fit_main_off,...
   1:1024,fit_main_on,1:1024,fit_plateau_off,...
   1:1024,fit_signal_off,1:1024,yfit,1:1024,ysignal)
xlim([600 1000])

%%
%figure(2)

hFig = figure(2);
xwidth=8.9;
ywidth=5;

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','centimeters')
set(hFig, 'Position', [0 0 xwidth ywidth])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

hold on


%h1=plot(x0,raw_off-fit_main_off-fit_plateau_off,'bx','Color',[.65 .65 1],'LineWidth',1)
h1=plot(x0,raw_off-fit_main_off-fit_plateau_off,'Color',[.65 .65 1],'LineWidth',1)

% 
% [nh, xh] = hist(EnergyIn,[-400:5:400]*1e3);
% %%plot(xh*1e-3,nh/max(nh)*.2175,'bs','LineWidth',3)
% 
%h3=plot(x0,raw_on-fit_main_on-fit_plateau_off,'rx','Color',[1 .65 .65],'LineWidth',1)
h3=plot(x0,raw_on-fit_main_on-fit_plateau_off,'Color',[1 .65 .65],'LineWidth',1)

h2=plot(x0,fit_signal_off,'LineWidth',1,'Color',[0 0 .75])
% 
% h4=plot(x0,model_signal_on,'LineWidth',1,'Color',[.75 0 0])
%     
%  [nh, xh] = hist(EnergyOut,[-400:8:400]*1e3);
%  h5=plot(xh*1e-3,nh/max(nh)*.1,'ko','LineWidth',1)
% 
% set(h5                         , ...
%   'Marker'          , 'o'         , ...
%   'MarkerSize'      , 2           , ...
%   'MarkerEdgeColor' , 'k'      , ...
%   'MarkerFaceColor' , 'k' );

h6=plot(x0(820:end),ysignal(820:end)-fit_main_on(820:end)-.8*fit_plateau_off(820:end),'LineWidth',1,'Color',[.75 0 0])

hx=xlabel('Energy deviation {\it\DeltaE} (keV)')
hy=ylabel('Charge density (arb. unit)')
%title('Adding Simulation Data')

axis([-120 120 -.005 .225])

box on
%enhance_plot;


%xwidth=600;
%ywidth=300;

%hL=legend('Laser Off','Spec. Fit','Laser On','Model','Simulation','Location','Northeast')
%hL=legend('Laser Off','Location','Northeast')
hL=legend('Laser Off','Spec. Fit Off','Laser On','Spec. Fit On','Location','Northwest')
set(hL, 'FontSize'   , 6           );
legend boxoff

 set( gca                       , ...
    'FontName'   , 'Arial' );
set([hx, hy], ...
    'FontName'   , 'Arial');
set(gca             , ...
    'FontSize'   , 6           );
set([hx,hy]  , ...
    'FontSize'   , 7          );

set(gca, ...
  'TickDir'     , 'out'     , ...
  'Box'         , 'on'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );

% set(gcf, 'renderer', 'painters');
% print -depsc2 test.eps;  %for plots
 
 export_fig test6.eps -painters -rgb

 
 
 %print -dtiffn test2.tiff; %for images
close all
%folder='V:\ARDB\E163\Data\130814\';
folder='~/Documents/MATLAB/';
filename='AllTheTeeth_Corrected.tsv';
eval(['load ',folder,filename]);

hFig = figure(2);
%xwidth=8.9;
%ywidth=5;
width=4;
height=2;

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
%set(hFig,'Units','centimeters')
set(hFig,'Units','inches')

%set(hFig, 'Position', [0 0 width height])

 set(gcf, 'Color', 'w');
set(hFig,'Units','points')

hold on

Energy=AllTheTeeth_Corrected(:,1);
y1c=AllTheTeeth_Corrected(:,2);

E=(Energy-824*ones(size(Energy)))*1.2;

col=[0 .6 .4];
Green=[0 .8 .2];

col1=[230 97 1]/255;
col2=[253 184 99]/255;
col3=[178 171 210]/255;
col4=[94 60 153]/255;

%plot(E,y1c./max(y1c)*1.15,'.','Color',Green,'MarkerSize',4)
%plot(E,y1c./max(y1c)*1.15,'.','Color','k','MarkerSize',4)
%h1=plot(E,y1c./max(y1c)*1.15,'.','Color','c','MarkerSize',4)


%load 'V:\ARDB\E163\Data\130409\run2216_f18-4_1.mat';
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
      h2=plot(xVar,spectra(ind,x)./Peak1amp(ind),'LineWidth',3,'Color',col2);
      
      hold on
      %hOff_fit=plot(xF,yfit,'b','LineWidth',2);
   end
end

h1=plot(E,y1c./max(y1c)*1.15,'.','Color',col3,'MarkerSize',4)


hold on
c=[14.69 .1 .7324 .2796];
fShift0=c(1)*exp(-(0-c(2))^c(3)/c(4));
xVar=xVar-fShift0*ones(1,length(x));

%plot(xVar,yfit,'Color',[1 .7 0],'LineWidth',.75);%,'LineStyle','--')
%plot(xVar,yfit,'Color','r','LineWidth',.75);%,'LineStyle','--')
%h3=plot(xVar,yfit,'Color',[1 .7 0],'LineWidth',.75);%,'LineStyle','--')
h3=plot(xVar,yfit,'Color',col4,'LineWidth',1);%,'LineStyle','--')

hx=xlabel('Energy deviation {\it\DeltaE} [keV]')
hy=ylabel('Charge density [a.u.]')
xlim([-750 200])
ylim([0 1])
%hL=legend('Simulation','Data','Spec. Fit')
hL=legend([h2,h1,h3],{'Data','Simulation','Spec. Fit'})

set(hL, 'FontSize'   , 10           );
legend boxoff

%  set( gca                       , ...
%     'FontName'   , 'Arial' );
% set([hx, hy], ...
%     'FontName'   , 'Arial');
set(gca             , ...
    'FontSize'   , 10           );
set([hx,hy]  , ...
    'FontSize'   , 10          );

set(gca,'Box'         , 'on');
%   'TickDir'     , 'out'     , ...
%   'Box'         , 'on'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'LineWidth'   , 1         );

% set(gcf, 'renderer', 'painters');
% print -depsc2 test.eps;  %for plots
 
% export_fig SpecAll.eps -painters -cmyk%-rgb%-cmyk



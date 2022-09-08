%% AnalyicalMode.m
%
%   Script to calculate broadening of electron enegry spectrum due to laser
%   acceleration when the electron pulse is longer than a small fraction of
%   the optical pulse length
%
% --------------------------------

clear all
close all
xwidth=3.25;             % Figure width
ywidth=1.5;             % Figure height

Amp=50;                % MeV/m

DETAILS=0;
SMALL=0;                % Whether to apply small amplitude broadening

C= 299792458;
wave=800e-9;
tau=wave/C*1e15;        % single cycle in fs (~2.7 fs)

wp = 8.23;              % leading edge HWHM width in keV
wm = 15.4;              % trailing edge HWHM width in keV


dE=.01;                 % Energy resolution

Elim=(wp+wm)+2*Amp;     % Energy limits

E1 = [-Elim:dE:Elim];   %

N= 2*Amp;               % Number of time slices to use in analysis.


wp = wp/sqrt(2*log(2)); % initial distribution leading edge RMS width in keV
wm = wm/sqrt(2*log(2)); % initial distribution trailing edge RMS width in keV


f1 = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*(wm/wp*c(3))^2)).*(x<=c(2)) ...
    + c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)).*(x>c(2));

%%

cols(1,:)=[255 205 135]/255; %light Orange
cols(2,:)=[166 160 220]/255;  %light purple
cols(3,:)=[205 66 0]/255;      %dark orange
cols(4,:)=[47 08 117]/255;     %dark purple

% newcols=morgenstemning(20);



%%

cmap=colormap(morgenstemning(34));
cshift=7;

hFig1 = figure(7);
%set(hFig1,'ActivePositionProperty','position')
%set(hFig1,'Units','centimeters')
set(hFig1,'Units','inches')
set(hFig1, 'Position', [0 0 xwidth ywidth])

set(gcf, 'Color', 'w');
set(hFig1,'Units','points')

hold on

hx1=xlabel('Time (fs)');
hy1=ylabel('Amplitude (arb. unit)');
hold on

% cosine

t1 = [0:.0001:1]*tau;
F1 = -cos(t1*(2*pi)/tau);        % 180 phase offset for display convenience
plot(t1,F1,'--k','LineWidth',1)

% ________________%
% inject N bunches in the time domain
t2 = [0:.025:.25]*tau;
F2 = -cos(t2*(2*pi)/tau);
%cm = colormap(autumn(length(t2)));
cm = cmap(cshift:cshift+10,:);


for i1 = 1:length(t2)
    stem(t2(i1),F2(i1),'LineWidth',1,'Color',cm(i1,:),'Marker','none')
end

% corresponding N pulses in the energy domain

hFig2 = figure(8);
%set(hFig2,'ActivePositionProperty','position')
%set(hFig2,'Units','centimeters')
set(hFig2,'Units','inches')

set(hFig2, 'Position', [0 0 xwidth ywidth])

set(gcf, 'Color', 'w');
set(hFig2,'Units','points')

hold on

hx2=xlabel('Energy Deviation, {\it\DeltaE} (keV)');
hy2=ylabel('Charge density (arb. unit)');


for i1 = 1:length(t2)
    s1 = f1([1,Amp*cos(2*pi*t2(i1)/tau),wp],E1);
    plot(E1,s1,'Color',cm(i1,:),'LineWidth',1)
end
%--------------------%


% inject N bunches in the time domain
figure (7)
hold on

t2 = [.5:-.025:.25]*tau;
F2 = -cos(t2*(2*pi)/tau);
%cm = colormap(summer(length(t2)));
cm = flipud(cmap(cshift+11:cshift+21,:));

for i1 = 1:length(t2)
    stem(t2(i1),F2(i1),'LineWidth',1,'Color',cm(i1,:),'Marker','none')
end

% corresponding N pulses in the energy domain
figure (8)

hold on

for i1 = 1:length(t2)
    s1 = f1([1,-Amp*cos(2*pi*t2(i1)/tau),wp],E1);
    plot(E1,s1,'Color',cm(i1,:),'LineWidth',1)
end

% inject N bunches in the time domain
figure (7)
hold on
t2 = [.5:.025:.75]*tau;
F2 = -cos(t2*(2*pi)/tau);
for i1 = 1:length(t2)
    stem(t2(i1),F2(i1),'LineWidth',1,'Color',cm(i1,:),'Marker','none')
end


% corresponding N pulses in the energy domain
figure (8)

hold on

for i1 = 1:length(t2)
    s1 = f1([1,-Amp*cos(2*pi*t2(i1)/tau),wp],E1);
    plot(E1,s1,'Color',cm(i1,:),'LineWidth',1)
end

% inject N bunches in the time domain
figure (7)
hold on
t2 = [1:-.025:.75]*tau;
F2 = -cos(t2*(2*pi)/tau);
%cm = colormap(autumn(length(t2)));
cm = cmap(cshift:cshift+10,:);
for i1 = 1:length(t2)
    stem(t2(i1),F2(i1),'LineWidth',1,'Color',cm(i1,:),'Marker','none')
end
axis([0 tau -1.01 1.01])

ax1=gca;
set(ax1,'YTick',-1:.5:1)
set(ax1,'XTick',0:.5:2.5)
set(ax1, 'Color','none','FontSize', 10 );
set([hx1,hy1], 'FontSize',10);

tightfig;
% set(gca, ...
%     'TickDir'     , 'out'     , ...
%     'Box'         , 'on'     , ...
%     'TickLength'  , [.02 .02] , ...
%     'LineWidth'   , 1         );

export_fig Anal1.eps -cmyk -r300 -painters


% corresponding N pulses in the energy domain
figure (8)

hold on

for i1 = 1:length(t2)
    s1 = f1([1,-Amp*cos(2*pi*t2(i1)/tau),wp],E1);
    plot(E1,s1,'Color',cm(i1,:),'LineWidth',1)
end
%--------------------%

%ylabel('Intensity')
%title('Energy Domain')

axis([-Elim Elim 0 1.01])

%enhance_plot
%xlabel('\Delta Energy [keV]','FontSize',18)
%set(gca,'xticklabel',{[]})

ax2=gca;
set(ax2,'YTick',0:.25:1)
set(ax2, 'Color','none','FontSize', 10 );
set([hx2,hy2], 'FontSize',10);

tightfig;


% set( gca                       , ...
%     'FontName'   , 'Arial' );
% set([hx2, hy2], ...
%     'FontName'   , 'Arial');
% set(gca             , ...
%     'FontSize'   , 6           );
% set([hx2,hy2]  , ...
%     'FontSize'   , 7          );
% set(gca, ...
%     'TickDir'     , 'out'     , ...
%     'Box'         , 'on'     , ...
%     'TickLength'  , [.02 .02] , ...
%     'LineWidth'   , 1         );

export_fig Anal2.eps -cmyk -r300 -painters

%%
% The resulting distribution - Spectrometer image
% -------------------- 4 ------------------------%

hFig = figure(9);

set(hFig,'Units','inches')

set(hFig, 'Position', [0 0 xwidth ywidth])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')

hold on

hx3=xlabel('Energy gain {\it\DeltaE} [keV]');
hy3=ylabel('Charge density [arb. unit]');
hold on

s2 = 0;
s2p=0;
s2p2=0;
s2p3=0;
s2n=0;

if Amp<100
    ii=N/4;    %25
    fi=3*N/4;  %75
elseif Amp<400 
    ii=3*N/8;  %75
    fi=5*N/8;  %125
else
    ii=7*N/16;
    fi=9*N/16;
end

for i1 = 0:N-1
    %for i1 = N/4:3*N/4%0:N-1
    s1 = f1([1,-Amp*cos(2*pi/N*i1),wp],E1);
    s2 = s2+s1;
% used this when i1=1:N-1 to test partial integration    
%     if i1>=N/4 && i1<3*N/4
%         s2p = s2p+s1;
%         if i1>=3*N/8 && i1<5*N/8
%             s2p2=s2p2+s1;
%             if i1>=15*N/32 && i1<17*N/32%7*N/16 && i1<9*N/16
%                 s2p3=s2p3+s1;    
%             end
%         end
%     else
%         s2n = s2n+s1;
%     end
 if i1>=ii && i1<fi
                s2p=s2p+s1;    
 end

    if (i1==0 || i1==N/2) && DETAILS
        hFig = figure(9);
        plot(E1,s1,':k','LineWidth',1)
        hold on
    end
end


s0 = f1([N,0,wp],E1);

w0=wp*sqrt(2*log(2));  %convert back to HWHM from RMS

hold on
h1=plot(E1,s0/N,'LineWidth',2,'Color',cols(2,:));
hold on

h2=plot(E1,s2/N,'LineWidth',2,'Color',cols(4,:));

%h4=plot(E1,s2n/N,'r--','LineWidth',2);

% this is to plot a partial positive sum
%h3=plot(E1,s2p/N,'k--','LineWidth',2);

%h5=plot(E1,s2n/N+s2p/N,'r:','LineWidth',2);

%h5=plot(E1,s2p2/N,'k-.','LineWidth',2);
%h6=plot(E1,s2p3/N,'k:','LineWidth',2);
%legend('Initial','Final')


hold on


[Vmax Pmax]=max(s2p);
HWP=interp1(s2p(Pmax:Pmax+Amp/dE/2),E1(Pmax:Pmax+Amp/dE/2),Vmax/2)
% if SMALL
%     [cout sfit]=FitSpectrum(s2,4,2*Amp/dE);
%     %[cout sfit]=FitSpectrum(s2p,4,Amp/dE);
%     w1=cout(4)*sqrt(log(2))*dE;
%     pP=E1(round(cout(2)));
%     aP=cout(1)/N;
% else
%     [cout sfit temp temp  y1 y2 temp temp coutstd]=FitSpectrum(s2,14,2*Amp/dE);
%     %[cout sfit temp temp  y1 y2 temp temp coutstd]=FitSpectrum(s2p,14,Amp/dE);
%     w1=cout(6)*sqrt(log(2))*dE;
%     pP=E1(round(cout(2)+cout(7)));
%     aP=cout(5)/N;
%     %pP1=((cout(2)+cout(7))/length(E1)*2*E1(end)+E1(1))
%     %pP_er=sqrt(coutstd(2)^2+coutstd(7)^2)/length(E1)*2*E1(end)
%     %w1_er=coutstd(6)*sqrt(log(2))*dE
%     if DETAILS
%         hold on
%         plot(E1,y1/N,'g')
%         
%         hold on
%         plot(E1,y2/N,'g')
%     end
% end
% 
% hold on
% plot(E1,sfit/N,'--k','LineWidth',2)
% 
% 
% hold on
% line([w0 pP+w1],[1/2 1/2],'LineWidth',1,'Color','r');
% hold on
% line([pP+w1 pP+w1],[9*aP/20 11/20],...
%     'LineWidth',1,'Color','r');
% hold on
% line([w0 w0],[9/20 11/20],...
%     'LineWidth',1,'Color','b');
% 
% DE=pP+w1-w0
% I1=trapz(s2(E1<0));
% I2=trapz(s2(E1>=0));
% I2/I1;


%xlim([-Elim/2 Elim/2])
%ylim([0 .5])
%legend('Initial','Final', 'Neg. Int.','Pos. Int.','Pos. 1/2 Int.','Pos. 1/4 Int',...
%    'Location','NorthWest')
%title('400 Mev Gradient')

%
% figure
% MaxE=0:200;
% Grad=50*sqrt((.0484*MaxE+1.0463).^2-1)-23;
% plot(MaxE,Grad)
% 
% 
% xlim([50 200])

%  ylim([0 .5])

% --------- load particle simulation and plot -------------- %

FOLDER='~/Dropbox/Research/Data/';
load([FOLDER,'trnsprtd_100prd_2401prtcls_NormDist_EsprdNominal_50KeVGain.mat']);

[nh, xh] = hist(EnergyOut,[-400:6:400]*1e3);
%h3=plot(xh*1e-3,nh/max(nh)*30/100,'kx','LineWidth',1)
h3=plot(xh*1e-3,nh/max(nh)*30/100,'LineStyle','none','Marker','o',...
    'Color',cols(3,:),'Linewidth',1,...
    'MarkerSize',2,'MarkerFaceColor',cols(1,:));
    
%hL=legend('Laser Off','Laser On','Simulation')

% set(hL, 'FontSize'   , 6           );
% legend boxoff
%
%  set( gca                       , ...
%     'FontName'   , 'Arial' );
% set([hx3, hy3], ...
%     'FontName'   , 'Arial');
% set(gca             , ...
%     'FontSize'   , 6           );
% set([hx3,hy3]  , ...
%     'FontSize'   , 7          );

axis([-Elim Elim 0 1.01])
ax2=gca;
set(ax2,'YTick',0:.25:1)
set(ax2, 'Color','none','FontSize', 10 );
%set([hx3,hy3,hL], 'FontSize',10);
set([hx3,hy3], 'FontSize',10);


tightfig;
% set(gca, ...
%   'TickDir'     , 'out'     , ...
%   'Box'         , 'on'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'LineWidth'   , 1         );
%
export_fig Anal3.eps -painters -cmyk -r300
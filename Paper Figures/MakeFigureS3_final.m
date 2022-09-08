%% inject N bunches to show time - energy domain effect of modulation
% -------------------- 3 ------------------------%
xwidth=8.9;
ywidth=4.5;

hFig1 = figure(7);
set(hFig1,'ActivePositionProperty','position')
set(hFig1,'Units','centimeters')
set(hFig1, 'Position', [0 0 xwidth ywidth])

set(gcf, 'Color', 'w');
set(hFig1,'Units','points')

hold on

hx1=xlabel('Time (ps)')
hy1=ylabel('Amplitude (arb. unit)')
hold on

% cosine
x1 = [0:.0001:1]*2.7;
m1 = -cos(x1*(2*pi)/2.7);
plot(x1,m1,'--k','LineWidth',1)
%________________%
% inject N bunches in the time domain
x2 = [0:.025:.25]*2.7;
m2 = -cos(x2*(2*pi)/2.7);
cm = colormap(autumn(length(x2)));
for i1 = 1:length(x2)
    stem(x2(i1),m2(i1),'LineWidth',1,'Color',cm(i1,:),'Marker','none')
end


%% corresponding N pulses in the energy domain

hFig2 = figure(8);
set(hFig2,'ActivePositionProperty','position')
set(hFig2,'Units','centimeters')
set(hFig2, 'Position', [0 0 xwidth ywidth])

set(gcf, 'Color', 'w');
set(hFig2,'Units','points')

hold on

hx2=xlabel('Energy Deviation {\it\DeltaE} (keV)')
hy2=ylabel('Charge density (arb. unit)')

wp = 8.23; % leading edge width
wm = 15.4; % trailing edge width
x1 = [-100:.001:100];

hold on
f1 = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*(wm/wp*c(3))^2)).*(x<=c(2)) ...
    + c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)).*(x>c(2));

for i1 = 1:length(x2)
    s1 = f1([1,50*cos(2*pi*x2(i1)/2.7),wp],x1);
    plot(x1,s1,'Color',cm(i1,:),'LineWidth',1)
end
%--------------------%

% inject N bunches in the time domain
figure (7)
hold on

x2 = [.5:-.025:.25]*2.7;
m2 = -cos(x2*(2*pi)/2.7);
cm = colormap(summer(length(x2)));
for i1 = 1:length(x2)
    stem(x2(i1),m2(i1),'LineWidth',1,'Color',cm(i1,:),'Marker','none')
end

% corresponding N pulses in the energy domain
figure (8)

hold on
f1 = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*(wm/wp*c(3))^2)).*(x<=c(2)) ...
    + c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)).*(x>c(2));

for i1 = 1:length(x2)
    s1 = f1([1,-50*cos(2*pi*x2(i1)/2.7),wp],x1);
    plot(x1,s1,'Color',cm(i1,:),'LineWidth',1)
end

% inject N bunches in the time domain
figure (7)
hold on
x2 = [.5:.025:.75]*2.7;
m2 = -cos(x2*(2*pi)/2.7);
cm = colormap(summer(length(x2)));
for i1 = 1:length(x2)
    stem(x2(i1),m2(i1),'LineWidth',1,'Color',cm(i1,:),'Marker','none')
end


% corresponding N pulses in the energy domain
figure (8)

hold on
f1 = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*(wm/wp*c(3))^2)).*(x<=c(2)) ...
    + c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)).*(x>c(2));

for i1 = 1:length(x2)
    s1 = f1([1,-50*cos(2*pi*x2(i1)/2.7),wp],x1);
    plot(x1,s1,'Color',cm(i1,:),'LineWidth',1)
end

% inject N bunches in the time domain
figure (7)
hold on
x2 = [1:-.025:.75]*2.7;
m2 = -cos(x2*(2*pi)/2.7);
cm = colormap(autumn(length(x2)));
for i1 = 1:length(x2)
    stem(x2(i1),m2(i1),'LineWidth',1,'Color',cm(i1,:),'Marker','none')
end
axis([0 2.7 -1.01 1.01])



 set( gca                       , ...
    'FontName'   , 'Arial' );
set([hx1, hy1], ...
    'FontName'   , 'Arial');
set(gca             , ...
    'FontSize'   , 6           );
set([hx1,hy1]  , ...
    'FontSize'   , 7          );

set(gca, ...
  'TickDir'     , 'out'     , ...
  'Box'         , 'on'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );

 export_fig Anal1.eps -painters -rgb
 
 
% corresponding N pulses in the energy domain
figure (8)

hold on
f1 = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*(wm/wp*c(3))^2)).*(x<=c(2)) ...
    + c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)).*(x>c(2));

for i1 = 1:length(x2)
    s1 = f1([1,-50*cos(2*pi*x2(i1)/2.7),wp],x1);
    plot(x1,s1,'Color',cm(i1,:),'LineWidth',1)
end
%--------------------%

%ylabel('Intensity')
%title('Energy Domain')
axis([-100 100 0 1.01])
%enhance_plot
%xlabel('\Delta Energy [keV]','FontSize',18)
%set(gca,'xticklabel',{[]}) 

 set( gca                       , ...
    'FontName'   , 'Arial' );
set([hx2, hy2], ...
    'FontName'   , 'Arial');
set(gca             , ...
    'FontSize'   , 6           );
set([hx2,hy2]  , ...
    'FontSize'   , 7          );
set(gca, ...
  'TickDir'     , 'out'     , ...
  'Box'         , 'on'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );

 export_fig Anal2.eps -painters -rgb
 
%% 
% The resulting distribution - Spectrometer image
% -------------------- 4 ------------------------%

hFig3 = figure(9);
set(hFig3,'ActivePositionProperty','position')
set(hFig3,'Units','centimeters')
set(hFig3, 'Position', [0 0 xwidth ywidth])

set(gcf, 'Color', 'w');
set(hFig3,'Units','points')

hold on

hx3=xlabel('Energy Deviation {\it\DeltaE} (keV)')
hy3=ylabel('Charge Density (arb. unit)')
hold on

N=100;
s2 = x1*0;
f1 = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*(wm/wp*c(3))^2)).*(x<=c(2)) ... 
    + c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)).*(x>c(2));

for i1 = 0:N-1
    s1 = f1([1,-50*cos(2*pi/N*i1),wp],x1);
    s2 = s2+s1;
end

s0 = f1([N,0,wp],x1);

hold on
h1=plot(x1,s0/100,'--b','LineWidth',1)
h2=plot(x1,s2/100,'r','LineWidth',1)
%legend('Initial','Final')
axis([-100 100 0 (N+1)/100])

% --------- load particle simulation and plot -------------- %

load('trnsprtd_100prd_2401prtcls_NormDist_EsprdNominal_50KeVGain.mat');

[nh, xh] = hist(EnergyOut,[-400:6:400]*1e3);
h3=plot(xh*1e-3,nh/max(nh)*33/100,'kx','LineWidth',1)
%title(['Energy gain: ' num2str(dE(i2)) 'keV'])
hL=legend('Laser Off','Laser On','Simulation')

set(hL, 'FontSize'   , 6           );
legend boxoff

 set( gca                       , ...
    'FontName'   , 'Arial' );
set([hx3, hy3], ...
    'FontName'   , 'Arial');
set(gca             , ...
    'FontSize'   , 6           );
set([hx3,hy3]  , ...
    'FontSize'   , 7          );

set(gca, ...
  'TickDir'     , 'out'     , ...
  'Box'         , 'on'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );

export_fig Anal3.eps -painters -rgb
%cols(0,:)=[1 1 1];  %use 'w'
cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple
%cols(5,:)=[0 0 0]; %use 'k'

folder='~/Dropbox/Research/Data/';

width=3.6;
height=2;
close all
hFig = figure(94);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')


%drift distance
%dl = 1.296; % 530 drift to screen 585
%dl = 0.703; % 560 drift to screen 585
dl=30.20988-28.76048;%-0.1306/2;

BeamE=60;
gamma=BeamE/.511;
beta=1;
c=299792458;
m=9.11e-31;
q=1.602e-19;
factor=.1*q/m/c/gamma;  %.1 because 1T=10kG; converts to 1/m2
%factor=.1;
%data is from QuadscanX_130224

filename='QuadscanX_130224.dat';
data=importdata([folder,filename]);

bees = data(:,1); % bvals in kG
s11 = data(:,2); % spot size
s22 = data(:,3); % spot size

s11 = s11/sqrt(2);
s22 = s22/sqrt(2);

% convert units
k = bees*factor;
kf=linspace(k(1),k(end),100);
s11 = s11*1e-6;
s22 = s22*1e-6;

% find fit
f=@(c) sum((c(1)*(k-c(2)).^2+c(3)-s11).^2);
[c, rm] = fminsearch(f,[1e-3 4*factor 1e-4]);

% calculate normalized emittance
emx = sqrt(c(1)*c(2))/dl^2*gamma;

h1a=plot(k,s22*1e6,'v','LineWidth',1,'Color',cols(3,:),...
    'MarkerSize',4,'MarkerFaceColor',cols(1,:));
hold on
h1b=plot(k,s11*1e6,'v','LineWidth',1,'Color',cols(4,:),...
    'MarkerSize',4,'MarkerFaceColor',cols(2,:));
hold on
h1f=plot(kf,(c(1)*(kf-c(2)).^2+c(3))*1e6,':','LineWidth',1,'Color',cols(4,:));

%data2 is from QuadscanY_130224
filename='QuadscanY_130224.dat';
data2=importdata([folder,filename]);
bees = data2(:,1); % bvals in kG
s11 = data2(:,2); % spot size
s22 = data2(:,3); % spot size

s11 = s11/sqrt(2);
s22 = s22/sqrt(2);

% convert units
k = bees*factor;
kf=linspace(k(1),k(end),100);
s11 = s11*1e-6;
s22 = s22*1e-6;

% find fit
f=@(c) sum((c(1)*(k-c(2)).^2+c(3)-s22).^2);
[c, rm] = fminsearch(f,[1e-3 4*factor 1e-4]);

% calculate normalized emittance
emy = sqrt(c(1)*c(2))/dl^2*gamma;

h2a=plot(k,s22*1e6,'o','LineWidth',1,'Color',cols(3,:),...
    'MarkerSize',4,'MarkerFaceColor',cols(1,:));
hold on
h2b=plot(k,s11*1e6,'o','LineWidth',1,'Color',cols(4,:),...
    'MarkerSize',4,'MarkerFaceColor',cols(2,:));
hold on
h2f=plot(kf,(c(1)*(kf-c(2)).^2+c(3))*1e6,':','LineWidth',1,'Color',cols(3,:));

hy=ylabel('RMS spotsize [\mu m]');
hx=xlabel('Focusing strength [1/m^2]');

ax1=gca;

%title(['emittance = ' num2str(emx) ' X ' num2str(emy)])
fprintf(['emittance = ' num2str(emx) ' X ' num2str(emy) '\n'])

axis tight

tightfig;
%hL=legend([h1a,h1b,h2a,h2b,h1f,h2f],{'Scan 1 X','Scan 1 Y','Scan 2 X','Scan 2 Y','Quad. fit X','Quad. fit Y'});
%set(hL,'Location','northeastoutside','Box','Off')
set(ax1,'FontSize', 10);
set([hx,hy],'FontSize', 10);

export_fig QuadScan.eps -cmyk -r300 -painters%-painters

%%
hFig = figure(95);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')

filename='PMQscan_121022.dat';
data=importdata([folder,filename]);

bees = data(:,1); % bvals in kG
s11 = data(:,2); % spot size
s22 = data(:,3); % spot size

s11 = s11*sqrt(2*log(2));
s22 = s22*sqrt(2*log(2));

% convert units
k = bees;
kf=linspace(k(1),k(end),100);
s11 = s11;
s22 = s22;

% find fit
%f=@(c) sum((c(1)*(k-c(2)).^2+c(3)-s11).^2);
%[c, rm] = fminsearch(f,[1 73 -4000]);
c=[1.2 -73 1170];
h1f=plot(kf,(c(1)*kf.^2+c(2)*kf+c(3)),':','LineWidth',1,'Color',cols(3,:));
hold on

c=[-.727 140];
h2f=plot(kf,(c(1)*kf+c(2)),':','LineWidth',1,'Color',cols(4,:));
hold on
h1a=plot(k,s22,'o','LineWidth',1,'Color',cols(3,:),...
    'MarkerSize',4,'MarkerFaceColor',cols(1,:));
hold on
h1b=plot(k,s11,'o','LineWidth',1,'Color',cols(4,:),...
    'MarkerSize',4,'MarkerFaceColor',cols(2,:));


hy=ylabel('RMS spotsize [\mu m]');
hx=xlabel('Distance from last PMQ [mm]');

ax1=gca;

axis tight

tightfig;
%hL=legend([h1a,h1b,h1f],{'Scan 1 X','Scan 1 Y','Quad. fit X'});
%set(hL,'Location','northeastoutside','Box','Off')
set(ax1,'FontSize', 10);
set([hx,hy],'FontSize', 10);

export_fig PMQscan.eps -cmyk -r300 -painters%-painters


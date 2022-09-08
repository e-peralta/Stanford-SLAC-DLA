%% Setup stuff
%set(0,'defaultfigureposition',[296 342 560 420]);
set(0,'defaultfigureposition',[10 10 560 420]);
SIG_OPTIONS=optimset('Display','off','TolFun',1e-05,'TolX',1e-05);
%folder='Y:\LEAP\Data_2013_FebRun\';

%savefolder=['V:\ARDB\E163\Data\',datestr(now,'yymmdd'),'\'];%Fluence1\'];
%if exist(savefolder,'dir')==0
%   mkdir(savefolder);
%end

col='b';         col2=[0 0 .7];
col3=[.5 .5 1];  col4=[0 .6 .4];
col5='c';

% col='r';         col2=[.7 0 0];
% col3=[1 .5 .5];  col4=[1 .8 0];
% col5='m';

%%  ----  Gradient vs Polarization  ------------

%load 'V:\ARDB\E163\Data\130510\Polarization\PolScan1.mat'  %weighted, HWHM fitted, 
%load 'V:\ARDB\E163\Data\130512\Polarization2\PolScan2.mat'  %weighted, HWHM fitted, 
load 'V:\ARDB\E163\Data\130513\Polarization1\PolScan2.mat'  %weighted, HWHM fitted, 
%load 'C:\Users\electron\Documents\AARD\LEAP Data\EdgarsData\LEAP Analysis\PolScan2.mat'  %weighted, HWHM fitted, 

%%
convF=1.2;
RUN=Scan.RUN;

 %Problem with 2284,   %Also 2283,2285,2276,2298,2297 had no signal and way below noise
%Ndel=[11]; %2284 only
indExc=RUN==2284;
% Low signal to noise
%indExc=(RUN==2284|RUN==2276|RUN==2282|RUN==2283|RUN==2285|RUN==2297|RUN==2298); 
%indExc=(RUN==2284|RUN==2276|RUN==2282|RUN==2285|RUN==2296|RUN==2297|RUN==2298); 

%indExc=(RUN==2296|RUN==2297|RUN==2298); 

% %Too many pts filtered
% indExc=(indExc|RUN==2278|RUN==2287|RUN==2289|RUN==2291|RUN==2294|RUN==2295); 

% %Ndel=[11];                     
Ndel=find(indExc);
Scan.RUN(Ndel)=[];   
Scan.Meas(Ndel,:,:)=[];
Scan.Meas_er(Ndel,:,:)=[];
Scan.ChiSqrd(Ndel,:,:)=[];
Scan.MeanE(Ndel)=[];
Scan.Jitter(Ndel)=[];
Scan.NoiseFloor(Ndel)=[];
Scan.Angle(Ndel)=[];


%
RUN=Scan.RUN;
wH=Scan.Meas(:,8,3)*convF;
wH_er=Scan.Meas_er(:,8,3)*convF;
ChiSqrd=Scan.ChiSqrd(:,8,3);
Angle=Scan.Angle;
%NF=Scan.NoiseFloor;
NF=Scan.NoiseFloor(:,3)*convF;

Penergy=Scan.Meas(:,9,3);
Pe_std=Scan.Meas_er(:,9,3);
PE(:,1) =Scan.Meas(:,9,1);
PE(:,2) =Scan.Meas(:,9,2);
PE(:,3) =Scan.Meas(:,9,3);
PE_er(:,1) =Scan.Meas_er(:,9,1);
PE_er(:,2) =Scan.Meas_er(:,9,2);
PE_er(:,3) =Scan.Meas_er(:,9,3);

Eshift(:,1)=Scan.Meas(:,8,1)*convF;
Eshift(:,2)=Scan.Meas(:,8,2)*convF;
Eshift(:,3)=Scan.Meas(:,8,3)*convF;
Eshift_er(:,1)=Scan.Meas_er(:,8,1)*convF;
Eshift_er(:,2)=Scan.Meas_er(:,8,2)*convF;
Eshift_er(:,3)=Scan.Meas_er(:,8,3)*convF;

PE= (-36.39*PE.^2 + 419.95*PE - 16.842*ones(size(PE)));

Pe_std= Pe_std.*(-36.39*2*Penergy + 419.95*ones(size(Penergy)));
Penergy= (-36.39*Penergy.^2 + 419.95*Penergy - 16.842*ones(size(Penergy)));

%{
for i=1:length(RUN)
   %if Scan.Meas(i,5,2)*convF>MaxE(i)
   if Scan.ChiSqrd(i,8,2)<ChiSqrd(i)
      fprintf('Single pass is better for run %d \n', RUN(i))
      wH(i)=Scan.Meas(i,8,2)*convF;
      wH_er(i)=Scan.Meas_er(i,8,2)*convF;
   end
   %if Scan.Meas(i,5,1)*convF>MaxE(i)
   if Scan.ChiSqrd(i,8,1)<ChiSqrd(i)
      fprintf('Raw is better for run %d \n', RUN(i))
      wH(i)=Scan.Meas(i,8,1)*convF;
      wH_er(i)=Scan.Meas_er(i,8,1)*convF;
   end
end
%}

% figure (16)
% subplot(3,1,1)
% plot(RUN,Angle,'o')
% axis tight
% ylabel('Angle (deg)')
% grid on
% enhance_plot;
% 
% subplot(3,1,2)
% %errorbar(RUN,Penergy,Pe_std,'o')
% % RUN3(:,1)=RUN;
% % RUN3(:,2)=RUN;
% % RUN3(:,3)=RUN;
% % errorbar(RUN3,PE,PE_er,'o')
% plot(RUN,PE,'o')
% axis tight
% %ylim([65 100])
% ylabel('Pulse Energy (uJ)')
% grid on
% enhance_plot;
% 
% subplot(3,1,3)
% plot(RUN,Scan.NoiseFloor, 'o')
% %xlabel('Incidence Angle \theta (deg)')
% ylabel('Noise Floor (keV)')
% xlabel('RUN #')
% axis tight
% grid on
% enhance_plot;
% 
% figure (17)
% plot(Angle,Scan.NoiseFloor(ind),'o')
% axis tight

indBreak=find(RUN>2286);
Nbreak=indBreak(1)-1;
%Nbreak=11;  % Original
%Nbreak=8; %After taking out 2284

Angle1=Angle(1:Nbreak);
wH1=wH(1:Nbreak);
wH1_er=wH_er(1:Nbreak);
RUN1=RUN(1:Nbreak);

Angle2=Angle(Nbreak+1:end);
wH2=wH(Nbreak+1:end);
wH2_er=wH_er(Nbreak+1:end);
RUN2=RUN(Nbreak+1:end);

[Angle, ind]=sort(Angle);
wH=wH(ind);
wH_er=wH_er(ind);
RUN=RUN(ind);

[Angle1, ind]=sort(Angle1);
wH1=wH1(ind);
wH1_er=wH1_er(ind);
RUN1=RUN1(ind);

[Angle2, ind]=sort(Angle2);
wH2=wH2(ind);
wH2_er=wH2_er(ind);
RUN2=RUN2(ind);
%
% 
% figure (2)
% subplot(3,1,1)
% errorbar(Angle1,wH1,wH1_er,'o');%,'o','LineWidth',2,'Color',col4)
% hold on
% errorbar(Angle2,wH2,wH2_er,'x');%,'o','LineWidth',2,'Color',col4)
% 
% ylabel('Energy Shift [keV]')
% axis tight
% ylim([0 30])
% grid on
% 
% enhance_plot;
% hold on
% 
% for i=1:length(Angle)
%    text(Angle(i),wH(i)+2,num2str(RUN(i),'%g'),'Rotation',-25)
%    hold on
% end
% 
% subplot(3,1,2)
% %errorbar(Angle,Penergy,Pe_std,'o')
% plot(Angle,PE,'o')
% axis tight
% grid on
% %xlabel('Polarization \phi (deg)')
% ylabel('Pulse Energy (uJ)')
% axis tight
% %ylim([0 30])
% hold on
% enhance_plot;
% for i=1:length(Angle)
%    text(Angle(i),Penergy(i)+2,num2str(RUN(i),'%g'),'Rotation',-25)
%    hold on
% end
% 
% subplot(3,1,3)
% plot(Angle,Scan.NoiseFloor, 'o')
% %xlabel('Incidence Angle \theta (deg)')
% ylabel('Noise Floor (keV)')
% xlabel('Polarization \phi (deg)')
% axis tight
% grid on
% hold on
% enhance_plot;

% 
% MaxE=5.1;  % Noise floor level
% MaxE_er=1.5;

Grad=50*sqrt((.0484*wH+1.0463).^2-1)-23;
Grad_er=50*.0484*(.0484*wH+1.0463)./sqrt((.0484*wH+1.0463).^2-1).*wH_er;

load V:\ARDB\E163\Data\130809\calibration_values
Grad = interp1(xcal,ycal,wH,'pchip');

var=Grad;
var_er=Grad_er;

Grad1=50*sqrt((.0484*wH1+1.0463).^2-1)-23;
Grad1_er=50*.0484*(.0484*wH1+1.0463)./sqrt((.0484*wH1+1.0463).^2-1).*wH1_er;
Grad2=50*sqrt((.0484*wH2+1.0463).^2-1)-23;
Grad2_er=50*.0484*(.0484*wH2+1.0463)./sqrt((.0484*wH2+1.0463).^2-1).*wH2_er;
Grad1 = interp1(xcal,ycal,wH1,'pchip');
Grad2 = interp1(xcal,ycal,wH2,'pchip');

% scale the 2nd data set to match the first
Grad = [Grad2*1.15; Grad1];
Grad_er = [Grad2_er*1.15; Grad1_er];

flatW=@(c,x) c(1)./Pe_std;
[Pavg,~,~,~,~,~,~]=lsqcurvefitstd(flatW,mean(Penergy),Angle,...
       Penergy./Pe_std,min(Penergy),max(Penergy),SIG_OPTIONS);
         
Grad=Grad.*sqrt(Penergy./Pavg);
Grad_er=Grad_er.*sqrt(Penergy./Pavg);

indExc=(RUN==2297|RUN==2298|RUN==2282|RUN==2285|RUN==2276); 

Grad(indExc)=[];
Grad_er(indExc)=[];
Angle(indExc)=[];
% %Ndel=[11];                     
Ndel=find(indExc);
Scan.RUN(Ndel)=[];   
Scan.Meas(Ndel,:,:)=[];
Scan.Meas_er(Ndel,:,:)=[];
Scan.ChiSqrd(Ndel,:,:)=[];
Scan.MeanE(Ndel)=[];
Scan.Jitter(Ndel)=[];
Scan.NoiseFloor(Ndel)=[];
Scan.Angle(Ndel)=[];

hFig = figure(3);
xwidth=8.9;
ywidth=6.5;
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','centimeters')
set(hFig, 'Position', [0 0 xwidth ywidth])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

hold on

% *********  Fits   *********
shift=15;
Angle=Angle+15*ones(length(Angle),1);
z=-180+shift:180+shift;

%{
%Fit used in submitted figure
cosF =@(c,x) c(1)*ones(size(x))+abs(c(2)*cos(pi*(x-c(3))./180));
yW=Grad_er.^1;
cosFW =@(c,x) (c(1)*ones(size(x))+abs(c(2)*cos(pi*(x-c(3))./180)))./yW;
[cout1,ChiSqX,~,~,~,~,~,cout1_std]=lsqcurvefitstd(cosFW,[30 30 0],Angle,...
 Grad./yW,[-50 1 -45],[50 150 45],SIG_OPTIONS);
cout1
cout1_std
%}

% new fit with no offset
cosF =@(c,x) c(1)*ones(size(x))+abs(c(2)*cos(pi*(x-c(3))./180));
cosF1 =@(c,x) 0*ones(size(x))+abs(c(1)*cos(pi*(x-c(2))./180));
yW=Grad_er.^1;
cosFW =@(c,x) (0*ones(size(x))+abs(c(1)*cos(pi*(x-c(2))./180)))./yW;
[cout1,ChiSqX,~,~,~,~,~,cout1_std]=lsqcurvefitstd(cosFW,[30 0],Angle,...
 Grad./yW,[1 -45],[150 45],SIG_OPTIONS);
cout1
cout1_std


% cosF =@(c,x) abs(c(1)*cos(pi*(x-c(2))./180));
% yW=Grad_er.^1;
% cosFW =@(c,x) abs(c(1)*cos(pi*(x-c(2))./180))./yW;
% [cout1,ChiSqX,~,~,~,~,~,cout1_std]=lsqcurvefitstd(cosFW,[30 0],Angle,...
%  Grad./yW,[1 -45],[150 45],SIG_OPTIONS);
% cout1
% cout1_std

yW=Grad1_er.^1;
cosFW =@(c,x) (c(1)*ones(size(x))+abs(c(2)*cos(pi*(x-c(3))./180)))./yW;
[cout2,ChiSqX,~,~,~,~,~,cout2_std]=lsqcurvefitstd(cosFW,[30 30 0],Angle1,...
 Grad1./yW,[-50 1 -45],[50 150 45],SIG_OPTIONS);
%cout2
yW=Grad2_er.^1;
cosFW =@(c,x) (c(1)*ones(size(x))+abs(c(2)*cos(pi*(x-c(3))./180)))./yW;
[cout3,ChiSqX,~,~,~,~,~,cout3_std]=lsqcurvefitstd(cosFW,[30 30 0],Angle2,...
 Grad2./yW,[-50 1 -45],[50 150 45],SIG_OPTIONS);
%cout3

% yW=Grad1_er.^1;
% cosFW =@(c,x) abs(c(1)*cos(pi*(x-c(2))./180))./yW;
% [cout2,ChiSqX,~,~,~,~,~,cout2_std]=lsqcurvefitstd(cosFW,[30 0],Angle1,...
%  Grad1./yW,[1 -45],[150 45],SIG_OPTIONS);
% %cout2
% yW=Grad2_er.^1;
% cosFW =@(c,x) abs(c(1)*cos(pi*(x-c(2))./180))./yW;
% [cout3,ChiSqX,~,~,~,~,~,cout3_std]=lsqcurvefitstd(cosFW,[30 0],Angle2,...
%  Grad2./yW,[1 -45],[150 45],SIG_OPTIONS);
% %cout3

z1=z(end/2:end);
z2=z(1:end/2+1);

%hold on

% hold on
% h1f=plot(z1,cosF(cout2,z1),'b')
% hold on
% plot(z2,cosF(cout2,z2),':b')
% 
% hold on
% h2f=plot(z2,cosF(cout3,z2),'r')
% hold on
% plot(z1,cosF(cout3,z1),':r')

%Noise floors
NF=mean(NF)+std(NF);
%NF=50*sqrt((.0484*NF+1.0463).^2-1)-23
NF = interp1(xcal,ycal,NF,'pchip');

% NF1=mean(NF(1:12))+std(NF(1:12));
% NF1=50*sqrt((.0484*NF1+1.0463).^2-1)-23;
% NF2=mean(NF(13:end))+std(NF(13:end));
% NF2=50*sqrt((.0484*NF2+1.0463).^2-1)-23;
% 
hold on
hN=plot(z,NF*ones(size(z)),'-.k','LineWidth',1)%

% hold on
% plot(z1,NF1*ones(size(z1)),'b-.')
% hold on
% plot(z2,NF2*ones(size(z2)),'r-.')

dGreen=[0 .6 .4];
hold on
h0f=plot(z,cosF1(cout1,z),'Color',dGreen,'LineWidth',1)
hold on
h0=errorbar(Angle,Grad,Grad_er,'ob')


%h1=errorbar(Efield*1e-9,var,dEfield*1e-9,var_er,'o');
%set(h0,'LineWidth',2,'MarkerFaceColor',[0 .7 1],'MarkerSize',6);
set(h0,'LineWidth',1,'MarkerFaceColor','c','MarkerSize',4);


% %Used in submitted plot
% xlim([-68 68])
% ylim([-10 80])
xlim([-90 90])
ylim([0 80])

% Etick=[0,.5,1,1.5,2,2.5,3,3.5]
set(gca,'XTick',[-90,-45,0,45,90])
%set(h(2),'Color','b'), set(h(3),'Color','b')
%set(h1,'Color',col)
hx=xlabel('Polarization Angle \phi (deg)')
hy=ylabel('Gradient (MeV/m)')
%legend([h1 h1f h2 h2f],{'Data1','Fit1','Data2','Fit2'})
box on

hL=legend([h0 h0f hN],{'Data','Fit','Noise Level'})
set(hL, 'FontSize', 6 ,'Location','NorthEastOutside');
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
 
 export_fig EDFig5aL.eps -painters -rgb
 

% for i=1:length(Angle)
%    text(Angle(i)+3,Grad(i)+1,num2str(RUN(i),'%g'),'Rotation',-25)
%    hold on
% end

%% ----  Gradient vs TIP  ------------
% % 
% %load 'V:\ARDB\E163\Data\130514\TIP2\TipScan1.mat'  %weighted, HWHM fitted, 
% load 'V:\ARDB\E163\Data\130516\TIP\TipScan1.mat'  %weighted, HWHM fitted, 
% load 'V:\ARDB\E163\Data\130516\TIP2\TipScan1.mat'  %weighted, HWHM fitted, 
% 
% convF=1.2;
% RUN=Scan.RUN;
% MaxE=Scan.Meas(:,5,3)*convF;
% MaxE_er=Scan.Meas_er(:,5,3)*convF;
% ChiSqrd=Scan.ChiSqrd(:,5,3);
% TIP=[9335:10:9425 9425 9435 9435:10:9465 9465 9465:10:9485 9485]; 
% Angle=(TIP-9325)*0.0005955;
% Angle=Angle';
% NF=Scan.NoiseFloor*convF;
% Penergy=Scan.Meas(:,9,3);
% Pe_std=Scan.Meas_er(:,9,3);
% 
% NF=Scan.NoiseFloor(:,3)*convF;
% 
% for i=1:length(RUN)
%    %if Scan.Meas(i,5,2)*convF>MaxE(i)
% 
%    if Scan.ChiSqrd(i,5,2)<ChiSqrd(i)
%       fprintf('Single pass is better for run %d \n', RUN(i))
%       MaxE(i)=Scan.Meas(i,5,2)*convF;
%       MaxE_er(i)=Scan.Meas_er(i,5,2)*convF;
%       Penergy(i)=Scan.Meas(i,9,2)*convF;
%       Pe_std(i)=Scan.Meas_er(i,9,2)*convF;
%       NF(i)=Scan.NoiseFloor(i,2)*convF;
%    end
%    %if Scan.Meas(i,5,1)*convF>MaxE(i)
%    if Scan.ChiSqrd(i,5,1)<ChiSqrd(i)
%       fprintf('Raw is better for run %d \n', RUN(i))
%       MaxE(i)=Scan.Meas(i,5,1)*convF;
%       MaxE_er(i)=Scan.Meas_er(i,5,1)*convF;
%       Penergy(i)=Scan.Meas(i,9,1)*convF;
%       Pe_std(i)=Scan.Meas_er(i,9,1)*convF;
%       NF(i)=Scan.NoiseFloor(i,1)*convF;
%    end
% end
% 
% Pe_std= Pe_std.*(-36.39*2*Penergy + 419.95*ones(size(Penergy)));
% Penergy= (-36.39*Penergy.^2 + 419.95*Penergy - 16.842*ones(size(Penergy)));
% 
% Ndel=[10];
% MaxE(Ndel)=[];   %Worst initial Noise floorclose all
% MaxE_er(Ndel)=[];
% RUN(Ndel)=[];
% Penergy(Ndel)=[];
% Pe_std(Ndel)=[];
% Angle(Ndel)=[];
% NF(Ndel)=[];
% ChiSqrd(Ndel)=[];
% 
% [Angle, ind]=sort(Angle);
% MaxE=MaxE(ind);
% MaxE_er=MaxE_er(ind);
% RUN=RUN(ind);
% Penergy=Penergy(ind);
% Pe_std=Pe_std(ind);
% NF=NF(ind);
% 
% figure (16)
% subplot(2,1,1)
% plot(RUN,Angle,'o')
% axis tight
% grid on
% subplot(2,1,2)
% errorbar(RUN,Penergy,Pe_std, 'o')
% axis tight
% grid on
% 
% % figure (17)
% % plot(Angle,Scan.NoiseFloor(ind),'o')
% % axis tight
% %
% 
% figure (2)
% cla
% subplot(4,1,1)
% errorbar(Angle,MaxE,MaxE_er,'o');%,'o','LineWidth',2,'Color',col4)
% %plot(Angle,MaxE,'o')
% %xlabel('Incidence Angle \theta (deg)')
% ylabel('Energy Shift [keV]')
% axis tight
% %ylim([0 30])
% grid on
% hold on
% enhance_plot;
% for i=1:length(Angle)
%    text(Angle(i),MaxE(i)+2,num2str(RUN(i),'%g'),'Rotation',-25)
%    hold on
% end
% 
% 
% subplot(4,1,2)
% errorbar(Angle,Penergy,Pe_std, 'o')
% %xlabel('Incidence Angle \theta (deg)')
% ylabel('Laser Pulse Energy (uJ)')
% axis tight
% grid on
% hold on
% 
% enhance_plot;
% for i=1:length(Angle)
%    text(Angle(i),Penergy(i)+.001,num2str(RUN(i),'%g'),'Rotation',-25)
%    hold on
% end
% 
% subplot(4,1,3)
% plot(Angle,NF, 'o')
% %xlabel('Incidence Angle \theta (deg)')
% ylabel('Noise Floor (keVJ)')
% axis tight
% grid on
% hold on
% 
% enhance_plot;
% for i=1:length(Angle)
%    text(Angle(i),NF(i)+.1,num2str(RUN(i),'%g'),'Rotation',-25)
%    hold on
% end
% 
% subplot(4,1,4)
% plot(Angle,ChiSqrd, 'o')
% %xlabel('Incidence Angle \theta (deg)')
% ylabel('Flat line Fit Chi^2')
% axis tight
% grid on
% hold on
% 
% enhance_plot;
% for i=1:length(Angle)
%    text(Angle(i),ChiSqrd(i)+1,num2str(RUN(i),'%g'),'Rotation',-25)
%    hold on
% end
% 
% %legend('Vac. Dist. Lead.','Glass Dist. Lead.','Glass Dist. Trail.',...
% %    'Location','NorthWest')
% 
% %
% % MaxE=5.1;  % Noise floor level
% % MaxE_er=1.5;
% Grad=50*sqrt((.0484*MaxE+1.0463).^2-1)-23;
% Grad_er=50*.0484*(.0484*MaxE+1.0463)./sqrt((.0484*MaxE+1.0463).^2-1).*MaxE_er;
% 
% flatW=@(c,x) c(1)./Pe_std;
% [Pavg,~,~,~,~,~,~]=lsqcurvefitstd(flatW,mean(Penergy),Angle,...
%        Penergy./Pe_std,min(Penergy),max(Penergy),SIG_OPTIONS);
%          
% Grad=Grad.*sqrt(Penergy./Pavg);
% Grad_er=Grad_er.*sqrt(Penergy./Pavg);
% 
% figure (3)
% hold on
% % errorbar(Angle,Grad2,Grad2_er,'or');%,'Color',col);
% % hold on
% 
% % for i=1:length(Angle)
% %    text(Angle(i),Grad(i)+2,num2str(RUN(i),'%g'),'Rotation',-25)
% %    hold on
% % end
% %grid on
% %axis tight
% %ylim([-10 80])
% % 
% 
% % x=0:.00001:.001;
% % N=1000;
% % a=1./(1+sin(x));
% % %G=a.*sin(2*pi*N./a)./N/pi./(1-a.^2);
% % G2=(1+sin(x))./(2+sin(x))./sin(x).*sin(2*pi*N*sin(x))./pi/N;
% 
% % % Fit to data
% theta=-.01:.0001:.1;
% yW=Grad_er.^1;
% 
% % %Does not assume angle unit is correct
% % GradF=@(c,x) c(1)*(1+sin((x-c(2))./c(3)))./(2+sin((x-c(2))./c(3)))./sin((x-c(2))./c(3)).*sin(2*pi*c(4)*sin((x-c(2))./c(3)))./2/pi/c(4);
% % GradFW=@(c,x) (c(1)*(1+sin((x-c(2))./c(3)))./(2+sin((x-c(2))./c(3)))./sin((x-c(2))./c(3)).*sin(2*pi*c(4)*sin((x-c(2))./c(3)))./2/pi/c(4))./yW;
% % [cout1,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(GradFW,[100 0 1 10],Angle,Grad./yW,...
% %   [1 -.1 .1 1],[300 .2 20 1000],SIG_OPTIONS);
% 
% % Assumes angle unit is correct
% GradF=@(c,x) c(1)*(1+sin(pi*(x-c(2))./180))./(2+sin(pi*(x-c(2))./180))./sin(pi*(x-c(2))./180).*sin(2*pi*c(3)*sin(pi*(x-c(2))./180))./2/pi/c(3);
% GradFW=@(c,x) (c(1)*(1+sin(pi*(x-c(2))./180))./(2+sin(pi*(x-c(2))./180))./sin(pi*(x-c(2))./180).*sin(2*pi*c(3)*sin(pi*(x-c(2))./180))./2/pi/c(3))./yW;
% [cout1,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(GradFW,[100 0 10],Angle,Grad./yW,...
%   [1 -.1 1],[300 .2 1000],SIG_OPTIONS);
% 
% 
% cout1
% cout_std
% %figure (19)
% %plot(theta,GradF([100 0 10],theta))
% % hold on
% % h2=plot(theta,GradF(cout1,theta));
%  
% %Noise floors
% NFm=mean(NF);
% NFs=std(NF);
% NF=NFm+NFs;
% NF=50*sqrt((.0484*NF+1.0463).^2-1)-23;
% 
% hold on
% plot(theta,NF*ones(size(theta)),'k-.')
% 
% hold on
% h2=plot(theta,GradF(cout1,theta));
% 
% h1=errorbar(Angle,Grad,Grad_er,'ob');%,'Color',col);
% hold on
% 
% xlim([0 .1])
% ylim([0 140])
% 
% xlabel('Incidence \theta (deg)')
% ylabel('Gradient (MeV/m)')
% legend([h1(1) h2],{'Data','Fit'})
% box on
% enhance_plot;
% %grid on

%%

%linFW =@(c,x) c(1).*sqrt(x-c(2))./yW;
%[cout1,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(linFW,[1 0],Penergy,...
%   var./yW,[0.001 -10],[1000 10],SIG_OPTIONS);
% linF =@(c,x) c(1).*sqrt(x);
% linFW =@(c,x) c(1).*sqrt(x)./yW;
% [cout1,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(linFW,1,Penergy,...
%     var./yW,0.001,1000,SIG_OPTIONS);


% %% Fit to data
% yW=var_er.^1;
% linF =@(c,x) c(1).*sqrt(x-c(2));
% linFW =@(c,x) c(1).*sqrt(x-c(2))./yW;
% [cout1,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(linFW,[1 0],Penergy,...
%    var./yW,[0.001 -10],[1000 10],SIG_OPTIONS);
% % linF =@(c,x) c(1).*sqrt(x);
% % linFW =@(c,x) c(1).*sqrt(x)./yW;
% % [cout1,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(linFW,1,Penergy,...
% %     var./yW,0.001,1000,SIG_OPTIONS);
% 
% hold on
% 
% PenergyF=0:.05:max(Penergy)+.25;
% gradF=linF(cout1,PenergyF);
% h2=plot(PenergyF,gradF,'-.','Color',col4)
% box on
% enhance_plot;
% legend([h1(1) h2],{'data', 'theory'})

% %%
% figure
% n=1.5;
% eps0=8.85e-12;
% tau=1.3e-12;
% c=3e8;
% w1=60e-6/2;
% w2=600e-6/2;
% Energy=(1:350)*1e-6;
% 
% 
% F=2*Penergy./pi/w1/w2*1e-6;
% Efield=sqrt(2*F./eps0/c/tau);
% dEfield=Pe_std*1e-6./sqrt(Penergy*1e-6)/sqrt(pi*eps0*c*tau*w1*w2);
% 
% yW=var_er.^1;
% linF =@(c,x) c(1)*ones(size(x))+c(2)*x;
% linFW =@(c,x) (c(1)*ones(size(x))+c(2)*x)./yW;
% [cout1,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(linFW,[0 1e-7],Efield,...
%    var./yW,[-100 0],[100 1],SIG_OPTIONS);
% EfieldF=0:max(Efield)/200:max(Efield)+max(Efield)/100;
% gradF=linF(cout1,EfieldF);
% 
% figure (4)
% h1=ploterr(Efield*1e-9,var,dEfield*1e-9,var_er,'.');
% hold on
% plot(EfieldF*1e-9,gradF,'r:')
% hold on
% plot(EfieldF*1e-9,(18+5.7)*ones(size(EfieldF)),'--k')%,
% %h1=ploterr(Efield,MaxE,dEfield,MaxE_er,'.');%,
% xlabel('Peak Longitudinal Electric Field E_z (GV/m)')
% ylabel('Gradient (MeV/m)')
% axis tight
% 
% 
% 
% % errorbar(fluence1,EmodH,EmodH_er,'ob','LineWidth',2);
% % yW=EmodH_er;
% % linFW =@(c,x) c(1).*sqrt(x-c(2))./yW;
% % [cout2,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(linFW,[1 0],fluence1',...
% %     EmodH./yW,[0.01 -1],[max(EmodH) 1],SIG_OPTIONS);
% % hold on
% %
% % errorbar(fluence,EmodL,EmodL_er,'or','LineWidth',2)
% % yW=EmodL_er;
% % linFW =@(c,x) c(1).*sqrt(x-c(2))./yW;
% % [cout3,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(linFW,[1 0],fluence',...
% %     EmodL./yW,[0.01 -1],[max(EmodL) 1],SIG_OPTIONS);
% % hold on
% %
% % fluenceF=0:.05:max(fluence)+.25;
% % gradF=linF(cout1,fluenceF);
% % plot(fluenceF,gradF,'Color',col4)
% % hold on
% %
% % gradF=linF(cout2,fluenceF);
% % plot(fluenceF,gradF,'b')
% % hold on
% %
% % gradF=linF(cout3,fluenceF);
% % plot(fluenceF,gradF,'r')
% %
% % xlabel('Laser Fluence [J/cm^2]')
% % ylabel('Energy Shift [keV]')
% % axis tight
% % ylim([0 120])
% % grid on
% %
% % %legend('Vac. Dist. Lead.','Glass Dist. Lead.','Glass Dist. Trail.',...
% % %    'Location','NorthWest')
% % enhance_plot;
% %%
% 
% P_on_m=Scan.LaserP_on_m(ind);
% P_on_std=Scan.LaserP_on_std(ind);
% 
% Pow= (-36.39*P_on_m.^2 + 419.95*P_on_m - 16.842*ones(size(P_on_m)))./.8;
% Pow_std= P_on_std.*(-36.39*2*P_on_std + 419.95*ones(size(P_on_m)))./.8;
% 
% figure
% subplot(2,1,1)
% stem(Energy0,Pow,'o','LineWidth',2,'Color',col4);
% xlabel('E_{expected} [uJ]')
% ylabel('E_{meas} [uJ] ')
% grid on
% for i=1:length(Energy0)
%    text(Energy0(i),Pow(i)+2,num2str(RUN(i),'%g'),'Rotation',90)
%    hold on
% end
% xlim([-5 Energy0(end)+5])
% set(gca,'YTick',0:50:500)
% 
% hold on
% plot(Energy0,Energy0)
% % xlim([0 50])
% % ylim([0 50])
% % set(gca,'YTick',0:10:50)
% 
% % xlim([400 430])
% % ylim([400 435])
% % set(gca,'YTick',400:5:435)
% 
% enhance_plot;
% 
% subplot(2,1,2)
% stem(Energy0,Pow_std,'o','LineWidth',2,'Color',col4);
% hold on
% ylim([0 9])
% xlim([-5 Energy0(end)+5])
% xlabel('E_{expected} [uJ]')
% ylabel('\delta E_{meas} [uJ]')
% grid on
% enhance_plot;
% 
% %% Energy gain ^2 vs pulse energy
% col4=[0 .6 .4];
% figure (2)
% cla
% 
% % Energy0(9)=[];
% % MaxE(9)=[];
% % MaxE_er(9)=[];
% 
% Energy=Energy0;
% var=wH.^2;
% var_er=2*wH.*wH_er;
% 
% errorbar(Energy,var,var_er,'o','LineWidth',2)%,'Color',col4)
% 
% Energy=Energy0(1:end-3);
% var=wH(1:end-3).^2;
% var_er=2*wH(1:end-3).*wH_er(1:end-3);
% 
% 
% yW=var_er;
% linF =@(c,z) c(1)*ones(size(z))+ c(2)*z;
% linFW =@(c,x) (c(1)*ones(size(x))+ c(2).*x)./yW;
% [cout1,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(linFW,[1 1],Energy',...
%    var./yW,[-20 0],[20 100000],SIG_OPTIONS);
% hold on
% 
% EnergyF=0:.05:max(Energy)+5;
% gradF=linF(cout1,EnergyF);
% plot(EnergyF,gradF,'Color',col4,'LineWidth',2)
% hold on
% 
% xlabel('Laser Pulse Energy [\mu J]')
% ylabel('Gradient [MeV/m]')
% ylabel('(Energy Difference [keV] ) ^2')
% axis tight
% ylim([0 14400])
% grid on
% 
% val=[45, 65,85,105,120];
% gradV=50*sqrt((.0484*val+1.0463).^2-1)-23;
% set(gca,'YTick',[0,45^2,65^2,85^2,105^2,120^2])
% %set(gca,'YTick',[0,35^2,65^2,90^2,110^2,125^2])
% set(gca,'YTickLabel',{'0','45^2','65^2','85^2','105^2','120^2'})
% 
% %legend('Vac. Dist. Lead.','Glass Dist. Lead.','Glass Dist. Trail.',...
% %    'Location','NorthWest')
% enhance_plot;
% 
% 
% % figure (12)
% %
% % xlabel('45^2    65^2  85^2    105^2   120^2')
% %
% % enhance_plot;
% %
% % figure(13)
% %
% % xlabel('130^2  180^2  230^2  279^2  316^2')
% % ylabel('Gradient [MeV/m]')
% %
% %   enhance_plot;
% 
% 
% %%
% %% Gradient^2 vs pulse energy
% col4=[0 .6 .4];
% figure (3)
% cla
% 
% % Energy0(9)=[];
% % MaxE(9)=[];
% % MaxE_er(9)=[];
% 
% Energy=Energy0;
% 
% Grad=50*sqrt((.0484*wH+1.0463).^2-1)-23;
% Grad_er=50*.0484*(.0484*wH+1.0463)./sqrt((.0484*wH+1.0463).^2-1).*wH_er;
% var=Grad;
% var_er=Grad_er;
% %
% var=Grad.^2;
% var_er=2*Grad.*Grad_er;
% 
% errorbar(Energy,var,var_er,'o','LineWidth',2)%,'Color',col4)
% %%
% % Energy=Energy0(1:end-3);
% % var=MaxE(1:end-3).^2;
% % var_er=2*MaxE(1:end-3).*MaxE_er(1:end-3);
% 
% yW=var_er.^0;
% linF =@(c,z) c(1)*ones(size(z))+ c(2)*z;
% linFW =@(c,x) (c(1)*ones(size(x))+ c(2).*x)./yW;
% [cout1,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(linFW,[1 1000],Energy',...
%    var./yW,[-200 0],[200 1000000],SIG_OPTIONS);
% hold on
% 
% EnergyF=0:.05:max(Energy)+5;
% gradF=linF(cout1,EnergyF);
% plot(EnergyF,gradF,'Color',col4,'LineWidth',2)
% hold on
% 
% xlabel('Laser Pulse Energy [\mu J]')
% ylabel('(Gradient [MeV/m] ) ^2')
% axis tight
% ylim([0 124000])
% grid on
% 
% %val=[45, 65,85,105,120];
% set(gca,'YTick',[0,50^2,100^2,150^2,200^2,250^2,300^2])
% 
% set(gca,'YTickLabel',{'0','50^2','100^2','150^2','200^2','250^2','300^2'})
% 
% %legend('Vac. Dist. Lead.','Glass Dist. Lead.','Glass Dist. Trail.',...
% %    'Location','NorthWest')
% enhance_plot;
% 
% 
% 
%%

% hold on
% errorbar(RUN1,EmodH,EmodH_er,'ob','LineWidth',2);
% hold on
%
% errorbar(RUN,EmodL,EmodL_er,'or','LineWidth',2)
%
% xlabel('Run Number ylim')
% ylabel('Gradient [MeV/m]')
% axis tight
% ylim([0 250])
% grid on
%
% %legend('Vac. Dist. Lead.','Glass Dist. Lead.','Glass Dist. Trail.',...
% %    'Location','NorthWest')
% enhance_plot;
% %

% %% Temporal Overlap
%
% convF=1;
% MaxE=Scan.Other(:,5);
% ind=(MaxE~=0);
% MaxE=Scan.Other(ind,5)*convF;
% EmodH=Scan.Other(ind,3)*convF;
% EmodL=Scan.Other(ind,4)*convF;
%
% atten=Scan.atten(ind);
% fluence=1.558*sin(0.066*atten-1.703)+1.482;
%
% [fluence ind]=sort(fluence);
% MaxE=MaxE(ind);
% EmodH=EmodH(ind);
% EmodL=EmodL(ind);
% col4=[0 .6 .4];
% figure (2)
% cla
% plot(fluence,MaxE,'o-','LineWidth',2,'Color',col4)
% hold on
% plot(fluence,EmodH,'ob-','LineWidth',2);
% hold on
%
% plot(fluence,EmodL,'or-','LineWidth',2)
% xlabel('Laser Fluence [J/cm^2]')
% ylabel('Optimal Temporal Overlap [ps]')
% axis tight
% %ylim([0 250])
% grid on
%
% legend('Vac. Dist. Lead.','Glass Dist. Lead.','Glass Dist. Trail.',...
%     'Location','NorthWest')
% enhance_plot;
%
% %%
% %% Energy Jitter
%
% convF=2.4;
% Epeak=Scan.Meas(:,2);
% %ind=(Epeak~=0);
% ind=1:length(Epeak);
% Epeak=Scan.Meas(ind,2)*convF;
% Epeak_er=Scan.Meas_er(ind,3)*convF;
% Epeak_off=Scan.Other(ind,2)*convF;
% atten=Scan.atten(ind);
% fluence=1.558*sin(0.066*atten-1.703)+1.482;
% [fluence ind]=sort(fluence);
% Epeak=Epeak(ind);
% Epeak_er=Epeak_er(ind);
% Epeak_off=Epeak_off(ind);
% col4=[0 .6 .4];
% figure (3)
% cla
% errorbar(fluence,Epeak,Epeak_er,'o-r','LineWidth',2)
% hold on
% plot(fluence,Epeak_off,'x-b','LineWidth',2)
%
% xlabel('Laser Fluence [J/cm^2]')
% ylabel('Glass Peak Position [keV]')
% axis tight
% %ylim([-100 100])
% ylim([1320 1480])
% title('Grating Gap=400nm')
% legend('Laser On','Laser Off')
% grid on
%
% enhance_plot;
% %%
% RUN=Scan.RUN;
% Epeak_er=Scan.Meas_er(:,3);
%
% figure (4)
% stem(RUN,Epeak_er,'or','LineWidth',2)
% xlabel('Run number')
% ylabel('Main Peak Jitter [pix]')
% axis tight
% ylim([min(Epeak_er) 5])
% title('Grating Gap=400nm')
% grid on
% enhance_plot;
%

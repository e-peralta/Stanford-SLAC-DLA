%% Setup stuff
%set(0,'defaultfigureposition',[296 342 560 420]);
set(0,'defaultfigureposition',[10 10 560 420]);
SIG_OPTIONS=optimset('Display','off','TolFun',1e-05,'TolX',1e-05);
%folder='Y:\LEAP\Data_2013_FebRun\';
savefolder=['V:\ARDB\E163\Data\',datestr(now,'yymmdd'),'\Fluence1\'];
if exist(savefolder,'dir')==0
   mkdir(savefolder);
end



%% Setup stuff
%set(0,'defaultfigureposition',[296 342 560 420]);
set(0,'defaultfigureposition',[10 10 560 420]);
SIG_OPTIONS=optimset('Display','off','TolFun',1e-05,'TolX',1e-05);
%folder='Y:\LEAP\Data_2013_FebRun\';
savefolder=['V:\ARDB\E163\Data\',datestr(now,'yymmdd'),'\'];%Fluence1\'];
if exist(savefolder,'dir')==0
   mkdir(savefolder);
end

col='b';         col2=[0 0 .7];
col3=[.5 .5 1];  col4=[0 .6 .4];
col5='c';


%% ----  Energy Gain vs Fluence plot

%load 'V:\ARDB\E163\Data\130423\Fluence1\FluenceScan1.mat'  %weighted, HWHM fitted, 
%with ChiSqd values and manual removal of tweaked portions and manual adjustment of bad fits
 

load ~/Dropbox/Research/Data/FluenceScan1.mat  %weighted, HWHM fitted, 
RUN=Scan.RUN;

convF=1.2;

RUN=Scan.RUN;
MaxE=Scan.Meas(:,5,3)*convF;
MaxE_er=Scan.Meas_er(:,5,3)*convF;
ChiSqrd=Scan.ChiSqrd(:,5,3);
Penergy=Scan.Meas(:,8,3);
Pe_std=Scan.Meas_er(:,8,3);
% Penergy=Scan.Pon_mean;
% Pe_std=Scan.Pon_std;


for i=1:length(RUN)
    
   %if Scan.Meas(i,5,2)*convF>MaxE(i)
   if Scan.ChiSqrd(i,5,2)<ChiSqrd(i)
      fprintf('Single pass is better for run %d \n', RUN(i))
      MaxE(i)=Scan.Meas(i,5,2)*convF;
      MaxE_er(i)=Scan.Meas_er(i,5,2)*convF;
      ChiSqrd(i)=Scan.ChiSqrd(i,5,2);
      Penergy(i)=Scan.Meas(i,8,2);
      Pe_std(i)=Scan.Meas_er(i,8,2);
   end
   %if Scan.Meas(i,5,1)*convF>MaxE(i)
   if Scan.ChiSqrd(i,5,1)<ChiSqrd(i)
      fprintf('Raw is better for run %d \n', RUN(i))
      MaxE(i)=Scan.Meas(i,5,1)*convF;
      MaxE_er(i)=Scan.Meas_er(i,5,1)*convF;
      Penergy(i)=Scan.Meas(i,8,1);
      Pe_std(i)=Scan.Meas_er(i,8,1);
   end
end
indExc=(RUN==2205|RUN==2209);  %%These had too few points after tweak removal
MaxE(indExc)=[];
MaxE_er(indExc)=[];
RUN(indExc)=[];
Penergy(indExc)=[];
Pe_std(indExc)=[];
% 
% % indExc=(RUN==2228);  %%These had too few points after tweak removal
% % MaxE(indExc)=1.12*1.2;
% % MaxE_er(indExc)=6.57*1.2;
% 
indExc=(RUN==2228|RUN==2229);  %%These had low signal to noise
MaxE(indExc)=[];
MaxE_er(indExc)=[];
RUN(indExc)=[];
Penergy(indExc)=[];
Pe_std(indExc)=[];

% % % These had more than half the pts removed in filtering
% indExc=(RUN==2203|RUN==2207|RUN==2208|RUN==2210|RUN==2211|RUN==2212|RUN==2214|RUN==2215|RUN==2220|RUN==2221|RUN==2228|RUN==2230|RUN==2231|RUN==2235|RUN==2239|RUN==2240);  %%These are low lying - no real reason to exclude
% MaxE(indExc)=[];
% MaxE_er(indExc)=[];
% RUN(indExc)=[];
% Penergy(indExc)=[];
% Pe_std(indExc)=[];

% Manually fixed (see 130812 data folder)
ind=RUN==2227;
MaxE(ind)=11.42*1.2;
MaxE_er(ind)=1.37*1.2;

% % Justified as not adjusting beam quality during run
indExc=(RUN==2204);
MaxE(indExc)=[];
MaxE_er(indExc)=[];
RUN(indExc)=[];
Penergy(indExc)=[];
Pe_std(indExc)=[];


Pe_std= .93*Pe_std.*(-36.39*2*Penergy + 419.95*ones(size(Penergy)));
Penergy= .93*(-36.39*Penergy.^2 + 419.95*Penergy - 16.842*ones(size(Penergy)));
%93% accounts for the losses between vacuum window and IP

indX1=RUN<2225;
Penergy1=Penergy(indX1);
Pe_std1=Pe_std(indX1);
MaxE1=MaxE(indX1);
MaxE1_er=MaxE_er(indX1);
RUN1=RUN(indX1);

Penergy2=Penergy(~indX1);
Pe_std2=Pe_std(~indX1);
MaxE2=MaxE(~indX1);
MaxE2_er=MaxE_er(~indX1);
RUN2=RUN(~indX1);

%Energy0= 219.729*sin(0.066*atten-1.703)+208.919; %uJ

%[fluence ind]=sort(fluence);
%Energy0=Energy0(ind);

[Penergy, ind]=sort(Penergy);
Pe_std=Pe_std(ind);
MaxE=MaxE(ind);
MaxE_er=MaxE_er(ind);
RUN=RUN(ind);

col4=[0 .6 .4];

% MaxE=5.1;  % Noise floor level
% MaxE_er=1.5;
Grad=50*sqrt((.0484*MaxE+1.0463).^2-1)-23;
Grad_er=50*.0484*(.0484*MaxE+1.0463)./sqrt((.0484*MaxE+1.0463).^2-1).*MaxE_er;
var=Grad;
var_er=Grad_er;


%load V:\ARDB\E163\Data\130809\calibration_values % load simulation data
load ~/Dropbox/Research/Data/calibration_values
var = interp1(xcal,ycal,MaxE,'pchip');

% var1=50*sqrt((.0484*MaxE1+1.0463).^2-1)-23;
% var1_er=50*.0484*(.0484*MaxE1+1.0463)./sqrt((.0484*MaxE1+1.0463).^2-1).*MaxE1_er;
% var2=50*sqrt((.0484*MaxE2+1.0463).^2-1)-23;
% var2_er=50*.0484*(.0484*MaxE2+1.0463)./sqrt((.0484*MaxE2+1.0463).^2-1).*MaxE2_er;
% 
% var1 = interp1(xcal,ycal,MaxE1,'pchip');
% var2 = interp1(xcal,ycal,MaxE2,'pchip');

%%
% figure (5)
% plot(MaxE,Grad,'xb')
% hold on
% x=1:max(MaxE);
% plot(x,2.48.*x+19.3*ones(size(x)),'r')
% xlabel(' Energy shift [keV] ')
% ylabel(' Gradient [MeV/m]')
% axis tight
% enhance_plot;
% Relationship is linear after 40keV (118.5 MeV/m):   G=2.48 DE +  19.3

%%
% var=Grad.^2;
% var_er=2*Grad.*Grad_er;
%
figure (3)
hold on
%errorbar(Penergy,var,var_er,'o','LineWidth',2,'Color','r');%col4)
%h1=ploterr(Penergy,var,Pe_std,var_er,'.');%,'Color',col);
h1=errorbar(Penergy,var,var_er,'.');%,'Color',col);

set(h1,'Color',col)

% h1=ploterr(Penergy1,var1,Pe_std1,var1_er,'.b');%,'Color',col);
% hold on
% h2=ploterr(Penergy2,var2,Pe_std2,var2_er,'.r');%,'Color',col);

%set(h(2),'Color','b'), set(h(3),'Color','b')
xlabel('Laser Pulse Energy [\mu J]')
ylabel('Gradient [MeV/m]')
%grid on
axis tight
enhance_plot;
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

%%
n=1.5;
eps0=8.85e-12;
tau=1.3e-12;
c=3e8;
w1=60e-6/2;
w2=600e-6/2;

%%
Energy=(0:350)*1e-6;

F=2*Penergy./pi/w1/w2*1e-6;

%%
Efield=sqrt(2*F./eps0/c/tau);
dEfield=Pe_std*1e-6./sqrt(Penergy*1e-6)/sqrt(pi*eps0*c*tau*w1*w2);

% F1=2*Penergy1./pi/w1/w2*1e-6;
% Efield1=sqrt(2*F1./eps0/c/tau);
% dEfield1=Pe_std1*1e-6./sqrt(Penergy1*1e-6)/sqrt(pi*eps0*c*tau*w1*w2);
% 
% F2=2*Penergy2./pi/w1/w2*1e-6;
% Efield2=sqrt(2*F2./eps0/c/tau);
% dEfield2=Pe_std2*1e-6./sqrt(Penergy2*1e-6)/sqrt(pi*eps0*c*tau*w1*w2);


yW=var_er.^1;
linF =@(c,x) c(1)*ones(size(x))+c(2)*x;
linFW =@(c,x) (c(1)*ones(size(x))+c(2)*x)./yW;
[cout1,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(linFW,[0 1e-7],Efield,...
   var./yW,[-100 0],[100 1],SIG_OPTIONS);
EfieldF=0:max(Efield)/200:max(Efield)+max(Efield)/100;
gradF=linF(cout1,EfieldF);

% yW=var1_er.^1;
% linFW =@(c,x) (c(1)*ones(size(x))+c(2)*x)./yW;
% [cout1,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(linFW,[0 1e-7],Efield1,...
%    var1./yW,[-100 0],[100 1],SIG_OPTIONS);
% gradF1=linF(cout1,EfieldF);
% 
% yW=var2_er.^1;
% linFW =@(c,x) (c(1)*ones(size(x))+c(2)*x)./yW;
% [cout1,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(linFW,[0 1e-7],Efield2,...
%    var2./yW,[-100 0],[100 1],SIG_OPTIONS);
% gradF2=linF(cout1,EfieldF);

figure (4)
h2=plot(EfieldF*1e-9,gradF,'b:','LineWidth',3)
hold on
% h2a=plot(EfieldF*1e-9,gradF1,'b:','LineWidth',3)
% hold on
% h2b=plot(EfieldF*1e-9,gradF2,'r:','LineWidth',3)
% hold on

h5=plot(EfieldF*1e-9,(18+5.7)*ones(size(EfieldF)),'--k','LineWidth',2)%,
%h1=ploterr(Efield,MaxE,dEfield,MaxE_er,'.');%,
hold on

xlabel('Peak Incident Electric Field E_0 (GV/m)')
ylabel('Gradient (MeV/m)')
axis tight

% Etick=[6.6, 25,50,75,100,125];
% Gtick=50*sqrt((.0484*Etick+1.0463).^2-1)-23;
% set(gca,'YTick',Gtick)
% set(gca,'YTickLabel',{'6.6','25','50','75','100','125'});%,'85^2','105^2','120^2'})


Etick=[0,.5,1,1.5,2,2.5,3,3.5]
set(gca,'XTick',Etick)
enhance_plot;
%set(gca,'XTickLabel',{'0','0.5','1','1.5','2.0','2.5','3.0','3.5'});

% % Additional Axes
% Jtick=[6.6,15,30,50,75,100,125];
% Gtick=50*sqrt((.0484*Jtick+1.0463).^2-1)-23;
% set(gca,'YTick',Gtick)
% set(gca,'YTickLabel',{'6.6','15','30','50','75','100','125'});%,'85^2','105^2','120^2'})
% ylabel('Measured Max. Energy Shift (keV)')
% set(gca,'YAxisLocation','right')
% % 
% Jtick=[10 50:50:200 300];
% Etick=2*sqrt(Jtick./eps0/c/tau/pi/w1/w2*1e-6)*1e-9;
% set(gca,'XTick',Etick)
% set(gca,'XTickLabel',{'0.01','0.05','0.1','0.15','0.2','0.3'});%,'85^2','105^2','120^2'})
% xlabel('Laser Pulse Energy (mJ)')
% set(gca,'XAxisLocation','top')
% % 
% Boundary from fabrication

% 400nm gap best
delta1=1;
fD1=.2145;
%fD1=.249;
eta1=1.295 ;

% 400nm gap worst
delta2=.75;
fD2=.1546;
eta2=1.198 ;

%%
%Energy=259.4;

F=2*Energy./pi/w1/w2;
Efield0=sqrt(2*F./eps0/c/tau);
G1=delta1*fD1*Efield0./n/eta1
G2=delta2*fD2*Efield0./n/eta2
%
dGreen=[0 .6 .4];
dYellow=[1 .8 0];

hold on
h3=plot(Efield0*1e-9,G1*1e-6,'-.','Color',dGreen,'LineWidth',2)
hold on
h4=plot(Efield0*1e-9,G2*1e-6,'-.','Color','r','LineWidth',2)
%h4=plot(Efield0*1e-9,G2*1e-6,'-.','Color',dGreen,'LineWidth',2)
hold on

h1=ploterr(Efield*1e-9,var,dEfield*1e-9,var_er,'o');
set(h1(1),'LineWidth',2,'MarkerFaceColor',[0 .7 1],'MarkerSize',4);
set(h1(2),'LineWidth',2);
set(h1(3),'LineWidth',2);%,'MarkerFaceColor',[0 .7 1],'MarkerSize',5);
%             set(h2,'MarkerFaceColor',[0 .7 1],'MarkerSize',5,'LineWidth',1.5)
%             set(h1,'MarkerFaceColor',[1 .7 0],'MarkerSize',5,'LineWidth',1.5)

% h1=ploterr(Efield1*1e-9,var1,dEfield1*1e-9,var1_er,'ob');
% set(h1(1),'LineWidth',2,'MarkerFaceColor',[0 .7 1],'MarkerSize',4);
% set(h1(2),'LineWidth',2);
% set(h1(3),'LineWidth',2);%,'MarkerFaceColor',[0 .7 1],'MarkerSize',5);
% 
% hold on
% h2=ploterr(Efield2*1e-9,var2,dEfield2*1e-9,var2_er,'or');
% set(h2(1),'LineWidth',2,'MarkerFaceColor',[1 .7 0],'MarkerSize',4);
% set(h2(2),'LineWidth',2);
% set(h2(3),'LineWidth',2);%,'

box on

% hold on
% h2=plot(Energy*1e6,G1*1e-6,'-.','Color',col4)
% hold on
% plot(Energy*1e6,G2*1e-6,'-.','Color',col4)
% box on
ylim([-15 335])
legend([h1(1) h2 h3 h4 h5],{'Data', 'Fit','Sim:best','Sim:worst','Noise level'})
enhance_plot;

hold on

% % Run labels
% for i=1:length(Efield)
%    text(Efield(i)*1e-9-.15,var(i)+10,num2str(RUN(i),'%g'),'Rotation',-25)
%    hold on
% end

% %% Cross-correlation Plot
% 
% figure
% folder='V:\ARDB\E163\Data\130423\Fluence1\';
% RUNNO=2216;
% post='_f18-4_1.mat';
% 
% %folder='V:\ARDB\E163\Data\130517\';
% %RUNNO=2232;
% %post='_f18-4_3.mat';
% 
% filename=['run',num2str(RUNNO),post]
% 
% eval(['load ',folder,filename]);
%       
%       regen=Data.regen;
%       xVar=Data.delay./2;
%       yVar=Data.MaxEshift*1.2-362.35*ones(size(xVar));
%       yVar_er=Data.MaxEshift_er*1.2;
%      
%    eventNum=1:length(regen);
%    omitList=Data.omitList;
%    regen(omitList)=[];
%    eventNum(omitList)=[];
%    xVar(omitList)=[];
%    yVar(omitList)=[];
%    yVar_er(omitList)=[];
%    
%    ind_on=(regen==1);
%    ind_off=(regen==0);
%    
%    event_on=eventNum(ind_on);
%    event_off=eventNum(ind_off);
%    xVar_on=xVar(ind_on);
%    xVar_off=xVar(ind_off);
%    yVar_on=yVar(ind_on);
%    yVar_off=yVar(ind_off);
%    yVar_on_er=yVar_er(ind_on);
%    yVar_off_er=yVar_er(ind_off);
%    
%    [xVar_on, ind_on]=sort(xVar_on);
%    yVar_on=yVar_on(ind_on);
%    event_on=event_on(ind_on);
%    [xVar_off, ind_off]=sort(xVar_off);
%    yVar_off=yVar_off(ind_off);
%    event_off=event_off(ind_off);
%    yW=yVar_off_er;
%    
%          flatW=@(c,x) c(1)./yW;
%          [myVar_off,~,~,~,~,~,~]=lsqcurvefitstd(flatW,mean(yVar_off),xVar_off,...
%             yVar_off./yW,min(yVar_off),max(yVar_off),SIG_OPTIONS);
%    stdyVar_off=rms(yVar_off-myVar_off*ones(size(xVar_off)));
%   
%          sech2 =@(c,x) c(1)*sech((x-c(2))./c(3)).^2+c(4);
%          % initial guess - laser on
%          c1 = max(yVar_on);
%          c2=(xVar(end)+xVar(1))/2;
%          c3=c2/4;
%          c4=5;
%          guess=[c1 c2 c3 c4];
%          
%          cmin=[0 min(xVar_on) 1 0];
%          cmax=[2*max(yVar_on) max(xVar_on) max(xVar_on) max(yVar_on)];
%          
%             pow=1;
%             yW=yVar_on_er.^pow;  %standard deviations
%             sech2W =@(c,x) (c(1)*sech((x-c(2))./c(3)).^2+c(4))./yW;          
%       
%             [cout,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(sech2W,guess,xVar_on,yVar_on./yW,cmin,cmax,SIG_OPTIONS);
%             
%         
%         width=cout(3)*2/1.1346
%         width_er=cout_std(3)*2/1.1346
%           
%         xxvar=linspace(xVar_on(1),xVar_on(end));
%          h1f=plot(xxvar,sech2(cout,xxvar),'LineWidth',2.5,'Color',[1 0 0])
%          
%          hold on; 
%          h2f=line([min(xVar) max(xVar)],[myVar_off myVar_off],...
%                'LineWidth',2.5,'Color','b');
%          h2f2=line([min(xVar) max(xVar)],[myVar_off+stdyVar_off myVar_off+stdyVar_off],...
%                'LineWidth',2.5,'Color','k','LineStyle','--');
%             hold on
%             h2=errorbar(xVar_off,yVar_off,yVar_off_er,'ob');
%             
%             hold on
%             h1=errorbar(xVar_on,yVar_on,yVar_on_er,'or');
%             
%             set(h2,'MarkerFaceColor',[0 .7 1],'MarkerSize',5,'LineWidth',1.5)
%             set(h1,'MarkerFaceColor',[1 .7 0],'MarkerSize',5,'LineWidth',1.5)
%             axis tight
%             xlabel('Laser Delay (ps)')
%            
%             ylabel('Max. Energy Shift (keV)')
%             xlim([6 12])
%             ylim([-15 85])
%             legend([h1 h1f h2 h2f],{'Laser On','Sech Fit','Laser Off','Fit'})
%             enhance_plot;
            
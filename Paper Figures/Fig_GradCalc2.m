% Made figures LengthOptim.eps and GradGainSim.eps with this script
clear all
close all

%cols(0,:)=[1 1 1];  %use 'w'
cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple
cols(5,:)=[0 0 0]; %use 'k'

c=299792458;
eps0=8.85e-12;

ps=1e-12;
fs=1e-15;
um=1e-6;

Fin=.85/(.01)^2;
n=1.548;

sigx=78.7*um;       %RMS spotsize (radius)
sigz0=309.7*um;
%Uin=.5*Fin*pi*sigx*sigz
Uin=3.26e-4;
Fin0=2*Uin/(pi*sigx*sigz0)

tau=1.25*ps;
%Fin0=1.8e3;
Efield0=sqrt(2*Fin0./eps0/c/tau)


%400nm structure
fA1=.131;
eta1=2.76;

fA2=.094;
eta2=2.82;


% 800nm structure
fA1=.0769;
eta1=3.561 ;

fA2=.0575;
eta2=3.472;


A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-9;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-9;

%%
% Equation for Fluence vs pulsewidth:

c1=957.0339;
c2=7.1512e6;
Fmeas=@(t) c1+c2*t.^.25;

%%
% ******* Display G & E values for a given L and taup
% clc
% 
% L=1000*um;
% 
% %Fin=2*Uin/(pi*sigx*L/2)
% Fin=2*Uin/(pi*sigx*sigz)
% 
% taup=2.66*ps;
% 
% Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./taup).^2;
% 
% A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-6;
% A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-6;
% 
% Grad1=A1*quad(Pulse,0,L/2/c)/L/sqrt(taup)
% Grad2=A2*quad(Pulse,0,L/2/c)/L/sqrt(taup)
% 
% E1=Grad1*L*1e3
% E2=Grad2*L*1e3


L=450*um;
tau=1.25*ps;

Fin=2*Uin/(pi*sigx*L/2);  

if Fin>Fmeas(tau)
    Fin=Fmeas(tau);
    sigz=2*Uin/(pi*sigx*Fin);
    %fprintf('j=%d,i=%d, Fluence exceeded -> adjusted \n',j,i)
else
    sigz=L/2;
end
A1=(eta1*fA1/n)*sqrt(2*Fin/eps0/c);
A2=(eta2*fA2/n)*sqrt(2*Fin/eps0/c);

Pulse=@(z,t) exp(-z.^2/sigz^2/2).*(sech(2*asech(1/sqrt(2))*t./tau).^2);
Grad1=A1*4*quad2d(Pulse,0,L/2,0,L/2/c)/L/(L/c)/sqrt(tau)
Grad2=A2*4*quad2d(Pulse,0,L/2,0,L/2/c)/L/(L/c)/sqrt(tau)
E1=Grad1*L
E2=Grad2*L


%%  Didn't quite get this to work to show graphically a good vs bad match

% T=(0:.05:10)*ps;
% Nt=length(T);
% del=T(round(length(T)/2));%1*ps;
% 
% figure
% plot(T*1e12,Pulse(T-del))
% 
% Ntaup=round(taup/(T(2)-T(1)));
% Taup=linspace(0,taup,Nt-1);
% 
% % Z=linspace(0,L,Nt-1);
% % for i=1:Nt
% % 
% % %z=linspace(0,L,Ntaup-1)
% % %for i=1:Ntaup
% %     PulseMap(i,:)=circshift(Pulse(T-del)',i-round(Nt/2)-round(Ntaup/2)+2);
% % end
% 
% Z=linspace(0,L,Ntaup-1);
% for i=1:Nt
% 
% %z=linspace(0,L,Ntaup-1)
% %for i=1:Ntaup
%     PulseMap(i,:)=circshift(Pulse(T-del)',i-round(Nt/2)-round(Ntaup/2)+2);
% end
% 
% figure
% imagesc(Z*1e6,Taup*1e12,PulseMap)
% %imagesc(z*1e6,1:Ntaup,PulseMap)


%% Gradient calculation with fixed fluence

close all

%width=4.5; % if alone
width=3.25; % if paired
height=2.5;

hFig = figure(95);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
%set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')



%L=[150 310 620 1240]*um;
%L=[33 100 260 620 1240]*um;
L=[4 100 200 400 1000 620]*um;


tau=logspace(-14,-11,1000);
Grad1=zeros(size(tau));
Grad2=zeros(size(tau));

%L=[100 250 550 1100]*um;

%L=1000*um;

%for j=1:length(L)
for j=1:1
for i=1:length(tau)

Fin=2*Uin/(pi*sigx*L(j)/2);  %VARIABLE

if Fin>Fmeas(tau(i))
    Fin=Fmeas(tau(i));
    sigz=2*Uin/(pi*sigx*Fin);
    %fprintf('j=%d,i=%d, Fluence exceeded -> adjusted \n',j,i)
else
    sigz=L(j)/2;
end
A1=(eta1*fA1/n)*sqrt(2*Fin/eps0/c);
A2=(eta2*fA2/n)*sqrt(2*Fin/eps0/c);

Pulse=@(z,t) exp(-z.^2/sigz^2/2).*(sech(2*asech(1/sqrt(2))*t./tau(i)).^2);
Grad1(i)=A1*4*quad2d(Pulse,0,L(j)/2,0,L(j)/2/c)/L(j)/(L(j)/c)/sqrt(tau(i));
Grad2(i)=A2*4*quad2d(Pulse,0,L(j)/2,0,L(j)/2/c)/L(j)/(L(j)/c)/sqrt(tau(i));
Sigz(i)=sigz;
end


% % % Uncomment if Longitudinal Pulsewidth
% if j==6
% semilogx(tau*1e12,Sigz*1e6,'--k','Linewidth',1);
% hold on
% else
%     semilogx(tau*1e12,Sigz*1e6,'Color',cols(6-j,:),'Linewidth',1);
% hold on
% end

 
% % Uncomment if Energy gain
% Grad1=Grad1*L(j)*1e9*1e-3;
% Grad2=Grad2*L(j)*1e9*1e-3;

if j==6
semilogx(tau*1e12,Grad1*1e-9,'--k','Linewidth',1);
hold on
semilogx(tau*1e12,Grad2*1e-9,'--k','Linewidth',1);
hold on
else
    semilogx(tau*1e12,Grad1*1e-9,'Color',cols(6-j,:),'Linewidth',1);
hold on
semilogx(tau*1e12,Grad2*1e-9,'Color',cols(6-j,:),'Linewidth',1);
hold on
end


end

[maxV ind]=max(Grad1);
maxG=maxV
pulse=tau(ind)
ax1=gca;
set(gca,'Color','none');

% ylim([0 180])
% set(ax1,'YTick',20:20:180)
% set(ax1,'XTickLabel',{'0.03','','','','','','','0.1','','0.3','','0.5','','','','','1.0','','2','3','4','5'})
% hy=ylabel('Energy gain, \Delta{\itE} [keV]');


%ylim([0 1.12])
%set(ax1,'YTick',[.1:.1:1.1])
%hy=ylabel('Average gradient, {\itG} [GeV/m]');

hy=ylabel('Damage length, {\itL}_{dam} [\mum]');

%xlim([.03 5])
set(ax1,'XTick',[.03:.01:.09 .1:.1:1 1.25 2:5])
set(ax1,'XTickLabel',{'0.03','','','','','','','0.1','','0.3','','0.5','','','','','1.0','','2','3','4','5'})
hx=xlabel('FWHM pulse width, \tau_p [ps]');

ax2=copyobj(ax1,hFig);
delete(get(ax2,'Children'));
set(ax2,'Color','None','Box','off','Ygrid','on','Xgrid','on');
set(ax2,'Xcolor',[.9 .9 .9],'Ycolor',[.9 .9 .9]);%,'XTickLabel',[],'YTickLabel',[]);
%get(ax3)    %Unsure how to remove labels from this third axes.. will have
%to do by hand
uistack(ax2,'bottom');

linkaxes([ax1,ax2],'xy')

%hL=legend([h1a,h2a,h3a,h4a],{'{\it L} = 0.15 mm','{\it L} = 0.31 mm','{\it L} = 0.62 mm','{\it L} = 1.24 mm'});
%set(hL,'Location','northeast','Box','Off')

set([ax1,ax2,hx,hy],'FontSize', 10);
%set([hx,hy],'FontSize', 10);

%export_fig GradientvTauvLength.eps -cmyk -r300 -painters%-painters

 
%% %% Plot of energy gain vs taup for different structure lengths
% 
% %To match 2sigz=L ->sigz=L/2
% 
% %close all
% 
% width=3.25;
% height=2.5;
% hFig = figure(97);
% 
% %set(gcf,'PaperPositionMode','auto')
% set(hFig,'ActivePositionProperty','position')
% set(hFig,'Units','inches')
% %set(hFig, 'Position', [0 0 width height])
% 
% set(hFig, 'Color', 'w');
% set(hFig,'Units','points')
% 
% 
% L=[150 310 620 1240]*um;
% tau=logspace(-14,-11,100);
% 
% L=620*um;
% tau=logspace(-14,-11,1000);
% 
% 
% scale=1e-3;
% j=1;
% %Fin=2*Uin/(pi*sigx*L(j)/2)
% Fin=2*Uin/(pi*sigx*sigz)
% for i=1:length(tau)
% Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
% A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
% A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;
% 
% Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% end
% 
% E1=Grad1*L(j);
% E2=Grad2*L(j);
% 
% h1a=semilogx(tau*1e12,E1,'Color',cols(3,:),'Linewidth',1);
% hold on
% h1b=semilogx(tau*1e12,E2,'Color',cols(3,:),'Linewidth',1);%,'LineStyle','--');
% 
% %%
% hold on
% 
% j=2;
% %Fin=2*Uin/(pi*sigx*L(j)/2)
% Fin=2*Uin/(pi*sigx*sigz)
% for i=1:length(tau)
% Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
% 
% A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
% A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;
% Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% end
% 
% E1=Grad1*L(j);
% E2=Grad2*L(j);
% 
% h2a=semilogx(tau*1e12,E1,'Color',cols(1,:),'Linewidth',1);
% hold on
% h2b=semilogx(tau*1e12,E2,'Color',cols(1,:),'Linewidth',1);%,'LineStyle','--');
% 
% hold on
% 
% j=3;
% %Fin=2*Uin/(pi*sigx*L(j)/2)
% Fin=2*Uin/(pi*sigx*sigz)
% 
% for i=1:length(tau)
% Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
% A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
% A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;
% Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% end
% 
% E1=Grad1*L(j);
% E2=Grad2*L(j);
% 
% h3a=semilogx(tau*1e12,E1,'Color',cols(4,:),'Linewidth',1);
% hold on
% h3b=semilogx(tau*1e12,E2,'Color',cols(4,:),'Linewidth',1);%,'LineStyle','--');
% 
% 
% hold on
% 
% j=4;
% Fin=2*Uin/(pi*sigx*L(j)/2)
% %Fin=2*Uin/(pi*sigx*sigz)
% 
% for i=1:length(tau)
% Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
% 
% A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
% A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;
% Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% end
% 
% E1=Grad1*L(j);
% E2=Grad2*L(j);
% 
% h4a=semilogx(tau*1e12,E1,'Color',cols(2,:),'Linewidth',1);
% hold on
% h4b=semilogx(tau*1e12,E2,'Color',cols(2,:),'Linewidth',1)%,'LineStyle','--');
% 
% uistack(h2b, 'top');
% uistack(h2a, 'top');
% uistack(h1b, 'top');
% uistack(h1a, 'top');
% 
% 
% % j=5;
% % Fin=2*Uin/(pi*sigx*L(j)/2)
% % %Fin=2*Uin/(pi*sigx*sigz)
% % 
% % for i=1:length(tau)
% % Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
% % 
% % A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
% % A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;
% % Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% % Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% % end
% % 
% % E1=Grad1*L(j);
% % E2=Grad2*L(j);
% % 
% % h4a=semilogx(tau*1e12,E1,'Color',cols(2,:),'Linewidth',1);
% % hold on
% % h4b=semilogx(tau*1e12,E2,'Color',cols(2,:),'Linewidth',1,'LineStyle','--');
% 
% 
% 
% ax1=gca;
% set(gca,'xscale','log','Color','none');
% 
% %axis tight
% ylim([0 225])
% xlim([.1 5])
% 
% hx=xlabel('FWHM pulse width, \tau_p [ps]');
% set(ax1,'XTick',[.1:.1:1 1.25 1.5 2:5]);
% set(ax1,'XTickLabel',{'0.1','0.2','0.3','','0.5','','','','','1.0','','1.5','2','3','4','5'});
% 
% hy=ylabel('Average energy gain, \Delta{\it E} [keV]');
% set(ax1,'YTick',[0:25:225]);
% %set(ax1,'YTickLabel',{'','0.1','','0.3','','0.5','','0.7','','0.9','','1.1
% %','','1.3','','1.5','','1.7'});
% 
% ax2=copyobj(ax1,hFig);
% delete(get(ax2,'Children'));
% set(ax2,'Color','None','Box','off','Ygrid','on','Xgrid','on');
% set(ax2,'Xcolor',[.9 .9 .9],'Ycolor',[.9 .9 .9]);%,'XTickLabel',[],'YTickLabel',[]);
% %get(ax3)    %Unsure how to remove labels from this third axes.. will have
% %to do by hand
% uistack(ax2,'bottom');
% 
% linkaxes([ax1,ax2],'xy')
% 
% hL=legend([h1a,h2a,h3a,h4a],{'{\it L} = 0.15 mm','{\it L} = 0.31 mm','{\it L} = 0.62 mm','{\it L} = 1.24 mm'});
% set(hL,'Location','northwest','Box','Off')
% 
% set([ax1,ax2],'FontSize', 10);
% set([hx,hy,hL],'FontSize', 10);
% 
% export_fig LengthOptimE.eps -cmyk -r300 -painters%-painters

%%  Plot of energy gain vs L for different taup

close all
clear tau L E1 E2 Emax Lmax Grad1 Grad2

cols2=[cols(4,:); cols(1,:);cols(5,:);cols(2,:);cols(3,:)];
width=4.5;
height=2.5;
hFig = figure(97);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')

L=(10:10:5000)*um;
tau=[.100 .25 .613 1.67 2.5]*ps;

%tau=logspace(-14,-11,50);

%tau=(50:50:5000)*fs;

%tau=logspace(-13.5,-11,100);


col='rbkycgm';

scale=1e-3;
for j=1:length(tau)
    Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(j)).^2;
    for i=1:length(L)
        
        
        if L(i)<620*um
            Fin=2*Uin/(pi*sigx*sigz);
        else
            Fin=2*Uin/(pi*sigx*L(i)/2);
        end
        
        A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
        A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;
        
        %E1(i)=A1*quad(Pulse,0,L(i)/2/c)/sqrt(tau(j));
        %E2(i)=A2*quad(Pulse,0,L(i)/2/c)/sqrt(tau(j));
        
        Grad1(i)=A1*quad(Pulse,0,L(i)/2/c)/L(i)/sqrt(tau(j));
        Grad2(i)=A2*quad(Pulse,0,L(i)/2/c)/L(i)/sqrt(tau(j));
        E1(i)=Grad1(i)*L(i);
        E2(i)=Grad2(i)*L(i);
        
    end
    [Emax indMax]=max(E1);
    Lmax(j)=L(indMax);
    
    h1a=plot(L*1e3,E1,'Color',cols2(j,:),'Linewidth',1);
    hold on
    h1b=plot(L*1e3,E2,'Color',cols2(j,:),'Linewidth',1);%,'LineStyle','--');
    hold on
    
end


for j=1:length(L)
%for j=3:3
for i=1:length(tau)

Fin=2*Uin/(pi*sigx*L(j)/2);  %VARIABLE

if Fin>Fmeas(tau(i))
    Fin=Fmeas(tau(i));
    sigz=2*Uin/(pi*sigx*Fin);
    %fprintf('j=%d,i=%d, Fluence exceeded -> adjusted \n',j,i)
else
    sigz=L(j)/2;
end
A1=(eta1*fA1/n)*sqrt(2*Fin/eps0/c);
A2=(eta2*fA2/n)*sqrt(2*Fin/eps0/c);

Pulse=@(z,t) exp(-z.^2/sigz^2/2).*(sech(2*asech(1/sqrt(2))*t./tau(i)).^2);
Grad1(i)=A1*4*quad2d(Pulse,0,L(j)/2,0,L(j)/2/c)/L(j)/(L(j)/c)/sqrt(tau(i));
Grad2(i)=A2*4*quad2d(Pulse,0,L(j)/2,0,L(j)/2/c)/L(j)/(L(j)/c)/sqrt(tau(i));
Sigz(i)=sigz;
end




ax1=gca;
set(gca,'Color','none');

xlim([0 1.2])
ylim([0 210])

hx=xlabel('Structure length, {\it L} [mm]');
set(ax1,'XTick',[.1:.1:1.2]);
%set(ax1,'XTickLabel',{'0.1','0.2','0.3','','0.5','','','','','1.0','','1.5','2','3','4','5'});

hy=ylabel('Average energy gain, \Delta{\it E} [keV]');
set(ax1,'YTick',[0:25:225]);
%set(ax1,'YTickLabel',{'','0.1','','0.3','','0.5','','0.7','','0.9','','1.1
%','','1.3','','1.5','','1.7'});

ax2=copyobj(ax1,hFig);
delete(get(ax2,'Children'));
set(ax2,'Color','None','Box','off','Ygrid','on','Xgrid','on');
set(ax2,'Xcolor',[.9 .9 .9],'Ycolor',[.9 .9 .9]);%,'XTickLabel',[],'YTickLabel',[]);
%get(ax3)    %Unsure how to remove labels from this third axes.. will have
%to do by hand
uistack(ax2,'bottom');

linkaxes([ax1,ax2],'xy')

%hL=legend([h1a,h2a,h3a,h4a],{'{\it L} = 0.15 mm','{\it L} = 0.31 mm','{\it L} = 0.62 mm','{\it L} = 1.24 mm'});
%set(hL,'Location','northwest','Box','Off')

set([ax1,ax2],'FontSize', 10);
set([hx,hy],'FontSize', 10);

export_fig LengthOptim.eps -cmyk -r300 -painters%-painters

% 
% %% Fluence variation with pulsewidth
% 
% close all
% tau=(.01:.01:5)*ps;
% %at 1um
% tau2=5*ps;
% tau1=.1*ps;
% 
% F2=3.3e4/2;
% F1=1.3e4/2;
% 
% c2a=(F2-F1)/(tau2^.25-tau1^.25);
% c1a=F1-c2a*tau1^.25;
% 
% Fa=c1a+c2a*tau.^.25;
% 
% %at .5um
% tau2=2*ps;
% tau1=.2*ps;
% 
% F2=1.6e4/2;
% F1=1e4/2;
% 
% c2b=(F2-F1)/(tau2^.25-tau1^.25);
% c1b=F1-c2b*tau1^.25;
% 
% Fb=c1b+c2b*tau.^.25;
% %F800=.72*Fb+.28*Fa;   %Gets what i need
% c1=.72*c1b+.28*c1a
% c2=.72*c2b+.28*c2a
% F800=c1+c2*tau.^.25;
% 
% 
% figure
% semilogx(tau,Fa,tau,Fb,tau,F800)

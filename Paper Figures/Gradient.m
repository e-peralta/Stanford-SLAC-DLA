clear all
close all

%cols(0,:)=[1 1 1];  %use 'w'
cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple
%cols(5,:)=[0 0 0]; %use 'k'

c=299792458;
eps0=8.85e-12;

ps=1e-12;
fs=1e-15;
um=1e-6;

Fin=.85/(.01)^2;
n=1.548;

sigx=78.7*um;       %RMS spotsize (radius)
sigz=309.7*um;
%Uin=.5*Fin*pi*sigx*sigz
Uin=3.19e-4;


fA1=.131;
eta1=2.76;

fA2=.094;
eta2=2.82;


A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-6;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-6;

%%
clc

L=600*um;

Fin=2*Uin/(pi*sigx*L/2)

taup=1.25*ps;

Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./taup).^2;
Gradh=A1*quad(Pulse,0,L/2/c)/L/sqrt(taup)
Gradl=A2*quad(Pulse,0,L/2/c)/L/sqrt(taup)



A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-6;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-6;

GradOpth=A1*quad(Pulse,0,L/2/c)/L/sqrt(taup)
GradOptl=A2*quad(Pulse,0,L/2/c)/L/sqrt(taup)

%%
T=(0:.05:10)*ps;
Nt=length(T);
del=T(round(length(T)/2));%1*ps;

figure
plot(T*1e12,Pulse(T-del))

Ntaup=round(taup/(T(2)-T(1)));
Taup=linspace(0,taup,Nt-1);

% Z=linspace(0,L,Nt-1);
% for i=1:Nt
% 
% %z=linspace(0,L,Ntaup-1)
% %for i=1:Ntaup
%     PulseMap(i,:)=circshift(Pulse(T-del)',i-round(Nt/2)-round(Ntaup/2)+2);
% end

Z=linspace(0,L,Ntaup-1);
for i=1:Nt

%z=linspace(0,L,Ntaup-1)
%for i=1:Ntaup
    PulseMap(i,:)=circshift(Pulse(T-del)',i-round(Nt/2)-round(Ntaup/2)+2);
end

figure
imagesc(Z*1e6,Taup*1e12,PulseMap)
%imagesc(z*1e6,1:Ntaup,PulseMap)



% %% Gradient calculation with fixed fluence
% 
% 
% close all
% 
% width=4;
% height=2.5;
% hFig = figure(95);
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
% L=[100 250 500 1000]*um;
% %L=[100 250 550 1100]*um;
% tau=logspace(-14,-11,50);
% 
% j=1;
% for i=1:length(tau)
% Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
% 
% Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% end
% 
% Grad=(Grad1+Grad2)/2;
% dGrad=(Grad2-Grad1)/2;
% h1=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
%         'Marker','o','Color',cols(3,:),'Linewidth',1,'MarkerSize',3,...
%     'MarkerFaceColor',cols(3,:));
% 
% 
% hold on
% 
% j=2;
% tau=logspace(-14,-11,51);
% for i=1:length(tau)
% Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
% 
% Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% end
% 
% Grad=(Grad1+Grad2)/2;
% dGrad=(Grad2-Grad1)/2;
% h2=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
%         'Marker','^', 'Color',cols(1,:),'Linewidth',1,...
%     'MarkerSize',3,'MarkerFaceColor',cols(1,:));
% 
% hold on
% 
% j=3;
% clear tau Grad1 Grad2 Grad dGrad
% tau=logspace(-14,-11,49);
% %tau=logspace(-14,-11,100);
% for i=1:length(tau)
% Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
% 
% Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% end
% 
% Grad=(Grad1+Grad2)/2;
% dGrad=(Grad2-Grad1)/2;
% h3=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
% 'Marker','s','Color',cols(4,:),'Linewidth',1,'MarkerSize',3,...
% 'MarkerFaceColor',cols(4,:));
% 
% 
% hold on
% 
% j=4;
% tau=logspace(-14,-11,50);
% for i=1:length(tau)
% Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
% 
% Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% end
% 
% Grad=(Grad1+Grad2)/2;
% dGrad=(Grad2-Grad1)/2;
% h4=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
%     'Marker','v','Color',cols(2,:),'Linewidth',1,...
%     'MarkerSize',3,'MarkerFaceColor',cols(2,:));
%     
% 
% ax1=gca;
% set(gca,'xscale','log','Color','none');
% 
% %axis tight
% ylim([0 900])
% xlim([.1 5])
% 
% set(ax1,'XTick',[.1:.1:1 1.5 2:5])
% set(ax1,'YTick',[100:100:800])
% set(ax1,'XTickLabel',{'0.1','0.2','0.3','','0.5','','','','','1.0','1.5','2','3','4','5'})
% 
% hx=xlabel('FWHM pulse width, \tau_p [ps]');
% hy=ylabel('Average gradient, G_0 [MeV/m]');
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
% 
% 
% hL=legend([h1,h2,h3,h4],{'{\it L} = 0.1 mm','{\it L} = 0.25 mm','{\it L} = 0.5 mm','{\it L} = 1 mm'});
% set(hL,'Location','northeast','Box','Off')
% set([ax1,ax2],'FontSize', 10);
% set([hx,hy,hL],'FontSize', 10);
% 
% %export_fig LengthOptim.eps -cmyk -r300 -painters%-painters


%% Gradient calculation with fixed pulse engery


%To match 2sigz=L ->sigz=L/2

close all

width=3.25;
height=2.5;
hFig = figure(96);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')


L=[120 240 480 960]*um;

j=1;
tau=logspace(-14,-11,60);
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;

Fin=2*Uin/(pi*sigx*L(j)/2);
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-9;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-9;

Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

Grad=(Grad1+Grad2)/2;
dGrad=(Grad2-Grad1)/2;
h1=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
        'Marker','o','Color',cols(3,:),'Linewidth',1,'MarkerSize',3,...
    'MarkerFaceColor',cols(3,:));
hold on

j=2;
tau=logspace(-14,-11,61);
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;

Fin=2*Uin/(pi*sigx*L(j)/2);
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-9;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-9;
Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

Grad=(Grad1+Grad2)/2;
dGrad=(Grad2-Grad1)/2;
h2=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
        'Marker','^', 'Color',cols(1,:),'Linewidth',1,...
    'MarkerSize',3,'MarkerFaceColor',cols(1,:));
hold on

j=3;
clear tau Grad1 Grad2 Grad dGrad
tau=logspace(-14,-11,59);
%tau=logspace(-14,-11,100);
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
Fin=2*Uin/(pi*sigx*L(j)/2);
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-9;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-9;
Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

Grad=(Grad1+Grad2)/2;
dGrad=(Grad2-Grad1)/2;
h3=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
'Marker','s','Color',cols(4,:),'Linewidth',1,'MarkerSize',3,...
'MarkerFaceColor',cols(4,:));

hold on

j=4;
tau=logspace(-14,-11,60);
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
Fin=2*Uin/(pi*sigx*L(j)/2);
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-9;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-9;
Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

Grad=(Grad1+Grad2)/2;
dGrad=(Grad2-Grad1)/2;
h4=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
    'Marker','v','Color',cols(2,:),'Linewidth',1,...
    'MarkerSize',2,'MarkerFaceColor',cols(2,:));
ax1=gca;
set(gca,'xscale','log','Color','none');

%axis tight
ylim([0 1.75])
xlim([.1 5])

hx=xlabel('FWHM pulse width, \tau_p [ps]');
set(ax1,'XTick',[.1:.1:1 1.5 2:5]);
set(ax1,'XTickLabel',{'0.1','0.2','0.3','','0.5','','','','','1.0','1.5','2','3','4','5'});

hy=ylabel('Average gradient, {\it G}_0 [GeV/m]');
set(ax1,'YTick',[0:.1:1.7]);
set(ax1,'YTickLabel',{'','0.1','','0.3','','0.5','','0.7','','0.9','','1.1','','1.3','','1.5','','1.7'});

ax2=copyobj(ax1,hFig);
delete(get(ax2,'Children'));
set(ax2,'Color','None','Box','off','Ygrid','on','Xgrid','on');
set(ax2,'Xcolor',[.9 .9 .9],'Ycolor',[.9 .9 .9]);%,'XTickLabel',[],'YTickLabel',[]);
%get(ax3)    %Unsure how to remove labels from this third axes.. will have
%to do by hand
uistack(ax2,'bottom');

linkaxes([ax1,ax2],'xy')

hL=legend([h1,h2,h3,h4],{'{\it L} = 0.12 mm','{\it L} = 0.24 mm','{\it L} = 0.48 mm','{\it L} = 0.96 mm'});
set(hL,'Location','northeast','Box','Off')

set([ax1,ax2],'FontSize', 10);
set([hx,hy,hL],'FontSize', 10);

export_fig LengthOptim.eps -cmyk -r300 -painters%-painters
%% Plot energy gain now

%To match 2sigz=L ->sigz=L/2

close all

width=3.25;
height=2.5;
hFig = figure(97);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')


L=[120 240 480 960]*um;


tau=logspace(-14,-11,100);

scale=1e-3;
j=1;
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;

Fin=2*Uin/(pi*sigx*L(j)/2);
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;

Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

E1=Grad1*L(j);
E2=Grad2*L(j);

h1a=semilogx(tau*1e12,E1,'Color',cols(3,:),'Linewidth',1);
hold on
h1b=semilogx(tau*1e12,E2,'Color',cols(3,:),'Linewidth',1,'LineStyle','--');

hold on

j=2;
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;

Fin=2*Uin/(pi*sigx*L(j)/2);
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;
Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

E1=Grad1*L(j);
E2=Grad2*L(j);

h2a=semilogx(tau*1e12,E1,'Color',cols(1,:),'Linewidth',1);
hold on
h2b=semilogx(tau*1e12,E2,'Color',cols(1,:),'Linewidth',1,'LineStyle','--');

hold on

j=3;
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
Fin=2*Uin/(pi*sigx*L(j)/2);
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;
Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

E1=Grad1*L(j);
E2=Grad2*L(j);

h3a=semilogx(tau*1e12,E1,'Color',cols(4,:),'Linewidth',1);
hold on
h3b=semilogx(tau*1e12,E2,'Color',cols(4,:),'Linewidth',1,'LineStyle','--');


hold on

j=4;
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
Fin=2*Uin/(pi*sigx*L(j)/2);
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;
Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

E1=Grad1*L(j);
E2=Grad2*L(j);

h4a=semilogx(tau*1e12,E1,'Color',cols(2,:),'Linewidth',1);
hold on
h4b=semilogx(tau*1e12,E2,'Color',cols(2,:),'Linewidth',1,'LineStyle','--');

ax1=gca;
set(gca,'xscale','log','Color','none');

%axis tight
ylim([0 225])
xlim([.1 5])

hx=xlabel('FWHM pulse width, \tau_p [ps]');
set(ax1,'XTick',[.1:.1:1 1.5 2:5]);
set(ax1,'XTickLabel',{'0.1','0.2','0.3','','0.5','','','','','1.0','1.5','2','3','4','5'});

hy=ylabel('Average energy gain, \Delta{\it E} [keV]');
set(ax1,'YTick',[0:25:225]);
%set(ax1,'YTickLabel',{'','0.1','','0.3','','0.5','','0.7','','0.9','','1.1
%','','1.3','','1.5','','1.7'});

%**************************

ax2=copyobj(ax1,hFig);
delete(get(ax2,'Children'));
set(ax2,'Color','None','Box','off','Ygrid','on','Xgrid','on');
set(ax2,'Xcolor',[.9 .9 .9],'Ycolor',[.9 .9 .9]);%,'XTickLabel',[],'YTickLabel',[]);
%get(ax3)    %Unsure how to remove labels from this third axes.. will have
%to do by hand
uistack(ax2,'bottom');

linkaxes([ax1,ax2],'xy')

hL=legend([h1a,h2a,h3a,h4a],{'{\it L} = 0.12 mm','{\it L} = 0.24 mm','{\it L} = 0.48 mm','{\it L} = 0.96 mm'});
set(hL,'Location','southeast','Box','Off')

set([ax1,ax2],'FontSize', 10);
set([hx,hy,hL],'FontSize', 10);

export_fig LengthOptimE.eps -cmyk -r300 -painters%-painters

%%

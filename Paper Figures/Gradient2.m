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
sigz=309.7*um;
%Uin=.5*Fin*pi*sigx*sigz
Uin=3.25e-4;
Fin0=2*Uin/(pi*sigx*sigz)


fA1=.131;
eta1=2.76;

fA2=.094;
eta2=2.82;


A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-9;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-9;

%%
clc

L=1000*um;

%Fin=2*Uin/(pi*sigx*L/2)
Fin=2*Uin/(pi*sigx*sigz)

taup=2.66*ps;

Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./taup).^2;

A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-6;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-6;

Grad1=A1*quad(Pulse,0,L/2/c)/L/sqrt(taup)
Grad2=A2*quad(Pulse,0,L/2/c)/L/sqrt(taup)

E1=Grad1*L*1e3
E2=Grad2*L*1e3

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



L=[150 310 620 1240]*um;
tau=logspace(-14,-11,60);
%L=[100 250 550 1100]*um;

L=1000*um;
tau=logspace(-14,-11,1000);

j=1;
Fin=2*Uin/(pi*sigx*sigz)       %FIXED
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-9;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-9;
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

h1a=semilogx(tau*1e12,Grad1,'Color',cols(3,:),'Linewidth',1);
hold on
h1b=semilogx(tau*1e12,Grad2,'Color',cols(3,:),'Linewidth',1);
hold on

grid on
%%
% Grad=(Grad1+Grad2)/2;
% dGrad=(Grad2-Grad1)/2;
% h1=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
%         'Marker','o','Color',cols(3,:),'Linewidth',1,'MarkerSize',3,...
%     'MarkerFaceColor',cols(3,:));


j=2;
Fin=2*Uin/(pi*sigx*sigz)       %FIXED
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-9;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-9;

tau=logspace(-14,-11,61);
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;

Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

h2a=semilogx(tau*1e12,Grad1,'Color',cols(1,:),'Linewidth',1);
hold on
h2b=semilogx(tau*1e12,Grad2,'Color',cols(1,:),'Linewidth',1);
hold on

% Grad=(Grad1+Grad2)/2;
% dGrad=(Grad2-Grad1)/2;
% h2=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
%         'Marker','^', 'Color',cols(1,:),'Linewidth',1,...
%     'MarkerSize',3,'MarkerFaceColor',cols(1,:));
% 
% hold on

% Fin=2*Uin/(pi*sigx*L(j)/2)     % VARIABLE
% for i=1:length(tau)
% Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
% 
% A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-9;
% A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-9;
% 
% Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
% end
% 
% Grad=(Grad1+Grad2)/2;
% dGrad=(Grad2-Grad1)/2;
% h1=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
%         'Marker','o','Color',cols(1,:),'Linewidth',1,'MarkerSize',3,...
%     'MarkerFaceColor',cols(1,:));
% hold on

j=3;
clear tau Grad1 Grad2 Grad dGrad
tau=logspace(-14,-11,59);
%tau=logspace(-14,-11,100);
Fin=2*Uin/(pi*sigx*sigz)       %FIXED
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-9;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-9;
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;

Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

h3a=semilogx(tau*1e12,Grad1,'Color',cols(4,:),'Linewidth',1);
hold on
h3b=semilogx(tau*1e12,Grad2,'Color',cols(4,:),'Linewidth',1);
hold on

% Grad=(Grad1+Grad2)/2;
% dGrad=(Grad2-Grad1)/2;
% h3=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
% 'Marker','s','Color',cols(4,:),'Linewidth',1,'MarkerSize',3,...
% 'MarkerFaceColor',cols(4,:));
% hold on

j=4;
tau=logspace(-14,-11,60);

% Fin=2*Uin/(pi*sigx*sigz)       %FIXED
% A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-9;
% A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-9;
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
% hold on
Fin=2*Uin/(pi*sigx*L(j)/2)       %VARIABLE
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*1e-9;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*1e-9;
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;

Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

h4a=semilogx(tau*1e12,Grad1,'Color',cols(2,:),'Linewidth',1);
hold on
h4b=semilogx(tau*1e12,Grad2,'Color',cols(2,:),'Linewidth',1);
hold on

% Grad=(Grad1+Grad2)/2;
% dGrad=(Grad2-Grad1)/2;
% h4=errorbar(tau*1e12,Grad,dGrad,'LineStyle','none',...
%     'Marker','v','Color',cols(2,:),'Linewidth',1,...
%     'MarkerSize',3,'MarkerFaceColor',cols(2,:));

ax1=gca;
%set(gca,'xscale','log','Color','none');
set(gca,'Color','none');

%axis tight
ylim([0 .7])
xlim([.1 5])

set(ax1,'XTick',[.1:.1:1 1.5 2:5])
set(ax1,'YTick',[.1:.1:.7])
set(ax1,'XTickLabel',{'0.1','0.2','0.3','','0.5','','','','','1.0','1.5','2','3','4','5'})

hx=xlabel('FWHM pulse width, \tau_p [ps]');
hy=ylabel('Average gradient, {\itG}_0 [MeV/m]');

ax2=copyobj(ax1,hFig);
delete(get(ax2,'Children'));
set(ax2,'Color','None','Box','off','Ygrid','on','Xgrid','on');
set(ax2,'Xcolor',[.9 .9 .9],'Ycolor',[.9 .9 .9]);%,'XTickLabel',[],'YTickLabel',[]);
%get(ax3)    %Unsure how to remove labels from this third axes.. will have
%to do by hand
uistack(ax2,'bottom');

linkaxes([ax1,ax2],'xy')


%hL=legend([h1,h2,h3,h4],{'{\it L} = 0.12 mm','{\it L} = 0.24 mm','{\it L} = 0.48 mm','{\it L} = 0.96 mm'});
hL=legend([h1a,h2a,h3a,h4a],{'{\it L} = 0.15 mm','{\it L} = 0.31 mm','{\it L} = 0.62 mm','{\it L} = 1.24 mm'});
set(hL,'Location','northeast','Box','Off')

set([ax1,ax2],'FontSize', 10);
set([hx,hy,hL],'FontSize', 10);

export_fig LengthOptim.eps -cmyk -r300 -painters%-painters


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

%close all

width=3.25;
height=2.5;
hFig = figure(97);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
%set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')


L=[150 310 620 1240]*um;
tau=logspace(-14,-11,100);

L=620*um;
tau=logspace(-14,-11,1000);





scale=1e-3;
j=1;
%Fin=2*Uin/(pi*sigx*L(j)/2)
Fin=2*Uin/(pi*sigx*sigz)
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;

Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

E1=Grad1*L(j);
E2=Grad2*L(j);

h1a=semilogx(tau*1e12,E1,'Color',cols(3,:),'Linewidth',1);
hold on
h1b=semilogx(tau*1e12,E2,'Color',cols(3,:),'Linewidth',1);%,'LineStyle','--');

%%
hold on

j=2;
%Fin=2*Uin/(pi*sigx*L(j)/2)
Fin=2*Uin/(pi*sigx*sigz)
for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;

A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;
Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

E1=Grad1*L(j);
E2=Grad2*L(j);

h2a=semilogx(tau*1e12,E1,'Color',cols(1,:),'Linewidth',1);
hold on
h2b=semilogx(tau*1e12,E2,'Color',cols(1,:),'Linewidth',1);%,'LineStyle','--');

hold on

j=3;
%Fin=2*Uin/(pi*sigx*L(j)/2)
Fin=2*Uin/(pi*sigx*sigz)

for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;
A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;
Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

E1=Grad1*L(j);
E2=Grad2*L(j);

h3a=semilogx(tau*1e12,E1,'Color',cols(4,:),'Linewidth',1);
hold on
h3b=semilogx(tau*1e12,E2,'Color',cols(4,:),'Linewidth',1);%,'LineStyle','--');


hold on

j=4;
Fin=2*Uin/(pi*sigx*L(j)/2)
%Fin=2*Uin/(pi*sigx*sigz)

for i=1:length(tau)
Pulse=@(t) sech(2*asech(1/sqrt(2))*(t)./tau(i)).^2;

A1=(2*eta1*fA1/n)*sqrt(2*Fin*c/eps0)*scale;
A2=(2*eta2*fA2/n)*sqrt(2*Fin*c/eps0)*scale;
Grad1(i)=A1*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
Grad2(i)=A2*quad(Pulse,0,L(j)/2/c)/L(j)/sqrt(tau(i));
end

E1=Grad1*L(j);
E2=Grad2*L(j);

h4a=semilogx(tau*1e12,E1,'Color',cols(2,:),'Linewidth',1);
hold on
h4b=semilogx(tau*1e12,E2,'Color',cols(2,:),'Linewidth',1)%,'LineStyle','--');

uistack(h2b, 'top');
uistack(h2a, 'top');
uistack(h1b, 'top');
uistack(h1a, 'top');


% j=5;
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
% h4b=semilogx(tau*1e12,E2,'Color',cols(2,:),'Linewidth',1,'LineStyle','--');



ax1=gca;
set(gca,'xscale','log','Color','none');

%axis tight
ylim([0 225])
xlim([.1 5])

hx=xlabel('FWHM pulse width, \tau_p [ps]');
set(ax1,'XTick',[.1:.1:1 1.25 1.5 2:5]);
set(ax1,'XTickLabel',{'0.1','0.2','0.3','','0.5','','','','','1.0','','1.5','2','3','4','5'});

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

hL=legend([h1a,h2a,h3a,h4a],{'{\it L} = 0.15 mm','{\it L} = 0.31 mm','{\it L} = 0.62 mm','{\it L} = 1.24 mm'});
set(hL,'Location','northwest','Box','Off')

set([ax1,ax2],'FontSize', 10);
set([hx,hy,hL],'FontSize', 10);

export_fig LengthOptimE.eps -cmyk -r300 -painters%-painters

%%  Optimal energy gain vs L

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
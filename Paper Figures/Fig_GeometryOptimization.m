close all
%clear all
clc
c       = 299792458;  %m/s
mu0     = 4*pi*1e-7;

folder='V:\ARDB\E163\Grating Structures\Lumerical\';
%folder='F:\Lumerical\';

%filename='Lum200g325wHscan';           % 200nm gap, 20nm radius, 5nm mesh%
%filename='Lum200g700hWscan';           % 200nm gap, 20nm radius, 5nm mesh

%Gap scan
%filename='Lum475w675hGscan';           % Ideal structure, Fig 2.3

%Width scans
%filename='Lum200g700hWscan';           % Ideal structure
%filename='Lum400g700hWscan';           % Target structure
%filename='Lum800g700hWscan';           % Plan B structure

%filename='Lum400g900hWscan';           % Target structure

%Height scans
%filename='Lum200g475wHscan';           % Ideal structure
%filename='Lum400g325wHscan';           % Target structure
%filename='Lum800g325wHscan';           % Plan B structure

%filename='Lum400g575wHscan';           % Target structure

%Misalignment scans
%filename='Lum200g670h450wSscan';
filename='Lum400g700h325wSscan';
%filename='Lum800g700h325wSscan';

%filename='Lum400g900h325wSscan';
%filename='Lum800g1075h325wSscan';


load([folder,filename]);


%%
lambda  = 800e-9;
betaf   = .9999637;
period  = x(end)-x(1);
tau     = lambda/c;
dt      = t(2)-t(1);
dy      = model.res;

Ntx=round(period/(betaf*c*dt));      % time steps to traverse 1 period
Ntopt=round(tau/dt);                 % time steps in one optical cyle

%Normalize fields:
trans=2*(1)/(1+model.n);
%Ex2=Ex2./.942./t;
%Ey2=Ey2./.942./t;
%Hz2=Hz2./.942./t;

Nx=length(x);
Ny=length(y);
Nt=length(t);
N=size(EX,3);
G=zeros(1,N);
D=zeros(1,N);

EmaxT=zeros(N,2);
EmaxB=zeros(N,2);
EmaxV=zeros(N,2);

Ex=zeros(Nx,Nt);
Ey=zeros(Nx,Nt);
Hz=zeros(Nx,Nt);

close all
PLOTS=0;
nplot=18;

height  = model.height;
gap     = model.gap;
width   = model.width; 

%param=model.height;
%param=model.width;
%param=model.gap
%param=25e-9:25e-9:1250e-9; %for filename='Lum200g475wHscan';   
param=model.shift;

for i=1:N
%for i=nplot:nplot

    %height  = model.height(i);
    %gap     = model.gap(i);
    %width   = model.width(i);  

    Xcent=round(Nx/2);
    Xedge=round(Nx/2-width/2/dy);
    
    Ex=squeeze(EX(:,:,i));
    Ey=squeeze(EY(:,:,i));
    Hz=squeeze(HZ(:,:,i));
    E_longt=flipud(Ex);
    E_trans=flipud(Ey);
    H_vert=flipud(Hz);
    
    Def=E_trans+betaf*c*mu0*H_vert;
    
    Emax=EMax(:,:,i);
    
    %EmaxV=max(max(Emax));
    
    
    num=1;              % number of optical cycles to average
    
    E_longt_e = zeros(Ntx,num*Ntopt);
    Def_e=zeros(Ntx,num*Ntopt);
    
    % Determine Emax value
    
    % Look at the bottom grating tooth
    Ytop=round(Ny/2-gap/dy/2);
    Ysub=round(Ny/2-gap/2/dy-height/dy);
    Emax1=0;
    for s=Ysub:Ytop+1
        if max(Emax(Xedge:Xedge+width/dy+1,s))>Emax1
            [Emax1, xmax1]=max(Emax(Xedge:Xedge+width/dy+1,s));
            MaxIndB=[Xedge+xmax1-1, s];
        end
    end
    Emax2=max(Emax(Xedge:Xedge+width/dy+1,MaxIndB(2)-1));
    Emax3=max(Emax(Xedge:Xedge+width/dy+1,MaxIndB(2)+1));
    
    if s==Ytop+1
        EmaxB(i,1)=mean([Emax1, Emax2]);
        EmaxB(i,2)=std([Emax1, Emax2]);
        indB=-1;
    else
        %Xmax=[xmax1, xmax2, xmax3];
        %Ymax=[MaxInd, MaxInd-1,MaxInd+1];
        [EmaxBot, indB]=sort([Emax2, Emax1, Emax3]);
        indB=indB-2;
        EmaxB(i,1)=mean([EmaxBot(2) EmaxBot(3)]);
        EmaxB(i,2)=std([EmaxBot(2) EmaxBot(3)]);
    end
    EmaxV(i,:)=EmaxB(i,:);
    indMax=indB;
    MaxInd=MaxIndB;
    
    % Look into top grating tooth too
    Ybot=round(Ny/2+gap/dy/2);
    Ysub=round(Ny/2+gap/2/dy+height/dy);
    Emax1=0;
    for s=Ybot-1:Ysub
        if max(Emax(Xedge:Xedge+width/dy+1,s))>Emax1
            [Emax1, xmax1]=max(Emax(Xedge:Xedge+width/dy+1,s));
            MaxIndT=[Xedge+xmax1-1, s];
        end
    end
    Emax2=max(Emax(Xedge:Xedge+width/dy+1,MaxIndT(2)-1));
    Emax3=max(Emax(Xedge:Xedge+width/dy+1,MaxIndT(2)+1));
    
    if s==Ybot-1
        EmaxT(i,1)=mean([Emax1, Emax3]);
        EmaxT(i,2)=std([Emax1, Emax3]);
        indT=1;
    else
        [EmaxTop, indT]=sort([Emax2, Emax1, Emax3]);
        indT=indT-2;
        EmaxT(i,1)=mean([EmaxTop(2) EmaxTop(3)]);
        EmaxT(i,2)=std([EmaxTop(2) EmaxTop(3)]);
    end
    
    maxTop=0;
    if EmaxT(i,1)>EmaxV(i,1)
        EmaxV(i,:)=EmaxT(i,:);
        indMax=indT;
        maxTop=1;
        MaxInd=MaxIndT;
        fprintf('Max is on top grating for step %d \n',i)
    end

    %% Plot of Maximum Field
    if PLOTS && i==nplot
        figure(10)
        subplot(2,2,[1 3])
                
        %pcolor(x*1e6,y*1e6,Emax')
        %axis equal
        %axis tight
        
        pcolor(Emax')
        %colormap(hot)
        ylabel('Transverse position [\Lambda]')
        xlabel('Long. pos. [\Lambda]')
        title('Max |E| Field [E_0]');
        shading interp
       
        hold on
        line([MaxInd(1) MaxInd(1)],[MaxInd(2)-.5 MaxInd(2)+.5],'LineWidth',3);
        hold on
        line([MaxInd(1)-1 MaxInd(1)+1],[MaxInd(2) MaxInd(2)],'LineWidth',3);
        
        %colorbar
        
        Cols='bbkrr';
        LinStyl=':---:';
        
        subplot(2,2,4)
        for k=-1:1:1
            plot(x*1e6,Emax(:,Ytop+k),'Color',Cols(k+3),'LineStyle',LinStyl(k+3))
            %plot(Y(y/4-5:3*y/4+5+1),Dam(MaxInd1(1)+k,y/4-5:3*y/4+5+1),'Color',Cols(k+3),'LineStyle',LinStyl(k+3))
            hold on
        end
        if maxTop==0
            hold on
            plot(x*1e6,Emax(:,MaxInd(2)),'--g')
            hold on
            plot(x*1e6,Emax(:,MaxInd(2)+indMax),':g')
            title(['Bottom tooth Emax = ',num2str(EmaxV(i,1)),'+-',num2str(EmaxV(i,2))])
        end
        
        axis tight
        ylabel('Maximum field [E_0]')
        xlabel('Position Along Grating Tooth [\Lambda]')
        %legend('Max-3nm','Max','Max+3nm');%,'Location','NorthwestOutside')
       
        
        Cols='rrkbb';
        subplot(2,2,2)
        for k=-1:1:1
            plot(x*1e6,Emax(:,Ytop+gap/dy+k),'Color',Cols(k+3),'LineStyle',LinStyl(k+3))
            hold on
        end
        if maxTop==1
            hold on
            plot(x*1e6,Emax(:,MaxInd(2)),'--g')
            hold on
            plot(x*1e6,Emax(:,MaxInd(2)+indMax),':g')
            title(['Top tooth Emax = ',num2str(EmaxV(i,1)),'+-',num2str(EmaxV(i,2))])
        
        end
        axis tight
        ylabel('Maximum field [E_0]')
        xlabel('Position Along Grating Tooth [\Lambda]')
        %legend('Max-3nm','Max','Max+3nm')
        %title(['Top tooth Surf, Emax = ',num2str(Emax2)])
        
        %% Other ramdon plots
        % figure
        % plot(E_longt(round(Nx/2),:))
        
        % figure
        % plot(y,Emax(round(Nx/2),:))
        %
        % %%
        % close all
        % figure
        % for k=1:length(t)
        %     plot(y*1e6,Ex2(round(Nx/2),:,k)')
        %     %pcolor(x*1e6,y*1e6,Dam(:,:,k)')
        %     drawnow;
        % end
        %
        % %% Plot Index map
        % pcolor(x*1e6,y*1e6,real(nx'))
        % shading interp
        
 
        %if PLOTS
        figure (2)
        subplot(2,2,1)
        %pcolor((1:Nx)/Nx,(1:Ntopt)/Ntopt,E_longt(:,end-Ntopt+1:end)')
        offset=0;%round(Ntopt/4)-1;
        pcolor((1:Nx)/Nx,(1:Ntopt)/Ntopt,E_longt(:,1+offset:Ntopt+offset)')
        ylabel('Time [\tau]')
        xlabel('Position Along Grating [\Lambda]')
        %title('Longitudinal field at channel center [E_{max}]');
        shading interp
        colorbar
        hold on
        axis xy
        line([0 1],[0 1],'LineWidth',2,'Color','w');
        
        subplot(2,2,3)
        %pcolor((1:Nx)/Nx,(1:Ntopt)/Ntopt,E_trans(:,1+offset:Ntopt+offset)'+betaf*c*H_vert(:,1+offset:Ntopt+offset)');
        pcolor((1:Nx)/Nx,(1:Ntopt)/Ntopt,Def(:,1+offset:Ntopt+offset)');
        ylabel('Time [\tau]')
        xlabel('Position Along Grating [\Lambda]')
        %title('Longitudinal field at channel center [E_{max}]');
        shading interp
        colorbar
        hold on
        axis xy
        line([0 1],[0 1],'LineWidth',2,'Color','w');
        
    end
    %% Transform onto electron frame
    for nt=1:Ntx
        
        xt = (nt-1)*(betaf*c)*dt;
        j = round(xt/dy);
        j = round((j/(Nx-1) - floor(j/(Nx-1)))*(Nx-1));
        %jj(nt)=j;
        %zz(nt)=z;
        
        %E_longt_e(nt,1:Ntopt) = circshift(E_longt(j+1,nt:nt+Ntopt-1)',nt-1)';
        %Def_e(nt,1:Ntopt) = circshift(Def(j+1,nt:nt+Ntopt-1)',nt-1)';
        E_longt_e(nt,1:Ntopt) = circshift(E_longt(j+1,nt:nt+Ntopt-1),1-nt)';
        Def_e(nt,1:Ntopt) = circshift(Def(j+1,nt:nt+Ntopt-1),1-nt)';
    end
    
    E_longt_e=fliplr(E_longt_e);
    Def_e=fliplr(Def_e);
    % E_trans_e=fliplr(E_trans_e);
    % B_vert_e=fliplr(B_vert_e);
    
    % Integrate and find max Gradient
    
    E_longt_av=sum(E_longt_e);
    E_longt_av=E_longt_av*(betaf*c*dt)/lambda;          %dx/lambda factor from integrated avg
    Def_av=sum(Def_e);
    Def_av=Def_av*(betaf*c*dt)/lambda;          %dx/lambda factor from integrated avg
    
    %E_longt_av=E_longt_av/fNORM;                       %normalize to incident field
    %G_lr=max(E_longt_av);                              %choose phase for max gradient
    %phi_max=find(E_longt_av==G_lr)
    %opt_phi=phi_max;                                   %comment this out once opt_phi found
    [G_lrmax, opt_phi]=max(E_longt_av);                     %choose phase for max gradient
    
    % sanity check
    %G_lr=sum(E_longt_e(:,opt_phi))*(beta*c*dt)/lambda
    %D_lr=sum(Def_e(:,opt_phi))*(beta*c*dt)/lambda
    G_lr=sum(E_longt_e)*(betaf*c*dt)/lambda;
    D_lr=sum(Def_e)*(betaf*c*dt)/lambda;
    
    format
    %lr_max=G_lr(opt_phi)
    G(i)=G_lr(opt_phi);
    format shorte
    %D_lr_max=D_lr(opt_phi)
    D(i)=D_lr(opt_phi);
    
    if PLOTS && i==nplot
        phi_max=opt_phi;%phi_max1;
        OptPhase=360*phi_max/length(E_longt_av);
        
        figure (2)
        subplot(2,2,2)
        pcolor((1:num*Ntopt)*2*180/Ntopt,(1:Ntx)/Ntx,E_longt_e)
        %title(['E_{longt,e} with \lambda =',num2str(lambda*1e9),'nm, \beta = ',num2str(beta)]);
        xlabel('Laser Phase [deg]')
        ylabel('Position Along Grating [\Lambda]')
        shading interp
        colorbar
        hold on
        axis xy
        line([OptPhase OptPhase],[0 1],'LineWidth',2,'Color','w');
        
        subplot(2,2,4)
        %pcolor((1:num*Ntopt)*2*180/Ntopt,(1:Ntx)/Ntx,E_trans_e+betaf*c*B_vert_e)
        pcolor((1:num*Ntopt)*2*180/Ntopt,(1:Ntx)/Ntx,Def_e)
        %title(['E_{longt,e} with \lambda =',num2str(lambda*1e9),'nm, \beta = ',num2str(beta)]);
        xlabel('Laser phase [deg]')
        ylabel('Position Along Grating [\Lambda]')
        shading interp
        colorbar
        hold on
        axis xy
        line([OptPhase OptPhase],[0 1],'LineWidth',2,'Color','w');
        
        %%
        figure (3)
        subplot(1,2,1)
        plot((1:num*Ntopt)*2*180/Ntopt,E_longt_av,'b','LineWidth',2);
        hold on
        plot((1:num*Ntopt)*2*180/Ntopt,Def_av,'r--','LineWidth',2);
        ylabel('Average Force [E_{max}]')
        xlabel('Laser Phase [deg]')
        grid on
        xlim([0 360])
        hold on
        line([OptPhase OptPhase],[-.25 .25],'LineWidth',2,'Color','k');
        ylim([-max(E_longt_av) max(E_longt_av)])
        
        subplot(1,2,2)
        
        plot((1:Ntx)/Ntx,E_longt_e(:,phi_max),'b','LineWidth',2);
        hold on
        plot((1:Ntx)/Ntx,Def_e(:,phi_max),'r--','LineWidth',2);
        xlabel('Position Along Grating [nm]')
        ylabel('Optimal Instantaneous Force [E_{max}]')
        grid on
        
    end
    
end
%%

cols(1,:)=[255 205 135]/255; %light Orange
cols(2,:)=[166 160 220]/255;  %light purple
cols(3,:)=[205 66 0]/255;      %dark orange
cols(4,:)=[47 08 117]/255;     %dark purple
close all

hFig=figure (1)
xwidth=3.5;%24;
ywidth=2;%5;

set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 xwidth ywidth])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

h1=line(param/lambda,G,'LineStyle','none','Marker','o',...
    'Color',cols(4,:),'Linewidth',1,...
    'MarkerSize',3,'MarkerFaceColor',cols(2,:));
ax1 = gca;
set(ax1,'Color','none','YColor',cols(4,:))
axis tight
%xlim([0 1.25])
%ylim([ 0 .52]);  %gap
%ylim([ 0 .42]);  % width
%ylim([.2 .38]);  %shift

ylimits = get(ax1,'YLim');
xlimits = get(ax1,'XLim');
yinc = (ylimits(2)-ylimits(1))/4;
set(ax1,'YTick',ylimits(1):yinc:ylimits(2));
%set(ax1,'XTick',0:.25:1.5);
%set(ax1,'XTick',[-.5:.25:.5]);

hy1=ylabel('Average gradient, {\it G}_0 [{\it E}_0]')
ax2 = axes('Position',get(ax1,'Position'),...
    'YAxisLocation','right',...
    'Color','none',...  %background color
    'YColor',cols(3,:));
h2=line(param/lambda,D*10,'LineStyle','none','Marker','o',...
    'Color',cols(3,:),'Linewidth',1,...
    'MarkerSize',3,'MarkerFaceColor',cols(1,:),...
    'Parent',ax2);

axis tight

set(ax2,'YLim',[-.2 .2])
set(ax2,'XLim',xlimits)
%set(ax2,'XLim',[0 1.5])
%ylim([-ymax ymax]);
ylimits = get(ax2,'YLim');
yinc = (ylimits(2)-ylimits(1))/4;
set(ax2,'YTick',ylimits(1):yinc:ylimits(2))
%set(ax2,'XTick',0:.25:1.5);
%set(ax2,'XTick',[-.5:.25:.5]);
hy2=ylabel('Avg. deflection, {\it D}_0 (x 10) [{\it E}_0]')
hx=xlabel('Gap Size [\lambda_g]')
%xlabel('Misalignment [\lambda_g]')
%grid on

ax3=copyobj(ax2,hFig);
delete(get(ax3,'Children'));
set(ax3,'Color','None','Box','off','Ygrid','on','Xgrid','on');
set(ax3,'Xcolor',[.9 .9 .9],'Ycolor',[.9 .9 .9],'XTickLabel',[],'YTickLabel',[]);
%get(ax3)    %Unsure how to remove labels from this third axes.. will have
%to do by hand
uistack(ax3,'bottom');
%linkaxes([ax1,ax2],'x')
set([ax1,ax2,ax3], 'FontSize', 10 );
set([hx,hy1,hy2], 'FontSize',10);

tightfig;
%export_fig G&DvsGap.eps -cmyk -r300 -painters



%close all
hFig=figure (2)
xwidth=3;%24;
ywidth=2.3;%5;

set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 xwidth ywidth])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')

scale=lambda;
%scale=1e-6;

h1=errorbar(param(2:end)/scale,EmaxV(2:end,1),EmaxV(2:end,2),...
    'LineStyle','none','Marker','o',...
    'Color',cols(4,:),'Linewidth',1,...
    'MarkerSize',2,'MarkerFaceColor',cols(2,:));
axis tight
ax1=gca;
set(ax1,'Color','none');
%set(ax1,'XTick',[0:.25:1.5]);   % gap (lambda)
%set(ax1,'XTick',[0:.25:1]);      % width (lambda)

%set(gca,'XTick',[-.5:.25:.5]);
hy=ylabel('{\it E}_{max} [{\it E}_0]');
hx=xlabel('Gap Size [\lambda_g]');
%hx=xlabel('Misalignment [\lambda_g]')

%ylim([1.7 3.3]);
%set(ax1,'YTick',[1.7:.4:3.3]);

ax3=copyobj(ax1,hFig);
delete(get(ax3,'Children'));
set(ax3,'Color','None','Box','off','Ygrid','on','Xgrid','on');
set(ax3,'Xcolor',[.9 .9 .9],'Ycolor',[.9 .9 .9]);%,'XTickLabel',[],'YTickLabel',[]);
% %get(ax3)    %Unsure how to remove labels from this third axes.. will have
% %to do by hand
uistack(ax3,'bottom');
% %linkaxes([ax1,ax2],'x')
set([ax1,ax3], 'FontSize', 10 );
set([hx,hy], 'FontSize',10);

linkaxes([ax1,ax3],'xy')
tightfig;

%export_fig EmaxvsGap.eps -cmyk -r300 -painters


% %%
% AF=G./EmaxV(:,1)';
% dAF=AF.*EmaxV(:,2)'./EmaxV(:,1)';
% 
% figure (3)
% errorbar(param*1e6,AF,dAF,'o')
% axis tight
% set(gca,'XTick',[0:.25:1.5]);
% %set(gca,'XTick',[-.5:.25:.5]);
% %ylim([ 0 .28]);
% %ylim([.08 .15]);
% ylabel('Acceleration factor, f_A')
% xlabel('Gap Size [\lambda_g]')
% %xlabel('Misalignment [\lambda_g]')
% grid on
% 
% %%
% DF=D./EmaxV(:,1)';
% dDF=DF.*EmaxV(:,2)'./EmaxV(:,1)';
% 
% figure (4)
% errorbar(param*1e6,10*DF,10*dDF,'o')
% axis tight
% set(gca,'XTick',[0:.25:1.5]);
% %set(gca,'XTick',[-.5:.25:.5]);
% %ylim([ 0 .28]);
% %ylim([.08 .15]);
% ylabel('Deflection factor, f_D')
% xlabel('Gap Size [\lambda_g]')
% %xlabel('Misalignment [\lambda_g]')
% grid on

%%
hFig=figure (5)
xwidth=3.5;%24;
ywidth=2;%5;

set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 xwidth ywidth])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

AF=G./EmaxV(:,1)';
dAF=AF.*EmaxV(:,2)'./EmaxV(:,1)';

DF=D./EmaxV(:,1)';
dDF=DF.*EmaxV(:,2)'./EmaxV(:,1)';

%errorbar(param*1e6,AF,dAF,'o')
scale=lambda;

h1=errorbar(param(2:end)/scale,AF(2:end),dAF(2:end),'LineStyle','none','Marker','o',...
    'Color',cols(4,:),'Linewidth',1,...
    'MarkerSize',2,'MarkerFaceColor',cols(2,:));
%h1=line(param/scale,G,'LineStyle','none','Marker','o',...
%    'Color',cols(4,:),'Linewidth',1,...
%    'MarkerSize',3,'MarkerFaceColor',cols(2,:));

ax1 = gca;
set(ax1,'Color','none','YColor',cols(4,:));
axis tight
%xlim([0 1.25])

%ylim([ 0 .28]); % gap
ylim([ 0 .24]); % width, 200 g
%ylim([ 0 .2]);  % width, 400g
%ylim([ 0 .12]);  % width, 800g

%ylim([.2 .38]);  %shift

ylimits = get(ax1,'YLim');
xlimits = get(ax1,'XLim');
yinc = (ylimits(2)-ylimits(1))/4;
set(ax1,'YTick',ylimits(1):yinc:ylimits(2));
%set(ax1,'XTick',0:.25:1.5);    %gap (lambda)
%set(ax1,'XTick',0:.25:1);       %width (lambda)
%set(ax1,'XTick',[-.5:.25:.5]);

hy1=ylabel('Acceleration factor, {\it f}_A');

ax2 = axes('Position',get(ax1,'Position'),...
    'YAxisLocation','right',...
    'Color','none',...  %background color
    'YColor',cols(3,:));

set(ax2,'YLim',[0 1])  %width 200
%set(ax2,'YLim',[0 .8])  %gap, width 400
%set(ax2,'YLim',[0 1.6])  %width 800

set(ax2,'XLim',xlimits)
%set(ax2,'XLim',[0 1.5])
%ylim([-ymax ymax]);
ylimits = get(ax2,'YLim');
yinc = (ylimits(2)-ylimits(1))/4;
set(ax2,'YTick',ylimits(1):yinc:ylimits(2))
%set(ax2,'XTick',0:.25:1.5); %gap (lambda)
%set(ax2,'XTick',0:.25:1);
%set(ax2,'XTick',[-.5:.25:.5]);
hy2=ylabel('Deflection factor, {\it -f}_D (x 100)');
%hx=xlabel('Gap Size [\lambda_g]');
hx=xlabel('Pillar width [\lambda_g]');

%xlabel('Misalignment [\lambda_g]')
%grid on

hold(ax2,'on')
h2=errorbar(param(2:end)/scale,-DF(2:end)*100,dDF(2:end)*100,'LineStyle','none','Marker','o',...
    'Color',cols(3,:),'Linewidth',1,...
    'MarkerSize',2,'MarkerFaceColor',cols(1,:),...
    'Parent',ax2);

%axis tight

ax3=copyobj(ax2,hFig);
delete(get(ax3,'Children'));
set(ax3,'Color','None','Box','off','Ygrid','on','Xgrid','on');
set(ax3,'Xcolor',[.9 .9 .9],'Ycolor',[.9 .9 .9]);%,'XTickLabel',[],'YTickLabel',[]);
%get(ax3)    %Unsure how to remove labels from this third axes.. will have
%to do by hand
uistack(ax3,'bottom');
%linkaxes([ax1,ax2],'x')
set([ax1,ax2,ax3], 'FontSize', 10 );
set([hx,hy1,hy2], 'FontSize',10);

tightfig;
%export_fig fa&fDvsGap.eps -cmyk -r300 -painters
%export_fig fa&fDvsWidth200.eps -cmyk -r300 -painters


%%

% AF1=AF;
% dAF1=dAF;
% DF1=DF;
% dDF1=dDF;

% AF2=AF;
% dAF2=dAF;
% DF2=DF;
% dDF2=dDF;

% AF3=AF;
% dAF3=dAF;
% DF3=DF;
% dDF3=dDF;

%%

% %%
% %close all
% hFig=figure (6)
% xwidth=3.4;%24;
% ywidth=2.5;%5;
% 
% set(hFig,'Units','inches')
% set(hFig, 'Position', [0 0 xwidth ywidth])
% 
% set(hFig, 'Color', 'w');
% set(hFig,'Units','points')
% 
% scale=1e-6;
% 
% h1=errorbar(param/scale,AF1,dAF1,...
%     'LineStyle','none','Marker','.',...
%     'Color',cols(1,:),'Linewidth',1,...
%     'MarkerSize',3);%,'MarkerFaceColor','k');
% hold on
% h1d=plot(param/scale,AF1,'k','Linewidth',1,...
%     'LineStyle','none','Marker','.',...
%     'MarkerSize',4);
% hold on
% h3=errorbar(param/scale,AF3,dAF3,...
%     'LineStyle','none','Marker','v',...
%     'Color',cols(3,:),'Linewidth',1,...
%     'MarkerSize',3,'MarkerFaceColor',cols(1,:));
% hold on
% h2=errorbar(param/scale,AF2,dAF2,...
%     'LineStyle','none','Marker','o',...
%     'Color',cols(4,:),'Linewidth',1,...
%     'MarkerSize',3,'MarkerFaceColor',cols(2,:));
% 
% axis tight
% ax1=gca;
% set(ax1,'Color','none');
% %set(ax1,'XTick',[0:.25:1.5]);   % gap (lambda)
% %set(ax1,'XTick',[0:.25:1]);      % width (lambda)
% %set(ax1,'XTick',[0:.1:1]);      % width (micron)
% %set(ax1,'XTick',[0:.2:1.25]);      % height (micron)
% set(ax1,'YTick',[0:.04:.24]);      % width (micron)
% set(ax1,'XTick',[-.4:.1:.4]);      % shift (micron)
% 
% %set(gca,'XTick',[-.5:.25:.5]);
% hy=ylabel('Acceleration factor, {\it f}_A');
% %hx=xlabel('Pillar width [\mum]');
% %hx=xlabel('Pillar height [\mum]');
% hx=xlabel('Misalignment [\mum]')
% 
% %ylim([1.7 3.3]);
% %set(ax1,'YTick',[1.7:.4:3.3]);
% 
% ax3=copyobj(ax1,hFig);
% delete(get(ax3,'Children'));
% set(ax3,'Color','None','Box','off','Ygrid','on','Xgrid','on');
% set(ax3,'Xcolor',[.9 .9 .9],'Ycolor',[.9 .9 .9]);%,'XTickLabel',[],'YTickLabel',[]);
% % %get(ax3)    %Unsure how to remove labels from this third axes.. will have
% % %to do by hand
% uistack(ax3,'bottom');
% % %linkaxes([ax1,ax2],'x')
% set([ax1,ax3], 'FontSize', 10 );
% set([hx,hy], 'FontSize',10);
% 
% linkaxes([ax1,ax3],'xy')
% tightfig;
% 
% export_fig fAvsShift.eps -cmyk -r300 -painters
% 
% %%
% hFig=figure (7)
% xwidth=3.4;%24;
% ywidth=2.5;%5;
% 
% set(hFig,'Units','inches')
% set(hFig, 'Position', [0 0 xwidth ywidth])
% 
% set(hFig, 'Color', 'w');
% set(hFig,'Units','points')
% 
% scale=1e-6;
% 
% h1=errorbar(param/scale,-DF1*100,dDF1*100,...
%     'LineStyle','none','Marker','.',...
%     'Color',cols(1,:),'Linewidth',1,...
%     'MarkerSize',3);%,'MarkerFaceColor','k');
% hold on
% h1d=plot(param/scale,-DF1*100,'k','Linewidth',1,...
%     'LineStyle','none','Marker','.',...
%     'MarkerSize',4);
% hold on
% h3=errorbar(param/scale,-DF3*100,dDF3*100,...
%     'LineStyle','none','Marker','v',...
%     'Color',cols(3,:),'Linewidth',1,...
%     'MarkerSize',3,'MarkerFaceColor',cols(1,:));
% hold on
% h2=errorbar(param/scale,-DF2*100,dDF2*100,...
%     'LineStyle','none','Marker','o',...
%     'Color',cols(4,:),'Linewidth',1,...
%     'MarkerSize',3,'MarkerFaceColor',cols(2,:));
% 
% axis tight
% 
% ax1=gca;
% set(ax1,'Color','none');
% %set(ax1,'XTick',[0:.25:1.5]);   % gap (lambda)
% %set(ax1,'XTick',[0:.25:1]);      % width (lambda)
% %set(ax1,'XTick',[0:.1:1]);      % width (micron)
% %set(ax1,'YTick',[-.2:.2:1.4]);      % width (micron)
% %set(ax1,'XTick',[0:.2:1.25]);      % height (micron)
% %set(ax1,'YTick',[-.8:.2:1.4]);      % width (micron)
% set(ax1,'XTick',[-.4:.1:.4]);      % shift (micron)
% set(ax1,'YTick',[.2:.1:.8]);      % width (micron)
% 
% %set(gca,'XTick',[-.5:.25:.5]);
% hy=ylabel('Deflection factor, {\it -f}_D (X100)');
% %hx=xlabel('Pillar width [\mum]');
% %hx=xlabel('Pillar height [\mum]');
% hx=xlabel('Misalignment [\mum]')
% 
% %ylim([1.7 3.3]);
% %set(ax1,'YTick',[1.7:.4:3.3]);
% 
% ax3=copyobj(ax1,hFig);
% delete(get(ax3,'Children'));
% set(ax3,'Color','None','Box','off','Ygrid','on','Xgrid','on');
% set(ax3,'Xcolor',[.9 .9 .9],'Ycolor',[.9 .9 .9]);%,'XTickLabel',[],'YTickLabel',[]);
% % %get(ax3)    %Unsure how to remove labels from this third axes.. will have
% % %to do by hand
% uistack(ax3,'bottom');
% % %linkaxes([ax1,ax2],'x')
% set([ax1,ax3], 'FontSize', 10 );
% set([hx,hy], 'FontSize',10);
% 
% linkaxes([ax1,ax3],'xy')
% tightfig;
% 
% export_fig fDvsShift.eps -cmyk -r300 -painters
% %
% %  save('fAfDWidthScan.mat','param','AF1','dAF1','AF2','dAF2','AF3','dAF3',...
% %      'DF1','dDF1','DF2','dDF2','DF3','dDF3');
% 
%%
clc
[fA, ind]=min(AF);
fA
%dAF(ind)
%DF(ind)
eta=EmaxV(ind,1);
eta
param(ind)*1e6
% % 
% %%
% AF(param==400e-9)
% dAF(param==400e-9)
%%
[val, ind]=max(EmaxV(:,1));
val+EmaxV(ind,2)
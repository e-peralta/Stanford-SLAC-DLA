close all
clear all

% CHECK
% figure
% imagesc(cdata)
% map=colormap;
%
% beam=rgb2ind(cdata,map);
% figure
% image(beam)
% colormap(hot)
%
img=imread('~/Dropbox/Research/Data/IRspotsizeIP.png'); %15
img2=imread('~/Dropbox/Research/Data/IRspotsizeIP_colorbar.png');
map=flipud(squeeze(im2double(img2)));

Vshift=0;


%cols(0,:)=[1 1 1];  %use 'w'
cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple
%cols(5,:)=[0 0 0]; %use 'k'


N=102;
col1=colormap(morgenstemning(N,'invert',1));
col3=colormap(morgenstemning(3*N,'invert',1));
col0=[col1(1:round(N/8),:);col3(round(3*N/8)+1:round(3*N*7/8),:);col1(round(7*N/8)+1:end,:)];
col0size=size(col0)

figure
image(img)
%map=colormap;
colormap(map)
colorbar
%shading interp

%
beam=rgb2ind(img,map);
beam=beam-min(min(beam));
beam=beam(:,1:48);
beam=beam';

figure
image(beam)
colormap(col0);

bad(1,:)= [36*ones(1,5) ...
    37 37 37 ...
    38*ones(1,4) ...
    39 39 39 39 ...
    40 40 40 40 ...
    41*ones(1,8) ...
    42*ones(1,4) ...
    43 43 45 ...
    18 21 21 23 26 26 29 30 33 33 35];
    
bad(2,:)= [159 162 226 229 253 ...
    156 165 187 ...
    175 193 199 283 ...
    173 184 208 310 ...
250 262 266 366 ...
145 152 164 219 224 245 273 302 ...
213 241 277 300 ...
96 158 141 ...
59 55 62 59 55 61 52 344 256 310 304];

for i=1:length(bad(1,:))
    beam(bad(1,i),bad(2,i))=0;
end


figure
image(beam)
colormap(col0);

    %%
[Nx Ny]=size(beam);
x=-115+(0:Nx-1)*230/47;
y=-600+(0:Ny-1)*1200/369;


VERBOSE=1;
PLOTEACHIMG=0;
% SIG_OPTIONS=optimset('Display','off',...
%     'TolFun',1e-05,...
%     'TolX',1e-05);
showplots=1;


CalX=1;% Sony, Before 121217
CalY=1;% 


limx=[x(1) x(end)];
limy=[y(1) y(end)];

cXi=0;
cYi=0;


img=beam;

ImageAnalyze_E; %gives c1,c2,hor,ver,sigx,sigy

% add ROI start point to centroid
cxni = cXi+c1(2);
cyni = cYi+c2(2);

width=7;
height=4;

figure (1)
hFig = gcf;
%set(hFig,'PaperPositionMode','auto')
%set(hFig,'ActivePositionProperty','position')
set(hFig,'ActivePositionProperty','outerposition')
%set(hFig,'Units','centimeters')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(hFig, 'PaperUnits','inches', 'PaperPosition',[0 0 width height])
%set(gcf, 'PaperPosition',[0 0 width height])
set(hFig, 'Color', 'w');
set(hFig,'Units','points')

if(VERBOSE)
    
    %subplot(3,3,[1 5])
    subplot(2,4,1:3)
    
    imagesc(y,x,img)
    axis xy
    %colormap(morgenstemning(256,'invert',1))
    colormap(col0)
    caxis([0 510])
    %caxis([1 280])
    %xlim(limx)
    %ylim(limy)
    set(gca,'XTickLabel',[],'YtickLabel',[]);
    
    %subplot(5,4,[13 18]);
    subplot(2,4,5:7);
    size(y)
    size(1:length(hor))
    
    h1=plot(y,f(c1,1:length(hor)),'LineWidth',2,'Color',cols(1,:));
    
    hold on;
    line([y(round(c1(2)-c1(3)*sqrt(2*log(2)))) y(round(c1(2)+c1(5)*sqrt(2*log(2))))],...
        [c1(4)+c1(1)/2 c1(4)+c1(1)/2], 'LineWidth',2,'Color',cols(3,:),'LineStyle',':');
    hold on
    h2=plot(y,hor,'Color',cols(4,:),'LineWidth',1);
    axis tight
    xlim(limy)
    
    %grid on;
    %title(sprintf('Line Out X. sx = %f',sigx))
    %title({sprintf('FWHM_x = %f [pix] %f [um]',sigx/CalX ,sigx),sprintf('Centroid at %f [pix] %f [um]',cxni,cxni*CalX)});
    ax1=gca;
    set(ax1,'YAxisLocation','right');
    set(ax1,'YtickLabel',[]);
    hx1=xlabel('{\it Z } [\mum]');
    %hL=legend([h2,h1],{'Projected Profile','Gauss Fit'});
    %set(hL, 'FontSize', 6);
    %legend boxoff
    
    
    %subplot(5,4,[3 11]);
    subplot(2,4,4);
    
    %sumy=sum(f(c2,1:length(ver)));
    plot(f(c2,1:length(ver)),x,'g','LineWidth',2,'Color',cols(1,:));
    hold on;
    %line([c2(2)-c2(3)*sqrt(2*log(2)) c2(2)+c2(5)*sqrt(2*log(2))],[c2(4)+c2(1)/2 c2(4)+c2(1)/2],'LineWidth',2,'Color','k');
    line([c2(4)+c2(1)/2 c2(4)+c2(1)/2],[x(round(c2(2)-c2(3)*sqrt(2*log(2)))) ...
        x(round(c2(2)+c2(5)*sqrt(2*log(2))))],'LineWidth',2,'Color',cols(3,:),'LineStyle',':');
    hold on
    plot(ver,x,'Color',cols(4,:),'LineWidth',1);
    axis tight
    ylim(limx)
    ax2=gca;
    set(ax2,'XtickLabel',[],'Xdir','reverse');
    hy1=ylabel('{\it X } [\mum]');
    set(ax2,'YAxisLocation','right');
    
    %grid on;
    
    %title(sprintf('Line Out Y. sy = %f',sigy))
    
    %title({sprintf('FWHM_y = %f [pix] %f
    %[um]',sigy/CalY,sigy),sprintf('Centroid at %f [pix] %f
    %[um]',cyni,cyni*CalY)})
    %ylabel('Y line out')
    
    
    [RN CN]=size(img);
    V=reshape(img,RN*CN,1);
    
    %subplot(5,4,[15 19]);
    subplot(2,4,8);
    [nb,xb]=hist(double(V),256);
    semilogy(xb,nb,'k.','Color',cols(3,:))
    %xlim([0,256]);
    axis tight
    %ylim([0 10000]);
    ax3=gca;
    set(ax3,'YAxisLocation','right');
    %grid on;
    hx2=xlabel('Pixel Value');
    hy2=ylabel('Pixel Count');
    set(ax3,'Ytick',[10 100 1000 10000],'TickLength'  , [.06 .06]);
    set(ax3,'Xtick',[170 340 510]);
    
    %set(ax3,'YGrid','off');
    
    set([hx1,hx2,hy1,hy2], 'FontSize' , 12 );
    
    set([ax1,ax2,ax3], 'FontSize', 10  );
end

%%
hold on

%tightfig;

dX=(x(2)-x(1));
dY=(y(2)-y(1));
FWHMx=(c1(3)+c1(5))*dX*sqrt(2*log(2));
FWHMy=(c2(3)+c2(5))*dY*sqrt(2*log(2));

if c1std(3)>c1std(5)
    FWHMerx=c1std(3)*dX*sqrt(2*log(2));
else
    FWHMerx=c1std(5)*dX*sqrt(2*log(2));
end
if c2std(3)>c2std(5)
    FWHMery=c2std(3)*dY*sqrt(2*log(2));
else
    FWHMery=c2std(5)*dY*sqrt(2*log(2));
end

scale=2*sqrt(2*log(2));
fprintf('FWHM_x [um] = %f +- %f\n ',FWHMx,FWHMerx);
fprintf('FWHM_y [um] = %f +- %f\n',FWHMy,FWHMery);
fprintf('OR\n sig_x [um] = %f +- %f\n ',FWHMx/scale,FWHMerx/scale);
fprintf('sig_y [um] = %f +- %f\n',FWHMy/scale,FWHMery/scale);

%spaceplots;

%export_fig IRspotsizeIP2.eps -cmyk -r300 -painters
%export_fig eBeamScreen.eps -rgb -r300 %-painters
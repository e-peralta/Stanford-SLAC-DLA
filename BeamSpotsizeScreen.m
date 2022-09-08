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
%img=imread('/Users/Milkshakes/Dropbox/Research/Data/ebeam_121025-2zoom.png');
img=imread('/Users/Milkshakes/Dropbox/Research/Data/ebeam_121025-100x100_16bit.png');
%img=imread('/Users/Milkshakes/Dropbox/Research/Data/ebeam_130222_80x50_8bit.png');


figure
image(img)
map=colormap;

% For spotsize measurement
beam=rgb2ind(img,map);
figure
image(beam)
colormap(hot)

% figure
% plot(sum(beam'))
%%

% Diverging Colorbrewer Colors (B&W and Colorblind safe)
col1=[230 97 1]/255;
col2=[253 184 99]/255;
col3=[178 171 210]/255;
col4=[94 60 153]/255;

[Ny Nx]=size(beam);

VERBOSE=1;
PLOTEACHIMG=0;
% SIG_OPTIONS=optimset('Display','off',...
%     'TolFun',1e-05,...
%     'TolX',1e-05);
showplots=0;


CalX=3.9218;% Sony, Before 121217
CalY=3.9218;% 

%CalX=7.7836;% Xybion, After 121217
%CalY=7.7836;% 

% x=linspace(2,148,Nx)*CalX; %For first file
% y=linspace(5,135,Ny)*CalY;
% limx=[30 120]*CalX;
% limy=[20 110]*CalY;

x=linspace(1,100,Nx)*CalX; %For first file
y=linspace(1,100,Ny)*CalY;


% x=linspace(1,80,Nx)*CalX;
% y=linspace(1,50,Ny)*CalY;
limx=[1 x(end)];
limy=[1 y(end)];

cXi=0;
cYi=0;

% filepath='V:\ARDB\E163\Data\140115\';
% img=imread([folder,'2.tif']);

%filepath=['V:\ARDB\E163\DATA\',datestr(now,'yymmdd'),'\'];

%[filename,filepath] = uigetfile({'*.tif;*.jpg;*.png'},'Choose file to process',filepath);

%img=imread([filepath filename]);


img=beam;

ImageAnalyze_E; %gives c1,c2,hor,ver,sigx,sigy

% add ROI start point to centroid
cxni = cXi+c1(2);
cyni = cYi+c2(2);

width=4.5;
height=4.5;

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
    
    %subplot(5,4,[1 10])
    subplot(3,3,[1 5])
    
    imagesc(x,y,img)
    axis xy
    colormap(morgenstemning(256,'invert',1))
    
    caxis([10 68])
    %caxis([1 280])
    xlim(limx)
    ylim(limy)
    set(gca,'XTickLabel',[],'YtickLabel',[]);
    
    %subplot(5,4,[13 18]);
    subplot(3,3,7:8);
    
    h1=plot(x,f(c1,1:length(hor)),'LineWidth',2,'Color',col2);
    
    hold on;
    line([x(round(c1(2)-c1(3)*sqrt(2*log(2)))) x(round(c1(2)+c1(5)*sqrt(2*log(2))))],...
        [c1(4)+c1(1)/2 c1(4)+c1(1)/2], 'LineWidth',2,'Color',col1,'LineStyle',':');
    hold on
    h2=plot(x,hor,'.','Color',col4);
    axis tight
    xlim(limx)
    
    %grid on;
    %title(sprintf('Line Out X. sx = %f',sigx))
    %title({sprintf('FWHM_x = %f [pix] %f [um]',sigx/CalX ,sigx),sprintf('Centroid at %f [pix] %f [um]',cxni,cxni*CalX)});
    ax1=gca;
    set(ax1,'YAxisLocation','right');
    set(ax1,'YtickLabel',[]);
    hx1=xlabel('{\it X } [\mum]');
    %hL=legend([h2,h1],{'Projected Profile','Gauss Fit'});
    %set(hL, 'FontSize', 6);
    %legend boxoff
    
    
    %subplot(5,4,[3 11]);
    subplot(3,3,[3 6]);
    
    %sumy=sum(f(c2,1:length(ver)));
    plot(f(c2,1:length(ver)),y,'g','LineWidth',2,'Color',col2);
    hold on;
    %line([c2(2)-c2(3)*sqrt(2*log(2)) c2(2)+c2(5)*sqrt(2*log(2))],[c2(4)+c2(1)/2 c2(4)+c2(1)/2],'LineWidth',2,'Color','k');
    line([c2(4)+c2(1)/2 c2(4)+c2(1)/2],[y(round(c2(2)-c2(3)*sqrt(2*log(2)))) ...
        y(round(c2(2)+c2(5)*sqrt(2*log(2))))],'LineWidth',2,'Color',col1,'LineStyle',':');
    hold on
    plot(ver,y,'.','Color',col4);
    axis tight
    ylim(limy)
    ax2=gca;
    set(ax2,'XtickLabel',[],'Xdir','reverse');
    hy1=ylabel('{\it Y } [\mum]');
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
    subplot(3,3,9);
    [nb,xb]=hist(double(V),256);
    semilogy(xb,nb,'k.','Color',col1)
    xlim([0,256]);
    axis tight
    ylim([0 10000]);
    ax3=gca;
    set(ax3,'YAxisLocation','right');
    %grid on;
    hx2=xlabel('Pixel Value');
    hy2=ylabel('Pixel Count');
    set(ax3,'Ytick',[10 100 1000 10000],'TickLength'  , [.06 .06]);
    %set(ax3,'YGrid','off');
    
    set([hx1,hx2,hy1,hy2], 'FontSize' , 12 );
    
    set([ax1,ax2,ax3], 'FontSize', 10  );
end

%%
hold on

tightfig;

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

%export_fig eBeamScreen.eps -cmyk -r300 %-painters
export_fig eBeamScreen.eps -rgb -r300 %-painters
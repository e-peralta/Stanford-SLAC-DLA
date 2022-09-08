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
img=imread('/Users/Milkshakes/Dropbox/Research/Data/SpatialOverlap_130222-OTR.png'); %15
img2=imread('/Users/Milkshakes/Dropbox/Research/Data/SpatialOverlap_130222-IR.png');

Vshift=0;

cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple


N=102;
col1=colormap(morgenstemning(N,'invert',1));
col3=colormap(morgenstemning(3*N,'invert',1));
col0=[col1(1:round(N/8),:);col3(round(3*N/8)+1:round(3*N*7/8),:);col1(round(7*N/8)+1:end,:)];
col0size=size(col0)

figure
image(img)
map=colormap;
shading interp
beam=rgb2ind(img,map);
beam=beam-min(min(beam));
[y x]=size(beam);
%x=339+x*4/15; for 0215 OTR
%y=172+y*20/76;%
x=275+(1:x)/4;
y=188+(1:y)*20/63;

% Fix the missing lines in OTR image
L=4; 
for i=1:L
    N=35;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
    N=92;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
    N=111;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
end

L=3; 
for i=1:L
    N=42;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
    N=48;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
    N=54;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
    N=61;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
    N=67;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
    N=73;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
    N=80;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
    N=86;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
    N=99;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
    N=105;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
    N=118;   beam(N+i,:)=((L+1-i)*beam(N,:)+i*beam(N+L+1,:))./(L+1);
end


figure
image(img2)
map=colormap;
beam2=rgb2ind(img2,map);
beam2=beam2-min(min(beam2));


[y2 x2]=size(beam2);
x2=257+(1:x2)*4/15;
y2=125+(1:y2)/3-Vshift;

pixCal=7.7836;% Xybion, After 121217

%lim=[276 374]*pixCal;
x=x*pixCal;
y=y*pixCal;
x2=x2*pixCal;
y2=y2*pixCal;

lim=[276 350]*pixCal;

width=4; %6.5" max
height=4;

%%
figure (1)
hFig = gcf;
%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
%set(hFig,'Units','centimeters')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

subplot(6,4,1:8)
imagesc(x,y,beam)
axis xy
%colormap(morgenstemning(64,'invert',1))
colormap(col0)
caxis([15 63])
caxis
xlim(lim)
set(gca,'XTickLabel',[]);
hy1=ylabel('{\it Y } [\mum]');
ax1=gca;

subplot(6,4,9:16)
imagesc(x2,y2,beam2)
axis xy
%colormap(morgenstemning(64,'invert',1))
colormap(col0)
caxis
caxis([15 63])
xlim(lim)
set(gca,'XTickLabel',[]);
hy2=ylabel('{\it Y } [\mum]');
ax2=gca;


trace=sum(beam)-min(sum(beam));
trace2=sum(beam2)-min(sum(beam2));
%%
subplot(6,4,17:24)
%h1=plot(x,trace/max(trace),'Color',col2);
%hold on
%h2=plot(x2,trace2/max(trace2),'Color',col4);

img=beam;

ImageAnalyze_E;

peak=x(round(c1(2)));
%figure

h1=plot(x,hor,'.','LineWidth',1,'Color',cols(1,:));
hold on
h1b=plot(x,f(c1,1:length(hor)),'LineWidth',1,'Color',cols(3,:));
hold on;
line([peak peak],[0 c1(1)+c1(4)], 'LineWidth',2,'Color',cols(3,:),'LineStyle',':');


dX=(x(2)-x(1));
FWHMx=(c1(3)+c1(5))*dX*sqrt(2*log(2))

img=beam2;
ImageAnalyze_E;
peak2=x2(round(c1(2)));

hold on

h2=plot(x2,hor,'.','LineWidth',1,'Color',cols(2,:));
hold on
h2b=plot(x2,f(c1,1:length(hor)),'LineWidth',1,'Color',cols(4,:));
hold on;
line([peak2 peak2],[0 c1(1)+c1(4)], 'LineWidth',2,'Color',cols(4,:),'LineStyle',':');

dX=(x2(2)-x2(1));
FWHMx2=(c1(3)+c1(5))*dX*sqrt(2*log(2))

xlim(lim)
ylim([5 40])
ax3=gca;

hx3=xlabel('{\it X } [\mum]');
hy3=ylabel('Intensity [a.u.]');

set([hx3,hy1,hy2,hy3], 'FontSize' , 10 );
set([ax1,ax2,ax3], 'FontSize', 10  );

hL=legend([h1b,h2b],{'OTR','IR'});
set(hL, 'FontSize', 10 ,'Box','Off' );


export_fig SpatialOverlap.eps -cmyk -r300 %-painters


% figure
% plot(sum(beam'),'r')
% figure
% [nb,xb]=imhist(beam);
% semilogy(xb,nb,'k')
% xlim([0,64]);
% 
% hold on
% [nb2,xb2]=imhist(beam2);
% semilogy(xb2,nb2,'b')
% xlim([0,64]);
% 

%%
% N=0;
% i=1;
% img2(1,:)=img(1,:);
% 
% img=beam;
% [y x]=size(img);
% 
% S=sum(beam');
% while(i<y-1)
%     if S(i+1)==S(i)
%     
%         %img2(i+1,:)=[];
%         N=N+1;
%         i
%     else
%         img2(i+1-N,:)=img(i+1,:);
%     end
%     i=i+1;
% end
% 
% % [y x]=size(img2);
% % % img2=img;
% % for i=2:y-1
% %     if sum(img2(i,:))<sum(img2(i-1,:))/2
% %         img2(i,:)=(img2(i-1,:)+img2(i+1,:))/2;
% %     end
% % end


%export_fig eBeamScreen.eps -cmyk -r300 %-painters
%export_fig eBeamScreen.eps -rgb -r300 %-painters
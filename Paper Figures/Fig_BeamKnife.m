clear all
close all

% Improved Morgenstemning colormap
N=102;
col1=colormap(morgenstemning(N,'invert',1));
col3=colormap(morgenstemning(3*N,'invert',1));
col0=[col1(1:round(N/8),:);col3(round(3*N/8)+1:round(3*N*7/8),:);col1(round(7*N/8)+1:end,:)];
col0size=size(col0)

gauss_f1 = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)).*(x<=c(2)) ...
    + c(1)*exp(-(x-c(2)).^2./(2*c(4)^2)).*(x>c(2));

% Discrete colors
%cols(0,:)=[1 1 1];  %use 'w'
cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple
%cols(5,:)=[0 0 0]; %use 'k'
%%
[xx yy]=FasterRasterProcess('KnifeEdgeX_130222.dat');

%figure (2)
hFig = gcf;

% width=12;
% height=6;
% set(hFig,'ActivePositionProperty','position')
% set(hFig,'Units','centimeters')
width=7;
height=3;
set(hFig,'ActivePositionProperty','outerposition')
set(hFig,'Units','inches')

set(hFig, 'Position', [0 0 width height])
set(hFig, 'PaperUnits','inches', 'PaperPosition',[0 0 width height])

set(gcf, 'Color', 'w');
set(hFig,'Units','points')

%figure (2)
subplot(2,4,1:2)
[cout cout_ci  xx yy yfit fe] = KnifeEdgeMeas(xx,yy);

h1x=plot(xx,fe*max(yy)/max(fe),'Color',cols(2,:),'LineWidth',2);

hold on
h2x=plot(xx,yy,'Color',cols(1,:),'LineWidth',1);

hold on
h3x=plot(xx,yfit,':','Color',cols(4,:),'LineWidth',2);

axis tight
xlim([50 250])
hx1=xlabel('{\it X} [\mum]');
hy1=ylabel('Integ. signal [a.u]');
ax1=gca;
%%
FWHMx=(cout(3)+cout(4))*sqrt(2*log(2));
if cout_ci(3)>cout_ci(4)
    FWHMerx=cout_ci(3)*sqrt(2*log(2));
else
    FWHMerx=cout_ci(4)*sqrt(2*log(2));
end
%drawnow

[xx2 yy2]=FasterRasterProcess('KnifeEdgeY_130222-2.dat');

xx2=xx2(.225*end:end);
yy2=yy2(.225*end:end);

%figure (2)
subplot(2,4,5:6)
[cout2 cout_ci2 xx2 yy2 yfit2 fe2] = KnifeEdgeMeas(xx2,yy2);

h1y=plot(xx2,fe2*max(yy2)/max(fe2),'Color',cols(2,:),'LineWidth',2);


%hold on
%h1a=plot(xx,circshift(awgn(fe*A/max(fe),20),-50),'Color',cols(4,:));

hold on
h2y=plot(xx2,yy2,'Color',cols(1,:),'LineWidth',1);

hold on
h3y=plot(xx2,yfit2,':','Color',cols(4,:),'LineWidth',2);

axis tight
xlim([60 range(xx2)])
hx2=xlabel('{\it Y} [\mum]');
hy2=ylabel('Integ. signal [a.u]');
ax2=gca;

FWHMy=(cout2(3)+cout2(4))*sqrt(2*log(2));
if cout_ci2(3)>cout_ci2(4)
    FWHMery=cout_ci2(3)*sqrt(2*log(2));
else
    FWHMery=cout_ci2(4)*sqrt(2*log(2));
end
%drawnow
%%
fex=gauss_f1([1,150,cout(3),cout(4)],1:300);
fey=gauss_f1([1,150,cout2(3),cout2(4)],1:300);
beam=fey'*fex;
subplot(2,4,[3 8])
imagesc(beam)
shading interp
axis xy
xlim([50 250])
ylim([100 200])

caxis([0 .98])
ax3=gca;
%axis equal
%colormap(morgenstemning(256,'invert',1))
colormap(col0)
%grid on
set(gca,'YAxisLocation','right');
ht3=title('Reconstructed profile');
hx3=xlabel('{\it X} [\mum]');
hy3=ylabel('{\it Y} [\mum]');
%%
%drawnow
%tightfig;
%drawnow


set([hx1,hx2,hy1,hy2,hx3,hy3,ht3], 'FontSize' , 10 );
    
set([ax1,ax2,ax3], 'FontSize', 10  );
%%
scale=2*sqrt(2*log(2));
fprintf('FWHM_x [um] = %f +- %f\n ',FWHMx,FWHMerx);
fprintf('FWHM_y [um] = %f +- %f\n',FWHMy,FWHMery);
fprintf('OR\n sig_x [um] = %f +- %f\n ',FWHMx/scale,FWHMerx/scale);
fprintf('sig_y [um] = %f +- %f\n',FWHMy/scale,FWHMery/scale);
 
 
export_fig eBeamKnife.eps -cmyk -r300 -painters

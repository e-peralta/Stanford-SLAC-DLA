function [cout cout_ci xx yy yfit fe] = KnifeEdgeMeas(xx,yy)
%% KnifeEdgeMeas.m

%
% Script calculates the beam FWHM from a knife edge scan
%
% RUN FasterRasterProcess first to generate vroll and data variables

% clear all
% [vroll data]=FasterRasterProcess('KnifeEdgeX_130222.dat');
% 
% close all
% figure
%plot(vroll,data)


%xx=vroll;
%yy=data;


% normalize Y
yy = yy - min(yy);

% normalize X
xx = xx - min(xx);

% sort 1
[xx xi] = sort(xx);
yy = yy(xi);

% find which direction is scanned
if mean(yy(1:10)) > mean(yy(end-10:end))
    stdy = mean(yy(1:10))/4;
    
    % shift plot for constant aligment
    xt = min(-xx);
    xmean = mean(xx);
    xx = -xx - xt;
    shifted = 1;
else
    stdy = mean(yy(end-10:end))/4;
    shifted = 0;
end

% sort #2
[xx xi] = sort(xx);
yy = yy(xi);

% smooth function
N = round(length(yy)/100);
ysm = conv(yy,ones(N,1)/N,'same');
xsm = conv(xx,ones(N,1)/N,'same');

% find outliers
rih = find(yy > ysm + stdy);
ril = find(yy < ysm - stdy);

xxo = xx([rih; ril]);
yyo = yy([rih; ril]);

% remove outliers
xx([rih; ril]) = [];
yy([rih; ril]) = [];
ysm(1:N/2) = []; ysm(end-N/2:end) = [];
xsm(1:N/2) = []; xsm(end-N/2:end) = [];

% if shifted
%     xx = xx-mean(xx)+xmean;
% end

% smooth function - round 2
N = round(length(yy)/100);
ysm = conv(yy,ones(N,1)/N,'same');
xsm = conv(xx,ones(N,1)/N,'same');
ysm(1:N/2) = []; ysm(end-N/2:end) = [];
xsm(1:N/2) = []; xsm(end-N/2:end) = [];

% To visualize the smoother data
% hold on
% plot(xsm,ysm,'r')

%initial guess
yoff = mean(ysm(1:N*2));
xoff = xsm(find(ysm > max(ysm)/2,1,'first'));
xs1 = xsm(find(ysm > max(ysm)*.1353,1,'first'));
xs2 = xsm(find(ysm > max(ysm)*.8647,1,'last'));
sig = abs(xs1+xs2)-2*xoff;
A = max(ysm);%/(0.5*sqrt(pi/2)*sig)
A=max(yy);

%define gaussian function with identical 
gauss_f1 = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)).*(x<=c(2)) ...
    + c(1)*exp(-(x-c(2)).^2./(2*c(4)^2)).*(x>c(2))+c(5);
%Define the error function as a gaussian integral
%erf_f2=@(c,x) cumsum((c(1)*exp(-(x-c(2)).^2./(2*c(3)^2))+c(5)).*(x<=c(2)) ...
%    + (c(1)*exp(-(x-c(2)).^2./(2*c(4)^2))+c(5)).*(x>c(2)))/length(x)*10;
erf_f2=@(c,x) cumsum(gauss_f1(c,x))/length(x)*10;
cin2=[A xoff sig sig yoff];
[cout,~,~,exitflag,~,~,~,cout_ci]=lsqcurvefitstd(erf_f2,cin2,xx,yy);


fe=gauss_f1(cout,xx);
yfit=erf_f2(cout,xx);

% h1=plot(xx,fe*A/max(fe),'Color',cols(2,:),'LineWidth',2);
% 
% 
% %hold on
% %h1a=plot(xx,circshift(awgn(fe*A/max(fe),20),-50),'Color',cols(4,:));
% 
% hold on
% h2=plot(xx,yy,'Color',cols(1,:));
% 
% hold on
% h3=plot(xx,yfit,':','Color',cols(4,:),'LineWidth',2);
% 
% axis tight

%hl=legend([h2,h3,h1],{'Scan Data','Erf Fit','Ext. Gaus'});
%set(hl, 'FontSize', 6  );

% hold on
% plot(xx,cumsum(fe)/length(xx)*10,'k:')

% FWHM=(cout(3)+cout(4))*sqrt(2*log(2))
% if cout_ci(3)>cout_ci(4)
%     FWHMer=cout_ci(3)*sqrt(2*log(2))
% else
%     FWHMer=cout_ci(4)*sqrt(2*log(2))
% end

%% Compare with Kens fit
%[cout2 cout_ci2 cout1 cout_ci1 shifted] = KnifeEdgeFit(vroll,data,1)
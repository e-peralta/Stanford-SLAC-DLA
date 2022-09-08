function [cout2 cout_ci2 cout1 cout_ci1 shifted] = KnifeEdgeFit(xx,yy,plot_on)
% ex: [cout cout_ci] = KnifeEdgeFit(pos,energy,[plot?])
% fits a set of data points to an erf function
% also returns the confidence interval
% cout2 = asummetric erf
% cout1 = symmetric erf

% kensoong@slac.stanford.edu 14 - 04 - 30

%%

if nargin < 3
    plot_on = 0;
end


% r1 = load('V:\ARDB\E163\Data\140616\knife_edge_01.dat');
% 
% %find which axis is scanned
% if std(r1(:,1)) > std(r1(:,2))
%     xi = 1;
% else
%     xi = 2;
% end
% 
% %define variables
% xx = r1(:,xi)/1e3; % x in mm
% yy = r1(:,5);


% normalize Y
yy = yy - min(yy);

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

% smooth function - round 2
N = round(length(yy)/100);
ysm = conv(yy,ones(N,1)/N,'same');
xsm = conv(xx,ones(N,1)/N,'same');
ysm(1:N/2) = []; ysm(end-N/2:end) = [];
xsm(1:N/2) = []; xsm(end-N/2:end) = [];

%{
hold on
plot(xx,yy,'r');
plot(xxo,yyo,'o');
plot(xsm,ysm)
%}

% init guess
yoff = mean(ysm(1:N*2));
xoff = xsm(find(ysm > max(ysm)/2,1,'first'));
xs1 = xsm(find(ysm > max(ysm)*.1353,1,'first'));
xs2 = xsm(find(ysm > max(ysm)*.8647,1,'last'));
sig = abs(xs1+xs2)-2*xoff;
A = max(ysm)/(0.5*sqrt(pi/2)*sig);

cin1 = [A sig xoff yoff];

% symmetric erf function
%xoff1 = xsm(find(ysm > max(ysm)/2,1,'first'));
%xoff2 = xsm(find(ysm < max(ysm)/2,1,'last'));
%lb = [0 (max(xsm)-min(xsm))*.01 xoff2 min(ysm)];
%ub = [max(ysm) max(xsm)-min(xsm) xoff1 max(ysm)];
erf_f1 = @(c,x) 0.5 * sqrt(pi/2) * c(1) * c(2) * erf(sqrt(2) * (x - c(3)) / c(2)) + c(4);
[cout1,~,~,~,~,~,~,cout_ci1]=lsqcurvefitstd(erf_f1,cin1,xsm,ysm);
%[cout1,~,~,~,~,~,~,cout_ci1]=lsqcurvefitstd(erf_f1,cin1,xsm,ysm,lb,ub);

% asymmetric erf function
lb = [0 (max(xsm)-min(xsm))*.01 cout1(3)-cout_ci1(3)*10 min(ysm) (max(xsm)-min(xsm))*.01];
ub = [max(ysm) max(xsm)-min(xsm) cout1(3)+cout_ci1(3)*10 max(ysm) max(xsm)-min(xsm)];
%lb = cout1-cout_ci1*50; lb(1) = cout1(1)*cout1(2)-cout_ci1(1)*cout1(2)*50;
%ub = cout1+cout_ci1*50; ub(1) = cout1(1)*cout1(2)+cout_ci1(1)*cout1(2)*50;

cin2 = [cout1(1)*cout1(2) cout1(2:end) cout1(2)];
erf_f2 = @(c,x) (0.5 * sqrt(pi/2) * c(1) * erf(sqrt(2) * (x - c(3)) / c(2)) + c(4)).*(x<c(3)) + ...
   (0.5 * sqrt(pi/2) * c(1) * erf(sqrt(2) * (x - c(3)) / c(5)) + c(4)).*(x>c(3)) ; 
[cout2,~,~,~,~,~,~,cout_ci2]=lsqcurvefitstd(erf_f2,cin2,xx,yy,lb,ub);

if shifted
    % restore original alignment
    cout1(3) = cout1(3)-mean(xx)+xmean;
    cout2(3) = cout2(3)-mean(xx)+xmean;
    xx = xx-mean(xx)+xmean;
end
    
%plot_on = 1;
if plot_on == 1
        
    figure()
    plot(xx,yy,'bx')
    hold on
    xplot = [min(xx):range(xx)/500:max(xx)];
    
    
    % symmetric
    plot(xplot,erf_f1(cout1,xplot),'g')
    %[ylb yub] = make_conf_intervals(erf_f1,cout1,cout_ci1*3,xx,1);
    %plot(xx,ylb,'g--')
    %plot(xx,yub,'g--')
    
    
    % asymmetric
    plot(xplot,erf_f2(cout2,xplot),'r')
  %  [ylb yub] = make_conf_intervals(erf_f2,cout2,cout_ci2*3,xx,1);
  %  plot(xx,ylb,'r--')
  %  plot(xx,yub,'r--')
    
    title(['\sigma = ' num2str((cout2(2)+cout2(5))/2) ' \pm ' num2str(mean([cout_ci2(2),cout_ci2(5)])) ' mm'])
    enhance_plot
    axis('tight')
end

% percent difference between symmetric and asymmetric
%((cout2(2)+cout2(5))-cout1(2)*2)/cout1(2)*100


function [cout cout_ci f] = double_lor(spec2,shift,cmin,cmax)

SIG_OPTIONS=optimset('Display','off','TolFun',1e-05,'TolX',1e-05);

% fitting function - double lorentzian with baseline
f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
    (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+...
    (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(6))).^2)).*(x < c(2)+c(8))+ ...
    (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(7))).^2)).*(x >= c(2)+c(8))+c(9);

[g1,g2]=max(spec2);  %faster, no averaging
g3=sum(spec2)/g1/sqrt(2*log(2));
g4=g3;
g5=g1;
g6=g4;
g7=g6;
g8=10;
g9=0;
g10=0;

guess=[g1,g2,g3,g4,g5,g6,g7,g8,g9];
if nargin<4
    cmin=[min(spec2) 1 10 10 min(spec2) 10 10 0 -2*abs(min(spec2))];
    cmax=[2*max(spec2) length(spec2) length(spec2) length(spec2) 2*max(spec2) length(spec2) length(spec2) 1024 2*max(spec2)];
end
% this is the fit routine, the last output is the standard deviation
[cout,~,~,~,~,~,~,cout_ci]=lsqcurvefitstd(f,guess,1:length(spec2),spec2,cmin,cmax,SIG_OPTIONS);

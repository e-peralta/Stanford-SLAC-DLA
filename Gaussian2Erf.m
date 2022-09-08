%% AnalyicalMode.m
%
%   Script to calculate broadening of electron enegry spectrum due to laser
%   acceleration when the electron pulse is longer than a small fraction of
%   the optical pulse length
%
% --------------------------------

%clear all
close all
xwidth=8.9;             % Figure width
ywidth=4.5;             % Figure height

Amp=50;                % MeV/m

DETAILS=0;
SMALL=0;                % Whether to apply small amplitude broadening

C= 299792458;
wave=800e-9;
tau=wave/C*1e15;        % single cycle in fs (~2.7 fs)

wp = 8.23;              % leading edge HWHM width in keV
wn = 15.4;              % trailing edge HWHM width in keV

figure
line([-wn wp],[.5 .5],'LineWidth',1,'Color','r');

dE=.01;                 % Energy resolution

Elim=(wp+wn)+2*Amp;     % Energy limits

E1 = [-Elim:dE:Elim];   %

N= 2*Amp;               % Number of time slices to use in analysis.


wp = wp/sqrt(2*log(2)); % initial distribution leading edge RMS width in keV
wn = wn/sqrt(2*log(2)); % initial distribution trailing edge RMS width in keV


f1 = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)).*(x<=c(2)) ...
    + c(1)*exp(-(x-c(2)).^2./(2*c(4)^2)).*(x>c(2));
f2=@(c,x) cumsum(c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)).*(x<=c(2)) ...
    + c(1)*exp(-(x-c(2)).^2./(2*c(4)^2)).*(x>c(2)))/length(x)*10;
f2b=@(c,x) cumsum(c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)).*(x<=c(2)) ...
    + c(1)*exp(-(x-c(2)).^2./(2*c(4)^2)).*(x>c(2)))/length(x)*10+c(5);

s1 = f1([1,0,wn,wp],E1);
%s2 = f2([1,0,wn,wp],E1);
s2 = f2b([1,0,wn,wp,0],E1);
s2=awgn(s2,30);

SIG_OPTIONS=optimset('Display','off','TolFun',1e-05,'TolX',1e-05);
cin2=[1,0,wn,wp,0];
[cout2,~,~,exitflag,~,~,~,cout_ci2]=lsqcurvefitstd(f2b,cin2,E1,s2,[],[],SIG_OPTIONS);

exitflag
%w0=wp*sqrt(2*log(2));  %convert back to HWHM from RMS
hold on
h1=plot(E1,s1,'b','LineWidth',3);
hold on
h1=plot(E1,s2,'b:','LineWidth',3);

%See difference in scan direction
hold on
h1=plot(E1,cumsum(s1)/length(E1)*10,'r','LineWidth',3);
% hold on
% h1=plot(E1,cumsum(fliplr(s1))/length(E1)*10,'r','LineWidth',3);


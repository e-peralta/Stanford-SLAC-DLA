function [cout, coutstd, yfit, ChiSq, ybkgd, ymain, ysignal, FOM, modeList]=FitSpectrum5c(ydata,modef,shift,tmod,ploton,list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function [cout yfit ChiSq ybkgd ymain ysignal FOM] = 
%                      FitSpectrum5(ydata,modef,shift,tmod,ploton)
%
%  This function fits the acquired E163 spectra.
%
% % %  Input arguments:
%
%    ydata = data to fit to
%    modef = selects one of the predefined fitting functions:
%       0: Single Landau
%       1: Single Asymmetric Gaussian
%       2: Single Asymmetric Lorentzian
%       3: Single Asymmetric Sech^2
%       4: Single Mixed Asymmetric 
%       5: Double Gaussian, fixed energy shift
%       6: Double Gaussian, fitting energy shift
%       7: Asymmetric Gaussian and small Gaussian, fixed energy shift
%       8: Asymmetric Gaussian and small Gaussian, fitting energy shift
%       9: Double Asymmetric Gaussian, fixed energy shift 
%      10: Double Asymetric Gaussian, fitting energy shift
%      11: Main Gaussian + Asymetric Gaussian, fixed energy shift
%      12: Main Gaussian + Asymetric Gaussian, fitting energy shift
%      13: Double Asymetric Gaussian, fitting energy shift, force both sigma- equal
%      14: Half gaussian with Assymetric gaussian, *** NEEDS WORK **
%      15: Double Assymetric Lorentzian, fixed energy shift 
%      16: Double Assymetric Lorentzian, fitting energy shift 
%      17: Triple Assymetric Lorentzian (double transmitted), fitting energy shifts  
%      18: Test Mode (used for tweaking algorithms) 
%    shift = expected shift between peaks (if modef > 4)
%    tmod = additional parameter used in FOM calculation (i.e. width of integration),
%             also expected width of a modulated transmission spectra (modef 17) 
%    ploton = displays data and fit when function is called 
%    list = whether or not to output the list of Fits available
%
% % % Output arguments
%
%   cout   = fit coefficients 
%   yfit   = fitted raw signal
%   ChiSq  =  the squared 2-norm of the residuals
%   ybkgd = subtracted background/baseline
%   ymain  = extracted main spectrum component
%   ysignal= extracted smaller component
%   FOM    = user defined FOM to calculate
%   modeList = list of available fits
%
%  121112, eperalta@slac.stanford.edu

clear guess
if nargin<6
    list=0;
    if nargin<5  
        ploton=0;
        if nargin<4  
            tmod=0;
            if nargin<3
                shift=0;
            end
        end
    end
end

if list
    modeList{1}='0: Landau';
    modeList{2}='1: Asym. Gaussian';
    modeList{3}='2: Asym. Lorentzian';
    modeList{4}='3: Asym. Hyperb Secant';
    modeList{5}='4: Mixed Asymmetric';
    modeList{6}='5: 1+1, fix sep';
    modeList{7}='6: 1+1, FIT sep';
    modeList{8}='7: 4+1, fix sep';
    modeList{9}='8: 4+1, FIT sep';
    modeList{10}='9: 4+4, fix sep';
    modeList{11}='10: 4+4, FIT sep';
    modeList{12}='11: 1+4, fix sep';
    modeList{13}='12: 1+4, FIT sep';
    modeList{14}='13: 10, sig+=sig-';
    modeList{15}='14: Needs work';
    modeList{16}='15: 3+3, fix sep';
    modeList{17}='16: 3+3, FIT sep';
    modeList{18}='17: 3+3+3, fix sep';
    modeList{19}='18: 3+3+3, FIT sep';
    modeList{20}='19: Test mode';
    cout=0;
    yfit=0;
    ChiSq=0;
    ybkgd=0;
    ymain=0;
    ysignal=0;
    FOM=0;
else
    modeList='Didnt ask for list';
    
[a b]=size(ydata);
if a>b
    ydata=ydata';
end
SIG_OPTIONS=optimset('Display','off','TolFun',1e-05,'TolX',1e-05);

%xdata=1:length(ydata);


z=ydata;
t=1:length(z);
Nthresh=30;  %threshold for large data analysis (truncates ends)
if length(z)<Nthresh
    buff1=0;
    buff2=0;
    avgback=0;
else
buff1=10;
buff2=10;
avgback=10;
end

xdata=1+buff1:length(ydata)-buff2;
xshift=buff1*ones(1,length(z));

ydata=ydata(xdata);

% %Original Background subtraction
% X1=min(xdata);
% X2=max(xdata);
% avgback=round(length(xdata)/20);
% Y1=mean(ydata(1:avgback));
% Y2=mean(ydata((length(xdata)-avgback):length(xdata)));
% ybkgd=((Y2-Y1)/(X2-X1)*(xdata-X1)+Y1);
% ydata=ydata-ybkgd;
% ybkgd=((Y2-Y1)/(X2-X1)*(t-xshift-X1)+Y1);

X1=min(xdata);
X2=max(xdata);
%avgback=round(length(xdata)/20);


Y1=mean(ydata(1:1+avgback));
Y2=mean(ydata((length(xdata)-avgback):length(xdata)));

minY=Y1;
if Y2<minY
    minY=Y2;
end
ybkgd=minY*ones(1,length(xdata));
ydata=ydata-ybkgd;
ybkgd=minY*ones(1,length(t));

%ydata2=[ydata,zeros(1,500)];
ydata2=ydata;

% gaussian
%f=inline('c(1)*exp(-(x-c(2)).^2/c(3)^2)+c(4)','c','x');
[g1,g2]=max(ydata);  %faster, no averaging
%g4=min(ydata);
%g3=(sum(ydata)-length(ydata)*g4)/g1/sqrt(2*log(2));

% [g1,g2]=peakfind(ydata,10);
 g3=sum(ydata)/g1/sqrt(2*log(2));

% Took this catch out EAP 12/12/12
%  if g3<0
%      fprintf('*** FitSpectrum error: data must have positive integral *** \n')
% %      cout=0;
% %      yfit=0;
% %      ChiSq=0;
% %      ybkgd=0;
% %      ymain=0;
% %      ysignal=0;
% %      FOM=0;
%  else
 


if modef==0
%Landau function  
f=@(c,x) c(1)*Landau((x-c(2))./c(3))+c(4);
guess=[g1,g2,g3,0];
cmin=[min(ydata) 1 1 -g1/10];
cmax=[2*max(ydata) length(ydata) length(ydata) g1/10];
elseif modef==1
%%asymmetric gaussian
f=@(c,x) (c(1)*exp(-.5*(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
         (c(1)*exp(-.5*(x-c(2)).^2/c(4)^2)).*(x >= c(2))+c(5);
g4=g3;
guess=[g1,g2,g3,g4,0];
cmin=[min(ydata) 1 1 1 -g1/10];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) g1/10];
elseif modef==2%14
% asymmetric lorentzian
f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
        (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+c(5);
g4=g3;
guess=[g1,g2,g3,g4,0];
cmin=[min(ydata) 1 1 1 -g1/10];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) g1/10];
elseif modef==3
% assymetric sech^2
f=@(c,x) (c(1)*sech((x-c(2))./c(3)).^2).*(x < c(2))+ ...
         (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2))+c(5);
g4=g3;
guess=[g1,g2,g3,g4,0];
cmin=[min(ydata) 1 1 1 -g1/10];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) g1/10];
elseif modef==4
% mixed asymetric 
f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
             (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2))+c(5);
g4=g3;
guess=[g1,g2,g3,g4,0];
cmin=[min(ydata) 1 1 1 -g1/10];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) g1/10];

elseif modef==5
%%double gaussian
f=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2)+ ...
          c(4)*exp(-(x-(c(2)+shift)).^2/c(5)^2);
f1=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
g4=g1/2;
g5=g3;
guess=[g1,g2,g3,g4,g5];
cmin=[min(ydata) 1 1 min(ydata)/2 1];
cmax=[2*max(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata)];

elseif modef==6
%%double gaussian with variable shift
f=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2)+ ...
          c(4)*exp(-(x-(c(2)+c(6))).^2/c(5)^2);
f1=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
g4=g1/2;
g5=g3;
g6=shift;
guess=[g1,g2,g3,g4,g5,g6];
cmin=[min(ydata) 1 1 min(ydata)/2 1 .75*shift];
cmax=[2*max(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) 1.25*shift];

elseif modef==7
%%asymetric gaussian and small gaussian    
f=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2))+...
          c(5)*exp(-(x-(c(2)+shift)).^2/c(6)^2);
f0=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
f1=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
g4=g3;
g5=g1/2;
g6=g4/2;
guess=[g1,g2,g3,g4,g5,g6];
cmin=[min(ydata) 1 1 1 min(ydata) 1];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata)];

elseif modef==8
%%asymetric gaussian and small gaussian with var shift
f=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2))+...
          c(5)*exp(-(x-(c(2)+c(7))).^2/c(6)^2);
f0=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
f1=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
g4=g3;
g5=g1/2;
g6=g4/2;
g7=shift;
guess=[g1,g2,g3,g4,g5,g6,g7];
cmin=[min(ydata) 1 1 1 min(ydata) 1 .75*shift];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) 1.25*shift];

elseif modef==9
%%double asymetric gaussian 
f=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2))+...
          (c(5)*exp(-(x-(c(2)+shift)).^2/c(6)^2)).*(x < c(2)+shift)+ ...
          (c(5)*exp(-(x-(c(2)+shift)).^2/c(7)^2)).*(x >= c(2)+shift);
f1=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
g4=g3;
g5=g1/2;
g6=g4/2;
g7=g6;
guess=[g1,g2,g3,g4,g5,g6,g7];
cmin=[min(ydata) 1 1 1 min(ydata) 1 1];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) length(ydata)];

elseif modef==10
%%double asymetric gaussian with variable shift
f=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2))+...
          (c(5)*exp(-(x-(c(2)+c(8))).^2/c(6)^2)).*(x < c(2)+c(8))+ ...
          (c(5)*exp(-(x-(c(2)+c(8))).^2/c(7)^2)).*(x >= c(2)+c(8));
f1=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
g4=g3;
g5=g1/2;
g6=g4/2;
g7=g6;
g8=shift;
guess=[g1,g2,g3,g4,g5,g6,g7,g8];
cmin=[min(ydata) 1 1 1 min(ydata) 1 1 .75*shift];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) length(ydata) 1.25*shift];

elseif modef==11
%%gaussian with asymetric gaussian with fixed shift
f=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2)+...
          (c(4)*exp(-(x-(c(2)+shift)).^2/c(5)^2)).*(x < c(2)+shift)+ ...
          (c(4)*exp(-(x-(c(2)+shift)).^2/c(6)^2)).*(x >= c(2)+shift);
f0=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
f1=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));

g4=g1/2;
g5=g3/2;
g6=g5;
guess=[g1,g2,g3,g4,g5,g6];
cmin=[min(ydata) 1 1 min(ydata) 1 1 ];
cmax=[2*max(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) length(ydata)];

elseif modef==12
%%gaussian with asymetric gaussian with variable shift
f=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2)+...
          (c(4)*exp(-(x-(c(2)+c(7))).^2/c(5)^2)).*(x < c(2)+c(7))+ ...
          (c(4)*exp(-(x-(c(2)+c(7))).^2/c(6)^2)).*(x >= c(2)+c(7));
f0=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
f1=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));

g4=g1/2;
g5=g3/2;
g6=g5;
g7=shift;
guess=[g1,g2,g3,g4,g5,g6,g7];
cmin=[min(ydata) 1 1 min(ydata) 1 1 .75*shift];
cmax=[2*max(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) length(ydata) 1.25*shift];

elseif modef==13
%%double asymetric gaussian with variable shift, both lowE sigmas equal
f=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2))+...
          (c(5)*exp(-(x-(c(2)+c(7))).^2/c(3)^2)).*(x < c(2)+c(7))+ ...
          (c(5)*exp(-(x-(c(2)+c(7))).^2/c(6)^2)).*(x >= c(2)+c(7));
f1=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
g4=g3;
g5=g1/2;
g6=g3;
g7=shift;
guess=[g1,g2,g3,g4,g5,g6,g7];
cmin=[min(ydata) 1 1 1 min(ydata) 1 .75*shift];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) 1.25*shift];

elseif modef==14    
%%half gaussian with assymetric gaussian, with variable shift, both lowE sigmas equal
f=@(c,x)  (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x >= c(2))+...
          (c(5)*exp(-(x-(c(2)+c(7))).^2/c(4)^2)).*(x < c(2)+c(7))+ ...
          (c(5)*exp(-(x-(c(2)+c(7))).^2/c(6)^2)).*(x >= c(2)+c(7));
f1=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
          (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
f0=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x >= c(2));
      
g4=g3;
g5=g1/2;
g6=g3;
g7=shift;
guess=[g1,g2,g3,g4,g5,g6,g7];
cmin=[min(ydata) 1 1 1 min(ydata) 1 .75*shift];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) 1.25*shift];

elseif modef==15
%%double asymmetric lorentzian, fixed shift
f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
        (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+...
        (c(5)*ones(size(x))./(1+((x-(c(2)+shift))./(0.5.*c(6))).^2)).*(x < c(2)+shift)+ ...
        (c(5)*ones(size(x))./(1+((x-(c(2)+shift))./(0.5.*c(7))).^2)).*(x >= c(2)+shift);

f1=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
          (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2));
g4=g3;
g5=g1/2;
g6=g4/2;
g7=g6;
guess=[g1,g2,g3,g4,g5,g6,g7];
cmin=[min(ydata) 1 1 1 min(ydata) 1 1];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) length(ydata)];
   
elseif modef==16
%double asymmetric lorentzian, fitting shift
f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
        (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+...
        (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(6))).^2)).*(x < c(2)+c(8))+ ...
        (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(7))).^2)).*(x >= c(2)+c(8));
f1=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
          (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2));
g4=g3;
g5=g1/2;
g6=g4/2;
g7=g6;
g8=shift;

guess=[g1,g2,g3,g4,g5,g6,g7,g8];
cmin=[min(ydata) 1 1 1 min(ydata) 1 1 .75*shift];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) length(ydata) 1.25*shift];
    
elseif modef==17
%triple asymmetric lorentzian, fixed  shift(two small ones close by)
%(c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+...
f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
        (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2))+...                
        (c(5)*ones(size(x))./(1+((x-(c(2)+shift-c(8)/2))./(0.5.*c(6))).^2)).*(x < c(2)+shift-c(8)/2)+ ...
        (c(5)*ones(size(x))./(1+((x-(c(2)+shift-c(8)/2))./(0.5.*c(7))).^2)).*(x >= c(2)+shift-c(8)/2)+...
        (c(11)*sech((x-(c(2)+shift+c(8)/2))./c(10)).^2).*(x >=c(2)+shift+c(8)/2);    
        %(c(11)*ones(size(x))./(1+((x-(c(2)+shift+c(8)/2))./(0.5.*c(9))).^2)).*(x < c(2)+shift+c(8)/2)+ ...        
f0=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
          (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2));
f1=@(c,x) (c(1)*ones(size(x))./(1+((x-(c(2)-c(5)/2))./(0.5.*c(3))).^2)).*(x < c(2)-c(5)/2)+ ...
          (c(1)*ones(size(x))./(1+((x-(c(2)-c(5)/2))./(0.5.*c(4))).^2)).*(x >= c(2)-c(5)/2)+ ...
          (c(8)*ones(size(x))./(1+((x-(c(2)+c(5)/2))./(0.5.*c(6))).^2)).*(x < c(2)+c(5)/2)+ ...
          (c(8)*sech((x-(c(2)+c(5)/2))./c(7)).^2).*(x >=c(2)+c(5)/2);        
          %(c(8)*ones(size(x))./(1+((x-(c(2)+c(5)/2))./(0.5.*c(7))).^2)).*(x >= c(2)+c(5)/2);    
g4=g3;
g5=g1/2;
g6=g4/2;
g7=g6;
g8=tmod;
g9=g7;
g10=g9;
g11=g5;

% 
% guess=[g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11];
% cmin=[min(ydata) 1 1 1 min(ydata) 1 1 0 1];
% cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) length(ydata) 100 length(ydata)];

guess=[g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11];
cmin=[0 1 1 1 0 10 10 0 10 10 .04];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) max(ydata)...
    length(ydata) length(ydata) 250 length(ydata) length(ydata) max(ydata)];

elseif modef==18
%triple asymmetric lorentzian, fitting shift (two small ones close by)
%(c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+...
f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
        (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2))+...        
        (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)-c(9)/2))./(0.5.*c(6))).^2)).*(x < c(2)+c(8)-c(9)/2)+ ...
        (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)-c(9)/2))./(0.5.*c(7))).^2)).*(x >= c(2)+c(8)-c(9)/2)+...
        (c(12)*ones(size(x))./(1+((x-(c(2)+c(8)+c(9)/2))./(0.5.*c(10))).^2)).*(x < c(2)+c(8)+c(9)/2)+ ...
        (c(12)*sech((x-(c(2)+c(8)+c(9)/2))./c(11)).^2).*(x >=c(2)+c(8)+c(9)/2);    
        %(c(12)*ones(size(x))./(1+((x-(c(2)+c(8)+c(9)/2))./(0.5.*c(11))).^2)).*(x >= c(2)+c(8)+c(9)/2);     
f0=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
          (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2));        
f1=@(c,x) (c(1)*ones(size(x))./(1+((x-(c(2)-c(5)/2))./(0.5.*c(3))).^2)).*(x < c(2)-c(5)/2)+ ...
          (c(1)*ones(size(x))./(1+((x-(c(2)-c(5)/2))./(0.5.*c(4))).^2)).*(x >= c(2)-c(5)/2)+ ...
          (c(8)*ones(size(x))./(1+((x-(c(2)+c(5)/2))./(0.5.*c(6))).^2)).*(x < c(2)+c(5)/2)+ ...
          (c(8)*sech((x-(c(2)+c(5)/2))./c(7)).^2).*(x >=c(2)+c(5)/2);        
          %(c(8)*ones(size(x))./(1+((x-(c(2)+c(5)/2))./(0.5.*c(7))).^2)).*(x >= c(2)+c(5)/2);    
g4=g3;
g5=g1/2;
g6=g4/2;
g7=g6;
g8=shift;
g9=tmod;
g10=g7;
g11=g10;
g12=g5;

guess=[g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12];
%cmin=[min(ydata) 1 1 1 min(ydata) 1 1 .75*shift 0 1 1];
cmin=[0 1 1 1 0 10 10 .75*shift 0 10 10 .04];
cmin=[0 1 1 1 0 10 10 .75*shift 0 7 7 .04];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) max(ydata)...
    length(ydata) length(ydata) 1.25*shift 250 length(ydata) length(ydata) max(ydata)];
elseif modef==19
% TEST MODE
% f=@(c,x) ((c(1)-c(9))*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)+c(9)*ones(size(x))).*(x < c(2))+ ...
%         (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+...
%         (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(6))).^2)).*(x < c(2)+c(8))+ ...
%         (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(7))).^2)).*(x >= c(2)+c(8));
% f0=@(c,x) ((c(1)-c(5))*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)+c(5)*ones(size(x))).*(x < c(2))+ ...
%           (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2));
% f1=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
%           (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2));

f=@(c,x) ((c(1)-c(9))*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)+c(9)*ones(size(x))).*(x < c(2))+ ...
        (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2))+...
        (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(6))).^2)).*(x < c(2)+c(8))+ ...
        (c(5)*sech((x-(c(2)+c(8)))./c(7)).^2).*(x >= c(2)+c(8));
     %(c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(7))).^2)).*(x >= c(2)+c(8));
f0=@(c,x) ((c(1)-c(5))*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)+c(5)*ones(size(x))).*(x < c(2))+ ...
          (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2));
f1=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
          (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2));          
       %(c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2));
g4=g3;
g5=g1/2;
g6=g4/2;
g7=g6;
g8=shift;
g9=.01;

guess=[g1,g2,g3,g4,g5,g6,g7,g8,g9];
%cmin=[min(ydata) 1 1 1 min(ydata) 1 1 .75*shift 0];
cmin=[min(ydata) 1 1 1 min(ydata) 1 1 .5*shift 0];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) length(ydata) 1.25*shift .02];
end


% xmin=xdata(1)+50;
% xmax=xdata(end)-50;

% xmin=xdata(1);
% xmax=xdata(end);
% y=ydata(xmin:xmax);

y=ydata2;
xmin=buff1;

% This is where the actual fitting takes place 
% (uses fminsearch if the optim. toolbox not available)
if exist('lsqcurvefit')==2
   %[c1, ChiSqX, ~, exitflag]=lsqcurvefit(f,guess,1:length(y),y,cmin,cmax,SIG_OPTIONS);
   
   %if modef==18
   %    [c1,ChiSqX,~,exitflag,~,~,~,c1std]=lsqcurvefitstd(@TripLorentz,guess,1:length(y),y,cmin,cmax,SIG_OPTIONS);       
   %else
       [c1,ChiSqX,~,exitflag,~,~,~,c1std]=lsqcurvefitstd(f,guess,1:length(y),y,cmin,cmax,SIG_OPTIONS);       
   %end
   c1(2)=c1(2)+xmin;
else
    fprintf('*** NOT USING OPTIM TOOLBOX .... NO IDEA if fits will be good*** \n')
    options = optimset('MaxFunEvals',1e6,'Display','off');
    k = @(c) sum((f(c,1:length(y))-y).^2)+sum(c<cmin)*1e9+sum(c>cmax)*1e9;
    %k = @(c) sum((f(c,1:length(y))-y).^2);
    c1=fminsearch(k,guess,options);
    ChiSqX=0;
end
%c1=abs(c1);    12/11/12 - Took this out (modef 17 cases could have c1(9)<0)

% if length(z)<Nthresh
%     xdata=linspace(xdata(1),xdata(end),100);
% end

if modef<1
    ymain=f(c1(1:4),xdata);
elseif modef<5
    ymain=f(c1(1:5),xdata);
elseif modef==5
    ymain=f1(c1(1:3),xdata);
    ysignal=f1([c1(4) c1(2)+shift c1(5)],xdata);    
elseif modef==6
    ymain=f1(c1(1:3),xdata);
    ysignal=f1([c1(4) c1(2)+c1(6) c1(5)],xdata);    
elseif modef==7
    ymain=f0(c1(1:4),xdata);
    ysignal=f1([c1(5) c1(2)+shift c1(6)],xdata);
elseif modef==8
    ymain=f0(c1(1:4),xdata);
    ysignal=f1([c1(5) c1(2)+c1(7) c1(6)],xdata);
elseif modef==9
    ymain=f1(c1(1:4),xdata);
    ysignal=f1([c1(5) c1(2)+shift c1(6) c1(7)],xdata);    
elseif modef==10
    ymain=f1(c1(1:4),xdata);
    ysignal=f1([c1(5) c1(2)+c1(8) c1(6) c1(7)],xdata);  
elseif modef==11
    ymain=f0(c1(1:3),xdata);
    ysignal=f1([c1(4) c1(2)+shift c1(5) c1(6)],xdata);  
elseif modef==12
    ymain=f0(c1(1:3),xdata);
    ysignal=f1([c1(4) c1(2)+c1(7) c1(5) c1(6)],xdata);  
elseif modef==13
    ymain=f1(c1(1:4),xdata);
    ysignal=f1([c1(5) c1(2)+c1(7) c1(3) c1(6)],xdata);  
elseif modef==14
    ymain=f0(c1(1:3),xdata);
    ysignal=f1([c1(5) c1(2)+c1(7) c1(4) c1(6)],xdata);  
elseif modef==15
    ymain=f1(c1(1:4),xdata);
    ysignal=f1([c1(5) c1(2)+shift c1(6) c1(7)],xdata);  
    FOM=0;
elseif modef==16
    ymain=f1(c1(1:4),xdata);
    ysignal=f1([c1(5) c1(2)+c1(8) c1(6) c1(7)],xdata);
    yint=ydata(round(c1(2))+shift:round(c1(2))+shift+tmod);
    FOM=sum(yint)/max(yint);   %FOM3
    
    %FOM=sum(y(c1(2)+shift:c1(2)+shift+tmod))/c1(5);   %FOM1
    %FOM=sum(y(c1(2)+shift-tmod:c1(2)+shift+tmod))/c1(5);  $FOM2
    % FOM=c1(5);
elseif modef==17
    ymain=f0(c1(1:4),xdata);
    ysignal=f1([c1(5) c1(2)+shift c1(6) c1(7) c1(8) c1(9) c1(10) c1(11)],xdata);
    yint=ydata(c1(2)+shift:c1(2)+shift+abs(tmod)+20);
    FOM=sum(yint)/max(yint);   %FOM3     
elseif modef==18
    ymain=f0(c1(1:4),xdata);
    ysignal=f1([c1(5) c1(2)+c1(8) c1(6) c1(7) c1(9) c1(10) c1(11) c1(12)],xdata);
%    yint=ydata(c1(2)+shift:c1(2)+shift+abs(tmod)+20);
%    FOM=sum(yint)/max(yint);   %FOM3    
elseif modef==19
    ymain=f0([c1(1:4) c1(9)],xdata);
    ysignal=f1([c1(5) c1(2)+c1(8) c1(6) c1(7)],xdata);
 %   FOM=exitflag;
end
if exist('FOM')==0
    FOM=0;
end
if exist('ysignal')==0
    ysignal=0;
end

for i=length(c1)+1:12
    c1(i)=0;
end

yfit= f(c1,xdata);
yfit=[zeros(1,buff1),yfit,zeros(1,buff2-1)];
ymain=[zeros(1,buff1),ymain,zeros(1,buff2-1)];
ysignal=[zeros(1,buff1),ysignal,zeros(1,buff2-1)];

%c1(2)=c1(2)+1;


if ploton
    %plot(ydata); hold on;
    %plot(xdata,f(c1,xdata),'r');
    x=xdata;
    y=ydata;
    %plot(x,y,'kx',x,y-yfit-ybkgd,'cx');
%     hold on
%     plot(x,yfit+ybkgd,'g',x,ymain+ybkgd,'r',x,ysignal+ybkgd,'b','LineWidth',2);
%     
    if modef<5
    plot(t,z,'kx');
    hold on
%     plot(x,ymain,'r',x,y-yfit,'b','LineWidth',2);
%     legend ('Raw data','Fit','Difference');
    %plot(x,ymain,'r',x,ybkgd,'b','LineWidth',2);
    plot(t,ymain+ybkgd,'r',t,ybkgd,'b','LineWidth',2);
    legend ('Raw data','Fit','Background');

    else
    plot(t,z,'kx',t,yfit+ybkgd,'g');
    hold on
    plot(t,ymain+ybkgd,'r',t,ysignal+ybkgd,'b','LineWidth',2);
%     hold on
%     plot(x,y-yfit,'cx');
    hold on
    plot(t,ybkgd,'c');
    legend ('Raw data','Fit','Main','Secondary');%,'Residuals');
    
    end
    axis tight
    
    
   switch modef
       case 0
           typestr='mode 0: Single Landau';
       case 1
           typestr='mode 1: Single Asymm Gaussian';
       case 2
           typestr='mode 2: Single Asymm. Lorentzian';
       case 3
           typestr='mode 3: Single Asymm. Sech^2';
       case 4
           typestr='mode 4: Single Mixed Asymmmetric';
       case 5
           typestr='mode 5: Double Gauss., fixed \Delta E';
       case 6
           typestr='mode 6: Double Gauss., fit \Delta E';
       case 7
           typestr='mode 7: Asymm. Gauss. + small Gauss., fixed \Delta E';
       case 8
           typestr='mode 8: Asymm. Gauss. + small Gauss., fit \Delta E';
       case 9
           typestr='mode 9: Double Asymm. Gauss., fixed \Delta E';
       case 10
           typestr='mode 10: Double Asymm. Gauss., fit \Delta E';
       case 11
           typestr='mode 11: Main Gauss. + Asymm. Gauss., fixed \Delta E';
       case 12
           typestr='mode 12: Main Gauss. + Asymm. Gauss., fit \Delta E';
       case 13
           typestr='mode 13: Double Asymm. Gauss., fit \Delta E, both \sigma- equal';
       case 14
           typestr='mode 14: Half gaussian with Assymetric gaussian, *** NEEDS WORK **';
       case 15
           typestr='mode 15: Double Asymm. Lorentzian, fixed \Delta E';
       case 16
           typestr='mode 16: Double Asymm. Lorentzian, fit \Delta E';
       case 17          
           typestr='mode 17: Triple Assym. Lorentz. (double trans.), fixed \Delta E';
       case 18          
           typestr='mode 18: Triple Assym. Lorentz. (double trans.), fit \Delta E';
       case 19          
           typestr='mode 19: TEST MODE';
   end
    
%     title({[typestr],...
%        ['FOM = ',num2str(FOM), '   ChiSqX =',num2str(ChiSq),'   shift =',num2str(cout(end))]});
    title(typestr);

    drawnow;
    hold off;
end

% if modef>4   %12/11/12 EAP , took this out since FOM now defined
% %FOM=(sum(ymain)-sum(ysignal))/1e8;    
% end
ChiSq=ChiSqX;
%ChiSq=sum((ydata-yfit).^2);%/(ydata.^2));
cout=c1;
coutstd=c1std;

end
end


function F = TripLorentz(c,x)
% F=(c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
%     (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+...
%     (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(6))).^2)).*(x < c(2)+c(8))+ ...
%     (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(7))).^2)).*(x >= c(2)+c(8))+...
%     (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)+c(9)))./(0.5.*c(7))).^2)).*(x < c(2)+c(8)+c(9))+ ...
%     (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)+c(9)))./(0.5.*c(10))).^2)).*(x >= c(2)+c(8)+c(9));

%Original Tripple lorentz fitting function
f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
        (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+...
        (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(6))).^2)).*(x < c(2)+c(8))+ ...
        (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(7))).^2)).*(x >= c(2)+c(8))+...
        (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)+c(9)))./(0.5.*c(7))).^2)).*(x < c(2)+c(8)+c(9))+ ...
        (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)+c(9)))./(0.5.*c(10))).^2)).*(x >= c(2)+c(8)+c(9));     
f0=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
          (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2));
f1=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
          (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+ ...
          (c(1)*ones(size(x))./(1+((x-(c(2)+c(5)))./(0.5.*c(4))).^2)).*(x < c(2)+c(5))+ ...
          (c(1)*ones(size(x))./(1+((x-(c(2)+c(5)))./(0.5.*c(6))).^2)).*(x >= c(2)+c(5));    
g4=g3;
g5=g1/2;
g6=g4/2;
g7=g6;
g8=shift;
g9=tmod;
g10=g6;

guess=[g1,g2,g3,g4,g5,g6,g7,g8,g9,g10];
cmin=[min(ydata) 1 1 1 min(ydata) 1 1 .75*shift 0 1];
cmax=[2*max(ydata) length(ydata) length(ydata) length(ydata) 2*max(ydata) length(ydata) length(ydata) 1.25*shift 100 length(ydata)];

end

function y=Landau(x)
%function y=Landau(x).m
%
% Calculates the landau function for collisional straggling:
% L(x)=Int[exp(-y.*(log(y)+x)).*sin(pi*y),y,0,infty]
%
% 12/05/22 - eperalta@slac.stanford.edu

% %% Generate Look-up table 
% Integ=@(x,y) exp(-y.*(log(y)+x)).*sin(pi*y);
% PLOTEACH=0;
% X=[-3.4:.1:-.5 -.48:.02:0 .1:.1:5 5.5:.5:100 105:10:195];
% Y=0:.01:350;
% cm=colormap(jet(length(X))); 
% 
% for i=1:length(X)
%     for j=1:length(Y)
%         f(j)=Integ(X(i),Y(j));
%         if isnan(f(j)) f(j)=0; end
%     end
%     L(i)=sum(f);
%     if PLOTEACH 
%         figure (99)
%         plot(Y,L./max(L),'Color',cm(i,:))
%         hold on
%     end
% end
% X=[-3.5 X 200];
% L=[0 L 0];
%%invert axis from energy loss to energy 
% X=-fliplr(X);
% L=fliplr(L)
% 
% plot(X,L)
% title('Generated Landau Function')
% save('V:\ARDB\E163\matlab\FittingFunctions\LandauTable.mat','X','L');

%% Look up Landau Table and interpolate values
%load('V:\ARDB\E163\matlab\FittingFunctions\LandauTable.mat');
folder='~/Documents/MATLAB/';
filename='LandauTable.mat';
eval(['load ',folder,filename]);

for i=1:length(x)
    if x(i) < -200 || x(i) > 3.5 
        y(i)=0;
    else
        y(i)=interp1(X,L,x(i));
    end
end

end
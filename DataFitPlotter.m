function [yfit, ymain, ysignal, ysig1, ysig2]=DataFitPlotter(cin,modef)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function [yfit, ymain, ysigna]=DataFitPlotter(cin,modef)
%
%  This function gives the fitting curves to acquired E163 spectra given
%  the extracted fitting parameters.
%
% % %  Input arguments:
%
%    cin   = extracted fit parameters (from FitSpectrum#.m)
%    modef = selects one of the predefined fitting functions:
%       0: Single Lorentzian
%       1: Single Gaussian
%       2: Single Landau
%       3: Single Asymmetric Lorentzian
%       4: Single Asymmetric Gaussian
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
%
% % % Output arguments
%
%   yfit   = fitted raw signal
%   ymain  = extracted main spectrum component
%   ysignal= extracted smaller component
%
%  130324, eperalta@slac.stanford.edu

c1=cin;

if modef==0
   %%simple lorentzian
   f=@(c,x) c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2);
elseif modef==1
   %%simple gaussian
   f=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
elseif modef==2%14
   %%Landau function
   % f=@(c,x) c(1)*Landau((x-c(2))./c(3));
   %%assymetric secant squared
   f=@(c,x) (c(1)*sech((x-c(2))./c(3)).^2).*(x < c(2))+ ...
      (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2));
elseif modef==3
   %%asymmetric lorentzian
   f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
      (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2));
elseif modef==4
   %%asymmetric gaussian
   f=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
elseif modef==5
   %%double gaussian
   f=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2)+ ...
      c(4)*exp(-(x-(c(2)+shift)).^2/c(5)^2);
   f1=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
elseif modef==6
   %%double gaussian with variable shift
   f=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2)+ ...
      c(4)*exp(-(x-(c(2)+c(6))).^2/c(5)^2);
   f1=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
elseif modef==7
   %%asymetric gaussian and small gaussian
   f=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2))+...
      c(5)*exp(-(x-(c(2)+shift)).^2/c(6)^2);
   f0=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
   f1=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
elseif modef==8
   %%asymetric gaussian and small gaussian with var shift
   f=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2))+...
      c(5)*exp(-(x-(c(2)+c(7))).^2/c(6)^2);
   f0=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
   f1=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
elseif modef==9
   %%double asymetric gaussian
   f=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2))+...
      (c(5)*exp(-(x-(c(2)+shift)).^2/c(6)^2)).*(x < c(2)+shift)+ ...
      (c(5)*exp(-(x-(c(2)+shift)).^2/c(7)^2)).*(x >= c(2)+shift);
   f1=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
elseif modef==10
   %%double asymetric gaussian with variable shift
   f=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2))+...
      (c(5)*exp(-(x-(c(2)+c(8))).^2/c(6)^2)).*(x < c(2)+c(8))+ ...
      (c(5)*exp(-(x-(c(2)+c(8))).^2/c(7)^2)).*(x >= c(2)+c(8));
   f1=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
elseif modef==11
   %%gaussian with asymetric gaussian with fixed shift
   f=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2)+...
      (c(4)*exp(-(x-(c(2)+shift)).^2/c(5)^2)).*(x < c(2)+shift)+ ...
      (c(4)*exp(-(x-(c(2)+shift)).^2/c(6)^2)).*(x >= c(2)+shift);
   f0=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
   f1=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
elseif modef==12
   %%gaussian with asymetric gaussian with variable shift
   f=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2)+...
      (c(4)*exp(-(x-(c(2)+c(7))).^2/c(5)^2)).*(x < c(2)+c(7))+ ...
      (c(4)*exp(-(x-(c(2)+c(7))).^2/c(6)^2)).*(x >= c(2)+c(7));
   f0=@(c,x) c(1)*exp(-(x-c(2)).^2/c(3)^2);
   f1=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
elseif modef==13
   %%double asymetric gaussian with variable shift, both lowE sigmas equal
   f=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2))+...
      (c(5)*exp(-(x-(c(2)+c(7))).^2/c(3)^2)).*(x < c(2)+c(7))+ ...
      (c(5)*exp(-(x-(c(2)+c(7))).^2/c(6)^2)).*(x >= c(2)+c(7));
   f1=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
elseif modef==14
   %%half gaussian with assymetric gaussian, with variable shift, both lowE sigmas equal
   f=@(c,x)  (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x >= c(2))+...
      (c(5)*exp(-(x-(c(2)+c(7))).^2/c(4)^2)).*(x < c(2)+c(7))+ ...
      (c(5)*exp(-(x-(c(2)+c(7))).^2/c(6)^2)).*(x >= c(2)+c(7));
   f1=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
      (c(1)*exp(-(x-c(2)).^2/c(4)^2)).*(x >= c(2));
   f0=@(c,x) (c(1)*exp(-(x-c(2)).^2/c(3)^2)).*(x >= c(2));
elseif modef==15
   %%double asymmetric lorentzian, fixed shift
   f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
      (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+...
      (c(5)*ones(size(x))./(1+((x-(c(2)+shift))./(0.5.*c(6))).^2)).*(x < c(2)+shift)+ ...
      (c(5)*ones(size(x))./(1+((x-(c(2)+shift))./(0.5.*c(7))).^2)).*(x >= c(2)+shift);
   f1=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
      (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2));
elseif modef==16
   %double asymmetric lorentzian, fitting shift
   f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
      (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+...
      (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(6))).^2)).*(x < c(2)+c(8))+ ...
      (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(7))).^2)).*(x >= c(2)+c(8));
   f1=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
      (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2));
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
   f1a=@(c,x) (c(1)*ones(size(x))./(1+((x-(c(2)-c(5)/2))./(0.5.*c(3))).^2)).*(x < c(2)-c(5)/2)+ ...
      (c(1)*ones(size(x))./(1+((x-(c(2)-c(5)/2))./(0.5.*c(4))).^2)).*(x >= c(2)-c(5)/2);
   f1b=@(c,x)(c(1)*ones(size(x))./(1+((x-(c(2)+c(5)/2))./(0.5.*c(3))).^2)).*(x < c(2)+c(5)/2)+ ...
      (c(1)*sech((x-(c(2)+c(5)/2))./c(4)).^2).*(x >=c(2)+c(5)/2);
   % 8->1
   % 6->3
   % 7->4
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
end

xdata=1:1024;

if modef<2
   ymain=f(c1(1:3),xdata);
elseif modef<5
   ymain=f(c1(1:4),xdata);
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
elseif modef==16
   ymain=f1(c1(1:4),xdata);
   ysignal=f1([c1(5) c1(2)+c1(8) c1(6) c1(7)],xdata);
elseif modef==17
   ymain=f0(c1(1:4),xdata);
   ysignal=f1([c1(5) c1(2)+shift c1(6) c1(7) c1(8) c1(9) c1(10) c1(11)],xdata);
elseif modef==18
   ymain=f0(c1(1:4),xdata);
   ysignal=f1([c1(5) c1(2)+c1(8) c1(6) c1(7) c1(9) c1(10) c1(11) c1(12)],xdata);
   ysig1=f1a([c1(5) c1(2)+c1(8) c1(6) c1(7) c1(9)],xdata);
   ysig2=f1b([c1(12) c1(2)+c1(8) c1(10) c1(11) c1(9)],xdata);
elseif modef==19
   ymain=f0([c1(1:4) c1(9)],xdata);
   ysignal=f1([c1(5) c1(2)+c1(8) c1(6) c1(7)],xdata);
end
if exist('ysignal')==0
   ysignal=0;
end

yfit= f(c1,xdata);

end

function [c1,resnorm,residual,exitflag,output,lambda,J,c1std] = lsqcurvefitstd_sep(varargin)
% this is a wrapper that calculates the confidence interval in addition to
% performing the lsqcurvefit.
% run this exactly as you would lsqcurvefit - see: help lsqcurvefit
% additional output: c1std = standard error for each fit parameters

[c1,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(varargin{:});
s2 = sum((residual).^2)/(length(residual)-length(c1));
vx = inv(J'*J)*s2;
c1std = full(real(diag(sqrt(vx))'));


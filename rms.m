function y=rms(x)
y = sqrt(mean(x .* conj(x)));

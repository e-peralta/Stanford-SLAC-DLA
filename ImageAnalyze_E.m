%  ImageAnalyze_E.m

%  Obtains horizontal and vertical traces through beam center to perform
%  ASSYMETRIC gaussian fit. 
%
%  130215-eperalta, used same intro as Eric Colby's ImageAnalyze.m script
%---------------------------------------------------------------------

%f=inline('c(1)*exp(-(x-c(2)).^2/c(3)^2)+c(4)','c','x');
f=inline('(c(1)*exp(-.5*(x-c(2)).^2/c(3)^2)+c(4)).*(x < c(2))+ (c(1)*exp(-.5*(x-c(2)).^2/c(5)^2)+c(4)).*(x >= c(2))','c','x');

SIG_OPTIONS=optimset('Display','off','TolFun',1e-05,'TolX',1e-05);
showplots=0;

se=sum(sum(img(1:2:end,:)'));
so=sum(sum(img(2:2:end,:)'));
%ver1=sum(img');
%[x i]=max(ver1);
%eo=mod(i,2);
if(se>=so)
    eo=0;
else
    eo=1;
end
%

%img=orientimage(img,Screen(ws));

hor=sum(img);
ver=sum(img');
hor=hor/length(ver);
ver=ver/length(hor);
[hv,hi]=max(hor);
[vv,vi]=max(img(:,hi));

% %interpolate even line outs
% if IsXyb==2 % is it a gigE cam?
%     %do nothing
% else % then it's a regular video cam
% if(Screen(ws).NROT==0 || Screen(ws).NROT==2)
%     for i=2:floor((length(ver)-1)/2-1)
%          ver(i*2+eo)=(ver(i*2-1+eo)+ver(i*2+1+eo))/2;
%     end
% else
%     for i=2:floor((length(hor)-1)/2)
%          hor(i*2+eo)=(hor(i*2-1+eo)+hor(i*2+1+eo))/2;
%     end   
% end
% end

    [g1,g2]=max(hor);
    g4=mean(hor(1:20));
    g3=(sum(hor)-length(hor)*g4)/g1/sqrt(2*log(2));
    g5=g3;
    %g3=30;
    guess=[g1,g2,g3,g4,g5];
    %[c1 ChiSqX]=lsqcurvefit(f,guess,1:length(hor),hor,[],[],SIG_OPTIONS);
    [c1, ~ , ~, ~, ~ , ~, ~, c1std]=lsqcurvefitstd(f,guess,1:length(hor),hor,[],[],SIG_OPTIONS);
    bga1=0;
    c1=abs(c1);
    [g1,g2]=max(ver);
    %g4=mean(ver(1:20));
    g4=mean(ver(1:min(20,length(ver))));
    g3=(sum(ver)-length(ver)*g4)/g1/sqrt(2*log(2));
    g5=g3;
    %g3=30;
    guess=[g1,g2,g3,g4,g5];
    %[c2 ChiSqY]=lsqcurvefit(f,guess,1:length(ver),ver,[],[],SIG_OPTIONS);
    [c2, ~ , ~, ~, ~ , ~, ~, c2std]=lsqcurvefitstd(f,guess,1:length(ver),ver,[],[],SIG_OPTIONS);
    bga1y=0;
    c2=abs(c2);


if showplots==1
fprintf('Gets here \n')
figure(3);
subplot(2,1,1);
plot(hor); hold on;
xov=[1:length(hor)];
plot(f(c1,xov)+bga1*xov,'r');
hold off;
subplot(2,1,2);
plot(ver); hold on;
yov=[1:length(ver)];
plot(f(c2,yov)+bga1y*yov,'r');
hold off;

end

%Oct 22 2010, added for use inside function where ver refers to command
vert=ver;
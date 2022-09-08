clear all
close all

%cols(0,:)=[1 1 1];  %use 'w'
cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple
%cols(5,:)=[0 0 0]; %use 'k'

folder='~/Dropbox/Research/Data/';

width=3.3;
height=2;
close all
hFig = figure(94);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')


filename='OscillatorSpectrum.txt';
data=importdata([folder,filename]);


wave = data(:,1); %
val = data(:,2); %


% convert units

w=linspace(wave(1),wave(end),100);

% find fit
%f=@(c) sum((c(1)*(k-c(2)).^2+c(3)-s11).^2);
%[c, rm] = fminsearch(f,[1e-3 4*factor 1e-4]);

h1a=plot(wave,val,'Color',cols(3,:),'Linewidth',1);
hold on
%h1f=plot(kf,(c(1)*(kf-c(2)).^2+c(3))*1e6,':','LineWidth',1,'Color',cols(4,:));


for i=1:4 %extend to 4 to test the other fitting algorithms (1 was best)
    [c1, c1std, yfit, ChiSq]=FitSpectrum(i,val,wave);%,0,0,1);
    h1f=plot(wave,yfit,'LineWidth',1,'Color',cols(2,:));%col4)
    hold on
    switch i
        case 1
            scale1=sqrt(2*log(2));
            scale2=sqrt(2*log(2));
            f=@(c,x) (c(1)*exp(-.5*(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
                (c(1)*exp(-.5*(x-c(2)).^2/c(4)^2)).*(x >= c(2))+c(5);
        case 2
            scale1=1/2;
            scale2=1/2;
            f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
                (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+c(5);
        case 3
            scale1=asech(sqrt(1/2));
            scale2=asech(sqrt(1/2));
            f=@(c,x) (c(1)*sech((x-c(2))./c(3)).^2).*(x < c(2))+ ...
                (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2))+c(5);
        case 4
            scale1=1/2;
            scale2=asech(sqrt(1/2));
            f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
                (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2))+c(5);
    end
    %h1f=plot(time,f(c1,time),':','LineWidth',1,'Color',cols(3,:));%col4)
    %hold on
    c1(3)*scale1
    c1(4)*scale2
    FWHM=(c1(3)*scale1+c1(4)*scale2);
    if c1std(3)*scale1>c1std(4)*scale2
        FWHMstd=c1std(3)*scale1;
    else
        FWHMstd=c1std(4)*scale2;
    end
    
    %line([c1(2)-c1(3)*scale1 c1(2)+c1(4)*scale2],...
    %    [c1(5)+c1(1)/2 c1(5)+c1(1)/2], 'LineWidth',2,'Color','r');
    fprintf('\nmode %g: FWHM= %g +- %g ps, Chi^2=%g\n',i,FWHM,FWHMstd,ChiSq)
end


hy=ylabel('Photon count [a.u.]');
hx=xlabel('Wavelength [nm]');

ax1=gca;

axis tight

tightfig;
%hL=legend([h1a,h1b,h2a,h2b,h1f,h2f],{'Scan 1 X','Scan 1 Y','Scan 2 X','Scan 2 Y','Quad. fit X','Quad. fit Y'});
%set(hL,'Location','northeastoutside','Box','Off')
set(ax1,'FontSize', 10);
set([hx,hy],'FontSize', 10);

%export_fig IRspectrum.eps -cmyk -r300 -painters%-painters

%%

width=3.5;
hFig = figure(95);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')


filename='Autocorrelation.txt';
data=importdata([folder,filename]);


time = data(:,1)*1e9/299792458*2; % from mm to ps

val = data(:,2)/1000; %
er=(data(:,3)-data(:,2))/1000;

% convert units

w=linspace(time(1),time(end),100);

% find fit
%f=@(c) sum((c(1)*(k-c(2)).^2+c(3)-s11).^2);
%[c, rm] = fminsearch(f,[1e-3 4*factor 1e-4]);



%h1a=plot(time,val,'Color',cols(3,:),'Linewidth',1);


for i=3:3 %extend to 4 to test the other fitting algorithms (1 was best)
    [c1, c1std, yfit, ChiSq]=FitSpectrum(i,val,time);%,0,0,1);
    h1f=plot(w-c1(2),yfit,'LineWidth',1,'Color',cols(2,:));%col4)
    hold on
    switch i
        case 1
            scale1=sqrt(2*log(2));
            scale2=sqrt(2*log(2));
            f=@(c,x) (c(1)*exp(-.5*(x-c(2)).^2/c(3)^2)).*(x < c(2))+ ...
                (c(1)*exp(-.5*(x-c(2)).^2/c(4)^2)).*(x >= c(2))+c(5);
        case 2
            scale1=1/2;
            scale2=1/2;
            f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
                (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2))+c(5);
        case 3
            scale1=asech(sqrt(1/2));
            scale2=asech(sqrt(1/2));
            f=@(c,x) (c(1)*sech((x-c(2))./c(3)).^2).*(x < c(2))+ ...
                (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2))+c(5);
        case 4
            scale1=1/2;
            scale2=asech(sqrt(1/2));
            f=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
                (c(1)*sech((x-c(2))./c(4)).^2).*(x >=c(2))+c(5);
    end
    %h1f=plot(time,f(c1,time),':','LineWidth',1,'Color',cols(3,:));%col4)
    %hold on
    c1(3)*scale1
    c1(4)*scale2
    FWHM=(c1(3)*scale1+c1(4)*scale2);
    if c1std(3)*scale1>c1std(4)*scale2
        FWHMstd=c1std(3)*scale1;
    else
        FWHMstd=c1std(4)*scale2;
    end
    
    %line([c1(2)-c1(3)*scale1 c1(2)+c1(4)*scale2],...
    %    [c1(5)+c1(1)/2 c1(5)+c1(1)/2], 'LineWidth',2,'Color','r');
    fprintf('\nmode %g: FWHM= %g +- %g ps, Chi^2=%g\n',i,FWHM,FWHMstd,ChiSq)
end

h1a=errorbar(time-c1(2),val,er,'o','Color',cols(3,:),'Linewidth',1,...
    'MarkerSize',4,'MarkerFaceColor',cols(1,:));

%hold on
%h1f=plot(kf,(c(1)*(kf-c(2)).^2+c(3))*1e6,':','LineWidth',1,'Color',cols(4,:));

hy=ylabel('Photon count [a.u.]');
hx=xlabel('time delay [ps]');

ax1=gca;

axis tight

tightfig;
hL=legend([h1a,h1f],{'Autocor.','Sech^2 fit'});
set(hL,'Location','northwest','Box','Off')
set(ax1,'FontSize', 10);
set([hx,hy,hL],'FontSize', 10);

%export_fig IRpulse.eps -cmyk -r300 -painters%-painters

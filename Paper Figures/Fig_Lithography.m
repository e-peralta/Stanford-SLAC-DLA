cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple
SIG_OPTIONS=optimset('Display','off','TolFun',1e-05,'TolX',1e-05);


width=4;
height=2.8;
close all
hFig = figure(94);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')

scale=800;

N=[21 06 07 10 11 14 15 19];
E=[95 100 105 110 115 120 125 150];
P=[798.1 788.0 783.3 785.3 787.2 789.1 781.7 790.7];
W=[476.9 447.6 437.0 407.1 407.2 408.7 386.5 350.0]/scale;
dW=[8.2 14.5 13.9 9.6 4.3 6.7 2.9 4.1]/scale;
%W=[476.9 447.6 437.0 407.1 407.2 408.7 386.5 350.0]./P;
%Wer=[8.2 14.5 13.9 9.6 4.3 6.7 2.9 4.1]./P;

guess=[-2.1 660]/scale;
cmin=[-5 0]/scale;
cmax=[0 1000]/scale;
%guess=[-003 .8];
%cmin=[-005 0];
%cmax=[0 1];
lin=@(c,x) (c(2)+x*c(1))./dW;
cout=lsqcurvefitstd(lin,guess,E,W./dW,cmin,cmax,SIG_OPTIONS);
%fitE=linspace(E(1),E(end),100);
fitE=95:.2:165;
fitW=cout(1)*fitE+cout(2);

%plot(E,W,'x');
%
h1=errorbar(E,W,dW,'o','Color',cols(3,:),'Linewidth',1,...
    'MarkerSize',4,'MarkerFaceColor',cols(1,:));
hold on
h1f=plot(fitE,fitW,':','Color',cols(1,:));
hold on
href=plot(fitE,400/scale*ones(size(fitE)),':k');
%href=plot(fitE,.5*ones(size(fitE)),':k');


%% A9 New Accel FOM

E2=[127.5:7.5:165];
W2=[445 379 327 312 291 249]/scale;
dW2=[8 15 7 15 12 7]/scale;

guess=[-4.8 1000]/scale;
cmin=[-10 0]/scale;
cmax=[0 2000]/scale;
lin=@(c,x) (c(2)+x*c(1))./dW2;
cout=lsqcurvefitstd(lin,guess,E2,W2./dW2,cmin,cmax,SIG_OPTIONS);
fitW2=cout(1)*fitE+cout(2);


E3=[105 112.5 120 127.5 135  135  142.5 150 157.7 165];
W3=[476 448   423 384    354 360  336  313 305   283]/scale;
dW3=[11  10    14   8     10   5   5   9   11    7  ]/scale;
%E3=[105 105 112.5 120 127.5 135  135  142.5 150 157.7 165];
%W3=[476 468 448   423 384    354 360  336  313 305   283];
%dW3=[11  7   10    14   8     10   5   5   9   11    7  ];

guess=[-3.3 810]/scale;
cmin=[-10 0]/scale;
cmax=[0 2000]/scale;
lin=@(c,x) (c(2)+x*c(1))./dW3;
cout=lsqcurvefitstd(lin,guess,E3,W3./dW3,cmin,cmax,SIG_OPTIONS);
fitW3=cout(1)*fitE+cout(2);


E4=[105 112.5 120 127.5 135  142.5 150 157.7 165];
W4=[457 427   402  380 363   341  319 306   288]/scale;
dW4=[13   9    3    5   16   6     8    5   13 ]/scale;
%E4=[105 112.5 120 127.5 127.5 135  142.5 150 157.7 165];
%W4=[457 427   402 384   380 363   341  319 306   288];
%dW4=[13   9    3    4    5   16   6     8    5   13 ];

guess=[-2.8 740]/scale;
cmin=[-10 0]/scale;
cmax=[0 2000]/scale;
lin=@(c,x) (c(2)+x*c(1))./dW4;
cout=lsqcurvefitstd(lin,guess,E4,W4./dW4,cmin,cmax,SIG_OPTIONS);
fitW4=cout(1)*fitE+cout(2);

hold on
%figure

%errorbar(E2,W2,dW2,'rx--')

h2=errorbar(E2,W2,dW2,'v','Color',cols(4,:),'Linewidth',1,...
    'MarkerSize',4,'MarkerFaceColor',cols(2,:));
hold on
h2f=plot(fitE,fitW2,':','Color',cols(2,:));


hold on
%errorbar(E3,W3,dW3,'kx--')
h3=errorbar(E3,W3,dW3,'p','Color',cols(4,:),'Linewidth',1,...
    'MarkerSize',4,'MarkerFaceColor',cols(4,:));
hold on
h3f=plot(fitE,fitW3,'--','Color',cols(4,:));

hold on
%errorbar(E4,W4,dW4,'bx--')
h4=errorbar(E4,W4,dW4,'^','Color',cols(4,:),'Linewidth',1,...
    'MarkerSize',4,'MarkerFaceColor',cols(2,:));
hold on
h4f=plot(fitE,fitW4,'--','Color',cols(2,:));

%%


axis tight
ylim([240 500]/scale)
hy=ylabel('Grating duty cycle');
hx=xlabel('Exposure dose [mJ/cm^2]');

ax1=gca;

tightfig;

%F<0 means focus moves AWAY from wafer
hL=legend([h1,h2,h3,h4],{'On Resist, F=0um','Etched, F=.2um','Etched, F=-.3um','Etched, F=-.8um'});
set(hL,'Location','southwest','Box','Off')
set(ax1,'FontSize', 10);
set([hx,hy,hL],'FontSize', 10);
set(ax1,'XTick',95:10:165)

export_fig Litho0.eps -cmyk -r300 -painters%-painters

%%







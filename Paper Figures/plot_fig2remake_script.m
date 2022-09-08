%%
load fig2remake_model.mat

%
subplot(1,2,1)
hold on

plot(x0,raw_off,'bx')
plot(x0,fit_full_off,'b','LineWidth',2)

plot(x0,raw_on,'rx')
plot(x0,fit_main_on+fit_plateau_off+model_signal_on,'r','LineWidth',2)

axis([675 975 -.01 .3])
title('Adding the plateau directly after modulation')
legend('laser on','model','laser off','model','Location','Best')

subplot(1,2,2)
hold on

plot(x0,raw_off,'bx')
plot(x0,fit_full_off,'b','LineWidth',2)

plot(x0,raw_on,'rx')
plot(x0,fit_main_on+model_plateau_signal_on,'r','LineWidth',2)

axis([675 975 -.01 .3])
title('Modulating the plateau with the signal')

%
figure(2)

subplot(1,2,1)
hold on

plot(x0,raw_off-fit_main_off-fit_plateau_off,'bx')
plot(x0,fit_signal_off,'b','LineWidth',2)

plot(x0,raw_on-fit_main_on-fit_plateau_off,'rx')
plot(x0,model_signal_on,'r','LineWidth',2)

axis([675 975 -.01 .275])

title('Adding the plateau directly after modulation')
legend('laser on','model','laser off','model','Location','Best')

%
subplot(1,2,2)
hold on

plot(x0,raw_off-fit_main_off,'bx')
plot(x0,fit_signal_off+fit_plateau_off,'b','LineWidth',2)

plot(x0,raw_on-fit_main_on,'rx')
plot(x0,model_plateau_signal_on,'r','LineWidth',2)


axis([675 975 -.01 .275])
title('Modulating the plateau with the signal')

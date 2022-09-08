%% load data
load('run2224_output.mat');
load('run2224_Spectra.mat');
c_off = cout_all(find(RS(1,:)==0),:);
c_on = cout_all(find(RS(1,:)==1),:);

%% find FWHM for dE = 0 to 400 pixels

g = [1 mean(c_on(:,2)) mean(c_on(:,3)) mean(c_on(:,4)) mean(c_off(:,5)./c_on(:,1)) mean(c_off(:,6)) mean(c_off(:,7)) mean(c_off(:,8))];

fwhms = [];
fwhmf = [];

i6L = [5:.25:200];
%i6L = [10];
ik6 = 0;
for i6 = i6L
    
    ik6 = ik6+1;
    dE = i6;
    x1 = [0:1:1600];
    N=200;
    s2 = x1*0;
    
    f1=@(c,x) (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(6))).^2)).*(x < c(2)+c(8))+ ...
        (c(5)*ones(size(x))./(1+((x-(c(2)+c(8)))./(0.5.*c(7))).^2)).*(x >= c(2)+c(8));
    
    g1 = g;
    for i1 = 0:N
        g1(8) = g(8)+dE*cos(2*pi/N*i1);
        s1 = f1(g1,x1);
        s2 = s2+s1;
    end
    s2 = s2/(N+1);
    s0 = f1(g,x1);
    
    [cout cout_ci f2] = double_lor(s2);
    
    f4=@(c,x) (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(3))).^2)).*(x < c(2))+ ...
        (c(1)*ones(size(x))./(1+((x-c(2))./(0.5.*c(4))).^2)).*(x >= c(2));
    
    nf = ceil(i6/15)+1;
    ms5 = s2-f4(cout,x1);
    ms2 = conv(ms5,ones(1,nf)/nf,'same');
    [a m0] = max(ms2);
    m1 = find(ms2(m0:end) < ms2(m0)/2,1,'first')-1;
    
    xhw = interp1(ms2(m0:end),x1(m0:end),a/2);
    fwhms(ik6) = xhw;
    fwhmf(ik6) = cout(7)/2+cout(2)+cout(8);
    
    %plot(x1,s2,x1,ms2,x1,f4(cout,x1),xhw,a/2,'o',[863.3266 863.3266],[0 .2])
    %hold on    
end


yref = f1(g,x1);
[yrm yri] = max(yref);
xref = interp1(yref(yri:end),x1(yri:end),yrm/2);

%
figure()
xhi = i6L*2.4;
plot(fwhms-xref,xhi)
xlabel('shift [pix]')
ylabel('Energy [GeV]')

figure()
xhi = i6L*2.4;
plot(xhi,fwhms-xref)
ylabel('shift [pix]')
xlabel('Energy [GeV]')

%%

xref400 = xref;
fwhms400 = fwhms;
xhi400 = xhi;

%% plot

%plot(fwhms100x-xref100,xhi100x,'.'); hold on
xp = [0:.001:180];
yp = interp1(fwhms100x-xref100,xhi100x,xp,'pchip');
plot(xp,yp)

%%
cout = [20.6649993852399,0.0202816723302518,-0.464135359959495,-21.6178770149860;]
A = cout(1);
B = cout(2);
C = cout(3);
D = cout(4);


yp2 = 1/B*(sqrt(((xp-D)/A).^2-1)+C);
hold on 
plot(xp,yp2,'r')

%% add labels
xlabel('shift [pix]')
ylabel('Energy Mod [MeV/m]')
grid on
enhance_plot
axis([0 130 -10 350])

%% more points is worse

plot(fwhms100-xref100,xhi100)
hold on 
plot(fwhms200-xref200,xhi200,'r')
plot(fwhms400-xref400,xhi400,'g')

xlabel('shift [pix]')
ylabel('Energy [GeV]')

%%
load calibration_values

xx = [0:.001:189];
yy = interp1(xcal,ycal,xx,'pchip');

plot(xx,yy)

title('Calibration')
xlabel('shift [pix]')
ylabel('Energy Mod [MeV/m]')
grid on
enhance_plot


axis([0 130 -10 350])
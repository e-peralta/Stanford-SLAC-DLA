%% FasterRasterProcess.m
function [vroll data] = FasterRasterProcess(filename)
% This script is used for post analysis of data acquired using the
% FasterRaster VI at the E163 experiment. Run the script and select the
% desired data file when prompted.
%
% Upgrades:
% 2011-10-12 kensoong@SLAC - original script
% 2011-10-19 eperalta@SLAC - added filtering and trend analysis
% 2012-05-14 behnamm@SLAC  - added auto detection of roll/step index
% 2012-05-14 eperalta@SLAC - added automatic file selection and detection
%                            to plot either a single save spectra, a screen
%                            capture, and a 1-D or 2-D scan with or without
%                            save spectra to analyze FOM after the fact

%% rasterscan processing

%clc
%clear all
%close all

%index guide: 1=x, 2=y, 3=tip, 4= rot.

FILTERON=1;
K_FILT=3; %change this to make the low pass filter stronger
FITON=0;

LARGEAPERTURE=0;

GATED=1;
SHIFT=230; %230 for gratings
DELTA=45;

GATED=0;
SHIFT=0; %230 for gratings


PLOTALLSPECTRA=0;
PAUSET=.5; %pause time for each spectrum plotted

SAVEPLOTS=1;


spectracatch=0;
% load from data file
%filepath=['V:\ARDB\E163\DATA\',datestr(now,'yymmdd'),'\'];
%filepath=['Y:\LEAP\Data_2012_OctRun','\'];

filepath='~/Dropbox/Research/Data/';

if ~exist('filename')
[filename,filepath] = uigetfile({'*.dat;*.mat'},'Choose file to process',filepath);
end

test=load([filepath,filename]);
if filename(end-2:end)=='mat'
    test=test.spectra;
end
filename = [filepath,filename];   
filename = filename(1:end-4);
sprintf('File chosen: %s',filename)


% %% I was planning to use this to fix the problem with the last column of some 
% of the data. 
% checkV=unique(test(:,2))
% %test(:,end);
% clear indmat
% for i=length(checkV)
% indmat(:,i)=(test(:,2)==checkV(i))
%     %test((test(:,2)==checkV(i)),end)=i-1;
% end
% %test(:,end)


%% Detect if selected data is a screen capture (square matrix) and plot

if size(test,1)==size(test,2) 
    
    
    % generate projection onto x/energy axis
    ROI=[1 1024 1 1024];
    
    figure (1)
     
    trace=mean(test(ROI(3):ROI(4),ROI(1):ROI(2)));
    trace=trace-mean(trace(1:10))*ones(1,length(trace));
    trace=trace/max(trace);
    trace=trace*size(test,1)*.6;
    trace=trace+.1*size(test,1)*ones(1,length(trace));
    
    %trace=800*trace/max(trace)+100*ones(1,length(trace));
    plot(trace)

    % create energy scale
    x_tpeak=654; %Get this number by hand from the previous plot
    %x_tpeak=0;
    x=1:length(test);
    xval2=1.7*(x-x_tpeak*ones(1,length(x)));  %1.7keV/px
    
    %re-scale the image for a nice plot
    minC=min(min(test));
    test=test-minC(ones(size(test)));
    maxC=max(max(test));
    test=test*255/maxC;

    figure (2)
    image(xval2,1:length(trace),test)
    hold on
    axis tight
    axis xy
    %axis equal
    caxis([min(min(test)) max(max(test))])
    ylim([100 length(trace)])
    plot(xval2,trace,'y','Color',[1 1 1],'LineWidth',3)
    xlabel('\Delta E [keV]')
    ylabel('Pixels/Counts[a.u]')
    
%    enhance_plot(0,0,-1);
    if SAVEPLOTS
    print(gcf, '-dpng',filename);
    end

else

%% Detect if selected data is a simple spectrum
if size(test,1)==1
    %%
%    figure(1)
    %[cout yfit ChiSq ybkgd ymain]=FitSpectrum3(test,16,SHIFT); 
    %[cout yfit ChiSq ybkgd ymain]=FitSpectrum3(fliplr(test),16,SHIFT,1); 

%     plot(test-ybkgd,'k')
%     hold on 
%     
%     plot(yfit,'r')
%     hold on
%     plot(ymain,'g:')
%     hold on 
%     plot(ysignal,'b:')
    
    %xtrans=cout(2)+cout(7);
    xtrans=0;
    ybkgd=zeros(1,length(test));
    x=1:length(test);
    xval2=1.7*(x-xtrans*ones(1,length(x)));  %1.7keV/px

    std1=std(test(1:400))
    std2=std(test(end-400:end))
    figure (1)
        
    plot(xval2,test-ybkgd,'LineWidth',2)
    axis tight
    xlabel('\Delta E [keV]')
    
    title(filename);
    grid on
    axis tight
%    enhance_plot;
    
    
    if SAVEPLOTS
    print(gcf, '-dpng',filename);
    end

else

    % ----- Finding the roll and step axis
% This hasn't yet been made to work with all types of scan yet, only with a
% 2D parameter scan
RasterData=test;

ind1=find(RasterData(:,end)-RasterData(1,end)~=0,1);
ref1=(RasterData(1,1:4)~=RasterData(2,1:4));
roll_index=find(ref1==1,1);
if RasterData(1,end)==RasterData(end,end)
    step_index=(roll_index==1)+1;
else
    ref2=(RasterData(1,1:4)~=RasterData(ind1,1:4));
    ref3=find(ref2==1,2);
    step_index=(ref3(find(ref3~=roll_index,1)));
end

clear ind1 ref1 ref2 ref3

% ------ separate inputs
vstep=test(:,step_index);
vroll=test(:,roll_index);
ncount=test(:,end);

nsteps = length(unique(ncount));
step = (vstep(end)-vstep(1))/nsteps;
roll0= vroll(1);
roll1 = vroll(find(ncount==0,1,'Last'));

switch step_index
    case 1
        step_coord='X';
    case 2
        step_coord='Y';
    case 3
        step_coord='TIP';
    case 4
        step_coord='ROT';
end

switch roll_index
    case 1
        roll_coord='X';
    case 2
        roll_coord='Y';
    case 3
        roll_coord='TIP';
    case 4
        roll_coord='ROT';
end


%% ---- Extract FOM if spectrums were saved
if size(test,2)>6
    buffer=1;
    bpf=.2; %baseline point fraction
    mrank=3; %median rank
    starti=1;
    %starti=1000;
    endi=size(test,1)
    %endi=1600
    spectracatch=1;
       
    for i=starti:endi
    %for i=29    
        clear xo yo
        spectrum(i,:)=test(i,5+buffer:end-1-buffer);
        
        if ~LARGEAPERTURE
        spectrum(i,:)=spectrum(i,:)./max(spectrum(i,:));   
        end
        specSort=sort(spectrum(i,:));
        baseline=mean(specSort(1:round(bpf*length(specSort))))*ones(1,length(specSort));
        spec=spectrum(i,:)-baseline;
        %spec2=medfilt1(spec,mrank);
        spec2=spec;
        [maxV ind]=max(spec2);        
        xo=[ind];
        yo=[maxV];
        
        if LARGEAPERTURE
            indlo=ind-SHIFT;
            indhi=ind+SHIFT;
            indmax=length(spec2);
            if indlo<1 indlo=1; end
            if indhi>indmax indhi=indmax; end
            Vlo=spec2(indlo);
            Vhi=spec2(indhi);
            if Vhi < Vlo
                ind=indlo;
                xo=[ind];
                yo=[Vlo];
            end
        end
        if GATED
            %use this part to detect large energy shifts
            indhi=ind+SHIFT+DELTA;
            if indhi>length(spec2)
                indhi=length(spec2);
                sprintf(['Integration window truncated at: i=',num2str(i),...
                     '     ',step_coord,'= ',num2str(test(i,step_index)),...
                  '     ',roll_coord,'= ',num2str(test(i,roll_index))]);
                %FOM(i)=0;
            end
            indlo=ind+SHIFT-DELTA;
            if indlo<1
                indlo=1;
                sprintf(['Integration window truncated at: i=',num2str(i),...
                     '     ',step_coord,'= ',num2str(test(i,step_index)),...
                  '     ',roll_coord,'= ',num2str(test(i,roll_index))]);
                %FOM(i)=0;
            end
            %FOM(i)=sum(spec2(indlo:indlo))*1e-5;
            xo=[xo ind+SHIFT];
            yo=[yo spec2(ind+SHIFT)];
            
        else
            %FOM(i)=sum(spec2)*1e-5;
        end
        
        if PLOTALLSPECTRA
        figure(2)
        %xshift is the center of main peak to plot spectra overlapped
        x=1:length(spectrum(i,:));
        xshift=x-ind*ones(1,length(x));  
        plot(xshift,spectrum(i,:),'k:')
        hold on
        %plot(xshift,specSort,'c:')
        %hold on
        %plot(xshift,spec,'g.-');
        %hold on
        plot(xshift,spec2,'r-');
        hold on        
        
        plot(xo-ind*[1 1],yo,'bx','MarkerSize',8,'LineWidth',2)
        xlim([-512 512])
%         titlestr=['Median rank= ',num2str(mrank),'    i = ',num2str(i),...
%                   '     ',step_coord,'= ',num2str(test(i,step_index)),...
%                   '     ',roll_coord,'= ',num2str(test(i,roll_index))];
        titlestr=['Median rank= ',num2str(mrank),'    i = ',num2str(i),...
                  '     ',step_coord,'= ',num2str(test(i,step_index)),...
                  '     ',roll_coord,'= ',num2str(test(i,roll_index)),...
                  '     Peak=',num2str(ind)];
        title(titlestr);
        drawnow
        pause(PAUSET)
                
        hold off
        
        end
%        test(i,5)=FOM(i);
    end
end  
clear i buffer bpf ind maxV baseline xo yo specSort spec spec2 PLOTALL
data=(test(:,5));

% if start ~= 1
%     figure (174)
%     plot(FOM(start:end))
% end


%% Single Scan data plot
% this section is to replace the jpg created by the VI with a useful one
if nsteps==1
    mrank=4;
    
    if SHIFT==0;
    ycaption=['FOM: Int ALL,    spect. med-rank=',num2str(mrank)];
    else
        %{
    ycaption=['FOM: Gated Int,    E_0=',num2str(SHIFT),'     \DeltaE=',...
                num2str(DELTA),'    spec. med-rank=',num2str(marank)];
    %}
    ycaption=['FOM: Gated Int,    E_0=',num2str(SHIFT),'     \DeltaE=',...
                num2str(DELTA),'    spec. med-rank=',num2str(0)];
    end

    

    switch roll_index
    case 1
        rolllabel='X';
        titlestr=['Y = ',num2str(test(1,2)),...
                  '   TIP = ',num2str(test(1,3)),...
                  '   TILT = ',num2str(test(1,4))];
    case 2
        rolllabel='Y';
        titlestr=['X = ',num2str(test(1,1)),...
                  '   TIP = ',num2str(test(1,3)),...
                  '   TILT = ',num2str(test(1,4))];
    case 3
        rolllabel='TIP';
        titlestr=['X = ',num2str(test(1,1)),...
                  '   Y = ',num2str(test(1,2)),...
                  '   TILT = ',num2str(test(1,4))];
    case 4
        rolllabel='TILT';
        titlestr=['X = ',num2str(test(1,1)),...
                  '   Y = ',num2str(test(1,2)),...
                  '   TIP = ',num2str(test(1,3))];
    end
 
    if PLOTALLSPECTRA
    figure (1)
    plot(vroll,data,'bx:')
    hold on
    %data2=medfilt1(data,mrank);
    data2=data;  %CHANGE THIS BACK
    plot(vroll,data2,'r')
    %axis tight
    if vroll(end)<vroll(1)
        xlim([vroll(end-1) vroll(2)])
    else
        xlim([vroll(2) vroll(end-1)])
    end
    xlabel(rolllabel)
    ylabel(ycaption)
    titlestr=[titlestr,'    FOM med-rank=',num2str(mrank)];
    title(titlestr)
    
    if SAVEPLOTS
    print(gcf, '-dpng',filename);
    end
    end
    
else


%% 3D line plot

dxp = 1; % interpolation mesh size
droll=roll1-roll0;
if droll<0
    dxp=-dxp;
end

%color map and marker map
cm=colormap(jet(nsteps)); 
mkr={'x','o','d','.','+','^','h','<','*','v','s','>','p'};

% figure (3)
% for i3=0:nsteps-1
%     hold on; 
%     plot3(vstep(ncount==i3),vroll(ncount==i3),data(ncount==i3),'x','Color',cm(i3+1,:),'Marker',mkr{mod(i3,13)+1})
% end
% view([1 -1 1])

cgrid = zeros(nsteps,length(roll0+dxp:dxp:roll0+droll));
xgrid = zeros(nsteps,1);

%% Overlapped traces plot
figure (4)

for i1 = 1:nsteps
    i1=i1-1;
    ind1 = ncount==i1;
    vstepi = vstep(ind1);
    vrolli = vroll(ind1);
    datai  = data(ind1);
    %% remove duplicate entries
    [vrolli,p] = sort(vrolli);
    datai = datai(p,:);
    h = diff(vrolli);
    vrolli(find(h==0)+1)=[];
    datai(find(h==0)+1)=[];
    vstepi(find(h==0)+1)=[];
    
    %% Filtering 
    if FILTERON
        alpha=abs(dxp)/K_FILT;
        % low pass filter
        F_datai = filter(alpha,[1 alpha-1],datai);
    end
    %% Gaussian fit    
    if FITON    
        center=(vrolli(end)+vrolli(1))/2;
        window=vrolli(end)-vrolli(1);
        peaks=1;
        type=1;   %1=gaussian, 2=lorentzian
        [results(i1+1,:),fiterr(i1+1)]=peakfit([vrolli F_datai],center,window,peaks,type);
    end
        
    if FILTERON
        hold on
        plot(vrolli,datai,'Color',cm(i1+1,:),'Marker',mkr{mod(i1,13)+1},'LineStyle','none')
        hold on
        plot(vrolli,F_datai,'Color',cm(i1+1,:),'Marker','none');%,'LineStyle','none')
    else
    hold on
    mkr{mod(i1,7)+1};
    plot(vrolli,datai,'Color',cm(i1+1,:),'Marker',mkr{mod(i1,13)+1});%,'LineStyle','none')
    end
    %% interpolate
    
    yi = roll0+dxp:dxp:roll0+droll; 
    zi = interp1(vrolli,datai,yi); 
    hold on
    %plot(vrolli,datai,'o',yi,zi,'-rx');
        
    %% save to matrix
    cgrid(i1+1,:)=zi;
    xgrid(i1+1)=vstepi(1);
end
axis tight


%% SCAN plot 
figure (5)
[xgrid,sorti]=sort(xgrid);cgrid=cgrid(sorti,:);

surf(yi,xgrid,cgrid,'EdgeColor','none')
colormap hot
%shading interp
%view([1 -1 1])
view([0 0 1])

%imagesc(yi,xgrid,cgrid);

colorbar
caxis([min(min(cgrid)) max(max(cgrid))])
h=colorbar();
ylabel(h,'FOM')
%zaxis([min(min(cgrid)) max(max(cgrid))])

%imagesc(cgrid)
%view([0 0 1])

axis('tight')
if roll_index==1
    xlabel('X');
elseif roll_index==2
    xlabel('Y');
elseif roll_index==3
    xlabel('Tip');
elseif roll_index==4
    xlabel('Rot');
end
    
if step_index==1
    ylabel('X');
elseif step_index==2
    ylabel('Y');
elseif step_index==3
    ylabel('Tip');
elseif step_index==4
    ylabel('Rot');
end

if spectracatch==1
    titlestr=['FOM: Gated Int,    E_0=',num2str(SHIFT),'     \DeltaE=',...
                num2str(DELTA),'    spec. med-rank=',num2str(mrank)];
    title(titlestr);
end

% %% TREND plots
% if FITON
%     figure
%     subplot(7,1,1:2)
%     plot(xgrid,results(:,2),'x')
%     ylabel('Peak Position')
%     subplot(7,1,3:4)
%     plot(xgrid,results(:,3),'x')
%     ylabel('Peak Height')
%     subplot(7,1,5:6)
%     plot(xgrid,results(:,4),'x')
%     ylabel('Peak Width')
%     subplot(7,1,7)
%     plot(xgrid,fiterr,'x')
%     ylabel('% RMS')
% end
% 
end
end
end


% %%
% figure
% 
% plot(data1,'b')
% hold on
% plot(data2,'r')
% grid on
%     axis tight
%     
%     enhance_plot;
%     
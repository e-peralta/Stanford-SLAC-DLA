% ebeam_profiler
% Takes a series of images from select screen and fits a gaussian to each
% image. Stop script with Ctrl+C.
% Would be useful to implement the GUI, but not vital.
% 11/18/2010 - KenSoong@SLAC.STANFORD.EDU
%
% Added a find centroid feature and now plots centroid value instead of
% pixels/intensity
% Centroid mislocation bug fixed. Now works properly with ROI selection.
% 02/21/2011 - KenSoong@SLAC.STANFORD.EDU
%
% Removed A&B software dependancies and made compatible with IMAQ toolbox
% 06/13/2011 - KenSoong
%
% Break out of infinite loop with the 'k' key!
% 08/10/2011 - KenSoong
%
% Now uses ImageAnalyze_E which implemented Asymmetric Gaussian fit,
% with a corrected gaussian expression (was missing factor of 1/2 in exp)
% and reports back the FWHM instead of sigma.
% 12/06/12 - eperalta

%UI updated
%07/27/2013 - Alex Kwiatkowski (akwiat@slac.stanford.edu)
close all
%---Choose screen and nimages---%
%SCR='Aero'
%SCR='PROF4250'
%SCR = 'ITRYAG2'; % change screen name manually!
%SCR = 'ITRYAG4'; % change screen name manually!

%SCR = 'SPYCAM'; % change screen name manually!
%SCR='IFEL D/S';
%SCR = 'IFEL U/S'; % Sony - temp replacement
%SCR = 'PROF0350'; % change screen name manually!
%SCR = 'PROF0585'; % change screen name manually!
SCR = 'PROF0790'; % change screen name manually!
%-------------------------------%

trigsource = 'Line2';
trigmode = 'On';
manualmode=1;
ROIcall=1;


% camera settings
nimages = 1; 
%rgain = 1000; % raw gain

%---Choose number of shots---%
nsteps = 20000; %20000

%----Startup parameters
clear MSX MSX0 MSY MSY0 SSX SSY BVAL Bset MCX MCY SX SY Sig Chi DSig ax cxn cyn bkout;

%global bkout
%bkout = 0;

PLOTEACHIMG=0;
DEBUG = 0;
VERBOSE = 1;
FASTFIT = 0;
%{
[output QD SCR FASTFIT DEBUG nimages nsteps VERBOSE PLOTEACHIMG ...
  QSMn QSMx usescp ROIcall]=quadscan_presets_gui;
%}

SIG_OPTIONS=optimset('Display','off',...
    'TolFun',1e-05,...
    'TolX',1e-05);

%----image acquisition settings
FPT=nimages*3; % frames per trigger (51 max, must be a multiple of 3)

bd=0; % Bit depth of camera 0=8-bit, 1=16-bit
%max gain is set by bit depth, if bd=8bit, max gain is 1023, 16bit-->511
MaxGainV=(-bd+2).*511.*ones(1,15);
%MaxGainV=MaxGainV.*0.5;
%MaxGainV=[511,511,511,511,511,511,511,511,511,511,511,511,511,511,511];

%----Quad Scan Acquisition Parameters
adaptor='ni';             %internal name for PCI-1409 card

%----Camera Screen data
clear Screen; Screens;                    %Load screen data from spreadsheet
ws=GetStructIndex(Screen,'Name',SCR);%Which screen is this
CalX=Screen(ws).XCal;               %microns per pixel - X direction (camera coords)
CalY=Screen(ws).YCal;               %microns per pixel - Y direction (camera coords)
ROIV=Screen(ws).ROIV;               %ROI [lowx, lowy, widthx, widthy]

tr=Screen(ws).NROT;                 %xy transposed?
xm=Screen(ws).IsXM;                 %x mirrored?
ym=Screen(ws).IsYM;                 %y mirrored?
IsXyb=Screen(ws).Xybn;              %Is this a Xybion? Trigger offset is different
imaqCH=Screen(ws).IMAQ;             %NI-IMAQ card input channel


if IsXyb==2 % Is it a gigE camera?
    %gigECh=Screen(ws).TVCh;
    gigECh=findcam(SCR);
    %gigECh=8;
else % then it's a regular video cam
    err=SetTVChannel(Screen(ws).TVCh);
    % temporary work around for ARDBW75/ARDBW57 TV channel discrepancies
    if strcmp(SCR,'IFEL U/S') || strcmp(SCR,'IFEL D/S')||strcmp(SCR,'ITRYAG2')
        [ret,compname]=system('hostname');
        if strcmp(compname(1:7),'ardbw75')
            imaqCH=2;
        end
    end    
end

%
%----Configure framegrabber
%
if IsXyb==2 %Is it a gigE camera?
    %{
    %axes_size=100*[0.5 1.16 4.51 3.01];
    %axes_size=[1 1 659 494];
    axes_size=[1 1 1034 779]; 
    %axes_size=[1 1 640 480]; 
    actxf=figure('Visible','off');
    h=actxcontrol('ActiveGIGE.ActiveGIGE.1',axes_size,actxf);
    h.set('Camera',gigECh);
    %if(SCR~='PROF1130')
    h.set('ExposureTimeAbs',80000);
    %end
    %h.set('GainRaw',511);
    gain=MaxGainV(gigECh+1);
    h.set('Format',0); %16-bit mono
    h.set('GainRaw',1023);
    h.set('OffsetX',ROIV(1));
    h.set('OffsetY',ROIV(2));
    h.set('SizeX',ROIV(3));
    h.set('SizeY',ROIV(4));
    h.set('AcquisitionMode','Continuous');
    h.set('TriggerSelector','AcquisitionStart');
    h.set('TriggerSource','Line1');
    %}
    vformat={'Mono12Packed','Mono16','Mono8','YUV422Packed'};
    vid = videoinput('gige', gigECh, vformat{3}); % choose a format
    src = getselectedsource(vid);
    vid.FramesPerTrigger = 1;
    src.ExposureTimeAbs = 80000;
    %set(vid.Source,'DigitalShift',4);   %ERC added 10-sep-2011
    %src.AllGainRaw = rgain;
    vid.ROIPosition = [ROIV(1) ROIV(2) ROIV(3) ROIV(4)];
    vid.TriggerRepeat = Inf;
    src.AcquisitionStartTriggerSource = trigsource; % correct line?
    triggerconfig(vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
    %src.AcquisitionStartTriggerMode = 'Off';
    src.AcquisitionStartTriggerMode = trigmode;
    %triggerconfig(vid, 'immediate');
    imaqmem(1500000000); % memeory for streaming - 1.5GB
else % then it must be a regular video
    
    vh=videoinput('ni',1); %Define handle to IMAQ card
    set(vh,'SelectedSourceName',['Channel ',num2str(imaqCH)]); %Select Input Channel
    %set(vh,'ROIPosition',ROIV);
    set(vh,'FramesPerTrigger',FPT);
    set(vh,'Timeout',30);
    %triggerconfig(vh,'hardware','risingedge','external0')
    triggerconfig(vh,'immediate')
    src = getselectedsource(vh);
	src.Brightness = 18;
    
end

% Continally probe the camera until it responds
if IsXyb==2
    failed=1;
    while failed==1
        start(vid); pause(0.5); % start the video acquisition
        [isrunning(vid), islogging(vid)];
        preview(vid); pause(0.5);
        stoppreview(vid); closepreview(vid);
        im1=peekdata(vid,1); flushdata(vid,'all');
        if isempty(im1)
            stop(vid); delete(vid);
            vformat={'Mono12Packed','Mono16','Mono8','YUV422Packed'};
            vid = videoinput('gige', gigECh, vformat{3}); % choose a format
            src = getselectedsource(vid);
            vid.FramesPerTrigger = 1;
            %src.ExposureTimeAbs = 80000;
            %src.AllGainRaw = rgain;
            vid.ROIPosition = [ROIV(1) ROIV(2) ROIV(3) ROIV(4)];
            vid.TriggerRepeat = Inf;
            src.AcquisitionStartTriggerSource = trigsource;
            triggerconfig(vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
            src.AcquisitionStartTriggerMode = trigmode;
            imaqmem(1500000000);
        else
            failed=0;
        end
    end
end

%%---------------------------------------------------------------------
%    Interactive ROI selection 
%%-------------------------------------------

if ROIcall
    npreviews=1;
    for il=1:npreviews
    
        if IsXyb==2 % ***** Is the camera a Basler gigE? ******
            %{
            h.Grab;
            M=transpose(h.GetRawData);
            if FPT > 1
                for i=2:3 %FPT    loop over FPT images
                    h.Grab;
                  img=transpose(h.GetRawData);
                  M=cat(4,M,img);
                end
            end
            %}
            M = peekdata(vid,1); flushdata(vid,'all');
            if FPT > 1
                for i=2:3 %loop over FPT images
                    img=peekdata(vid,1);  flushdata(vid,'all');
                    M=cat(4,M,img);
                end
            end
            flushdata(vid,'all');
            %stop(vid)
        else       % *****  camera is a regular video ******
        %bd=0;
        %set(vh,'FramesPerTrigger',3);
        %triggerconfig(vh,'hardware','risingedge','external0')
            start(vh) %Make card active
            pause(1) %wait for acquisition
            M=getdata(vh); %Get frame(s) from buffer%
        end
    
        MSUM1=squeeze(sum(sum(sum(M(:,:,1,1),1),2),4));
        MSUM2=squeeze(sum(sum(sum(M(:,:,1,2),1),2),4));
        MSUM3=squeeze(sum(sum(sum(M(:,:,1,3),1),2),4));
        [dummy igood]=max([MSUM1 MSUM2 MSUM3]);
        
        img=double(squeeze(M(:,:,1,igood))); 
     
        % keep highest pixel values across preview scan 
        if il==1
            img0=img;
            figure;
        else
            for ii=1:length(img0(:,1))
                for jj=1:length(img0(1,:))
                    if(img(ii,jj)>img0(ii,jj))
                        img0(ii,jj)=img(ii,jj);
                    end
                end
            end
        end
    
        imagesc(img0);
        if(PLOTEACHIMG)drawnow;end
    
    end
    pause(.5);
    title(['Select Window of interest, then double click']);

    % if IsXyb~=2
    % set(vh,'FramesPerTrigger',FTP);
    % triggerconfig(vh,'hardware','risingedge','external0')
    % end;
    [dummy Box]=imcrop; 
    %QuadsforROI(wq,Box);
    %Box = [1 1 770 770]
    cXi=2*round(Box(1)/2);
    cYi=2*round(Box(2)/2);
    cXf=cXi+2*round(Box(3)/2)-1;
    cYf=cYi+2*round(Box(4)/2)-1;
end

if cXi<1, cXi=1; end
if cYi<1, cYi=1; end

%
% Get the background image
%
if(DEBUG==0)
    if(IsXyb==2) %camera is a Basler gigE
        flushdata(vid,'all');
        
        if manualmode==0
            try
                pause(1);
                lcaPut('TA03:MISC:1037:UVBMSTOP',1);   %Block e-beam for bkg image
                pause(1) %wait for acquisition
            catch
                inputdlg('Please insert the beam stop then hit <CR>');
            end
        else
            inputdlg('Please insert the beam stop then hit <CR>\n');
        end
        
        %{
        h.Grab;
        M=transpose(h.GetRawData);
        if FPT > 1
            for i=2:FPT % loop over FPT images
                h.Grab;
                img=transpose(h.GetRawData);
                M=cat(4,M,img);
            end
        end
        % average FPT images to make bkg
        bkg=sum(double(M(cYi:cYf,cXi:cXf,1,1:3:FPT)),4)/(FPT/3);        
        %bkg=sum(double(M(:,:,1,1:3:FPT)),4)/(FPT/3);       
        %}
        %start(vid)
        M=peekdata(vid,1); flushdata(vid,'all');
        if FPT > 1
            for i=2:FPT % loop over FPT images
                img=peekdata(vid,1); flushdata(vid,'all');
                M=cat(4,M,img);
            end
        end
        % average FPT images to make bkg
        bkg=sum(double(M(cYi:cYf,cXi:cXf,1,1:3:FPT)),4)/(FPT/3);
        
        if manualmode==0
            try
                pause(1)
                lcaPut('TA03:MISC:1037:UVBMSTOP',2);
                pause(1)
            catch
                inputdlg('Please remove the beam stop then hit <CR>');
                
            end
        else
            inputdlg('Please remove the beam stop then hit <CR>');
            
        end
        
        
    else %camera is a regular video camera: IsXyb = 0 or 1
        pause(0.3);
        if manualmode==0
            try
                pause(1);
                lcaPut('TA03:MISC:1037:UVBMSTOP',1);   %Block e-beam for bkg image
                pause(1) %wait for acquisition
            catch
                inputdlg('Please insert the UV beam stop then hit <CR>');
                
            end
        else
            inputdlg('Please insert the UV beam stop then hit <CR>');
            
        end
        
        pause(1) %wait for acquisition
        start(vh);
        M=getdata(vh); %Get frame(s) from buffer
        bkg=sum(double(M(cYi:cYf,cXi:cXf,1,1:3:FPT)),4)/(FPT/3);
        %bkg=sum(double(M(:,:,1,1:3:FPT)),4)/(FPT/3);
        if manualmode==0
            try
                pause(1)
                lcaPut('TA03:MISC:1037:UVBMSTOP',2);
                pause(1)
            catch
                inputdlg('Please remove the UV beam stop then hit <CR>');
                
            end
        else
            inputdlg('Please remove the UV beam stop then hit <CR>');
            
        end
    end
elseif(DEBUG==1)
    %  start(vh);
    %  bkg=getdata(vh);
%     start(vh);
%     M=getdata(vh); %Get frame(s) from buffer
    if(IsXyb==2) %camera is a Basler gigE
        %{
        h.Grab;
        M=transpose(h.GetRawData);
        if FPT > 1
            for i=2:FPT % loop over FPT images
                h.Grab;
                img=transpose(h.GetRawData);
                M=cat(4,M,img);
            end
        end
        %}
        %start(vid)
        M=peekdata(vid,1); flushdata(vid,'all');
        if FPT > 1
            for i=2:FPT % loop over FPT images
                img=peekdata(vid,1); flushdata(vid,'all');
                M=cat(4,M,img);
            end
        end
    else %camera is a regular video camera: IsXyb = 0 or 1
        start(vh);
        M=getdata(vh); %Get frame(s) from buffer
    end
    % average FPT images to make bkg
    bkg=sum(double(M(cYi:cYf,cXi:cXf,1,1:3:FPT)),4)/(FPT/3);      
else
    bkg=zeros(cYf-cYi+1,cXf-cXi+1);
    %bkg=zeros(size(ROIV(3),ROIV(4)));
end

if(strcmp(SCR,'TESTSCRN'))
    bkg=12*rand(size(bkg));  %5% noise
end

% open initial plot
scrsz=get(0,'ScreenSize');
%figure('Position',[25 25 1200 1000]);
scrsz=[1 1 .98*scrsz(3) .98*scrsz(4)];
h=figure('NumberTitle','off','Name','Ebeam Profiler','Position',scrsz);
user_clicked_stop = 0;
user_clicked_quit = 0;
buttonhandle = uicontrol(h, 'Position', [500 40 100 40], 'String', 'Stop', 'Callback', 'user_clicked_stop = 1;');
buttonhandle2 = uicontrol(h, 'Position', [700 40 100 40], 'String', 'Quit', 'Callback', 'user_clicked_quit = 1;');

%%%%%
%%Modified by Alex Kwiatkowski (akwiat@slac.stanford.edu) 7/27/2013
%Buttons added to Stop or Quit the program
%Clicking the buttons activate the 'user_clicked_stop' or
%'user_clicked_quit' flags, which are checked by the acquisition and analysis while loop.
%Control pauses after the while loop until the 'user_clicked_quit' flag is
%set by the Quit button

%Commented out the calls to the lurker function (believe they are no longer
%necessary)

%%%%%

%figure('Position',scrsz);

%-----Start acquisition-------%

il = 0;
while il<nsteps && ~user_clicked_stop && ~user_clicked_quit
  il = il+1;
  %lurker
  
  BVAL = [1:il];
  NEWSPOT=1;      %1 = dual axis spot size plot
  if IsXyb==2 % ***** Is the camera a Basler gigE? ******
        %{
        h.Grab;
        M=transpose(h.GetRawData);
        if FPT > 1
            for i=2:FPT % loop over FPT images
                h.Grab;
                img=transpose(h.GetRawData);
                M=cat(4,M,img);
            end
        end
        %}
        M=peekdata(vid,1); flushdata(vid,'all');
        if FPT > 1
            for i=2:FPT % loop over FPT images
                img=peekdata(vid,1); flushdata(vid,'all');
                M=cat(4,M,img);
            end
        end
    else       % *****  camera is a regular video ******
        bd=0;
        start(vh) %Make card active
        pause(1) %wait for acquisition
        M=getdata(vh); %Get frame(s) from buffer%
    end
    imn=0;
    for jl=1:3:FPT   %Loop over images stored in PCI-1405 image buffer
        MSUM1=squeeze(sum(sum(sum(M(cYi:cYf,cXi:cXf,1,jl),1),2),4));
        MSUM2=squeeze(sum(sum(sum(M(cYi:cYf,cXi:cXf,1,jl+1),1),2),4));
        MSUM3=squeeze(sum(sum(sum(M(cYi:cYf,cXi:cXf,1,jl+2),1),2),4));
        [dummy igood]=max([MSUM1 MSUM2 MSUM3]);
        
        img=double(squeeze(M(cYi:cYf,cXi:cXf,1,jl+(igood-1))))-double(bkg); %Background subtraction
        imn=imn+1;
        showplots=0;
        % This is what does the gaussian fitting:
        ImageAnalyze_E; %gives c1,c2,hor,ver,sigx,sigy
                
        if(VERBOSE)
            subplot(4,4,1:2);imagesc(img)
            %title({'Background Subtracted & Rotated Image'})
            
            title({[SCR,' Image with ROI [',num2str(round(Box)),' ]']})
            %colorbar
            [RN CN]=size(img);
            V=reshape(img,RN*CN,1);
            if(PLOTEACHIMG)
                drawnow;
            end
            
            % add ROI start point to centroid
            cxni = cXi+c1(2);
            cyni = cYi+c2(2);
            
            subplot(4,4,3:4);
            plot(hor);
            hold on;
            sumx=sum(f(c1,1:length(hor)));
            plot(f(c1,1:length(hor)),'r');
            hold on
            line([c1(2)-c1(3)*sqrt(2*log(2)) c1(2)+c1(5)*sqrt(2*log(2))],[c1(4)+c1(1)/2 c1(4)+c1(1)/2], 'LineWidth',2,'Color','k');
            
            
            grid on;
            %title(sprintf('Line Out X. sx = %f',sigx))
            title({sprintf('FWHM_x = %f [pix] %f [um]',sigx/CalX ,sigx),sprintf('Centroid at %f [pix] %f [um]',cxni,cxni*CalX)});
            ax1 = gca;
            set(ax1,'YAxisLocation','right');
            ylabel('X line out')
            axis tight

            if(PLOTEACHIMG)
                drawnow;
            end
            hold off;
            
            
            subplot(4,4,7:8);
            plot(ver);
            hold on;
            sumy=sum(f(c2,1:length(ver)));
            plot(f(c2,1:length(ver)),'g');
            hold on
            line([c2(2)-c2(3)*sqrt(2*log(2)) c2(2)+c2(5)*sqrt(2*log(2))],[c2(4)+c2(1)/2 c2(4)+c2(1)/2],'LineWidth',2,'Color','k');
            
            grid on;
            
            %title(sprintf('Line Out Y. sy = %f',sigy))
                        
            title({sprintf('FWHM_y = %f [pix] %f [um]',sigy/CalY,sigy),sprintf('Centroid at %f [pix] %f [um]',cyni,cyni*CalY)});
            ax1 = gca;
            set(ax1,'YAxisLocation','right');
            ylabel('Y line out')
            axis tight
            
            if(PLOTEACHIMG)
                drawnow;
            end
            hold off;
            
        end
        
        SX(il,imn)=sigx;
        SY(il,imn)=sigy;
        ChiX(il,imn)=ChiSqX;
        ChiY(il,imn)=ChiSqY;
        
        cxn(il)=cxni; %cXi+c1(2);
        cyn(il)=cyni; %cYi+c2(2);
        
        if ~NEWSPOT
            subplot(4,4,9:16);
            hold on;
            plot(il,sigx,'r.',il,sigy,'g.')
        end
        
    end      %LOOP OVER IMAGES
    
     subplot(4,4,5:6);

     [ax,h1,h2]=plotyy(BVAL,cxn,BVAL,cyn);
     set(h1,'Marker','.','LineStyle','none','Color',[1 0 0]);
     set(h2,'Marker','.','LineStyle','none','Color',[0 0.5 0]);
     set(ax(1),'YColor',[1 0 0]);
     set(ax(2),'YColor',[0 0.5 0]);
     hold(ax(1),'on')
     hold(ax(2),'on')
     ylabel('Centroid Position [um]')
     if strcmp(SCR,'PROF0790')
     title({['E Jitter = ',num2str(std(cxn)),'[pix],  ',num2str(std(cxn)*CalX),'[um],    ',num2str(std(cxn)*CalX/5187.8),'[MeV]']});
     else
         title({['E Jitter = ',num2str(std(cxn)),'[pix],  ',num2str(std(cxn)*CalX),'[um]']});
     end
     grid on;
     % added code below to autoscale axis - EAP 12/04/12
     if length(cxn)>10
        minlim=round(min(cxn(end-10:end))-5);
        maxlim=round(max(cxn(end-10:end))+5);
        del=(maxlim-minlim)/4;
        ylim(ax(1),[minlim-2,maxlim+2])
        set(ax(1),'YTick',[minlim:del:maxlim])
        
        minlim=round(min(cyn(end-10:end))-5);
        maxlim=round(max(cyn(end-10:end))+5);
        del=(maxlim-minlim)/4;
        ylim(ax(2),[minlim-2,maxlim+2])
        set(ax(2),'YTick',[minlim:del:maxlim])
     end
     if(PLOTEACHIMG)
       drawnow;
     end
    
    
    if (NEWSPOT)
        subplot(4,4,9:16);
        [ax,h1,h2]=plotyy(BVAL,SX',BVAL,SY');
        set(h1,'Marker','.','LineStyle','none','Color',[1 0 0]);
        set(h2,'Marker','.','LineStyle','none','Color',[0 0.5 0]);
        set(ax(1),'YColor',[1 0 0]);
        set(ax(2),'YColor',[0 0.5 0]);
        set(get(ax(1),'Ylabel'),'String','FWHM_x [um]') 
        set(get(ax(2),'Ylabel'),'String','FWHM_y [um]') 
        hold(ax(1),'on')
        hold(ax(2),'on')
               
        %ylim(ax(1),[min(min(SX),5),max(max(SX),200)])
        %ylim(ax(2),[min(min(SY),5),max(max(SY),200)])
        
        % added this code below to have a continues update of the y-axis 
        % which keeps the last 10 measured values in range - EAP
       if length(SX)>10
        minlim=round(min(SX(end-10:end))-5);
        maxlim=round(max(SX(end-10:end))+5);
        del=(maxlim-minlim)/4;
        ylim(ax(1),[minlim-2,maxlim+2])
        set(ax(1),'YTick',[minlim:del:maxlim])
        
        minlim=round(min(SY(end-10:end))-5);
        maxlim=round(max(SY(end-10:end))+5);
        del=(maxlim-minlim)/4;
        ylim(ax(2),[minlim-2,maxlim+2])
        set(ax(2),'YTick',[minlim:del:maxlim])
       end
       drawnow;
    end
    
    %cut out points where mean>2sigma
    MSX0(il)=mean(SX(il,:));
    SSX(il)=std(SX(il,:));
    UseX=find(abs(SX(il,:)-MSX0(il))<=1*SSX(il));
    MCX(il)=mean(ChiX(il,UseX));
    MSX(il)=mean(SX(il,UseX));
    SSX(il)=std(SX(il,UseX));
    
    MSY0(il)=mean(SY(il,:));
    SSY(il)=std(SY(il,:));
    UseY=find(abs(SY(il,:)-MSY0(il))<=1*SSY(il));
    MCY(il)=mean(ChiY(il,UseY));
    MSY(il)=mean(SY(il,UseY));
    SSY(il)=std(SY(il,UseY));
    
    
    
    if(VERBOSE)
        if ~NEWSPOT
            subplot(4,4,9:16),
            plot(il,MSX0(il),'ro',il,MSY0(il),'go',...
                il,MSX0(il)+2*SSX(il),'r+',il,MSY0(il)+2*SSY(il),'g+',...
                il,MSX0(il)-2*SSX(il),'r+',il,MSY0(il)-2*SSY(il),'g+',...
                il,MSX(il),'k*',il,MSY(il),'k*')
        else
            axes(ax(1))
            grid on
            hold on
            plot(BVAL,MSX0,'ro',...
                BVAL,MSX0+2*SSX,'r+',...
                BVAL,MSX0-2*SSX,'r+',...
                BVAL,MSX,'r*-')
            hold off
            axes(ax(2))
            grid on
            hold on
            plot(BVAL,MSY0,'go',...
                BVAL,MSY0+2*SSY,'g+',...
                BVAL,MSY0-2*SSY,'g+',...
                BVAL,MSY,'g*-')
            hold off
            drawnow
        end
    end
end

set(h,'Name', 'Stopped.');

%subplot(4,4,9:16)
hold on
title({['<X> [um] = ',num2str(mean(cxn)*CalX),' +- ',num2str(std(cxn)*CalX),',     <Y> [um] = ',num2str(mean(cyn)*CalY),' +- ',num2str(std(cyn)*CalY),',      <FWHM_x> [um] = ',num2str(mean(MSX)),' +- ',num2str(std(MSX)),'     <FWHM_y> [um] = ',num2str(mean(MSY)),' +- ',num2str(std(MSY))]});


% fprintf('Centroid [pixels] at %f %f\n',mean(cxn),mean(cyn));
% fprintf('Centroid [microns] at %f %f\n',mean(cxn)*CalX,mean(cyn)*CalY);
% fprintf('<Sx> [um] = %f+-%f\n',mean(MSX),std(MSX));
% fprintf('<Sy> [um] = %f+-%f\n',mean(MSY),std(MSY));
% fprintf('<X> [pixels] = %f+-%f\n',mean(cxn),std(cxn));
% fprintf('<Y> [pixels] = %f+-%f\n',mean(cyn),std(cyn));

%
if IsXyb==2
    flushdata(vid,'all')
    stop(vid)
    delete(vid)
end

   savefolder=['V:\ARDB\E163\Data\',datestr(now,'yymmdd'),'\'];
    if exist(savefolder,'dir')==0
        mkdir(savefolder);
    end
   filename='Eprofile_1';
   ImgNum=0;
   numD=1;
   
   while exist([savefolder,filename,'.png'],'file')~=0
      ImgNum=ImgNum+1;
      if ImgNum>10
          numD=2;
      end
      filename=[filename(1:end-numD),num2str(ImgNum,'%g')];
   end
   
   %enhance_plot;
   print(gcf, '-dpng',[savefolder,filename]);
   
fprintf('\n%s:\n',filename);
fprintf('<X> [um] = %f+-%f\n',mean(cxn)*CalX,std(cxn)*CalX);
fprintf('<Y> [um] = %f+-%f\n',mean(cyn)*CalY,std(cyn)*CalY);
fprintf('<FWHM_x> [um] = %f+-%f\n',mean(MSX),std(MSX));
fprintf('<FWHM_y> [um] = %f+-%f\n',mean(MSY),std(MSY));
if strcmp(SCR,'PROF0790')
    fprintf('<dE> [MeV] = %f\n',std(cxn)*CalX/5187.8);    
end
  

%Waiting for the 'user_clicked_quit' flag to be set to 1 by 
%clicking on the Quit button (see button initialization above) 
while ~user_clicked_quit
    pause(.1)
end

lh = findobj('title','lurker')
lh = gcf
close(lh);

 close all  
%% print into today's folder
% printdefault = 1 ;
% if printdefault == 1
%     d1 = datevec(date);
%     fdr=sprintf('v:\\ardb\\E163\\Data\\%02d%02d%02d\\',d1(1)-2000,d1(2),d1(3));
%     runno=input('Which profile scan number is this?');
%     eval(['print -dbitmap ' fdr 'EProfile_',num2str(runno),'.bmp;']);
% end


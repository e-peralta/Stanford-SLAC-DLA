% Script QUADSCAN
%
% Matlab script to conduct quad scans at the NLCTA
%
%------------------------------------------------------------
% List of Changes and Bug Fixes
%------------------------------------------------------------
%
% 090203 england@slac.stanford.edu
% Modified from quadscan.m (ecolby@slac.stanford.edu)
% to include use of gigE cameras.
%
% 090213 cmcg@slac.stanford.edu
% Modified to invoke quadscan_gui when scan is done
% and save data to struct indexed by runnum
%
% 090609 cmcg@slac.stanford.edu
% Modified to run ExtractTwissThickLens.m
%
% 091103 cmcg@slac.stanford.edu
% Modified to use ExtractTwissThickLens_Weighted.m
% This code weights each point in the quad scan
% and computes transfer matrix from quad to screen including inner quads
%
% 100324 jnelson
%  added "NEWSPOT" support - changes how the bottom spot size plot is drawn
%  moved display of alphax etc to within the if statement for save/don't
%  save because don't save was generating errors.
%
% 100707 eperalta
% added option to interactively select a region of interest from the screen
%
% 110613 kensoong
% removed all third party A&B software dependancies 
% quadscan's video acquisition now run using Matlab's IMAQ toolbox
%
% 111015 kensoong
% added error handling for occasional slow camera acquisition rate
%
% Dependencies:
%   Requires Image Acquisition Toolbox
%   Requires AIDA server and Java libraries
%   No longer requires ActiveGIGE (A&B Software)
%
%------------------------------------------------------------

close all;  % get rid of all prior figures
clear MSX MSX0 MSY MSY0 SSX SSY BVAL Bset MCX MCY SX SY Sig Chi DSig ax;
SIG_OPTIONS=optimset('Display','off',...
    'TolFun',1e-05,...
    'TolX',1e-05);

% % ----------------------------------------------------%
% %          These values are now set from GUI          %
% %          but left here for the time being...        %
% % ----------------------------------------------------%
% 
% VERBOSE=1;      %1=display all graphs; 0=display sigmas only
% PLOTEACHIMG=1;  %1=plots each image for a given quad setting, 0=only plot once per quad setting
% DEBUG=0;        %1=Disables magnet, 2=Disable magnets and laser stopper commands
% PLOTPROPTWISS=0;%1=plots propagated twiss values in selected region of beamline
% FASTFIT=0;      %For Gaussian fitting: 1=Use fast POLYFIT-based fitting; 0=use LSQNONLIN-based fitting
NEWSPOT=1;      %1 = dual axis spot size plot
% 
% if (~exist('scanallquads'))
%     %----Quad scan setup
%     SCR='PROF1550';
%     %SCR='PROF4250';
%     %SCR='PROF1120';
%     %SCR='PROF0790';
%     if(DEBUG==0)
%         QD ='QUAD0530';    %'QUAD1920'
%     else
%         QD ='TESTQUAD';         %Name of quad to use for quad scan
%     end
% end
% 
% %Number of quad K values to use in scan
% nsteps=10;
% FPT=51; % number of frames per trigger (51 max, must be a multiple of 3)


%----------------------------------------------------%
%     Calls the gui to set quad scan settings     %
%----------------------------------------------------%

brokenflipper=1;

[output QD SCR FASTFIT DEBUG nimages nsteps VERBOSE PLOTEACHIMG ...
  QSMn QSMx usescp ROIcall npreviews]=quadscan_presets_gui;
tic;
FPT=nimages*3; % frames per trigger (51 max, must be a multiple of 3)
wtime = 5*FPT/10+1; % time to wait if read from buffer fails

bd=0; % Bit depth of camera 0=8-bit, 1=16-bit
%max gain is set by bit depth, if bd=8bit, max gain is 1023, 16bit-->511
MaxGainV=(-bd+2).*511.*ones(1,15);
%MaxGainV=[511,511,511,511,511,511,511,511,511,511,511,511,511,511,511];

%----------------------------------------------%
%  Get run number for saved data files
%----------------------------------------------%
if(DEBUG==0)
    GetNewRunNumber;
    %quad scan run (qrun) and disp scan run (drun) are both saved in same file
    qrun=qrun+1;
    save(filerun,'qrun','drun','prun','brun','pmqrun');
    FILENAME = ['QS_MATLAB_',num2str(qrun),'_',QD(5:8)];   %Filename to save data
    sf=fopen([DEFAULTPATH,FILENAME,'.dat'],'w+t');
else
    DEFAULTPATH='V:\ARDB\E163\Data\TEST\';
    ROOTNAME=DEFAULTPATH;
    FILENAME='debug.dat';
    qrun=1;
end

%
%----Quad Scan Acquisition Parameters
adaptor='ni';             %internal name for PCI-1409 card
%----Read machine configuration data
clear Quad; Quads;                   %Load quad data from spreadsheet
wq=GetStructIndex(Quad,'Name',QD);   %Which quad is this
% Magnet parameters
quad =Quad(wq).Numb;       %Number of quad to use in quad scan
%now set this from gui output
%lowb =Quad(wq).QSMn;       %lower bound for quad scan (BDES)
%highb=Quad(wq).QSMx;       %upper bound for quad scan (BDES)
lowb=QSMn;
highb=QSMx;
l=Quad(wq).Leff;           %effective quad length
GetCurrentQuadSettings;     %Get BACT values for all quads prior to scan (uses AIDA)

%save current quad settings in case of an AIDA crash
initquadval = Quad(wq).BACT;
save runinitial.mat quad initquadval

%Checks upstream quads to see if any are on between the scanning quad and
%the screen
CheckUSQuads(QD,SCR,Quad);

clear Screen; 
Screens;                    %Load screen data from spreadsheet
ws=GetStructIndex(Screen,'Name',SCR);%Which screen is this
CalX=Screen(ws).XCal;               %microns per pixel - X direction (camera coords)
CalY=Screen(ws).YCal;               %microns per pixel - Y direction (camera coords)
ROIV=Screen(ws).ROIV;               %ROI [lowx, lowy, widthx, widthy]
d=Screen(ws).ZPos-Quad(wq).ZPos;    %compute quad-to-screen distance
dthick=Screen(ws).ZPos-Quad(wq).ZPos-Quad(wq).Leff/2;
tr=Screen(ws).NROT;                 %xy transposed?
xm=Screen(ws).IsXM;                 %x mirrored?
ym=Screen(ws).IsYM;                 %y mirrored?
IsXyb=Screen(ws).Xybn;              %Is this a Xybion? Trigger offset is different
imaqCH=Screen(ws).IMAQ;             %NI-IMAQ card input channel
%{
if(IsXyb==1)
    igood=3;
else
    igood=1;
end
%}
aidainit;               %Setup AIDA
BeamE = get_E;          %Get chicane BACT value to determine beam energy
%
if(BeamE==-99)
    disp('Warning-- chicane BACT reading failed, using E=60 MeV')
    BeamE=60;
end
%
fprintf('\n *******************************************************')
fprintf('\n Performing quad scan with quad %s and screen %s',QD,SCR)
fprintf('\n Separation distance: d=%4.2f m',d)
fprintf('\n Quad scan range: %6.2f to %6.2f',lowb,highb)
fprintf('\n Saving data to file %s%s',DEFAULTPATH,FILENAME)
fprintf('\n *******************************************************\n')
%
delb=(highb-lowb)/(nsteps-1);
%
%----Set proper TV channel on demodulator
%
if IsXyb==2 % Is it a gigE camera?
    
    %!@#$%^&*!#$%!^&!$%^@&*!(#&*@^#@##$%^#$@#%$%$^$Z
    % If the camera list gets messed up run the following command: reset_camlist
    %!@#$%^&*!#$%!^&!$%^@&*!(#&*@^#@##$%^#$@#%$%$^$Z
    
    %reset_camlist;
    gigECh=findcam(SCR);
    %gigECh=3;
else % then it's a regular video cam
    err=SetTVChannel(Screen(ws).TVCh);
    % temporary work around for ARDBW75/ARDBW57 TV channel discrepancies
    if strcmp(SCR,'IFEL U/S') || strcmp(SCR,'IFEL D/S')
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
    src.AllGainRaw = 1023;
    vid.ROIPosition = [ROIV(1) ROIV(2) ROIV(3) ROIV(4)];
    vid.TriggerRepeat = Inf;
    src.AcquisitionStartTriggerSource = 'Line1'; % correct line?
    triggerconfig(vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
    %src.AcquisitionStartTriggerMode = 'On';
    src.AcquisitionStartTriggerMode = 'Off';
    imaqmem(1500000000); % memeory for streaming - 1.5GB
else % then it must be a regular video
    
    vh=videoinput('ni',1); %Define handle to IMAQ card
    set(vh,'SelectedSourceName',['Channel ',num2str(imaqCH)]); %Select Input Channel
    %set(vh,'ROIPosition',ROIV);
    set(vh,'FramesPerTrigger',FPT);
    set(vh,'Timeout',20)
    %triggerconfig(vh,'hardware','risingedge','external0')
    triggerconfig(vh,'immediate')
    
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
            src.ExposureTimeAbs = 80000;
            src.AllGainRaw = 1023;
            vid.ROIPosition = [ROIV(1) ROIV(2) ROIV(3) ROIV(4)];
            vid.TriggerRepeat = Inf;
            src.AcquisitionStartTriggerSource = 'Line1';
            triggerconfig(vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
            src.AcquisitionStartTriggerMode = 'On';
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
    npreviews=3;
    
    for il=1:npreviews
        BVAL(il)=lowb+(il-1)*(highb-lowb)/(npreviews-1);
        %if(DEBUG==0)
        if(DEBUG==0)&&Quad(wq).BDES~=BVAL(il)
           set_quad(quad,BVAL(il));
        end
    
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
            % capture from camera loop
            captured=0;
            while captured==0
                try % implemented as a work around for slow frame capture
                    captured=1;
                    start(vh) %Make card active
                    pause(1) %wait for acquisition
                    M=getdata(vh); %Get frame(s) from buffer%.
                catch
                    stop(vh)
                    captured=0;
                    fprintf('Failed to read from buffer\n');
                    pause(wtime);
                end
            end
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
        title(['Quad BVALs plotted = ',num2str(BVAL)]);
        %if(PLOTEACHIMG)drawnow;end
    
    end
    pause(.5);
    %title({['Quad BVALs plotted = ',num2str(BVAL)];['Select Window of interest, then double click']});
    title(['Select Window of interest, then double click']);

    % if IsXyb~=2
    % set(vh,'FramesPerTrigger',FTP);
    % triggerconfig(vh,'hardware','risingedge','external0')
    % end
      clear BVAL;
    [dummy Box]=imcrop; 
    QuadsforROI(wq,Box);
    %cXi=1+2*round(ROIV(1)/2);
    %cYi=1+2*round(ROIV(2)/2);
    cXi=ceil(2*(Box(1)/2));
    cYi=ceil(2*(Box(2)/2));
    cXf=cXi+2*round(Box(3)/2)-1;
    cYf=cYi+2*round(Box(4)/2)-1;
else
    ROIV=QuadsforROI(wq); 
    %cXi=1+2*round(ROIV(1)/2);
    %cYi=1+2*round(ROIV(2)/2);
    cXi=ceil(2*(Box(1)/2));
    cYi=ceil(2*(Box(2)/2));
    cXf=cXi+2*round(ROIV(3)/2)-1;
    cYf=cYi+2*round(ROIV(4)/2)-1;
end

%%---------------------------------------------------------------------
%    End ROI selection 
%%-------------------------------------------


%------------------------------------------------------
% Get the background image
%------------------------------------------------------
if(DEBUG==0)
    if(IsXyb==2) %camera is a Basler gigE
        flushdata(vid,'all');
        try
            %fail %force fail
            pause(1);
            lcaPut('TA03:MISC:1037:UVBMSTOP',1);   %Block e-beam for bkg image
            pause(1) %wait for acquisition
        catch
            %fprintf('Please insert the UV beam stop then hit <CR>\n')
            %pause
            inputdlg('Please insert the UV beam stop then hit <CR>');
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
        try
            %fail %force fail
            pause(1)
            lcaPut('TA03:MISC:1037:UVBMSTOP',2);
            pause(1)
        catch
            %fprintf('Please remove the UV beam stop then hit <CR>\n')
            %pause
            inputdlg('Please remove the UV beam stop then hit <CR>');
        end
    else %camera is a regular video camera: IsXyb = 0 or 1
        pause(0.3);
        if brokenflipper==0
            lcaPut('TA03:MISC:1037:UVBMSTOP',1);   %Block e-beam for bkg image
        else
            inputdlg('Please insert the UV beam stop then hit <CR>');
        end
        % capture from camera loop
        captured=0;
        while captured==0
            try % implemented as a work around for slow frame capture
                captured=1;
                start(vh);
                pause(1) %wait for acquisition
                M=getdata(vh); %Get frame(s) from buffer%.
            catch
                stop(vh)
                captured=0;
                fprintf('Failed to read from buffer\n');
                pause(wtime);
            end
        end
        bkg=sum(double(M(cYi:cYf,cXi:cXf,1,1:3:FPT)),4)/(FPT/3);
        %bkg=sum(double(M(:,:,1,1:3:FPT)),4)/(FPT/3);
        if brokenflipper==0
            lcaPut('TA03:MISC:1037:UVBMSTOP',2);
        else
            inputdlg('Please remove the UV beam stop then hit <CR>');
        end
    end
elseif(DEBUG==1)

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
        
        % capture from camera loop
        captured=0;
        while captured==0
            try % implemented as a work around for slow frame capture
                captured=1;
                start(vh);
                pause(1)
                M=getdata(vh); %Get frame(s) from buffer%.
            catch
                stop(vh)
                captured=0;
                fprintf('Failed to read from buffer\n');
                pause(wtime);
            end
        end
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


figure('Position',[50 50 1200 1000]);
if(VERBOSE)
    bkg_oriented=orientimage(bkg,Screen(ws));
    subplot(4,4,1:2);
    imagesc(bkg_oriented);
    %colorbar;
    title('Averaged Background Image')
    subplot(4,4,9:16);hold on;
    title('Fitted Spot Sizes vs. Quad BDES')
    xlabel('BDES');ylabel('Spot Sizes [micron]');
    grid;legend('\sigma_x','\sigma_y');
end


%------------------------------------------------------
%   Begin Quad Scan
%------------------------------------------------------

for il=1:nsteps
    %BVAL(il)=lowb+(il-1)*delb;
    BVAL(il)=highb-(il-1)*delb;
    if(VERBOSE)
        fprintf('Changing QUAD value to %8.4f, and taking images\n',BVAL(il))
    end
    %if(DEBUG==0) 
    %if(DEBUG==0)&&il>1 %misses the first quad setting if ROI select is off
    
    if(DEBUG==0)
        set_quad(quad,BVAL(il));
    end
    % Image Aquisition for Baslers
    if IsXyb==2     
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
    %Image Aquisition for cameras on Cable TV network    
    else       
        bd=0;
        % capture from camera loop
        captured=0;
        while captured==0
            try % implemented as a work around for slow frame capture
                captured=1;
                start(vh) %Make card active
                pause(1) %wait for acquisition
                M=getdata(vh); %Get frame(s) from buffer%.
            catch
                stop(vh)
                captured=0;
                fprintf('Failed to read from buffer\n');
                pause(wtime);
            end
        end
    end
    imn=0;
    %Loop over images stored in PCI-1405 image buffer
    for jl=1:3:FPT   
        MSUM1=squeeze(sum(sum(sum(M(cYi:cYf,cXi:cXf,1,jl),1),2),4));
        MSUM2=squeeze(sum(sum(sum(M(cYi:cYf,cXi:cXf,1,jl+1),1),2),4));
        MSUM3=squeeze(sum(sum(sum(M(cYi:cYf,cXi:cXf,1,jl+2),1),2),4));
%         MSUM1=squeeze(sum(sum(sum(M(:,:,1,jl),1),2),4));
%         MSUM2=squeeze(sum(sum(sum(M(:,:,1,jl+1),1),2),4));
%         MSUM3=squeeze(sum(sum(sum(M(:,:,1,jl+2),1),2),4));
        [dummy igood]=max([MSUM1 MSUM2 MSUM3]);
        
        img=double(squeeze(M(cYi:cYf,cXi:cXf,1,jl+(igood-1))))-double(bkg); %Background subtraction
        imn=imn+1;
        ov=VERBOSE;
        VERBOSE=0;
        ImageAnalyze;
        VERBOSE=ov;
        
        if(VERBOSE)
            subplot(4,4,1:2);imagesc(img)
            title('Background Subtracted & Rotated Image')
            %colorbar
            [RN CN]=size(img);
            V=reshape(img,RN*CN,1);
            if(PLOTEACHIMG)
                drawnow;
            end
            subplot(4,4,5:6);
            [nb,xb]=hist(double(V),256);
            semilogy(xb,nb,'k.')
            if bd==1 % 16-bit
                xlim([0,4096]);
            else % 8-bit
                xlim([0,256]);
            end
            title('Pixel Intensity Histogram')
            grid on;
            if(PLOTEACHIMG)
                drawnow;
            end
            
            subplot(4,4,3:4);
            plot(hor);
            hold on;
            plot(f(c1,1:length(hor)),'r');
            grid on;
            title('Line Out X');
            if(PLOTEACHIMG)
                drawnow;
            end
            hold off;
            
            subplot(4,4,7:8);
            plot(ver);
            hold on;
            plot(f(c2,1:length(ver)),'g');
            grid on;
            title('Line Out Y');
            if(PLOTEACHIMG)
                drawnow;
            end
            hold off;
        end
        SX(il,imn)=sigx;
        SY(il,imn)=sigy;
        ChiX(il,imn)=ChiSqX;
        ChiY(il,imn)=ChiSqY;
        
        if ~NEWSPOT
            subplot(4,4,9:16);
            hold on;
            plot(BVAL(il),sigx,'r.',BVAL(il),sigy,'g.')
        end
        
    end      %END LOOP OVER IMAGES
    if (NEWSPOT)
        subplot(4,4,9:16);
        [ax,h1,h2]=plotyy(BVAL,SX',BVAL,SY');
        set(h1,'Marker','.','LineStyle','none','Color',[1 0 0]);
        set(h2,'Marker','.','LineStyle','none','Color',[0 0.5 0]);
        set(ax(1),'YColor',[1 0 0]);
        set(ax(2),'YColor',[0 0.5 0]);
        hold(ax(1),'on')
        hold(ax(2),'on')
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
    %
    if(VERBOSE)
        if ~NEWSPOT
            subplot(4,4,9:16),
            plot(BVAL(il),MSX0(il),'ro',BVAL(il),MSY0(il),'go',...
                BVAL(il),MSX0(il)+2*SSX(il),'r+',BVAL(il),MSY0(il)+2*SSY(il),'g+',...
                BVAL(il),MSX0(il)-2*SSX(il),'r+',BVAL(il),MSY0(il)-2*SSY(il),'g+',...
                BVAL(il),MSX(il),'k*',BVAL(il),MSY(il),'k*')
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
    
    if(DEBUG==0)
        fprintf(sf,'%12.4f %12.4f %12.4f %12.4e %12.4e \n',BVAL(il),MSX(il),MSY(il),MCX(il),MCY(il));
    end
end   %END LOOP OVER QUAD K VALUES
%
if(DEBUG==0)
    set_quad(quad,initquadval);
end
% Stop video and release resources
if IsXyb==2 % is it a gigE camera?
    %{
    h.set('Acquire',false);
    close(actxf); %shuts down the activex control
    %}
    stop(vid);
    flushdata(vid,'all');
else % then it's a regular video cam
    stop(vh)       %stop video card
    flushdata(vh,'all')  %release resources
    %imaqreset;     %disconnect and release all image acq objects
end

%
% Write list of current quad settings (prior to quad scan) as footer
% in data file.
%
for i=1:length(Quad)
    if(DEBUG==0)
        fprintf(sf,'%s %10.4f\n',char(Quad(i).Name),Quad(i).BACT);
    end
    display([char(Quad(i).Name) ' = ' num2str(Quad(i).BACT)]);
end

if(DEBUG==0)fclose(sf);end

toc; % print out elapsed time

%   Extract Twiss Parameters
Bset=BVAL;
ChiCut=5;
GuessX=[0 0 0];
GuessY=[0 0 0];
%ExtractTwiss3;   % Version used until 9/19/09
try
    ExtractTwiss4;% New version has per-point sigma weighting
catch
    display(sprintf('ExtractTwiss4 FAILED!!!'));
end
fignumthinlens=gcf;
Sig=[MSX;MSY]'.*1e-6;  Chi=[SSX;SSY]'.*1e-6;   DSig=[SSX',SSY'].*1e-6;
%twissthicklens=fExtractTwissThickLens(Bset',Sig,Chi,l,dthick)
twissthicklens=fExtractTwissThickLensWeighted(Bset',Sig,DSig,l,dthick,Quad,wq,usescp);

%--------------------------------------------
%     Option to cut bad points
%--------------------------------------------
[output croplist]=quadscan_croppoints_gui;
if length(croplist)>0
  k=1;
  for j=1:length(Bset)
    if ~sum(j==croplist)
      Bset2(k)=Bset(j);
      Sig2(k,:)=Sig(j,:);
      DSig2(k,:)=DSig(j,:);
      k=k+1;
    end
  end
  twissthicklens=fExtractTwissThickLensWeighted(Bset2',Sig2,DSig2,l,dthick,Quad,wq,usescp);
end


fignumthicklens=gcf;
FinalQuadSummary;

if(DEBUG==0)
    eval(['cd ',DEFAULTPATH]);
    print(fignumthinlens, '-djpeg',[FILENAME '_thinlens']);
    print(fignumthicklens, '-djpeg',[FILENAME '_thicklens']);
end


% run gui to input values and determine whether to save
[output irisopen collopen xfitqlty yfitqlty comments saveflag savefilename thicklens]=quadscan_gui;
if (saveflag)
    if(exist(savefilename)==2)
        load(savefilename);
        globalrunnum=length(quadscandata)+1;
    else
        quadscandata=struct([]);
        globalrunnum=1;
    end
    
    quadscandata(globalrunnum).date=clock;
    quadscandata(globalrunnum).runnum=qrun;
    quadscandata(globalrunnum).screen=SCR;
    quadscandata(globalrunnum).calx=CalX;
    quadscandata(globalrunnum).caly=CalY;
    quadscandata(globalrunnum).roiv=ROIV;
    quadscandata(globalrunnum).szpos=Screen(ws).ZPos;
    quadscandata(globalrunnum).xmir=xm;
    quadscandata(globalrunnum).ymir=ym;
    quadscandata(globalrunnum).trans=tr;
    quadscandata(globalrunnum).quad=quad;
    quadscandata(globalrunnum).lowb=lowb;
    quadscandata(globalrunnum).highb=highb;
    quadscandata(globalrunnum).qzpos=Quad(wq).ZPos;
    for i=1:(length(Quad)-1)
        quadscandata(globalrunnum).quadnames(i)=Quad(i).Name;
        quadscandata(globalrunnum).quadvals(i)=Quad(i).BACT;
    end
    quadscandata(globalrunnum).bval=BVAL;
    quadscandata(globalrunnum).sigx=SX;
    quadscandata(globalrunnum).sigy=SY;
    quadscandata(globalrunnum).msx=MSX;
    quadscandata(globalrunnum).msy=MSY;
    quadscandata(globalrunnum).mcx=MCX;
    quadscandata(globalrunnum).mcy=MCY;
    quadscandata(globalrunnum).ssx=SSX;
    quadscandata(globalrunnum).ssy=SSY;
    if thicklens
        quadscandata(globalrunnum).fitmethod=thicklens;
        quadscandata(globalrunnum).alphax=twissthicklens(1,1);
        quadscandata(globalrunnum).alphay=twissthicklens(1,2);
        quadscandata(globalrunnum).betax=twissthicklens(2,1);
        quadscandata(globalrunnum).betay=twissthicklens(2,2);
        quadscandata(globalrunnum).gammax=twissthicklens(3,1);
        quadscandata(globalrunnum).gammay=twissthicklens(3,2);
        quadscandata(globalrunnum).emitx=twissthicklens(5,1);
        quadscandata(globalrunnum).emity=twissthicklens(5,2);
        
        %errors to fit, only used for extracttwiss_thicklens_weighted
        quadscandata(globalrunnum).alphax_err=twissthicklens(1,3);
        quadscandata(globalrunnum).alphay_err=twissthicklens(1,4);
        quadscandata(globalrunnum).betax_err=twissthicklens(2,3);
        quadscandata(globalrunnum).betay_err=twissthicklens(2,4);
        quadscandata(globalrunnum).gammax_err=twissthicklens(3,3);
        quadscandata(globalrunnum).gammay_err=twissthicklens(3,4);
        quadscandata(globalrunnum).emitx_err=twissthicklens(5,3);
        quadscandata(globalrunnum).emity_err=twissthicklens(5,4);
        
    else
        quadscandata(globalrunnum).fitmethod=0;
        quadscandata(globalrunnum).alphax=alpha0x;
        quadscandata(globalrunnum).alphay=alpha0y;
        quadscandata(globalrunnum).salphax=Salpha0x;
        quadscandata(globalrunnum).salphay=Salpha0y;
        quadscandata(globalrunnum).betax=beta0x;
        quadscandata(globalrunnum).betay=beta0y;
        quadscandata(globalrunnum).sbetax=Sbeta0x;
        quadscandata(globalrunnum).sbetay=Sbeta0y;
        quadscandata(globalrunnum).gammax=gamma0x;
        quadscandata(globalrunnum).gammay=gamma0y;
        quadscandata(globalrunnum).emitx=emitx;
        quadscandata(globalrunnum).emity=emity;
        quadscandata(globalrunnum).semitx=Semitx;
        quadscandata(globalrunnum).semity=Semity;
    end
    quadscandata(globalrunnum).irisopen=irisopen;
    quadscandata(globalrunnum).collopen=collopen;
    quadscandata(globalrunnum).xfitqlty=xfitqlty;
    quadscandata(globalrunnum).yfitqlty=yfitqlty;
    quadscandata(globalrunnum).comments=comments;
    %{
    %addition made on 04/01/09
    if(quad==1030 || quad==1070 || quad==1110 || quad==1130)
        twissout=propagatetwiss_from_struct(quadscandata,globalrunnum,10,PLOTPROPTWISS);
        quadscandata(globalrunnum).alphax1130=twissout(1,end);
        quadscandata(globalrunnum).alphay1130=twissout(4,end);
        quadscandata(globalrunnum).betax1130=twissout(2,end);
        quadscandata(globalrunnum).betay1130=twissout(5,end);
        quadscandata(globalrunnum).emitx1130=twissout(3,end);
        quadscandata(globalrunnum).emity1130=twissout(6,end);
    elseif(quad==4210 || quad==4220 || quad==4260 || quad==4270 || quad==4280)
        twissout=propagatetwiss_from_struct(quadscandata,globalrunnum,20,PLOTPROPTWISS);
        quadscandata(globalrunnum).alphax1130=twissout(1,1);
        quadscandata(globalrunnum).alphay1130=twissout(4,1);
        quadscandata(globalrunnum).betax1130=twissout(2,1);
        quadscandata(globalrunnum).betay1130=twissout(5,1);
        quadscandata(globalrunnum).emitx1130=twissout(3,1);
        quadscandata(globalrunnum).emity1130=twissout(6,1);
    end
    %end addition
    %}
    
    %
    % JLN 24Mar2010
    % moved the next 6 lines up to avoid the error if you don't save data
    %
    display(['alpha_x=' num2str(quadscandata(globalrunnum).alphax)]);
    display(['alpha_y=' num2str(quadscandata(globalrunnum).alphay)]);
    display(['beta_x=' num2str(quadscandata(globalrunnum).betax)]);
    display(['beta_y=' num2str(quadscandata(globalrunnum).betay)]);
    display(['emit_xn=' num2str(quadscandata(globalrunnum).emitx)]);
    display(['emit_yn=' num2str(quadscandata(globalrunnum).emity)]);
    
    save(savefilename, 'quadscandata');
    display(['Data saved to ', savefilename, ':quadscandata(', num2str(globalrunnum), ')']);
    
else
    display('Data not Saved to quadscandata struct');
end

delete runinitial.mat

fprintf('Global Run Number: %d, Local Run Number: %d, Quad: %s range: [%0.2f %0.2f]\n',globalrunnum,qrun,QD,lowb,highb);
if yfitqlty==0
    fprintf('| %d | %d | %s | %s | %0.2f to %0.2f | X | %0.2f |\n',globalrunnum,qrun,QD,SCR,lowb,highb,quadscandata(globalrunnum).emitx*1e6);
elseif xfitqlty==0
    fprintf('| %d | %d | %s | %s | %0.2f to %0.2f | Y | %0.2f |\n',globalrunnum,qrun,QD,SCR,lowb,highb,quadscandata(globalrunnum).emity*1e6);
else
    fprintf('| %d | %d | %s | %s | %0.2f to %0.2f | X/Y | %0.2f/%0.2f |\n',globalrunnum,qrun,QD,SCR,lowb,highb,quadscandata(globalrunnum).emitx*1e6,quadscandata(globalrunnum).emity*1e6);
end

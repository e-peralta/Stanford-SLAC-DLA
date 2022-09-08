function varargout = DataBrowser(varargin)
% DATABROWSER MATLAB code for DataBrowser.fig
%
% This GUI allows the user to browse through a specified E163Vision data
% set or a series of measured spectra from CameraOnly and to perform any
% fitting routing from FitSpectrum. The fits can be visually inspected for
% accuracy and fit parameteres can be extracted and plotted against
% experiment parameters to look for correlations.
%
% 06/06/2012 - eperalta@slac.stanford.edu
%
% Please don't overwrite original, resave as different version

% Last Modified by GUIDE v2.5 22-Apr-2013 23:26:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
   'gui_Singleton',  gui_Singleton, ...
   'gui_OpeningFcn', @DataBrowser_OpeningFcn, ...
   'gui_OutputFcn',  @DataBrowser_OutputFcn, ...
   'gui_LayoutFcn',  [] , ...
   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
   [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
   gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%% --- Executes just before DataBrowser is made visible.

function DataBrowser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function initializes the GUI parameters to produce a 'ready' look

% Choose default command line output for DataBrowser
handles.output = hObject;

%handles.folder='Y:\LEAP\Data_2013_FebRun\';
handles.folder='Y:\LEAP\UCLA_2012_DecRun\';
savefolder=['V:\ARDB\E163\Data\',datestr(now,'yymmdd'),'\'];
if exist(savefolder,'dir')==0
   mkdir(savefolder);
end

handles.dfltFigPos=get(0,'defaultfigureposition');
%handles.dfltFigPos=[296 342 560 420] on Newton
handles.savefolder=savefolder;
set(handles.Slider,'Value',1);
RUNNO='1602_Copy';
set(handles.RUNNO,'String',RUNNO,'Value',1);
set(handles.PlaySpeed,'Value',1);
handles.PAUSE_T=0;
set(handles.PlayRange,'Value',1);
handles.range=1;
set(handles.HoldPlot,'Value',0);
handles.holdPlot=0;
set(handles.LPFilter,'Value',1);
handles.lpFilter=1;
set(handles.ImgFilter,'Value',0);
handles.imgFilter=0;
set(handles.XrayFilter,'Value',1);
handles.xrayFilter=1;

set(handles.LPStrength,'String','.1');
handles.lpStrength=.1;
set(handles.ShowLPF,'Value',1);
handles.showLPF=1;

screen={load('V:\ARDB\E163\Data\2012\120515\OneShot5.dat')};
%ROI=[1 1024 1 700];
ROI=[1 1024 1 1024];
set(handles.Xmin,'String',num2str(ROI(1)),'Value',1);
set(handles.Xmax,'String',num2str(ROI(2)),'Value',1);
set(handles.Ymin,'String',num2str(ROI(3)),'Value',1);
set(handles.Ymax,'String',num2str(ROI(4)),'Value',1);
handles.ROI=ROI;
set(handles.EventNum,'String',num2str(1,'%g'));

[~, ~, ~, ~, ~, ~, ~, ModeList]=FitSpectrum5c(0,0,0,0,0,1);
set(handles.Fit1ModeList,'String',ModeList,'Value',1);
handles.fit1Mode=20;
set(handles.Fit1ModeList,'Value',handles.fit1Mode);
set(handles.MakeFit1,'Value',1);
handles.makeFit1=1;

set(handles.PlotData,'Value',1);
handles.plotData=1;
set(handles.PlotFit,'Value',1);
handles.plotFit=1;
set(handles.PlotSec,'Value',1);
handles.plotSec=1;
set(handles.PlotMain,'Value',0);
handles.plotMain=0;
set(handles.PlotWidths,'Value',1);
handles.plotWidths=1;

set(handles.PeakSep,'String','272');
handles.peakSep=272;

set(handles.Fit2ModeList,'String',ModeList,'Value',1);
handles.fit2Mode=5;
set(handles.Fit2ModeList,'Value',handles.fit2Mode);
set(handles.MakeFit2,'Value',1);
handles.makeFit2=1;

handles.sMin=1;
handles.sMax=-1;
handles.SYrange=1;
set(handles.AutoStripe,'Value',1);
handles.autoStripe=1;

set(handles.Debug,'Value',0);
handles.debug=0;
set(handles.ShowMean,'Value',1)
handles.showMean=1;
set(handles.ShowStanDev,'Value',1)
handles.showStanDev=1;

handles.normMain=1;
set(handles.NormMain,'Value',1)
handles.centMain=0;
set(handles.CenterMain,'Value',0)

handles.maxAmp=-1;
set(handles.MaxAmp,'String','def');
handles.minAmp=-1;
set(handles.MinAmp,'String','def');
handles.maxSpec=-1;
set(handles.MaxSpec,'String','def');
handles.minSpec=-1;
set(handles.MinSpec,'String',num2str(handles.minSpec,'%g'));

handles.delayMin=-1;
set(handles.DelayMin,'String','def');
handles.delayMax=-1;
set(handles.DelayMax,'String','def');
handles.maxFOM=-1;
set(handles.MaxFOM,'String','def');
handles.minFOM=-1;
set(handles.MinFOM,'String','def');

handles.maxPix=-1;
set(handles.MaxPix,'String','def');
handles.minPix=-1;
set(handles.MinPix,'String','def');

handles.imag=screen;
handles.regen=0;
handles.delay=0;
handles.list=[];
handles.fixList=[];
handles.omitList=[];
handles.reasonList=[];

varList={'Peak1amp','Peak1pos','Peak1wH','Peak1wL','MaxEshift','Extra',...
   'ChiSqdE','Peak2amp','PeakSep','Peak2wH','Peak2wL','SprdAmp','SprdPos',...
   'SprdwH','SprdwL','ChiSqdS','clockT','charge','regenP','delay','laserP'};
varList=varList';
labelList={'Main peak amplitude [a.u.]',...         %1
   'Main peak position [pix]',...
   'Main peak high energy width [pix]',...
   'Main peak low energy width [pix]',...
   'Maximum Energy shift [pix]',...                %5
   'Extra (custom) [a.u.]',...
   'Spectrum Fit Chi^2 [a.u.]',...
   'Transmitted peak amplitude [a.u.]',...
   'Peak separation [pix]',...
   'Transmitted peak high energy width [pix]',...  $10
   'Transmitted peak low energy width [pix]',...
   'Spread amplitude [pix]',...
   'Spread peak position [pix]',...
   'Spread width high side [pix]',...
   'Spread width low side [pix]',...               %15
   'Spread Fit Chi^2 [a.u.]',...
   'Clock Time [a.u.]',...
   'Charge [a.u.]',...
   'Regen Power [a.u.]',...
   'Laser Delay [ps]',...                          %20
   'Laser Power [a.u.]'};
handles.varList=varList;
handles.labelList=labelList;
handles.fomChoice=5;
handles.varChoice=20;
set(handles.FOMChoice,'String',varList,'Value',handles.fomChoice);
set(handles.VarChoice,'String',varList,'Value',handles.varChoice);

set(handles.Record,'Value',0);
handles.record=0;

set(handles.SaveAll,'Value',0);
handles.saveAll=0;
handles.specOnly=0;

set(handles.CorFit,'Value',0);
handles.corFit=0;
set(handles.ErrorBar,'Value',1);
handles.errBar=1;
set(handles.ErrOnly,'Value',0);
handles.errOnly=0;
set(handles.OnEvents,'Value',1);
handles.onEvents=1;
set(handles.OffEvents,'Value',1);
handles.offEvents=1;

handles.useOffRef=0;
set(handles.UseOffRef,'Value',0);

handles.checkROI=1;
set(handles.CheckROI,'Value',1);

handles.MaxEshift=0;
handles.MaxEshift_er=0;
handles.makeAvgOnOff=0;

handles.correct=0;
set(handles.CheckROI,'Value',0);

handles.EmaxSet=0;

UpdateDisp(hObject, eventdata, handles,1);

handles.loadParams=0;
fprintf('Welcome to Edgars DataBrowser,\n LOAD your desired run number to start');
set(handles.Status,'String','Load a run number to start\n');
guidata(hObject, handles);

% UIWAIT makes DataBrowser wait for user response (see UIRESUME)
% uiwait(handles.DataBrowser);

% --- Outputs from this function are returned to the command line.
function varargout = DataBrowser_OutputFcn(~, ~, handles)
varargout{1} = handles.output;

function GetRUNNO_Callback(hObject, ~, handles)
folder=handles.folder;
[filename,folder] = uigetfile('*.dat;','Choose RUN data file',folder);
handles.folder=folder;
handles.filename=filename;

RUNNO=filename(4:end-4);
set(handles.RUNNO,'String',RUNNO,'Value',1);

file=LEAPFile([folder filename]);
numEvents=GetNumEvents(file);

sstring=['Run ',RUNNO,' has ',num2str(numEvents),' events, press LOAD to continue'];
set(handles.Status,'String',sstring);

guidata(hObject,handles);

function RUNNO_Callback(hObject, ~, handles)
RUNNO=get(hObject,'String');
filename=['run',RUNNO,'.dat'];
handles.filename=filename;
folder=handles.folder;

if exist([folder,filename])==0
   %folder='S:\Data';
   %fprintf('\n %s not in Y drive, switched to ard-leap07 D drive \n',['run',RUNNO]);
   %sstring=[folder,filename ' not found'];
   sstring=[filename ' not found, use Look up function?'];
   fprintf('\nFile %s not found in default folder,\n maybe use Look up function?\n',filename) 
else
   file=LEAPFile([folder filename]);
   numEvents=GetNumEvents(file);
   handles.numEvents=numEvents;
   sstring=['Run ',RUNNO,' has ',num2str(numEvents),' events, press LOAD to continue'];
end
set(handles.Status,'String',sstring);
guidata(hObject,handles);

function LoadSpec_Callback(hObject, eventdata, handles)
folder='V:\ARDB\E163\Data\130224\';
[filename,folder] = uigetfile('*.dat;','Choose SPECTRA data file',folder);
set(handles.RUNNO,'String','Spectra Only');

spectra=load([folder,filename]);
numEvents=length(spectra(:,1));

%filename = filename(1:end-4);
iden=ones(length(spectra(1,:)),1);

imag=cell(1,numEvents);

for i=1:numEvents
   regen(i)=0;
   imag{i}=iden*spectra(i,:);
   clockT(i)=i;
   delay(i)=0;
end

handles.imag=imag;
handles.regen=regen;
handles.delay=delay;
handles.clockT=clockT;
handles.filename=filename;
handles.numEvents=numEvents;

handles.autoStripe=0;
set(handles.AutoStripe,'Value',handles.autoStripe);
handles.fit1Mode=4;
set(handles.Fit1ModeList,'Value',handles.fit1Mode);

handles.specOnly=1;
UpdateDisp(hObject, eventdata, handles,1);

set(handles.Slider,'Max',numEvents)
set(handles.Slider,'SliderStep',[1/(numEvents-1) .1]); %[arrow-step bar-step] in % change
set(handles.EventNum,'String',num2str(1,'%g'));

sstring=[filename,' has been loaded'];
set(handles.Status,'String',sstring);

handles.list=[];
handles.fixList=[];
handles.omitList=[];
handles.varChoice=17;
set(handles.VarChoice,'Value',handles.varChoice);

set(handles.EventList,'String','');
set(handles.FixList,'String','');
set(handles.OmitList,'String','');

set(handles.RUNNO,'String',filename);
guidata(hObject,handles);


function Load_Button_Callback(hObject, eventdata, handles)
folder=handles.folder;
filename=handles.filename;
DBmode=handles.debug;

file=LEAPFile([folder filename]);
numEvents=GetNumEvents(file);

params=GetParameters(file);
handles.comment=params.Additional_comments;
fprintf('\nFile: %s\nComments: %s\n',filename,params.Additional_comments)

if DBmode
   numEvents=21; %FOR DEBUGGING
   fprintf('\nDEBUG MODE : NumEvents -> %g , Loading... \n',numEvents)
else
   numEvents=GetNumEvents(file);
   fprintf('\nLoading %g events, This can take several minutes... \n',numEvents)
end

bkgnd=GetBackground(file);
handles.bkgnd=bkgnd;
eventParams=GetEventParams(file);

%[width, height] = GetImageDims(file);
%ROI=[1 width 1 height];

% pre-allocate vectors
regen=zeros(1,numEvents);
delay=zeros(1,numEvents);
imag=cell(1,numEvents);
% presT=zeros(1,numEvents);
clockT=zeros(1,numEvents);
charge=zeros(1,numEvents);
regenP=zeros(1,numEvents);
% ebeam=zeros(1,numEvents);
laserP=zeros(1,numEvents);
% laserX=zeros(1,numEvents);
% laserY=zeros(1,numEvents);
% laserT=zeros(1,numEvents);


badOnes=[];

for i=1:numEvents
   %for i=numEvents-10:numEvents
   try
      event = GetEvent(file,i);
      regen(i)=event.Regen_Status;
      %delay(i)=event.Delay_Stage;
      delay(i)=(event.Delay_Stage+1.5687)*30.1328+1;
      imag{i}=event.img-bkgnd;
      %     presT(i)=event.Present_T;
      if i==1
         clockT0=event.Wall_Clock_Time;
      end
      clockT(i)=event.Wall_Clock_Time-clockT0;
      charge(i)=event.Bunch_Charge;
      regenP(i)=event.Regen_Power;
      %     ebeam(i)=event.Ebeam_Mon;
      laserP(i)=event.Laser_Power;
      %     laserX(i)=event.Laser_Position_x;
      %     laserY(i)=event.Laser_Position_y;;
      %     laserT(i)=event.Laser_Timing;
      %     if mod(i,50)==0
      %         fprintf('\n size of image: %d x %d \n',size(imag{i},1),size(imag{i},2))
      %     end
      if mod(i,100)==0
         fprintf('  %d events read so far \n',i)
      end
   catch
      fprintf('  Could not read event # %d \n',i)
      badOnes=[badOnes i];
      beep;
      %continue
   end
end
regen(badOnes)=[];
delay(badOnes)=[];
imag(badOnes)=[];
% presT(badOnes)=[];
clockT(badOnes)=[];
charge(badOnes)=[];
regenP(badOnes)=[];
% ebeam(badOnes)=[];
laserP(badOnes)=[];
% laserX(badOnes)=[];
% laserY(badOnes)=[];
% laserT(badOnes)=[];

numEvents=length(regen);

beep;
fprintf('\n...Finished loading events \n')
sstring=['Finished loading events, you can now Browse, tweak fits, or calculate FOM'];
set(handles.Status,'String',sstring);

CloseLEAPFile(file);

handles.numEvents=numEvents;
handles.imag=imag;
%handles.ROI=ROI;
handles.regen=regen;
handles.delay=delay;
%handles.presT=presT;
handles.clockT=clockT;
handles.charge=charge;
handles.regenP=regenP;
% handles.ebeam=ebeam;
handles.laserP=laserP;
% handles.laserX=laserX;
% handles.laserY=laserY;
% handles.laserT=laserT;
handles.list=[];
handles.fixList=[];
handles.omitList=[];
handles.reasonList=[];
set(handles.EventList,'String','');
set(handles.FixList,'String','');
set(handles.OmitList,'String','');

% GUIparamFile=[filename(1:end-4),'_GUIparam.mat'];
%
% if exist([folder,GUIparamFile],'file')~=0
%     eval(['load ',[folder,GUIparamFile]]);
%     handles.ROI=ROI;
%     handles.lpFilter=lpfilter;
%     handles.fit1Mode=modeF1;
%     handles.peakSep=peakSep;
%     handles.fit2Mode=modeF2;
%     handles.autoStripe=autoStripe;
%     handles.SYrange=SYrange;
%     handles.maxAmp=maxAmp;
%     handles.list=list;
%     fprintf('\n Found saved GUI parameters \n')
% end

UpdateDisp(hObject, eventdata, handles,1);

set(handles.Slider,'Max',numEvents)
set(handles.Slider,'SliderStep',[1/(numEvents-1) .1]); %[arrow-step bar-step] in % change
set(handles.EventNum,'String',num2str(1,'%g'));

handles.specOnly=0;
handles.varChoice=20;
set(handles.VarChoice,'Value',handles.varChoice);

guidata(hObject,handles);

function CalcFOM_Callback(hObject, eventdata, handles)
handles.loadParams=0;
fprintf('\n Calculating FOM for ALL data points...\n')
imag=handles.imag;
regen=handles.regen;
ROI=handles.ROI;
modeF1=handles.fit1Mode;
modeF2=handles.fit2Mode;
sMin=handles.sMin;
sMax=handles.sMax;
lpFilter=handles.lpFilter;
lpStrength=handles.lpStrength;
xrayFilter=handles.xrayFilter;
autoStripe=handles.autoStripe;
peakSep=handles.peakSep;
normMain=handles.normMain;
varList=handles.varList;

numEvents=length(imag);

for i=1:16
   vName=char(varList{i});
   eval([vName,'=zeros(1,numEvents);']);
   if strcmp(vName,'ChiSqdE')==0 && strcmp(vName,'ChiSqdS')==0
      eval([vName,'_er=zeros(1,numEvents);']);
   end
end

x=(ROI(1):ROI(2));
y=(ROI(3):ROI(4));
spectra=zeros(numEvents,length(x));
SYrange=handles.SYrange;
switch SYrange
   case 1
      sYmin=1;
      sYmax=1024;
   case 2
      sYmin=ROI(3);
      sYmax=ROI(4);
end

ys=sYmin:sYmax;
wS=150;
spread=zeros(numEvents,length(ys));

Ybkgd=zeros(numEvents,3);
Cout1=zeros(numEvents,12);
%Ybkgd1=zeros(size(spectra));
Yfit1=zeros(size(spectra));
Ymain1=zeros(size(spectra));
Ysignal1=zeros(size(spectra));
Cout2=zeros(numEvents,12);
%Ybkgd2=zeros(size(spread));
Yfit2=zeros(size(spread));
%PeakAmp=zeros(1,numEvents);
%PeakPos=zeros(1,numEvents);
SpreadW=zeros(numEvents,2);

F=[1 2 1.2011 1.1346];
%modeF           1       5        10        15         20
modeFtoFitSep=  [0 0 0 0 0 6 0 7 0 8 0 7 7 7 0 0 8 0 8  8];
modeFtoPeak1wH= [3 3 4 4 4 3 3 4 4 4 4 3 3 4 3 4 4 4 4  4];
modeFtoPeak1wL= [3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3  3];
modeFtoPeak2wH= [3 3 4 4 4 5 5 6 6 7 7 6 6 6 6 7 7 9 11 7];
modeFtoPeak2wL= [3 3 3 3 3 5 5 6 6 6 6 5 5 3 4 6 6 6 6  6];
modeFtoPeak2amp=[1 1 1 1 1 4 4 5 5 5 5 4 4 5 5 5 5 5 5  5];
modeFtoFactor=F([2 3 4 2 3 3 3 3 3 3 3 3 3 3 1 2 2 2 2  2]);
%modeFtoFactor=F([2 3 4 2 4 3 3 3 3 3 3 3 3 3 1 2 2 2 2  2]);
[b,a]=butter(8,lpStrength,'low'); %%%change the second parameter to change the degree of LP-FIltering. ([0,1], 0 for filtering everything, 1 for no filtering at all)
c=[14.69 .1 .7324 .2796];
fShift0=c(1)*exp(-(lpStrength-c(2))^c(3)/c(4));

modeF1o=modeF1;
%Ntest=77;

for i=1:numEvents
   %for i=Ntest:Ntest
   %for i=Ntest-1:Ntest+1
   
   % as acquired, the screen image is actually the transpose
   screen=imag{i}(y(1):y(end),x(1):x(end));
   spectra0=mean(screen);
   
   %re-scale the image for a nice plot
   spectraF=filter(b,a,spectra0);
   % This initial fit is used for various normalizations below
   [cout, ~, ~, ybkgd, ~, ~, ~,~,cout_ci]=FitSpectrum5c(spectraF,modeF1-1,peakSep);
   peakAmp0=cout(1);
   %peakPos0=cout(2);
   ybkgd0=ybkgd(1);
   
   Peak1amp(i)=peakAmp0;
   Peak1amp_er(i)=cout_ci(1);
   %PeakPos(i)=peakPos0;
   Ybkgd(i,1)=ybkgd0;
   
   screenD=imag{i};
   screenD=screenD-ybkgd0*ones(size(screenD));
   
   if xrayFilter==1
      screenD(screenD>2*peakAmp0)=0;
   end
   screen=screenD(y(1):y(end),x(1):x(end));
   
   spectra0=mean(screen);
   spectraF=filter(b,a,spectra0);
   spectra(i,:)=spectra0;
   
   if normMain==1
      spectra0=spectra0/peakAmp0;
      spectraF=spectraF/peakAmp0;
   end
   
   if lpFilter==1
      Spectra=spectraF;
      xvar=x-fShift0*ones(1,length(x));
   else
      Spectra=spectra0;
      xvar=x;
   end
   
   if handles.useOffRef==1
      if regen(i)==0
         %modeF1=modeF1o;
         modeF1=20;
      else
         %modeF1=modeF1o+1;
         modeF1=19;
      end
   end
   
   %used w=warning('query','last') to find identifier
   warning('off','MATLAB:colon:nonIntegerIndex');
   %[cout , ~, ChiSq , ~, ~, ~, FOMf]=FitSpectrum5b(spectra(i,:),modeF1-1,peakSep);
   [cout,yfit, ChiSq, ybkgd, ymain, ysignal, ~,~,cout_ci]=FitSpectrum5c(Spectra,modeF1-1,peakSep);
   if modeFtoFitSep(modeF1)==0
      fitSep=peakSep;
      fitSep_er=0;
   else
      fitSep=cout(modeFtoFitSep(modeF1));
      fitSep_er=cout_ci(modeFtoFitSep(modeF1));
      %          if regen(i)==0 && handles.useOffRef==1
      %              fprintf('gets in \n')
      %              peakSep=fitSep;
      %          end
   end
   
   Ybkgd(i,2)=ybkgd(1);
   Cout1(i,:)=cout';
   Yfit1(i,:)=yfit;
   if modeF1>5
      Ymain1(i,:)=ymain;
      Ysignal1(i,:)=ysignal;
   end
   
   %     Peak1wH(i)=cout(modeFtoPeak1wH(modeF1))/modeFtoFactor(modeF1);
   %     Peak1wH_er(i)=cout_ci(modeFtoPeak1wH(modeF1))/modeFtoFactor(modeF1);
   %     Peak1wL(i)=cout(modeFtoPeak1wL(modeF1))/modeFtoFactor(modeF1);
   %     Peak1wL_er(i)=cout_ci(modeFtoPeak1wL(modeF1))/modeFtoFactor(modeF1);
   Peak1wH(i)=cout(modeFtoPeak1wH(modeF1))/F(4);
   Peak1wH_er(i)=cout_ci(modeFtoPeak1wH(modeF1))/F(4);
   Peak1wL(i)=cout(modeFtoPeak1wL(modeF1))/F(2);
   Peak1wL_er(i)=cout_ci(modeFtoPeak1wL(modeF1))/F(2);
   Peak1pos(i)=cout(2);
   Peak1pos_er(i)=cout_ci(2);
   Peak2wH(i)=cout(modeFtoPeak2wH(modeF1))/modeFtoFactor(modeF1);
   Peak2wH_er(i)=cout_ci(modeFtoPeak2wH(modeF1))/modeFtoFactor(modeF1);
   Peak2wL(i)=cout(modeFtoPeak2wL(modeF1))/modeFtoFactor(modeF1);
   Peak2wL_er(i)=cout_ci(modeFtoPeak2wL(modeF1))/modeFtoFactor(modeF1);
   Peak2amp(i)=cout(modeFtoPeak2amp(modeF1));
   Peak2amp_er(i)=cout_ci(modeFtoPeak2amp(modeF1));
   PeakSep(i)=fitSep;
   PeakSep_er(i)=fitSep_er;
   
   ChiSqdE(i)=ChiSq;
   %Extra(i)=FOMf;
   if length(cout)>8
      Extra(i)=cout(9);
      Extra_er(i)=cout_ci(9);
   end
   
   if modeF1==19
      %Peak2amp(i)=(cout(5)+cout(12))/2;
      %Peak2amp_er(i)=sqrt(cout_ci(5)^2+cout_ci(12)^2)/4;
      %Peak2wH(i)=cout(11)/F(4)+cout(9)/2;
      %Peak2wH_er(i)=sqrt((cout_ci(11)/F(4))^2+cout_ci(9)^2/4);
      Peak2wL(i)=Peak2wL(i)+cout(9)/2;
      Peak2wL_er(i)=sqrt(Peak2wL_er(i)^2+cout_ci(9)^2/4);
      MaxEshift(i)=cout(8)+cout(9)/2+cout(11)/F(4);
      MaxEshift_er(i)=sqrt(cout_ci(8)^2+cout_ci(9)^2/4+(cout_ci(11)/F(4))^2);
      %MaxEshift(i)=cout(8)+cout(9)/2+cout(11)*1.522;
      %MaxEshift_er(i)=sqrt(cout_ci(8)^2+cout_ci(9)^2/4+(cout_ci(11)*1.522)^2);
      
      % Low fluence mode%%
      Peak2amp(i)=cout(12);
      Peak2amp_er(i)=cout_ci(12);
      Peak2wH(i)=cout(11)/F(4);
      Peak2wH_er(i)=cout_ci(11)/F(4);
      
   elseif modeF1==20
      MaxEshift(i)=cout(8)+cout(7)/F(4);
      MaxEshift_er(i)=sqrt(cout_ci(8)^2+(cout_ci(7)/F(4))^2);
      %MaxEshift(i)=cout(8)+cout(7)*1.522;
      %MaxEshift_er(i)=sqrt(cout_ci(8)^2+(cout_ci(7)*1.522)^2);
   end
   
   if modeF1==18
      Peak2wH(i)=Peak2wH(i)+cout(8)/2;
      Peak2wL(i)=Peak2wL(i)+cout(8)/2;
   end
   
   if autoStripe==1
      sMin=xvar(1)+round(cout(2)+fitSep)-wS/2;
      sMax=xvar(1)+round(cout(2)+fitSep)+wS/2;
   else
      sMin=xvar(1)+round(cout(2)+peakSep)-wS/2;
      sMax=xvar(1)+round(cout(2)+peakSep)+wS/2;
   end
   SpreadW(i,:)=[sMin,sMax];
   %strip=imag{i}(sYmin:sYmax,sMin:sMax);
   strip=screenD(sYmin:sYmax,sMin:sMax);
   spread0=mean(strip,2);
   spreadF=filter(b,a,spread0);
   spread(i,:)=spread0;
   if exist('Ntest','var')==1
      sMin
      spread0(end-5:end)
   end
   
   if lpFilter==1
      Spread=spreadF;
   else
      Spread=spread0;
   end
   
   %[cout, ~, ChiSq , ~, ~, ~, ~]=FitSpectrum5b(spread(i,:),modeF2-1);
   [cout,yfit, ChiSq, ybkgd, ~, ~, ~,~,cout_ci]=FitSpectrum5c(Spread,modeF2-1);
   
   Cout2(i,:)=cout;
   Ybkgd(i,3)=ybkgd(1);
   Yfit2(i,:)=yfit;
   
   SprdPos(i)=cout(2);
   SprdPos_er(i)=cout_ci(2);
   SprdwH(i)=cout(modeFtoPeak1wH(modeF2))/modeFtoFactor(modeF2);
   SprdwH_er(i)=cout_ci(modeFtoPeak1wH(modeF2))/modeFtoFactor(modeF2);
   SprdwL(i)=cout(modeFtoPeak1wL(modeF2))/modeFtoFactor(modeF2);
   SprdwL_er(i)=cout_ci(modeFtoPeak1wL(modeF2))/modeFtoFactor(modeF2);
   SprdAmp(i)=cout(1);
   SprdAmp_er(i)=cout_ci(1);
   ChiSqdS(i)=ChiSq;
end

handles.SpreadW=SpreadW;
handles.spectra=spectra;
handles.spread=spread;

for i=1:16
   vName=char(varList{i});
   eval(['handles.',vName,'=',vName,';']);
   if strcmp(vName,'ChiSqdE')==0 && strcmp(vName,'ChiSqdS')==0
      eval(['handles.',vName,'_er=',vName,'_er;']);
   end
end
handles.Ybkgd=Ybkgd;
handles.Cout1=Cout1;
% Cout1(Ntest-1:Ntest+1,:)'
handles.Yfit1=Yfit1;
handles.Ymain1=Ymain1;
handles.Ysignal1=Ysignal1;
handles.Cout2=Cout2;
handles.Yfit2=Yfit2;
%handles.PeakAmp=PeakAmp;
%handles.PeakPos=PeakPos;
omitList=[];
handles.omitList=[];
handles.reasonList=[];
set(handles.OmitList,'String',num2str(omitList,'%g'));

fprintf('\n ...Finished calculating FOM \n')
beep;
sstring='Finished calculating FOM, you can now plot a correlation';
set(handles.Status,'String',sstring);
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function figTitle=UpdateDisp(hObject, eventdata, handles,eventNum,printFig)
%handles = guidata(hObject);

set(0,'defaultfigureposition',handles.dfltFigPos);
%clear all -handles
if nargin<5
   printFig=0;
   figTitle='NoSaveRequested';
end

PAUSE_T=handles.PAUSE_T;
if eventNum==0
   clear imag regen ROI
   imag{1}=handles.bkgnd;
   regen=0;
   delay=0;
   ROI=[1 1024 1 1024];
   makeFit1=0;
   autoStripe=0;
   makeFit2=0;
else
   imag=handles.imag(eventNum);
   regen=handles.regen(eventNum);
   delay=handles.delay(eventNum);
   ROI=handles.ROI;
   makeFit1=handles.makeFit1;
   autoStripe=handles.autoStripe;
   makeFit2=handles.makeFit2;
end

record=handles.record;
if record==1
   mov = avifile('Sample2.avi','fps',2,'quality',100,...
      'compression','None');
end

modeF1=handles.fit1Mode;
modeF2=handles.fit2Mode;

specOnly=handles.specOnly;
holdPlot=handles.holdPlot;
sMin=handles.sMin;
sMax=handles.sMax;
lpFilter=handles.lpFilter;
showLPF=handles.showLPF;
imgFilter=handles.imgFilter;
xrayFilter=handles.xrayFilter;
lpStrength=handles.lpStrength;
minPix=handles.minPix;
maxPix=handles.maxPix;
minSpec=handles.minSpec;
maxSpec=handles.maxSpec;
checkROI=handles.checkROI;

SYrange=handles.SYrange;

plotData=handles.plotData;
plotSec=handles.plotSec;
plotMain=handles.plotMain;
plotFit=handles.plotFit;
plotWidths=handles.plotWidths;

peakSep=handles.peakSep;

normMain=handles.normMain;
centMain=handles.centMain;

AVGSCRNS=handles.makeAvgOnOff;
CORRECT=handles.correct;

startNum=eventNum(1);
numEvents=length(imag);

x=(ROI(1):ROI(2));
y=(ROI(3):ROI(4));
%spectra=zeros(1,ROI(2)-ROI(1)+1);

switch SYrange
   case 1
      sYmin=1;
      sYmax=size(imag{1},1);%1024;
   case 2
      sYmin=ROI(3);
      sYmax=ROI(4);
end

ys=sYmin:sYmax;
%spread=zeros(numEvents,length(ys));

if minSpec== -1
   minSpec=min(x);
   set(handles.MinSpec,'String',num2str(minSpec,'%g'));
end

if maxSpec== -1
   maxSpec=max(x);
   set(handles.MaxSpec,'String',num2str(maxSpec,'%g'));
end

if minPix== -1
   minPix=250;
   %minPix=ROI(3);
   set(handles.MinPix,'String',num2str(minPix,'%g'));
end
if maxPix== -1
   maxPix=1024;
   %maxPix=ROI(4);
   set(handles.MaxPix,'String',num2str(maxPix,'%g'));
end

% xD=minSpec:maxSpec;
% yD=minPix:maxPix;
xD=1:size(imag{1},2);%1024;
yD=1:size(imag{1},1);%xD;

F=[1 2 1/sqrt(log(2)) 1/asech(1/sqrt(2))];  %[none Lorentz Gauss Sech2]

%modeF           1       5        10        15         20
modeFtoFitSep=  [0 0 0 0 0 6 0 7 0 8 0 7 7 7 0 0 8 0 8  8];
modeFtoPeak1wH= [3 3 4 4 4 3 3 4 4 4 4 3 3 4 3 4 4 4 4  4];
modeFtoPeak1wL= [3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3  3];
modeFtoPeak2wH= [3 3 4 4 4 5 5 6 6 7 7 6 6 6 6 7 7 9 11 7];
modeFtoPeak2wL= [3 3 3 3 3 5 5 6 6 6 6 5 5 3 4 6 6 6 6  6];
modeFtoPeak2amp=[1 1 1 1 1 4 4 5 5 5 5 4 4 5 5 5 5 5 5  5];
modeFtoFactor=F([2 3 4 2 3 3 3 3 3 3 3 3 3 3 1 2 2 2 2 2]);

% screen=zeros(length(y),length(x));
% spectra0=mean(screen);
% spectraF=spectra0;
% screenD=zeros(length(yD),length(xD));
scrnAvgOn=zeros(length(yD),length(xD));
scrnAvgOff=zeros(length(yD),length(xD));
numScrnOn=0;
numScrnOff=0;
wS=150;
% strip=zeros(length(ys),wS);
% spread=mean(strip,2);

modeF1o=modeF1;

% MaxEshift=handles.MaxEshift;
% MaxEshift_er=handles.MaxEshift_er;
% MaxEshift(eventNum)
% MaxEshift_er(eventNum)

%MaxEshift=zeros(numEvents,1);
%MaxEshift_er=zeros(numEvents,1);

[b,a]=butter(8,lpStrength,'low'); %%%change the second parameter to change the degree of LP-FIltering. ([0,1], 0 for filtering everything, 1 for no filtering at all)
c=[14.69 .1 .7324 .2796];
fShift0=c(1)*exp(-(lpStrength-c(2))^c(3)/c(4));

for i=1:numEvents
   maxAmp=handles.maxAmp;
   minAmp=handles.minAmp;
   
   set(handles.EventNum,'String',num2str(eventNum(i),'%g'));
   set(handles.DelayDisp,'String',num2str(delay(i),'%3.2f'));
   if regen(i)==1
      set(handles.LaserStat,'String','ON');
      set(handles.LaserStat,'BackgroundColor',[1 0 0]);
      set(handles.LaserStat,'ForegroundColor',[0 0 0]);
   else
      set(handles.LaserStat,'String','off');
      set(handles.LaserStat,'BackgroundColor',[0 0 1]);
      set(handles.LaserStat,'ForegroundColor',[1 1 1]);
   end
   if eventNum(i)==0
      set(handles.Slider,'Value',1);
   else
      set(handles.Slider,'Value',eventNum(i));
   end
   
   %screen=imag{i}(ROI(3):ROI(4),ROI(1):ROI(2));
   screen=imag{i}(y(1):y(end),x(1):x(end));
   spectra0=mean(screen);
   
   %re-scale the image for a nice plot
   
   spectraF=filter(b,a,spectra0);
   % This initial fit is used for various normalizations below
   [cout, ~, ~, ybkgd, ~, ~, ~,~]=FitSpectrum5c(spectraF,modeF1-1,peakSep);
   peakAmp0=cout(1);
   peakPos0=cout(2);
   if i==1
      posRef=peakPos0;
   end
   
   %screenD=imag{i}(minPix:maxPix,ROI(1):ROI(2));
   %screenD=imag{i}(minPix:maxPix,minSpec:maxSpec);
   %screenD=imag{i}(yD(1):yD(end),xD(1):xD(end));
   screenD=imag{i};
   screenD=screenD-ybkgd(1)*ones(size(screenD));
   %screen=screen-ybkgd(1)*ones(size(screen));
   
   if xrayFilter==1
      %screen(screen>2*peakAmp0)=0;
      screenD(screenD>2*peakAmp0)=0;
   end
   
   if AVGSCRNS==0
      screen=screenD(y(1):y(end),x(1):x(end));
      
      spectra0=mean(screen);
      spectraF=filter(b,a,spectra0);
      
      if normMain==1
         spectra0=spectra0/peakAmp0;
         spectraF=spectraF/peakAmp0;
      end
      if centMain==1
         xVar=x-peakPos0*ones(1,length(x));
         minSpecV=minSpec-peakPos0;
         maxSpecV=maxSpec-peakPos0;
      else
         xVar=x;
         minSpecV=minSpec;
         maxSpecV=maxSpec;
      end
      xvar=xVar;
      
      if lpFilter==1
         xvar=xvar-fShift0*ones(1,length(x));
      end
      if showLPF==1
         spectra=spectraF;
         xVar=xVar-fShift0*ones(1,length(x));
      else
         spectra=spectra0;
      end
      
      if minAmp== -1
         %minAmp=min(spectra);
         minAmp=-.01;
         set(handles.MinAmp,'String','def');
      end
      if maxAmp== -1
         maxAmp=max(spectra);
         set(handles.MaxAmp,'String','def');
      end
      
      if regen(i)==0
         col='b';         col2=[0 0 .7];
         col3=[.5 .5 1];  col4=[0 .6 .4];
         col5='c';
      else
         col='r';         col2=[.7 0 0];
         col3=[1 .5 .5];  col4=[1 .8 0];
         col5='m';
         %col=[0 .9 0]; %For Byer
      end
      
      if printFig==0
         axes(handles.Spectra_fig)
      end
  
      if holdPlot
         hold on;
         lnWd=1;
      else
         cla;
         lnWd=2;
      end
      
      
%  %% TEMPORARY - Remember to remove
%       folder='V:\ARDB\E163\Data\130430\';
%       filename='JoshScatterData.mat';
%       eval(['load ',folder,filename]);
% 
%       filename='JoshScatterData2.mat';
%       eval(['load ',folder,filename]);
% 
%       
%       E=(Energy-60*ones(size(Energy)))*1000/1.2+875*ones(size(Energy));
% %       y1a=Spec400_500;
% %       y1b=Spec400_1000;
% %       y2a=Spec800_500;
% %       y2b=Spec800_1000;
%       y1c=Spec400_500_2;
% 
%       plot(E,y1c./max(y1c)*1.15,'k')
%       hold on
% 
% %%
      if plotData==1 && makeFit1==0
         plot(xVar,spectra,col,'LineWidth',lnWd)
      end
      
      if makeFit1
         %used w=warning('query','last') to find identifier
         warning('off','MATLAB:colon:nonIntegerIndex');
         
         if handles.useOffRef==1
            if regen(i)==0
               %modeF1=17;
               modeF1=20;
               %modeF1=modeF1o;
            else
               %modeF1=16;
               modeF1=19;
               %modeF1=modeF1o+1;
               
            end
         end
         
         if lpFilter==1
            spectra=spectraF;
         else
            spectra=spectra0;
         end
         [cout, yfit, ChiSq, ybkgd, ymain, ysignal, ~,~,cout_ci]=FitSpectrum5c(spectra,modeF1-1,peakSep);
         %         ChiSq=sum(((spectra(50:end)-ybkgd(50:end)-yfit(50:end))./spectra(50:end)).^2);
         %         ChiSq=ChiSq./length(spectra(50:end));
         Cout1(eventNum(i),:)=cout';
         Peak1amp=cout(1);
         Peak1pos=cout(2)+xvar(1);
         %         Peak1wH=cout(modeFtoPeak1wH(modeF1))/modeFtoFactor(modeF1);
         %         Peak1wL=cout(modeFtoPeak1wL(modeF1))/modeFtoFactor(modeF1);
         Peak1wH=cout(modeFtoPeak1wH(modeF1))/F(4);
         Peak1wL=cout(modeFtoPeak1wL(modeF1))/F(2);
         Peak2wH=cout(modeFtoPeak2wH(modeF1))/modeFtoFactor(modeF1);
         Peak2wL=cout(modeFtoPeak2wL(modeF1))/modeFtoFactor(modeF1);
         Peak2amp=cout(modeFtoPeak2amp(modeF1));
         if modeF1==19
            Peak2amp=(cout(5)+cout(12))/2;
            Peak2wH=cout(11)/F(4)+cout(9)/2;
            Peak2wL=Peak2wL+cout(9)/2;
            MaxEshift=cout(8)+cout(9)/2+cout(11)/F(4);
            MaxEshift_er=sqrt(cout_ci(8)^2+cout_ci(9)^2/4+(cout_ci(11)/F(4))^2);
            
            %MaxEshift(eventNum(i))=cout(8)+cout(9)/2+cout(11)*1.522; %'2sigma'
         elseif modeF1==20
            %MaxEshift=cout(8)+cout(7)/F(2);
            %MaxEshift_er=sqrt(cout_ci(8)^2+(cout_ci(7)/F(2))^2);
            %MaxEshift(eventNum(i))=cout(8)+cout(7)*3.474;  %'2sigma'
            %MaxEshift(eventNum(i))=cout(8)+cout(7);  %'2sigma'  Lorentzian
            %MaxEshift(eventNum(i))=cout(8)+cout(7)*1.522;  %sech2
            %MaxEshift=cout(8)+cout(7)*1.522;  %sech2 2sigma
            %MaxEshift_er=sqrt(cout_ci(8)^2+(cout_ci(7)*1.522)^2);
            MaxEshift=cout(8)+cout(7)/F(4);  %sech2 2sigma
            MaxEshift_er=sqrt(cout_ci(8)^2+(cout_ci(7)/F(4))^2);
         end
         
         if modeF1==18
            Peak2wH=Peak2wH+cout(8)/2;
            Peak2wL=Peak2wL+cout(8)/2;
         end
         
         if modeFtoFitSep(modeF1)==0
            fitSep=peakSep;
            set(handles.FitSep,'String','N/A');
         else
            fitSep=cout(modeFtoFitSep(modeF1));
            set(handles.FitSep,'String',num2str(fitSep,'%3.2f'));
            %             if regen(i)==0 && handles.useOffRef==1
            %                   fprintf('gets in \n')
            %                 peakSep=fitSep;
            %             end
         end
         
         if showLPF==1
            spectra=spectraF;
         else
            spectra=spectra0;
         end
         if plotData==1
            plot(xVar,spectra-ybkgd,'Color',col,'LineWidth',lnWd)
         end
         if plotFit
            hold on
            plot(xvar,yfit,'Color',col5,'LineWidth',2,'LineStyle',':')
         end
         
         if plotWidths==1
            %point=[Peak1pos+fitSep ybkgd(round(Peak1pos+fitSep))];
            p=[Peak1pos 0];
            line([p(1)-Peak1wL p(1)],[p(2)+Peak1amp/2 p(2)+Peak1amp/2],...
               'LineWidth',2,'Color',col4);
            hold on
            line([p(1) p(1)+Peak1wH],[p(2)+Peak1amp/2 p(2)+Peak1amp/2],...
               'LineWidth',2,'Color',col4);
            hold on
            line([p(1) p(1)],[p(2)+9*Peak1amp/20 p(2)+11*Peak1amp/20],...
               'LineWidth',2,'Color',col4);
            
            if modeF1>5
               %                 point=[Peak1pos+fitSep 0];
               %                 hold on
               %                 line([point(1)-Peak2wL point(1)],[point(2)+Peak2amp/2 point(2)+Peak2amp/2], 'LineWidth',2,'Color',col4);
               %                 hold on
               %                 line([point(1) point(1)+Peak2wH],[point(2)+Peak2amp/2 point(2)+Peak2amp/2], 'LineWidth',2,'Color',col4);
               %                 hold on
               %                 line([point(1) point(1)],[point(2)+Peak2amp/2-Peak1amp/30 point(2)+Peak2amp/2+Peak1amp/30], 'LineWidth',2,'Color',col4);
               hold on
               %                     line([p(1)+MaxEshift(eventNum(i)) p(1)+ MaxEshift(eventNum(i))],...
               %                         [0 Peak1amp/4], 'LineWidth',2,'Color',col4);
               line([p(1)+MaxEshift p(1)+ MaxEshift], [0 Peak1amp/4], 'LineWidth',2,'Color',col4);
            end
         end
         set(handles.PeakE,'String',num2str(Peak1pos,'%3.2f'));
         set(handles.ChiSqE,'String',num2str(ChiSq,'%5.2e'));
         if specOnly==0
            if plotSec
               hold on
               plot(xvar,ysignal,'Color',col2,'LineWidth',2,'LineStyle','--')
            end
            if plotMain
               hold on
               plot(xvar,ymain,'Color',col3,'LineWidth',2,'LineStyle','-.')
            end
            
            if modeF1>=19
               set(handles.WidthEh,'String',num2str(MaxEshift,'%3.2f'));
            else
               set(handles.WidthEh,'String',num2str(Peak2wH,'%3.2f'));
            end
            set(handles.WidthEl,'String',num2str(Peak2wL,'%3.2f'));
            
            
            if autoStripe==1
               sMin=xvar(1)+round(cout(2)+fitSep)-wS/2;
               sMax=xvar(1)+round(cout(2)+fitSep)+wS/2;
            else
               sMin=xvar(1)+round(cout(2)+peakSep)-wS/2;
               sMax=xvar(1)+round(cout(2)+peakSep)+wS/2;
            end
            set(handles.Smin,'String',num2str(round(sMin),'%g'));
            set(handles.Smax,'String',num2str(round(sMax),'%g'));
            
            if checkROI==1
               hold on
               line([sMin sMin],[minAmp maxAmp], 'LineStyle','-.','Color',[.2 .2 .2]);
               hold on
               line([sMax sMax],[minAmp maxAmp], 'LineStyle','-.','Color',[.2 .2 .2]);
            end
         end
      end
      
      axis tight
      %set(gca,'XGrid','on')
      grid on
      xlim([minSpecV maxSpecV])
      
      if holdPlot ==1 && normMain==1
         maxAmp=1;
      end
      ylim([minAmp maxAmp])
      
      if printFig==1
         continue;
      end
      
      if centMain==1
         sMin=sMin+peakPos0;
         sMax=sMax+peakPos0;
      end
      
      %strip=imag{i}(sYmin:sYmax,sMin:sMax);
      strip=screenD(sYmin:sYmax,sMin:sMax);
      spread0=mean(strip,2);
      spreadF=filter(b,a,spread0);
      
      axes(handles.Spread_fig)
      if holdPlot
         hold on;
      else
         cla;
      end
      
      if plotData==1 && makeFit2==0
         if showLPF==1
            spread=spreadF;
         else
            spread=spread0;
         end
         plot(spread',ys,col,'LineWidth',lnWd)
      end
      
      if makeFit2
         if lpFilter==1
            spread=spreadF;
         else
            spread=spread0;
         end
         maxAmp=max(spread);
         warning('off','MATLAB:colon:nonIntegerIndex');                          %used w=warning('query','last') to find identifier
         %     axes(handles.Screen_fig)
         %     [cout, yfit,ChiSq, ybkgd, ~, ~, ~]=FitSpectrum5b(spread(i,:),modeF2-1,0,0,1);
         %     axes(handles.Spread_fig)
         [cout, yfit,ChiSq, ybkgd, ~, ~, ~]=FitSpectrum5c(spread,modeF2-1);
         SprdPos=cout(2);
         SprdwH=cout(modeFtoPeak2wH(modeF2))/modeFtoFactor(modeF2);
         SprdwL=cout(modeFtoPeak2wL(modeF2))/modeFtoFactor(modeF2);
         Sheight=cout(modeFtoPeak2amp(modeF2));
         
         if plotData==1
            if showLPF==1
               spread=spreadF;
            else
               spread=spread0;
            end
            plot(spread'-ybkgd,ys,'Color',col,'LineWidth',lnWd)
            %plot(spread',ys,'Color',col,'LineWidth',lnWd)
         end
         if plotFit
            hold on
            plot(yfit,ys,'Color',col5,'LineWidth',2,'LineStyle',':');
         end
         if plotSec && plotFit==0
            hold on
            plot(yfit,ys,'Color',col2,'LineWidth',2,'LineStyle','--')
         end
         if plotWidths==1
            p=[SprdPos+1 0];
            hold on
            line([p(2)+Sheight/2 p(2)+Sheight/2],[p(1) p(1)+SprdwH], 'LineWidth',2,'Color',col4);
            hold on
            line([p(2)+Sheight/2 p(2)+Sheight/2],[p(1)-SprdwL p(1)], 'LineWidth',2,'Color',col4);
            hold on
            line([p(2)+9*Sheight/20 p(2)+11*Sheight/20],[p(1) p(1)], 'LineWidth',2,'Color',col4);
         end
         if checkROI==1
            hold on
            line([-1 maxAmp],[y(1) y(1)], 'LineStyle','-.','Color',[.2 .2 .2]);
            hold on
            line([-1 maxAmp],[y(end) y(end)], 'LineStyle','-.','Color',[.2 .2 .2]);
         end
         set(handles.PeakS,'String',num2str(round(cout(2)),'%g'));
         set(handles.WidthS,'String',num2str(SprdwL+SprdwH,'%3.2f'));
         set(handles.ChiSqS,'String',num2str(ChiSq,'%5.2e'));
      end
      %xlim([-2 28])
      xlim([-1 maxAmp])
      %xlim([-1 30])
      ylim([minPix maxPix])
      
      set(gca,'YGrid','on')
      set(gca,'XAxisLocation','Top')
      
      %     %Used this to check amount of xray hits filtered
      %     figure (2)
      %     ind_XR=screenD>2*peakAmp;
      %     sum(sum(ind_XR))
      %     imagesc(screenD.*ind_XR)
      %     axis xy
      %     axis tight
      %     title(['2 x Amp: ',num2str(sum(sum(ind_XR)),'%g'),' pixels above threshold'])
      
   end
   screenD(screenD<0)=0;
   screenD=screenD/peakAmp0;
   if imgFilter==1
      screenD=medfilt2(screenD,[5 5]);
   end
   if regen(i)==1
      scrnAvgOn=scrnAvgOn+circshift(screenD,[0 round(peakPos0-posRef)]);
      numScrnOn=numScrnOn+1;
   end
   if regen(i)==0
      scrnAvgOff=scrnAvgOff+circshift(screenD,[0 round(peakPos0-posRef)]);
      numScrnOff=numScrnOff+1;
   end
   
   if printFig==0 && AVGSCRNS==0
      
      %         [maxV ,ind1]=max(screenD);
      %         [maxV, ind2]=max(maxV);
      %         [minV ,ind1]=min(screenD);
      %         [minV, ind2]=min(minV);
      axes(handles.Screen_fig)
      imagesc(xD,yD,screenD)
      caxis([0 .7])
      %caxis
      %caxis('auto')
      %caxis('manual')
      
      %jet
      %colormap('default')
      %cmap=colormap
      % colormap(jet)
      % caxis([minV maxV])
      
      if checkROI==1
         hold on
         line([x(1) x(end)],[y(1) y(1)], 'LineWidth',1,'Color','w','LineStyle',':');
         hold on
         line([x(1) x(end)],[y(end) y(end)], 'LineWidth',1,'Color','w','LineStyle',':');
         hold on
         line([sMin sMin],[sYmin sYmax], 'LineWidth',1,'Color','w','LineStyle',':');
         hold on
         line([sMax sMax],[sYmin sYmax], 'LineWidth',1,'Color','w','LineStyle',':');
      end
      
      xlim([minSpec maxSpec])
      ylim([minPix maxPix])
      axis xy
      set(gca, 'YTick', []);
      %colorbar('off')
   end
   
   drawnow
   
   if record==1
      M=getframe(gcf,[80 30 380 625]);
      mov = addframe(mov,M);
   end
   pause(PAUSE_T)
   
end
%Cout1(eventNum,:)'
% MaxEshift(eventNum)
% MaxEshift_er(eventNum)

if CORRECT==1 && length(eventNum)==1
   handles.MaxEshift(eventNum)=MaxEshift;
   handles.MaxEshift_er(eventNum)=MaxEshift_er;
end

if AVGSCRNS==1 && numEvents>1
   scrnAvgOn=scrnAvgOn/numScrnOn;
   scrnAvgOff=scrnAvgOff/numScrnOff;
   savefolder=handles.savefolder;
   filename=handles.filename(1:end-4);
   save([savefolder,filename,'_ScrnAvgs'],'eventNum','scrnAvgOff','scrnAvgOn');
   
   figure (97)
   imagesc(xD,yD,scrnAvgOn)
   hold on
   scrnAvgOn=filter(b,a,scrnAvgOn);
   scrnAvgOn=medfilt2(scrnAvgOn,[15,15]);
   contour(xD,yD,scrnAvgOn,[.05 .10 .15 .20 .45],...
      'LineColor','w','LineWidth',2);
   %contour(xD,yD,medfilt2(scrnAvgOn,[10,10]),10)
   %contour(xD,yD,filter(b,a,scrnAvgOn),10)
   %medfilt2(screenD,[5 5]);
   axis xy
   %axis equal
   colormap(hot)
   caxis([.05 .45])
   set(gca,'XGrid','on')
   set(gca,'XColor','w')
   title ('Averaged Spectrum Screen - Laser ON')
   
   figure (98)
   imagesc(xD,yD,scrnAvgOff)
   hold on
   scrnAvgOff=filter(b,a,scrnAvgOff);
   scrnAvgOff=medfilt2(scrnAvgOff,[15,15]);
   contour(xD,yD,scrnAvgOff,[.05 .10 .15 .20 .45],...
      'LineColor','w','LineWidth',2);
   %contour(xD,yD,medfilt2(scrnAvgOff,[10,10]),10)
   %contour(xD,yD,filter(b,a,scrnAvgOn),10)
   axis xy
   %axis equal
   colormap(hot)
   caxis([.05 .45])
   set(gca,'XGrid','on')
   set(gca,'XColor','w')
   
   title ('Averaged Spectrum Screen - Laser off')
end

if record==1
   mov = close(mov);
end

if printFig==1
   sstring='Spectrum plot saved';
   xlabel('Energy [pix]');
   ylabel('Counts [a.u.]');
   if makeFit1==0
      legend('Laser state1','Laser state2');
   end
   enhance_plot;
end
guidata(hObject,handles);

function figTitle=UpdateDisp2(hObject, eventdata, handles,eventNum,printFig)

%handles = guidata(hObject);
set(0,'defaultfigureposition',handles.dfltFigPos);
%clear all -handles
if nargin<5
   printFig=0;
   figTitle='NoSaveRequested';
end

numEvents=length(eventNum);
SpreadW=handles.Data.SpreadW(eventNum,:);
sMin=round(SpreadW(1));
sMax=round(SpreadW(2));

spectra=handles.Data.spectra(eventNum,:);
spread=handles.Data.spread(eventNum,:);

ROI=handles.Data.ROI;
x=(ROI(1):ROI(2));
y=(ROI(3):ROI(4));

imag1=ones(1024,1024);
imag2=zeros(1024,1024);
iden1=ones(1024,1);
imag2(1:1024,x(1):x(end))=iden1*spectra;
iden2=ones(1,1024-sMin+1);
imag1(1:1024,sMin:1024)=spread'*iden2./max(spread);
imag=imag1.*imag2;

regen=handles.Data.regen(eventNum);
delay=handles.Data.delay(eventNum);

peak1amp=handles.Data.Peak1amp(eventNum);
peak1pos=handles.Data.Peak1pos(eventNum);
ybkgd=handles.Data.Ybkgd(eventNum,:);

makeFit1=handles.makeFit1;
autoStripe=handles.autoStripe;
makeFit2=handles.makeFit2;

modeF1=handles.fit1Mode;
modeF2=handles.fit2Mode;

specOnly=handles.specOnly;
holdPlot=handles.holdPlot;
lpFilter=handles.lpFilter;
showLPF=handles.showLPF;
imgFilter=handles.imgFilter;
xrayFilter=handles.xrayFilter;
lpStrength=handles.lpStrength;
minPix=handles.minPix;
maxPix=handles.maxPix;
minSpec=handles.minSpec;
maxSpec=handles.maxSpec;
checkROI=handles.checkROI;

SYrange=handles.SYrange;

plotData=handles.plotData;
plotSec=handles.plotSec;
plotMain=handles.plotMain;
plotFit=handles.plotFit;
plotWidths=handles.plotWidths;

peakSep=handles.peakSep;

normMain=handles.normMain;
centMain=handles.centMain;

CORRECT=handles.correct;

startNum=eventNum(1);

%spectra=zeros(1,ROI(2)-ROI(1)+1);

switch SYrange
   case 1
      sYmin=1;
      sYmax=1024;
   case 2
      sYmin=ROI(3);
      sYmax=ROI(4);
end

ys=sYmin:sYmax;
%spread=zeros(numEvents,length(ys));

if minSpec== -1
   minSpec=min(x);
   set(handles.MinSpec,'String',num2str(minSpec,'%g'));
end

if maxSpec== -1
   maxSpec=max(x);
   set(handles.MaxSpec,'String',num2str(maxSpec,'%g'));
end

if minPix== -1
   minPix=250;
   %minPix=ROI(3);
   set(handles.MinPix,'String',num2str(minPix,'%g'));
end
if maxPix== -1
   maxPix=1024;
   %maxPix=ROI(4);
   set(handles.MaxPix,'String',num2str(maxPix,'%g'));
end

% xD=minSpec:maxSpec;
% yD=minPix:maxPix;
xD=1:1024;
yD=xD;

F=[1 2 1/sqrt(log(2)) 1/asech(1/sqrt(2))];  %[none Lorentz Gauss Sech2]

%modeF           1       5        10        15         20
modeFtoFitSep=  [0 0 0 0 0 6 0 7 0 8 0 7 7 7 0 0 8 0 8  8];
modeFtoPeak1wH= [3 3 4 4 4 3 3 4 4 4 4 3 3 4 3 4 4 4 4  4];
modeFtoPeak1wL= [3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3  3];
modeFtoPeak2wH= [3 3 4 4 4 5 5 6 6 7 7 6 6 6 6 7 7 9 11 7];
modeFtoPeak2wL= [3 3 3 3 3 5 5 6 6 6 6 5 5 3 4 6 6 6 6  6];
modeFtoPeak2amp=[1 1 1 1 1 4 4 5 5 5 5 4 4 5 5 5 5 5 5  5];
modeFtoFactor=F([2 3 4 2 3 3 3 3 3 3 3 3 3 3 1 2 2 2 2 2]);

% screen=zeros(length(y),length(x));
% spectra0=mean(screen);
% spectraF=spectra0;
% screenD=zeros(length(yD),length(xD));

wS=150;
% strip=zeros(length(ys),wS);
% spread=mean(strip,2);

modeF1o=modeF1;

% MaxEshift=handles.MaxEshift;
% MaxEshift_er=handles.MaxEshift_er;
% MaxEshift(eventNum)
% MaxEshift_er(eventNum)

%MaxEshift=zeros(numEvents,1);
%MaxEshift_er=zeros(numEvents,1);

[b,a]=butter(8,lpStrength,'low'); %%%change the second parameter to change the degree of LP-FIltering. ([0,1], 0 for filtering everything, 1 for no filtering at all)
c=[14.69 .1 .7324 .2796];
fShift0=c(1)*exp(-(lpStrength-c(2))^c(3)/c(4));

for i=1:numEvents
   maxAmp=handles.maxAmp;
   minAmp=handles.minAmp;
   
   set(handles.EventNum,'String',num2str(eventNum(i),'%g'));
   set(handles.DelayDisp,'String',num2str(delay(i),'%3.2f'));
   if regen(i)==1
      set(handles.LaserStat,'String','ON');
      set(handles.LaserStat,'BackgroundColor',[1 0 0]);
      set(handles.LaserStat,'ForegroundColor',[0 0 0]);
   else
      set(handles.LaserStat,'String','off');
      set(handles.LaserStat,'BackgroundColor',[0 0 1]);
      set(handles.LaserStat,'ForegroundColor',[1 1 1]);
   end
   if eventNum(i)==0
      set(handles.Slider,'Value',1);
   else
      set(handles.Slider,'Value',eventNum(i));
   end
   
   
   peakAmp0=peak1amp;
   peakPos0=peak1pos;
   if i==1
      posRef=peakPos0;
   end
   
   %screenD=imag{i}(minPix:maxPix,ROI(1):ROI(2));
   %screenD=imag{i}(minPix:maxPix,minSpec:maxSpec);
   %screenD=imag{i}(yD(1):yD(end),xD(1):xD(end));
   screenD=imag;
   screenD=screenD-ybkgd(1)*ones(size(screenD));
   %screen=screen-ybkgd(1)*ones(size(screen));
   
   if xrayFilter==1
      %screen(screen>2*peakAmp0)=0;
      screenD(screenD>2*peakAmp0)=0;
   end
   
   
   screen=screenD(y(1):y(end),x(1):x(end));
   
   spectra0=spectra;
   spectraF=filter(b,a,spectra0);
   
   if normMain==1
      spectra0=spectra0/peakAmp0;
      spectraF=spectraF/peakAmp0;
   end
   if centMain==1
      xVar=x-peakPos0*ones(1,length(x));
      minSpecV=minSpec-peakPos0;
      maxSpecV=maxSpec-peakPos0;
   else
      xVar=x;
      minSpecV=minSpec;
      maxSpecV=maxSpec;
   end
   xvar=xVar;
   
   if lpFilter==1
      xvar=xvar-fShift0*ones(1,length(x));
   end
   if showLPF==1
      spectra=spectraF;
      xVar=xVar-fShift0*ones(1,length(x));
   else
      spectra=spectra0;
   end
   
   if minAmp== -1
      %minAmp=min(spectra);
      minAmp=-.01;
      set(handles.MinAmp,'String','def');
   end
   if maxAmp== -1
      maxAmp=max(spectra);
      set(handles.MaxAmp,'String','def');
   end
   
   if regen(i)==0
      col='b';         col2=[0 0 .7];
      col3=[.5 .5 1];  col4=[0 .6 .4];
      col5='c';
   else
      col='r';         col2=[.7 0 0];
      col3=[1 .5 .5];  col4=[1 .8 0];
      col5='m';
      %col=[0 .9 0]; %For Byer
   end
   
   if printFig==0
      axes(handles.Spectra_fig)
   end
   
   if holdPlot
      hold on;
      lnWd=1;
   else
      cla;
      lnWd=2;
   end
   
   if plotData==1 && makeFit1==0
      plot(xVar,spectra,col,'LineWidth',lnWd)
   end
   
   if makeFit1
      %used w=warning('query','last') to find identifier
      warning('off','MATLAB:colon:nonIntegerIndex');
      
      if handles.useOffRef==1
         if regen(i)==0
            %modeF1=17;
            modeF1=20;
            %modeF1=modeF1o;
         else
            %modeF1=16;
            modeF1=19;
            %modeF1=modeF1o+1;
            
         end
      end
      
      if lpFilter==1
         spectra=spectraF;
      else
         spectra=spectra0;
      end
      [cout, yfit, ChiSq, ybkgd, ymain, ysignal, ~,~,cout_ci]=FitSpectrum5c(spectra,modeF1-1,peakSep);
      %         ChiSq=sum(((spectra(50:end)-ybkgd(50:end)-yfit(50:end))./spectra(50:end)).^2);
      %         ChiSq=ChiSq./length(spectra(50:end));
      Cout1(eventNum(i),:)=cout';
      Peak1amp=cout(1);
      Peak1pos=cout(2)+xvar(1);
      %         Peak1wH=cout(modeFtoPeak1wH(modeF1))/modeFtoFactor(modeF1);
      %         Peak1wL=cout(modeFtoPeak1wL(modeF1))/modeFtoFactor(modeF1);
      Peak1wH=cout(modeFtoPeak1wH(modeF1))/F(4);
      Peak1wL=cout(modeFtoPeak1wL(modeF1))/F(2);
      Peak2wH=cout(modeFtoPeak2wH(modeF1))/modeFtoFactor(modeF1);
      Peak2wL=cout(modeFtoPeak2wL(modeF1))/modeFtoFactor(modeF1);
      Peak2amp=cout(modeFtoPeak2amp(modeF1));
      if modeF1==19
         Peak2amp=(cout(5)+cout(12))/2;
         Peak2wH=cout(11)/F(4)+cout(9)/2;
         Peak2wL=Peak2wL+cout(9)/2;
         MaxEshift=cout(8)+cout(9)/2+cout(11)/F(4);
         MaxEshift_er=sqrt(cout_ci(8)^2+cout_ci(9)^2/4+(cout_ci(11)/F(4))^2);
         
         %MaxEshift(eventNum(i))=cout(8)+cout(9)/2+cout(11)*1.522; %'2sigma'
      elseif modeF1==20
         %MaxEshift(eventNum(i))=cout(8)+cout(7)/F(2);
         %MaxEshift(eventNum(i))=cout(8)+cout(7)*3.474;  %'2sigma'
         %MaxEshift(eventNum(i))=cout(8)+cout(7);  %'2sigma'  Lorentzian
         %MaxEshift(eventNum(i))=cout(8)+cout(7)*1.522;  %sech2
         %MaxEshift=cout(8)+cout(7)*1.522;  %sech2
         %MaxEshift_er=sqrt(cout_ci(8)^2+(cout_ci(7)*1.522)^2);
         MaxEshift=cout(8)+cout(7)/F(4);  %sech2
         MaxEshift_er=sqrt(cout_ci(8)^2+(cout_ci(7)/F(4))^2);
      end
      
      if handles.EmaxSet~=0
         MaxEshift=handles.EmaxSet;
         MaxEshift_er=2;
         
      end
      
      if modeF1==18
         Peak2wH=Peak2wH+cout(8)/2;
         Peak2wL=Peak2wL+cout(8)/2;
      end
      
      if modeFtoFitSep(modeF1)==0
         fitSep=peakSep;
         set(handles.FitSep,'String','N/A');
      else
         fitSep=cout(modeFtoFitSep(modeF1));
         set(handles.FitSep,'String',num2str(fitSep,'%3.2f'));
         %             if regen(i)==0 && handles.useOffRef==1
         %                   fprintf('gets in \n')
         %                 peakSep=fitSep;
         %             end
      end
      
      if showLPF==1
         spectra=spectraF;
      else
         spectra=spectra0;
      end
      if plotData==1
         plot(xVar,spectra-ybkgd,'Color',col,'LineWidth',lnWd)
      end
      if plotFit
         hold on
         plot(xvar,yfit,'Color',col5,'LineWidth',2,'LineStyle',':')
      end
      
      if plotWidths==1
         %point=[Peak1pos+fitSep ybkgd(round(Peak1pos+fitSep))];
         p=[Peak1pos 0];
         line([p(1)-Peak1wL p(1)],[p(2)+Peak1amp/2 p(2)+Peak1amp/2],...
            'LineWidth',2,'Color',col4);
         hold on
         line([p(1) p(1)+Peak1wH],[p(2)+Peak1amp/2 p(2)+Peak1amp/2],...
            'LineWidth',2,'Color',col4);
         hold on
         line([p(1) p(1)],[p(2)+9*Peak1amp/20 p(2)+11*Peak1amp/20],...
            'LineWidth',2,'Color',col4);
         
         if modeF1>5
            %                 point=[Peak1pos+fitSep 0];
            %                 hold on
            %                 line([point(1)-Peak2wL point(1)],[point(2)+Peak2amp/2 point(2)+Peak2amp/2], 'LineWidth',2,'Color',col4);
            %                 hold on
            %                 line([point(1) point(1)+Peak2wH],[point(2)+Peak2amp/2 point(2)+Peak2amp/2], 'LineWidth',2,'Color',col4);
            %                 hold on
            %                 line([point(1) point(1)],[point(2)+Peak2amp/2-Peak1amp/30 point(2)+Peak2amp/2+Peak1amp/30], 'LineWidth',2,'Color',col4);
            hold on
            %                     line([p(1)+MaxEshift(eventNum(i)) p(1)+ MaxEshift(eventNum(i))],...
            %                         [0 Peak1amp/4], 'LineWidth',2,'Color',col4);
            line([p(1)+MaxEshift p(1)+ MaxEshift], [0 Peak1amp/4], 'LineWidth',2,'Color',col4);
         end
      end
      set(handles.PeakE,'String',num2str(Peak1pos,'%3.2f'));
      set(handles.ChiSqE,'String',num2str(ChiSq,'%5.2e'));
      if specOnly==0
         if plotSec
            hold on
            plot(xvar,ysignal,'Color',col2,'LineWidth',2,'LineStyle','--')
         end
         if plotMain
            hold on
            plot(xvar,ymain,'Color',col3,'LineWidth',2,'LineStyle','-.')
         end
         
         if modeF1>=19
            set(handles.WidthEh,'String',num2str(MaxEshift,'%3.2f'));
         else
            set(handles.WidthEh,'String',num2str(Peak2wH,'%3.2f'));
         end
         
         set(handles.WidthEl,'String',num2str(Peak2wL,'%3.2f'));
         
         
         if autoStripe==1
            sMin=xvar(1)+round(cout(2)+fitSep)-wS/2;
            sMax=xvar(1)+round(cout(2)+fitSep)+wS/2;
         else
            sMin=xvar(1)+round(cout(2)+peakSep)-wS/2;
            sMax=xvar(1)+round(cout(2)+peakSep)+wS/2;
         end
         set(handles.Smin,'String',num2str(round(sMin),'%g'));
         set(handles.Smax,'String',num2str(round(sMax),'%g'));
         
         if checkROI==1
            hold on
            line([sMin sMin],[minAmp maxAmp], 'LineStyle','-.','Color',[.2 .2 .2]);
            hold on
            line([sMax sMax],[minAmp maxAmp], 'LineStyle','-.','Color',[.2 .2 .2]);
         end
      end
   end
   
   axis tight
   %set(gca,'XGrid','on')
   grid on
   xlim([minSpecV maxSpecV])
   
   if holdPlot ==1 && normMain==1
      maxAmp=1;
   end
   ylim([minAmp maxAmp])
   
   if printFig==1
      continue;
   end
   
   if centMain==1
      sMin=sMin+peakPos0;
      sMax=sMax+peakPos0;
   end
   
   spread0=spread;
   spreadF=filter(b,a,spread0);
   
   axes(handles.Spread_fig)
   if holdPlot
      hold on;
   else
      cla;
   end
   
   if plotData==1 && makeFit2==0
      if showLPF==1
         spread=spreadF;
      else
         spread=spread0;
      end
      plot(spread',ys,col,'LineWidth',lnWd)
   end
   
   if makeFit2
      if lpFilter==1
         spread=spreadF;
      else
         spread=spread0;
      end
      maxAmp=max(spread);
      warning('off','MATLAB:colon:nonIntegerIndex');                          %used w=warning('query','last') to find identifier
      %     axes(handles.Screen_fig)
      %     [cout, yfit,ChiSq, ybkgd, ~, ~, ~]=FitSpectrum5b(spread(i,:),modeF2-1,0,0,1);
      %     axes(handles.Spread_fig)
      [cout, yfit,ChiSq, ybkgd, ~, ~, ~]=FitSpectrum5c(spread,modeF2-1);
      SprdPos=cout(2);
      SprdwH=cout(modeFtoPeak2wH(modeF2))/modeFtoFactor(modeF2);
      SprdwL=cout(modeFtoPeak2wL(modeF2))/modeFtoFactor(modeF2);
      Sheight=cout(modeFtoPeak2amp(modeF2));
      
      if plotData==1
         if showLPF==1
            spread=spreadF;
         else
            spread=spread0;
         end
         plot(spread-ybkgd,ys,'Color',col,'LineWidth',lnWd)
         %plot(spread',ys,'Color',col,'LineWidth',lnWd)
      end
      if plotFit
         hold on
         plot(yfit,ys,'Color',col5,'LineWidth',2,'LineStyle',':');
      end
      if plotSec && plotFit==0
         hold on
         plot(yfit,ys,'Color',col2,'LineWidth',2,'LineStyle','--')
      end
      if plotWidths==1
         p=[SprdPos+1 0];
         hold on
         line([p(2)+Sheight/2 p(2)+Sheight/2],[p(1) p(1)+SprdwH], 'LineWidth',2,'Color',col4);
         hold on
         line([p(2)+Sheight/2 p(2)+Sheight/2],[p(1)-SprdwL p(1)], 'LineWidth',2,'Color',col4);
         hold on
         line([p(2)+9*Sheight/20 p(2)+11*Sheight/20],[p(1) p(1)], 'LineWidth',2,'Color',col4);
      end
      if checkROI==1
         hold on
         line([-1 maxAmp],[y(1) y(1)], 'LineStyle','-.','Color',[.2 .2 .2]);
         hold on
         line([-1 maxAmp],[y(end) y(end)], 'LineStyle','-.','Color',[.2 .2 .2]);
      end
      set(handles.PeakS,'String',num2str(round(cout(2)),'%g'));
      set(handles.WidthS,'String',num2str(SprdwL+SprdwH,'%3.2f'));
      set(handles.ChiSqS,'String',num2str(ChiSq,'%5.2e'));
   end
   %xlim([-2 28])
   xlim([-1 maxAmp])
   %xlim([-1 30])
   ylim([minPix maxPix])
   
   set(gca,'YGrid','on')
   set(gca,'XAxisLocation','Top')
   
   %     %Used this to check amount of xray hits filtered
   %     figure (2)
   %     ind_XR=screenD>2*peakAmp;
   %     sum(sum(ind_XR))
   %     imagesc(screenD.*ind_XR)
   %     axis xy
   %     axis tight
   %     title(['2 x Amp: ',num2str(sum(sum(ind_XR)),'%g'),' pixels above threshold'])
   
   screenD(screenD<0)=0;
   screenD=screenD/peakAmp0;
   if imgFilter==1
      screenD=medfilt2(screenD,[5 5]);
   end
   
   if printFig==0
      
      %         [maxV ,ind1]=max(screenD);
      %         [maxV, ind2]=max(maxV);
      %         [minV ,ind1]=min(screenD);
      %         [minV, ind2]=min(minV);
      axes(handles.Screen_fig)
      imagesc(xD,yD,screenD)
      caxis([0 .7])
      %caxis
      %caxis('auto')
      %caxis('manual')
      
      %jet
      %colormap('default')
      %cmap=colormap
      % colormap(jet)
      % caxis([minV maxV])
      
      if checkROI==1
         hold on
         line([x(1) x(end)],[y(1) y(1)], 'LineWidth',1,'Color','w','LineStyle',':');
         hold on
         line([x(1) x(end)],[y(end) y(end)], 'LineWidth',1,'Color','w','LineStyle',':');
         hold on
         line([sMin sMin],[sYmin sYmax], 'LineWidth',1,'Color','w','LineStyle',':');
         hold on
         line([sMax sMax],[sYmin sYmax], 'LineWidth',1,'Color','w','LineStyle',':');
      end
      
      xlim([minSpec maxSpec])
      ylim([minPix maxPix])
      axis xy
      set(gca, 'YTick', []);
      %colorbar('off')
   end
   
   drawnow
   
end
%Cout1(eventNum,:)'
% MaxEshift(eventNum)
% MaxEshift_er(eventNum)

if CORRECT==1 && length(eventNum)==1
   handles.Data.MaxEshift(eventNum)=MaxEshift;
   handles.Data.MaxEshift_er(eventNum)=MaxEshift_er;
end

if printFig==1
   sstring='Spectrum plot saved';
   xlabel('Energy [pix]');
   ylabel('Counts [a.u.]');
   if makeFit1==0
      legend('Laser state1','Laser state2');
   end
   enhance_plot;
end
guidata(hObject,handles);


function [output_txt]=cursorfnc(~,event_obj)
set(0,'defaultfigureposition',[-500 -500 560 420])
VERBOSE=0;

if VERBOSE==1, fprintf('--start of cursorfnc\n'); end

if waitforbuttonpress==0
   if VERBOSE==1, fprintf('    inside wait4buttonpress\n'); end
   hs=findobj('type','figure');
   if length(hs)>1
      close(hs(2:end))
   end
   
   pos=get(event_obj,'Position');
   xdata = get(event_obj,'DataIndex');
   handle=get(event_obj,'Target');
   hAxes = get(get(event_obj,'Target'),'Parent');
   eventInfo = getappdata(hAxes,'EventInfo');
   handles=eventInfo.handles;
   
   clear eventNum
   if handles.onEvents==1
      hOn=eventInfo.hOn;
      eventNum_on=eventInfo.On;
      if handle==hOn
         eventNum=eventNum_on(xdata);
      end
   end
   if handles.offEvents==1
      hOff=eventInfo.hOff;
      eventNum_off=eventInfo.Off;
      if handle==hOff
         eventNum=eventNum_off(xdata);
      end
   end
   
   if exist('eventNum','var')==0
      return
   end
   
   if VERBOSE==1, fprintf('    Event %g selected\n',eventNum); end
   % I couldn't call Update_Disp within this cursor function for some reason
   % so from here on it's basically Update_Disp again with minor modifications.
   
   source=handles.source;
   loadParams=handles.loadParams;
   varList=handles.varList;
   
   eval(['SpreadW=',source,'SpreadW(eventNum,:);']);
   sMin=round(SpreadW(1));
   sMax=round(SpreadW(2));
   eval(['Cout1=',source,'Cout1(eventNum,:);']);
   eval(['Cout2=',source,'Cout2(eventNum,:);']);
   
   if loadParams==0
      imag=handles.imag{eventNum};
      yfit1=handles.Yfit1(eventNum,:);
      ymain1=handles.Ymain1(eventNum,:);
      ysignal1=handles.Ysignal1(eventNum,:);
      chiSqd1=handles.ChiSqdE(eventNum);
      yfit2=handles.Yfit2(eventNum,:);
      chiSqd2=handles.ChiSqdS(eventNum);
      
      %peakAmp=handles.PeakAmp(eventNum);
      %peakPos=handles.PeakPos(eventNum);
      if VERBOSE==1, fprintf('    fit data read\n'); end
   else
      spectra=handles.Data.spectra(eventNum,:);
      spread=handles.Data.spread(eventNum,:);
      
      ROI=handles.Data.ROI;
      x=(ROI(1):ROI(2));
      y=(ROI(3):ROI(4));
      
      %imag1=zeros(1024,1024);
      %imag2=zeros(1024,1024);
      %iden1=ones(length(y),1);
      %imag2(y(1):y(end),x(1):x(end))=iden1*spectra;
      %iden2=ones(1,sMax-sMin+1);
      %imag1(1:1024,sMin:sMax)=spread'*iden2;
      %imag=imag1.*imag2;
      %imag=imag*1.5*max(spectra(sMin:sMax))/max(max(imag));
      %imag1(y(1):y(end),sMin:sMax)=zeros(length(y),sMax-sMin+1);
      %imag2(y(1):y(end),sMin:sMax)=zeros(length(y),sMax-sMin+1);
      %imag=imag+imag1+imag2;
      
      imag1=ones(1024,1024);
      imag2=zeros(1024,1024);
      iden1=ones(1024,1);
      imag2(1:1024,x(1):x(end))=iden1*spectra;
      iden2=ones(1,1024-sMin+1);
      %iden2=ones(1,1024-sMin);
      
      
      imag1(1:1024,sMin:1024)=spread'*iden2./max(spread);
      
      imag=imag1.*imag2;
      
      %         set(0,'defaultfigureposition',[296 342 560 420])
      %         figure
      %         imagesc(imag2)
      %         axis xy
      
      %peakAmp=1;
      %peakPos=0;
      if VERBOSE==1, fprintf('    using saved data\n'); end
   end
   
   eval(['regen=',source,'regen(eventNum);']);
   eval(['delay=',source,'delay(eventNum);']);
   eval(['ROI=',source,'ROI;']);
   eval(['lpFilter=',source,'lpFilter;']);
   %eval(['MaxEshift=',source,'MaxEshift(eventNum);']);
   
   eval(['peakSep=',source,'peakSep;']);
   %eval(['PeakSep=',source,'PeakSep(eventNum);']);
   eval(['peak1amp=',source,'Peak1amp(eventNum);']);
   eval(['peak1pos=',source,'Peak1pos(eventNum);']);
   eval(['ybkgd=',source,'Ybkgd(eventNum,:);']);
   
   for i=3:16
      vName=char(varList{i});
      eval([vName,'=',source,vName,'(eventNum);']);
   end
   
   eval(['modeF1=',source,'fit1Mode;']);
   eval(['modeF2=',source,'fit2Mode;']);
   
   
   %peakSep=handles.peakSep;
   makeFit1=handles.makeFit1;
   autoStripe=handles.autoStripe;
   makeFit2=handles.makeFit2;
   
   holdPlot=handles.holdPlot;
   
   imgFilter=handles.imgFilter;
   lpStrength=handles.lpStrength;
   showLPF=handles.showLPF;
   xrayFilter=handles.xrayFilter;
   minPix=handles.minPix;
   maxPix=handles.maxPix;
   maxAmp=handles.maxAmp;
   minAmp=handles.minAmp;
   minSpec=handles.minSpec;
   maxSpec=handles.maxSpec;
   
   SYrange=handles.SYrange;
   
   plotData=handles.plotData;
   plotSec=handles.plotSec;
   plotMain=handles.plotMain;
   plotFit=handles.plotFit;
   plotWidths=handles.plotWidths;
   checkROI=handles.checkROI;
   
   normMain=handles.normMain;
   centMain=handles.centMain;
   
   if VERBOSE==1, fprintf('    parameters loaded\n'); end
   x=(ROI(1):ROI(2));
   y=(ROI(3):ROI(4));
   xD=1:1024;
   yD=xD;
   
   switch SYrange
      case 1
         sYmin=1;
         sYmax=1024;
      case 2
         sYmin=ROI(3);
         sYmax=ROI(4);
   end
   
   ys=sYmin:sYmax;
   
   if minSpec== -1
      minSpec=min(x);
      set(handles.MinSpec,'String',num2str(minSpec,'%g'));
   end
   
   if maxSpec== -1
      maxSpec=max(x);
      set(handles.MaxSpec,'String',num2str(maxSpec,'%g'));
   end
   
   if minPix== -1
      minPix=250;
      set(handles.MinPix,'String',num2str(minPix,'%g'));
   end
   if maxPix== -1
      maxPix=1024;
      set(handles.MaxPix,'String',num2str(maxPix,'%g'));
   end
   
   F=[1 2 1.2013 1.1346];
   %modeF           1       5        10        15         20
   modeFtoFitSep=  [0 0 0 0 0 6 0 7 0 8 0 7 7 7 0 0 8 0 8  8];
   modeFtoPeak1wH= [3 3 4 4 4 3 3 4 4 4 4 3 3 4 3 4 4 4 4  4];
   modeFtoPeak1wL= [3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3  3];
   modeFtoPeak2wH= [3 3 4 4 4 5 5 6 6 7 7 6 6 6 6 7 7 9 10 7];
   modeFtoPeak2wL= [3 3 3 3 3 5 5 6 6 6 6 5 5 3 4 6 6 6 6  6];
   modeFtoPeak2amp=[1 1 1 1 1 4 4 5 5 5 5 4 4 5 5 5 5 5 5  5];
   modeFtoFactor=F([2 3 4 2 3 3 3 3 3 3 3 3 3 3 1 2 2 2 2  2]);
   
   [b,a]=butter(8,lpStrength,'low'); %%%change the second parameter to change the degree of LP-FIltering. ([0,1], 0 for filtering everything, 1 for no filtering at all)
   c=[14.69 .1 .7324 .2796];
   fShift0=c(1)*exp(-(lpStrength-c(2))^c(3)/c(4));
   
   set(handles.EventNum,'String',num2str(eventNum,'%g'));
   set(handles.Slider,'Value',eventNum);
   set(handles.DelayDisp,'String',num2str(delay,'%3.2f'));
   
   if regen==1
      set(handles.LaserStat,'String','ON');
      set(handles.LaserStat,'BackgroundColor',[1 0 0]);
      set(handles.LaserStat,'ForegroundColor',[0 0 0]);
   else
      set(handles.LaserStat,'String','off');
      set(handles.LaserStat,'BackgroundColor',[0 0 1]);
      set(handles.LaserStat,'ForegroundColor',[1 1 1]);
   end
   
   peakAmp0=peak1amp;
   peakPos0=peak1pos;
   
   %if loadParams==0
   
   screenD=imag;
   screenD=screenD-ybkgd(1)*ones(size(screenD));
   
   if xrayFilter==1
      %screen(screen>2*peakAmp0)=0;
      screenD(screenD>2*peakAmp0)=0;
   end
   screen=screenD(y(1):y(end),x(1):x(end));
   if loadParams==0
      spectra0=mean(screen);
   else
      spectra0=spectra;
   end
   
   spectraF=filter(b,a,spectra0);
   
   if normMain==1
      spectra0=spectra0/peakAmp0;
      spectraF=spectraF/peakAmp0;
   end
   
   if centMain==1
      xVar=x-peakPos0*ones(1,length(x));
      minSpecV=minSpec-peakPos0;
      maxSpecV=maxSpec-peakPos0;
   else
      xVar=x;
      minSpecV=minSpec;
      maxSpecV=maxSpec;
   end
   xvar=xVar;
   
   if lpFilter==1
      xvar=xvar-fShift0*ones(1,length(x));
   end
   if showLPF==1
      spectra=spectraF;
      xVar=xVar-fShift0*ones(1,length(x));
   else
      spectra=spectra0;
   end
   
   if minAmp== -1
      %minAmp=min(spectra);
      minAmp=-.01;
      set(handles.MinAmp,'String','def');
   end
   if maxAmp== -1
      maxAmp=max(spectra);
      set(handles.MaxAmp,'String','def');
   end
   
   if regen==0
      col='b';         col2=[0 0 .7];
      col3=[.5 .5 1];  col4=[0 .6 .4];
      col5='c';
   else
      col='r';         col2=[.7 0 0];
      col3=[1 .5 .5];  col4=[1 .8 0];
      col5='m';
      %col=[0 .9 0]; %For Byer
   end
   
   if VERBOSE==1, fprintf('    starting plots\n'); end
   
   axes(handles.Spectra_fig)
   
   if holdPlot
      hold on;
      lnWd=1;
   else
      cla;
      lnWd=2;
   end
   
   if plotData==1 && makeFit1==0
      plot(xVar,spectra,col,'LineWidth',lnWd)
   end
   
   if makeFit1
      %used w=warning('query','last') to find identifier
      warning('off','MATLAB:colon:nonIntegerIndex');
      
      if handles.useOffRef==1
         if regen==0
            %modeF1=17;
            modeF1=20;
            %modeF1=modeF1o;
         else
            %modeF1=16;
            modeF1=19;
            %modeF1=modeF1o+1;
            
         end
      end
      
      ybkgd1=ybkgd(2)*ones(size(spectra));
      %Peak1amp=Cout(1);
      %Peak1pos=Cout(2)+xvar(1);
      if loadParams==0
         cout=Cout1;
         yfit=yfit1;
         ymain=ymain1;
         ysignal=ysignal1;
         Peak1amp=cout(1);
         Peak1pos=cout(2)+xvar(1);
         Peak1wH=cout(modeFtoPeak1wH(modeF1))/F(4);
         Peak1wL=cout(modeFtoPeak1wL(modeF1))/F(2);
         %     Peak1wH=cout(modeFtoPeak1wH(modeF1))/modeFtoFactor(modeF1);
         %     Peak1wL=cout(modeFtoPeak1wL(modeF1))/modeFtoFactor(modeF1);
         Peak2wH=cout(modeFtoPeak2wH(modeF1))/modeFtoFactor(modeF1);
         Peak2wL=cout(modeFtoPeak2wL(modeF1))/modeFtoFactor(modeF1);
         Peak2amp=cout(modeFtoPeak2amp(modeF1));
         if modeF1==19
            Peak2amp=(cout(5)+cout(12))/2;
            Peak2wH=cout(11)/F(4)+cout(9)/2;
            Peak2wL=Peak2wL+cout(9)/2;
            MaxEshift=cout(8)+cout(9)/2+cout(11)/F(4);
            %MaxEshift=cout(8)+cout(9)/2+cout(11)*1.522; %'2sigma'
         elseif modeF1==20
            %MaxEshift(eventNum(i))=cout(8)+cout(7)/F(2);
            %MaxEshift=cout(8)+cout(7)*3.474;  %'2sigma'
            MaxEshift=cout(8)+cout(7)/F(4);
         end
         
         if modeF1==18
            Peak2wH=Peak2wH+cout(8)/2;
            Peak2wL=Peak2wL+cout(8)/2;
         end
         
         if modeFtoFitSep(modeF1)==0
            fitSep=peakSep;
            set(handles.FitSep,'String','N/A');
         else
            fitSep=cout(modeFtoFitSep(modeF1));
            %         if modeF1==18
            %            fitSep=fitSep+cout(9)/2;
            %         end
            set(handles.FitSep,'String',num2str(fitSep,'%3.2f'));
            if regen==0 && handles.useOffRef==1
               peakSep=fitSep;
            end
         end
      else
         if modeF1==19
            [yfit, ymain, ysignal,~,ysig2]=DataFitPlotter(Cout1,modeF1-1);
         else
            [yfit, ymain, ysignal]=DataFitPlotter(Cout1,modeF1-1);
         end
         
      end
      
      if plotData==1
         plot(xVar,spectra-ybkgd1,'Color',col,'LineWidth',lnWd)
      end
      if VERBOSE==1, fprintf('    spectra plot check !\n'); end
      plotFit=1;
      if plotFit
         hold on
         plot(xvar,yfit,'Color',col5,'LineWidth',2,'LineStyle',':')
      end
      if plotWidths==1
         %point=[peak1pos 0];
         %point=[peak1pos-xVar(1)-1 0];
         %point=[peak1pos-xvar(1)-1 0];
         
         %if loadParams==1
         point=[peak1pos-fShift0 0];
         %end
         peak1amp=max(spectra-ybkgd1);%%
         line([point(1)-Peak1wL point(1)],[point(2)+peak1amp/2 point(2)+peak1amp/2], 'LineWidth',2,'Color',col4);
         hold on
         line([point(1) point(1)+Peak1wH],[point(2)+peak1amp/2 point(2)+peak1amp/2], 'LineWidth',2,'Color',col4);
         hold on
         line([point(1) point(1)],[point(2)+9*peak1amp/20 point(2)+11*peak1amp/20], 'LineWidth',2,'Color',col4);
         
         if modeF1>5
            %             point=[Peak1pos+fitSep 0];
            %             hold on
            %             line([point(1)-Peak2wL point(1)],[point(2)+Peak2amp/2 point(2)+Peak2amp/2], 'LineWidth',2,'Color',col4);
            %             hold on
            %             line([point(1) point(1)+Peak2wH],[point(2)+Peak2amp/2 point(2)+Peak2amp/2], 'LineWidth',2,'Color',col4);
            %             hold on
            %             line([point(1) point(1)],[point(2)+Peak2amp/2-Peak1amp/30 point(2)+Peak2amp/2+Peak1amp/30], 'LineWidth',2,'Color',col4);
            hold on
            line([MaxEshift+point(1) MaxEshift+point(1)],[0 peak1amp/4], 'LineWidth',2,'Color',col4);
         end
      end
      if VERBOSE==1, fprintf('    spectra markers added\n'); end
      if plotSec
         hold on
         plot(xvar,ysignal,'Color',col2,'LineWidth',2,'LineStyle','--')
         if modeF1==19
            hold on
            plot(xvar,ysig2,'Color',col5,'LineWidth',2,'LineStyle',':')
         end
         
      end
      if plotMain
         hold on
         plot(xvar,ymain,'Color',col3,'LineWidth',2,'LineStyle','-.')
      end
      
      set(handles.PeakE,'String',num2str(peak1pos,'%3.2f'));
      set(handles.WidthEh,'String',num2str(Peak2wH,'%3.2f'));
      set(handles.WidthEl,'String',num2str(Peak2wL,'%3.2f'));
      %set(handles.ChiSqE,'String',num2str(chiSqd1,'%5.2e'));
      set(handles.ChiSqE,'String',num2str(ChiSqdE,'%5.2e'));
      
      if autoStripe
         %sMin=xvar(1)+round(cout(2)+ peakSep)-50;
         %sMax=xvar(1)+round(cout(2)+ peakSep)+50;
         set(handles.Smin,'String',num2str(round(sMin),'%g'));
         set(handles.Smax,'String',num2str(round(sMax),'%g'));
      end
      
      if checkROI==1
         hold on
         line([sMin sMin],[minAmp maxAmp], 'LineStyle','-.','Color',[.2 .2 .2]);
         hold on
         line([sMax sMax],[minAmp maxAmp], 'LineStyle','-.','Color',[.2 .2 .2]);
      end
   end
   
   axis tight
   set(gca,'XGrid','on')
   %set(gca,'XAxisLocation','Top')
   %xlim([minSpec maxSpec])
   xlim([minSpecV maxSpecV])
   ylim([minAmp maxAmp])
   
   if loadParams==0
      strip=screenD(sYmin:sYmax,sMin:sMax);
      spread0=mean(strip,2);
   else
      spread0=spread;
   end
   spreadF=filter(b,a,spread0);
   
   axes(handles.Spread_fig)
   if holdPlot
      hold on;
   else
      cla;
   end
   
   if showLPF==1
      spread=spreadF;
   else
      spread=spread0;
   end
   maxAmp=max(spread);
   
   if plotData==1 && makeFit2==0
      plot(spread',ys,col,'LineWidth',lnWd)
   end
   
   if makeFit2
      warning('off','MATLAB:colon:nonIntegerIndex');                          %used w=warning('query','last') to find identifier
      ybkgd2=ybkgd(3)*ones(size(spread));
      cout=Cout2;
      
      if loadParams==0
         yfit=yfit2;
         SprdPos=cout(2);
         SprdwH=cout(modeFtoPeak2wH(modeF2))/modeFtoFactor(modeF2);
         SprdwL=cout(modeFtoPeak2wL(modeF2))/modeFtoFactor(modeF2);
         Sheight=cout(modeFtoPeak2amp(modeF2));
      else
         [yfit]=DataFitPlotter(Cout2,modeF2-1);
      end
      
      if plotData==1
         %plot(spread'-ybkgd2,ys,'Color',col,'LineWidth',lnWd)
         plot(spread-ybkgd2,ys,'Color',col,'LineWidth',lnWd)
      end
      if VERBOSE==1, fprintf('    spread plot check!\n'); end
      
      if plotFit==1 %&& loadParams==0
         hold on
         plot(yfit,ys,'Color',col5,'LineWidth',2,'LineStyle',':');
      elseif plotSec==1 %&& loadParams==0
         hold on
         plot(yfit,ys,'Color',col2,'LineWidth',2,'LineStyle','--')
      end
      if plotWidths==1
         point=[SprdPos+1 0];
         Sheight=max(spread-ybkgd2);
         maxAmp=Sheight;
         hold on
         line([point(2)+Sheight/2 point(2)+Sheight/2],[point(1) point(1)+SprdwH], 'LineWidth',2,'Color',col4);
         hold on
         line([point(2)+Sheight/2 point(2)+Sheight/2],[point(1)-SprdwL point(1)], 'LineWidth',2,'Color',col4);
         hold on
         line([point(2)+9*Sheight/20 point(2)+11*Sheight/20],[point(1) point(1)], 'LineWidth',2,'Color',col4);
      end
      if checkROI==1
         hold on
         line([-1 maxAmp],[y(1) y(1)], 'LineStyle','-.','Color',[.2 .2 .2]);
         hold on
         line([-1 maxAmp],[y(end) y(end)], 'LineStyle','-.','Color',[.2 .2 .2]);
      end
      set(handles.PeakS,'String',num2str(round(cout(2)),'%g'));
      set(handles.WidthS,'String',num2str(SprdwL+SprdwH,'%3.2f'));
      set(handles.ChiSqS,'String',num2str(ChiSqdS,'%5.2e'));
   end
   axis tight
   %xlim([-2 28])
   xlim([-1 maxAmp])
   %xlim([-1 30])
   ylim([minPix maxPix])
   
   set(gca,'YGrid','on')
   set(gca,'XAxisLocation','Top')
   
   %if printFig==0
   axes(handles.Screen_fig)
   screenD(screenD<0)=0;
   screenD=screenD/peakAmp0;
   %         [maxV ,ind1]=max(screenD);
   %         [maxV, ind2]=max(maxV);
   %         [minV ,ind1]=min(screenD);
   %         [minV, ind2]=min(minV);
   if imgFilter==1
      screenD=medfilt2(screenD,[5 5]);
   end
   
   imagesc(xD,yD,screenD)
   caxis([0 .7])
   % colormap(jet)
   % caxis([minV maxV])
   
   if checkROI==1
      hold on
      line([x(1) x(end)],[y(1) y(1)], 'LineWidth',1,'Color','w','LineStyle',':');
      hold on
      line([x(1) x(end)],[y(end) y(end)], 'LineWidth',1,'Color','w','LineStyle',':');
      hold on
      line([sMin sMin],[sYmin sYmax], 'LineWidth',1,'Color','w','LineStyle',':');
      hold on
      line([sMax sMax],[sYmin sYmax], 'LineWidth',1,'Color','w','LineStyle',':');
   end
   
   xlim([minSpec maxSpec])
   ylim([minPix maxPix])
   axis xy
   set(gca, 'YTick', []);
   %colorbar('off')
   
   %end
   %drawnow
   % Couldn't get the cursor window to display properly but I don't need to
   output_txt={['event: ', num2str(eventNum,'%g')]};
   
end
if VERBOSE==1, fprintf('-- end of cursorfnc\n'); end
guidata(hObject,handles);

function PlotCor_Callback(hObject, eventdata, handles)
handles=guidata(hObject);
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function figTitle=UpdateCor(hObject, eventdata, handles,printFig)
%clear xVar yVar xLav yLab v2 v1
set(0,'defaultfigureposition',handles.dfltFigPos)
SIG_OPTIONS=optimset('Display','off','TolFun',1e-05,'TolX',1e-05);

if nargin<4
   printFig=0;
   figTitle='NoSaveRequested';
end

fomChoice=handles.fomChoice;
varChoice=handles.varChoice;

varList=handles.varList;
labelList=handles.labelList;
errBar=handles.errBar;
set(handles.ErrorBar,'Value',errBar);
errOnly=handles.errOnly;
corFit=handles.corFit;
set(handles.CorFit,'Value',corFit);

onEvents=handles.onEvents;
offEvents=handles.offEvents;

yLab=char(labelList{fomChoice});
v2=char(varList{fomChoice});
xLab=char(labelList{varChoice});
v1=char(varList{varChoice});

xCorVar={'Peak1wH','Peak1wL','MaxEshift','Peak2wH','Peak2wL','SprdwH','SprdwL'};
setCor=0;
for i=1:length(xCorVar)
   if strcmp(v2,xCorVar(i))==1
      setCor=setCor+1;
   end
   if setCor>0 && strcmp(v1,'delay')==1
      %corFit=1;
      %set(handles.CorFit,'Value',1);
   end
end
if strcmp(v2,'ChiSqdE')==1 || strcmp(v2,'ChiSqdS')==1
   errBar=0;
   set(handles.ErrorBar,'Value',0);
end

loadParams=handles.loadParams;
if loadParams==1
   source='handles.Data.';
else
   source='handles.';
end
handles.source=source;

if (fomChoice < 17 || varChoice < 17) && (isfield(handles,'Peak1pos')==0&&loadParams==0)
   fprintf('\n Please Calculate FOMs first! \n')
   sstring='Please Calculate FOMs first!';
   beep;
else
   
   eval(['regen=',source,'regen;']);
   eval(['xVar=',source,v1,';']);
   
   eventNum=1:length(regen);
   omitList=handles.omitList;
   regen(omitList)=[];
   eventNum(omitList)=[];
   xVar(omitList)=[];
   
   if errOnly==0
      eval(['yVar=',source,v2,';']);
   else
      eval(['yVar=',source,v2,'_er;']);
      yLab=[yLab,' Error'];
      v2=[v2,'Error'];
   end
   yVar(omitList)=[];
   if errBar==1
      eval(['yVar_er=',source,v2,'_er;']);
      yVar_er(omitList)=[];
   end
   v1FiltCond=varChoice<17 && strcmp(v1,'ChiSqdE')==0 && strcmp(v1,'ChiSqdS')==0;
   if v1FiltCond
      eval(['xVar_er=',source,v1,'_er;']);
      xVar_er(omitList)=[];
   end
   
%    % REMOVE THESE LINES
%     xVar=xVar./2;
%     yVar=yVar*1.2-360;
   
   ind_on=(regen==1);
   ind_off=(regen==0);
   
   event_on=eventNum(ind_on);
   event_off=eventNum(ind_off);
   xVar_on=xVar(ind_on);
   xVar_off=xVar(ind_off);
   yVar_on=yVar(ind_on);
   yVar_off=yVar(ind_off);
   if errBar==1
      yVar_on_er=yVar_er(ind_on);
      yVar_off_er=yVar_er(ind_off);
   end
   if v1FiltCond==1
      xVar_on_er=xVar_er(ind_on);
      xVar_off_er=xVar_er(ind_off);
   end
   
   [xVar_on, ind_on]=sort(xVar_on);
   yVar_on=yVar_on(ind_on);
   event_on=event_on(ind_on);
   [xVar_off, ind_off]=sort(xVar_off);
   yVar_off=yVar_off(ind_off);
   event_off=event_off(ind_off);
   
   if errBar==1
      yVar_on_er=yVar_on_er(ind_on);
      yVar_off_er=yVar_off_er(ind_off);
   end
   if v1FiltCond==1
      xVar_on_er=xVar_er(ind_on);
      xVar_off_er=xVar_er(ind_off);
   end
   
   Xmax=handles.delayMax;
   Xmin=handles.delayMin;
   Ymax=handles.maxFOM;
   Ymin=handles.minFOM;
   
   if Xmin== -1
      Xmin=min(xVar);
      %set(handles.DelayMin,'String','def');
      set(handles.DelayMin,'String',num2str(Xmin,'%4.2f'));
   end
   
   if Xmax== -1
      Xmax=max(xVar);
      %set(handles.DelayMax,'String','def');
      set(handles.DelayMax,'String',num2str(Xmax,'%4.2f'));
   end
   
   if Ymin== -1
      if onEvents-offEvents==1
         Ymin=min(yVar_on);
      elseif onEvents-offEvents==-1
         Ymin=min(yVar_off);
      else
         Ymin=min(yVar);
      end
      %set(handles.MinFOM,'String','def');
      set(handles.MinFOM,'String',num2str(Ymin,'%4.2f'));
   end
   
   if Ymax== -1
      if onEvents-offEvents==1
         Ymax=max(yVar_on);
      elseif onEvents-offEvents==-1
         Ymax=max(yVar_off);
      else
         Ymax=max(yVar);
      end
      %set(handles.MaxFOM,'String','def');
      set(handles.MaxFOM,'String',num2str(Ymax,'%4.2f'));
   end
   
   if ~printFig
      axes(handles.Correlation_fig)
   end
   
   cla
   
   ROIdepAvg=0;    % average only events seen on correlation plot
   
   if ROIdepAvg==0
      myVar_on=mean(yVar_on);
      myVar_off=mean(yVar_off);
      stdyVar_on=std(yVar_on);
      stdyVar_off=std(yVar_off);
      if v1FiltCond==1
         mxVar_on=mean(xVar_on);
         mxVar_off=mean(xVar_off);
      end
   else
      
      InRangeOn=[];
      InRangeOff=[];
      
      for i=1:length(xVar_on)
         if xVar_on(i)>=Xmin
            if xVar_on(i)<=Xmax
               InRangeOn=[InRangeOn yVar_on(i)];
            end
         end
      end
      for i=1:length(xVar_off)
         if xVar_off(i)>=Xmin
            if xVar_off(i)<=Xmax
               InRangeOff=[InRangeOff yVar_off(i)];
            end
         end
      end
      
      myVar_on=mean(InRangeOn);
      myVar_off=mean(InRangeOff);
      stdyVar_on=std(InRangeOn);
      stdyVar_off=std(InRangeOff);
   end
   
   if errBar==1
      try
         yW=yVar_off_er;
         flatW=@(c,x) c(1)./yW;
         [myVar_off,~,~,~,~,~,~]=lsqcurvefitstd(flatW,myVar_off,xVar_off,...
            yVar_off./yW,min(yVar_off),max(yVar_off),SIG_OPTIONS);
         %y_off./yW,myVar_off/10,myVar_off*5,SIG_OPTIONS);
         stdyVar_off=rms(yVar_off-myVar_off*ones(size(xVar_off)));
      catch
         fprintf('y Weighted average not possible for laser off points\n')
      end
      try
         yW=yVar_on_er;
         flatW=@(c,x) c(1)./yW;
         [myVar_on,~,~,~,~,~,~]=lsqcurvefitstd(flatW,myVar_on,xVar_on,...
            yVar_on./yW,min(yVar_on),max(yVar_on),SIG_OPTIONS);
         stdyVar_on=rms(yVar_on-myVar_on*ones(size(xVar_on)));
      catch
         fprintf('y Weighted average not possible for laser ON points\n')
      end
   end
   
   myVar_on
   stdyVar_on
   
   myVar_off
   stdyVar_off
   if v1FiltCond==1
      try
         [y_off, ind_off]=sort(yVar_off);
         x_off=xVar_off(ind_off);
         x_off_er=xVar_off_er(ind_off);
         
         xW=x_off_er;
         flatW=@(c,x) c(1)./xW;
         [mxVar_off,~,~,~,~,~,~]=lsqcurvefitstd(flatW,mxVar_off,y_off,...
            x_off./xW,min(x_off),max(x_off),SIG_OPTIONS);
         stdxVar_off=rms(x_off-mxVar_off*ones(size(x_off)));
      catch
         fprintf('x Weighted average not possible for laser off points\n')
      end
      try
         [y_on, ind_on]=sort(yVar_on);
         x_on=xVar_on(ind_on);
         x_on_er=xVar_on_er(ind_on);
         
         xW=x_on_er;
         flatW=@(c,x) c(1)./xW;
         [mxVar_on,~,~,~,~,~,~]=lsqcurvefitstd(flatW,mxVar_on,y_on,...
            x_on./xW,min(x_on),max(x_on),SIG_OPTIONS);
         stdxVar_on=rms(x_on-mxVar_on*ones(size(x_on)));
      catch
         fprintf('x Weighted average not possible for laser ON points\n')
      end
   end
   
   if printFig ==0
      set(handles.CorMeanOn,'String',num2str(myVar_on,'%3.2f'));
      set(handles.CorMeanOff,'String',num2str(myVar_off,'%3.2f'));
      set(handles.CorStdOn,'String',num2str(stdyVar_on,'%3.2f'));
      set(handles.CorStdOff,'String',num2str(stdyVar_off,'%3.2f'));
   end
   
   if handles.showMean==1
      if corFit==0 && onEvents==1
         line([Xmin Xmax],[myVar_on myVar_on], 'LineWidth',2,'Color','r');
         hold on
      end
      %if errBar==0 && offEvents==1;
      if offEvents==1
         line([Xmin Xmax],[myVar_off myVar_off], 'LineWidth',2,'Color','b');
         hold on
      end
      
      if v1FiltCond==1
         if onEvents==1
            line([mxVar_on mxVar_on],[Ymin Ymax], 'LineWidth',2,'Color','r');
            hold on
         end
         if offEvents==1
            line([mxVar_off mxVar_off],[Ymin Ymax], 'LineWidth',2,'Color','b');
            hold on
         end
      end
   end
   
   if corFit==1
      try
         % define gaussian fit
         gauss = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*c(3)^2))+c(4);
         sech2 =@(c,x) c(1)*sech((x-c(2))./c(3)).^2+c(4);
         func=sech2;
         % initial guess - laser on
         c1 = max(yVar_on);
         c2=(Xmax+Xmin)/2;
         c3=c2/4;
         c4=5;
         guess=[c1 c2 c3 c4];
         
         cmin=[0 min(xVar_on) 1 0];
         cmax=[2*max(yVar_on) max(xVar_on) max(xVar_on) max(yVar_on)];
         
         if errBar==0
            [cout,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(func,guess,xVar_on,yVar_on,...
               cmin,cmax,SIG_OPTIONS);
         else
            pow=0;
            pow2=0;
            %yW=yVar_on_er.^pow;  %standard deviations
            Nel=length(yVar_on);
            nCut=round(Nel/5);
            yW=[yVar_on_er(1:2*nCut).^pow yVar_on_er(2*nCut+1:3*nCut).^pow2 yVar_on_er(3*nCut+1:end).^pow];
            xVar_on(2*nCut)
            xVar_on(3*nCut+1)
            gaussW = @(c,x) (c(1)*exp(-(x-c(2)).^2./(2*c(3)^2))+c(4))./yW;
            sech2W =@(c,x) (c(1)*sech((x-c(2))./c(3)).^2+c(4))./yW;
            %         sech2W2 =@(c,x) (c(1)*sech((x-c(2))./c(3)).^2+c(4))./(yW.^c(5));
            %         guess=[guess .5];
            %         cmin=[cmin 0];
            %         cmax=[cmax 1];
            
            funcW=sech2W;
            [cout,ChiSqX,~,~,~,~,~,cout_std]=lsqcurvefitstd(funcW,guess,xVar_on,yVar_on./yW,cmin,cmax,SIG_OPTIONS);
            
            %         weight=cout(5)
            %         cout=cout(1:4);
         end
          cout(3)*2/1.1346
          cout_std(3)*2/1.1346
         %sum((funcW(cout,xVar_on)-yVar_on./yW).^2) %original ChiSqX
         %ChiSqX=sum((func(cout,xVar_on)-yVar_on).^2);       % ChiSqX i care about
         
         ChiSqX=sum(((func(cout,xVar_on)-yVar_on)./yVar_on_er).^2);       % ChiSqX i care about
         ChiSqX=ChiSqX/(length(xVar_on)-length(cout)-1);
         
         %plot(xVar_on,func(cout,xVar_on),'r','LineWidth',2)
         xxvar=linspace(xVar_on(1),xVar_on(end));
         plot(xxvar,func(cout,xxvar),'r','LineWidth',2)
         
         hold on;
         
         set(handles.CorMeanOn,'String',num2str(cout(1),'%3.2f'));
         set(handles.CorStdOn,'String',num2str(cout_std(1),'%3.2f'));
         set(handles.CorFitPos,'String',num2str(cout(2),'%3.2f'));
         set(handles.CorChiSq,'String',num2str(ChiSqX,'%g'));
      catch
         fprintf('Weighted correlation fit not possible for laser ON points\n')
      end
   else
      set(handles.CorFitPos,'String','N/A');
      set(handles.CorChiSq,'String','N/A');
   end
   
   grid on
   if handles.showStanDev==1
      lineT={':','-.','--'};
      Nsigma=3;
      if corFit==0 && onEvents==1
         for i=1:Nsigma
            line([Xmin Xmax],[myVar_on+i*stdyVar_on myVar_on+i*stdyVar_on],...
               'LineWidth',2,'Color','r','LineStyle',char(lineT(i)));
            hold on
            line([Xmin Xmax],[myVar_on-i*stdyVar_on myVar_on-i*stdyVar_on],...
               'LineWidth',2,'Color','r','LineStyle',char(lineT(i)));
            hold on
         end
      end
      if offEvents==1
         for i=1:Nsigma
            line([Xmin Xmax],[myVar_off+i*stdyVar_off myVar_off+i*stdyVar_off],...
               'LineWidth',2,'Color','b','LineStyle',char(lineT(i)));
            hold on
            line([Xmin Xmax],[myVar_off-i*stdyVar_off myVar_off-i*stdyVar_off],...
               'LineWidth',2,'Color','b','LineStyle',char(lineT(i)));
            hold on
         end
      end
      
      if v1FiltCond==1
         grid off
         if onEvents==1
            for i=1:Nsigma
               line([mxVar_on+i*stdxVar_on mxVar_on+i*stdxVar_on],[Ymin Ymax],...
                  'LineWidth',2,'Color','r','LineStyle',char(lineT(i)));
               hold on
               line([mxVar_on-i*stdxVar_on mxVar_on-i*stdxVar_on],[Ymin Ymax],...
                  'LineWidth',2,'Color','r','LineStyle',char(lineT(i)));
               hold on
            end
         end
         if offEvents==1
            for i=1:Nsigma
               line([mxVar_off+i*stdxVar_off mxVar_off+i*stdxVar_off],[Ymin Ymax],...
                  'LineWidth',2,'Color','b','LineStyle',char(lineT(i)));
               hold on
               line([mxVar_off-i*stdxVar_off mxVar_off-i*stdxVar_off],[Ymin Ymax],...
                  'LineWidth',2,'Color','b','LineStyle',char(lineT(i)));
               hold on
            end
         end
      end
      
   end
   
   eventInfo.handles=handles;
   if onEvents==1;
      if errBar==0
         hOn=plot(xVar_on,yVar_on,'or','LineWidth',2);
      else
         hOn=errorbar(xVar_on,yVar_on,yVar_on_er,'or');
      end
      hold on
      eventInfo.On=event_on;
      eventInfo.hOn=hOn;
   end
   if offEvents==1
      if errBar==0
         hOff=plot(xVar_off,yVar_off,'xb','LineWidth',2);
      else
         hOff=errorbar(xVar_off,yVar_off,yVar_off_er,'xb');
      end
      eventInfo.Off=event_off;
      eventInfo.hOff=hOff;
   end
   
   setappdata(gca,'EventInfo',eventInfo)
   
   xlim([Xmin Xmax])
   ylim([Ymin Ymax])
   
   figTitle=['_',v2,'_v_',v1];
   if printFig
      xlabel(xLab);
      ylabel(yLab);
%       if onEvents-offEvents==0;
%          legend([hOff hOn],{'Laser Off', 'Laser On'})
%       elseif onEvents-offEvents==1;
%          legend(hOn,'Laser On')
%       else
%          legend(hOff,'Laser Off')
%       end
      
      
      %     if regen(1)==0
      %         legend('Laser ON','Laser Off');
      %     else
      %         legend('Laser Off','Laser ON');
      %     end
      if corFit==1
         %         title(['Fit max=',num2str(cout(1),'%3.2f'),' \pm ',num2str(cout_std(1),'%3.2f')...
         %             ' @ ',num2str(cout(2),'%3.2f')]);
         title(['Fit max=',num2str(cout(1),'%3.2f'),' \pm ',num2str(cout_std(1),'%3.2f')...
            ' @ ',num2str(cout(2),'%3.2f'),',   ChiSq=',num2str(ChiSqX,'%3.2f'),'    omit ',num2str(length(omitList),'%g'),' pts.']);
      else
         title(['\Delta = ',num2str(myVar_on-myVar_off,'%3.2f'),' \pm ',num2str(sqrt(stdyVar_on^2+stdyVar_off^2),'%3.2f'),',     ',num2str(length(omitList),'%g'),' pts. omitted']);
      end
      
   else
      dcr= datacursormode();%handles.Correlation_fig);
      set(dcr,'Enable','on','DisplayStyle','window','UpdateFcn',@cursorfnc);
      set(dcr,'SnapToDataVertex','off')
      sstring='Data correlation plot ready, you can now inspect individual events';
      set(handles.Status,'String',sstring);
      a=findobj(gcf);
      % b=findobj(handles.Correlation_fig)
      set(a,'HitTest','off')
      if onEvents==1, set(hOn,'HitTest','on'); end
      if offEvents==1, set(hOff,'HitTest','on'); end
      
      %aT=get(a,'Tag');
      %aH=get(a,'HitTest');
      %aT2=get(a,'Type');
      %ahan=[aT, aH, aT2];
      % hLine=findall(gcf,'type','line');
      % hPlots=findall(gcf,'type','plot')
      
      
      
   end
   % folder=handles.savefolder;
   % filename=handles.filename;
   % filename=[folder,filename(1:end-4),figTitle,'.txt'];
   % save(filename,'xVar','yVar','regen','-ASCII');
   
end
guidata(hObject,handles);

function ApplyFiters_Callback(hObject, eventdata, handles)

SIG_OPTIONS=optimset('Display','off','TolFun',1e-05,'TolX',1e-05);
%
% omitList=[];
% reasonList=[];
omitList=handles.omitList;
reasonList=handles.reasonList;

loadParams=handles.loadParams;
if loadParams==1
   source='handles.Data.';
else
   source='handles.';
end

eval(['regen=',source,'regen;']);

%regen=handles.regen;
eventNum=1:length(regen);

regen(omitList)=[];
eventNum(omitList)=[];

fomChoice=handles.fomChoice;
varChoice=handles.varChoice;

varList=handles.varList;
labelList=handles.labelList;
errBar=handles.errBar;
errOnly=handles.errOnly;

yLab=char(labelList{fomChoice});
v2=char(varList{fomChoice});
xLab=char(labelList{varChoice});
v1=char(varList{varChoice});


filterVars=[2 1 8];

varChoice=20;
xLab=char(labelList{varChoice});
tit2=char(varList{varChoice});

passN=0;

omitList0=omitList;

%eval(['xVar=handles.',char(varList{varChoice}),';']);
eval(['xVar=',source,char(varList{varChoice}),';']);

xVar(omitList)=[];

for i=1:length(filterVars)
   
   fomChoice=filterVars(i);
   vName=char(varList{fomChoice});
   
   %     eval(['yVar=handles.',char(varList{fomChoice}),';']);
   %     eval(['yVar_er=handles.',char(varList{fomChoice}),'_er;']);
   eval(['yVar=',source,char(varList{fomChoice}),';']);
   eval(['yVar_er=',source,char(varList{fomChoice}),'_er;']);
   
   yVar(omitList0)=[];
   yVar_er(omitList0)=[];
   
   myVar=mean(yVar);
   stdyVar=std(yVar);
   
   yW=yVar_er;
   flatW=@(c,x) c(1)./yW;
   try
      [myVar,~,~,~,~,~,~]=lsqcurvefitstd(flatW,myVar,xVar,...
         yVar./yW,min(yVar),max(yVar),SIG_OPTIONS);
      stdyVar=rms(yVar-myVar*ones(size(xVar)));
   catch
      fprintf('could not do weighted average\n');
   end
   
   JitThresh=3;
   
   
   reason{1}=vName;
   
   if strcmp(vName,'Peak1pos')==1
      Nbef=length(omitList);
      indExc=(yVar<myVar-2*JitThresh)|(yVar>myVar+2*JitThresh);
      omitList=[omitList eventNum(indExc)];
      for s=1:sum(indExc)
         reasonList=[reasonList reason];
      end
      [omitList ind]=unique(omitList);
      reasonList=reasonList(ind);
      Naft=length(omitList);
      if Naft>Nbef
         fprintf('   %g pts filtered out due to Peak1pos\n', Naft-Nbef);
      end
      %fprintf('    Peak1pos filter applied \n');
   end
   
   Nsig=2;
   if strcmp(vName,'Peak1amp')==1
      Nbef=length(omitList);
      indExc=(yVar<myVar-Nsig*stdyVar)|(yVar>myVar+Nsig*stdyVar);
      omitList=[omitList eventNum(indExc)];
      for s=1:sum(indExc)
         reasonList=[reasonList reason];
      end
      [omitList ind]=unique(omitList);
      reasonList=reasonList(ind);
      Naft=length(omitList);
      if Naft>Nbef
         fprintf('   %g pts filtered out due to Peak1amp\n', Naft-Nbef);
      end
      %fprintf('    Peak1amp filter applied \n');
   end
   
%    TampThresh=.2;
%    if strcmp(vName,'Peak2amp')==1
%        indExc=(yVar<TampThresh);
%        omitList=[omitList eventNum(indExc)];
%        for s=1:sum(indExc)
%            reasonList=[reasonList reason];
%        end
%        [omitList ind]=unique(omitList);
%        reasonList=reasonList(ind);
%        fprintf('    Peak2amp filter applied \n');
%    end
   
   ind_on=(regen==1);
   ind_off=(regen==0);
   
   event_on=eventNum(ind_on);
   event_off=eventNum(ind_off);
   xVar_on=xVar(ind_on);
   xVar_off=xVar(ind_off);
   yVar_on=yVar(ind_on);
   yVar_off=yVar(ind_off);
   
   yVar_on_er=yVar_er(ind_on);
   yVar_off_er=yVar_er(ind_off);
   
   
   [xVar_on, ind_on]=sort(xVar_on);
   yVar_on=yVar_on(ind_on);
   event_on=event_on(ind_on);
   [xVar_off, ind_off]=sort(xVar_off);
   yVar_off=yVar_off(ind_off);
   event_off=event_off(ind_off);
   
   yVar_on_er=yVar_on_er(ind_on);
   yVar_off_er=yVar_off_er(ind_off);
   
   
   myVar_on=mean(yVar_on);
   myVar_off=mean(yVar_off);
   stdyVar_on=std(yVar_on);
   stdyVar_off=std(yVar_off);
   
   yW=yVar_off_er;
   flatW=@(c,x) c(1)./yW;
   try
      [myVar_off,~,~,~,~,~,~]=lsqcurvefitstd(flatW,myVar_off,xVar_off,...
         yVar_off./yW,min(yVar_off),max(yVar_off),SIG_OPTIONS);
      %y_off./yW,myVar_off/10,myVar_off*5,SIG_OPTIONS);
      stdyVar_off=rms(yVar_off-myVar_off*ones(size(xVar_off)));
   catch
      fprintf('could not do weighted average\n');
   end
%    
%    if strcmp(vName,'Peak1pos')==1
%       Nbef=length(omitList);
%       indExc=(yVar_off<myVar_off-2*JitThresh)|(yVar_off>myVar_off+2*JitThresh);
%       omitList=[omitList eventNum(indExc)];
%       for s=1:sum(indExc)
%          reasonList=[reasonList reason];
%       end
%       [omitList ind]=unique(omitList);
%       reasonList=reasonList(ind);
%       Naft=length(omitList);
%       if Naft>Nbef
%          fprintf('   %g pts filtered out due to Peak1pos\n', Naft-Nbef);
%       end
%       %fprintf('    Peak1pos filter applied \n');
%    end
%    
%     if strcmp(vName,'Peak1pos')==1
%       Nbef=length(omitList);
%       indExc=(yVar_on<myVar_on-2*JitThresh)|(yVar_on>myVar_on+2*JitThresh);
%       omitList=[omitList eventNum(indExc)];
%       for s=1:sum(indExc)
%          reasonList=[reasonList reason];
%       end
%       [omitList ind]=unique(omitList);
%       reasonList=reasonList(ind);
%       Naft=length(omitList);
%       if Naft>Nbef
%          fprintf('   %g pts filtered out due to Peak1pos\n', Naft-Nbef);
%       end
%       %fprintf('    Peak1pos filter applied \n');
%     end
%    
   
   if strcmp(vName,'Peak2amp')==1
      Nbef=length(omitList);
      indExc=(yVar_off<myVar_off-Nsig*stdyVar_off);
      omitList=[omitList event_off(indExc)];
      for s=1:sum(indExc)
         reasonList=[reasonList reason];
      end
      [omitList ind]=unique(omitList);
      reasonList=reasonList(ind);
      Naft=length(omitList);
      if Naft>Nbef
         fprintf('   %g OFF pts filtered out due to Peak2amp\n', Naft-Nbef);
      end
      %fprintf('    Peak2amp-offOnly filter applied \n');
   end
   
   
   
end
handles.omitList=omitList;
handles.reasonList=reasonList;
UpdateCor(hObject, eventdata, handles);

slist=[];
for i=1:length(omitList)
   slist=[slist,', ',num2str(omitList(i))];
end
set(handles.OmitList,'String',slist);
guidata(hObject,handles);

% --- Executes on button press in Play_Button.
function Play_Button_Callback(hObject, eventdata, handles)
numEvents=handles.numEvents;
range=handles.range;
ERR=0;
switch range
   case 1
      startNum=int16(str2double(get(handles.EventNum,'String')));
      endNum=startNum+10;
      if endNum>numEvents
         endNum=numEvents;
      end
      list=(startNum:endNum);
   case 2
      list=1:numEvents;
   case 3
      list=handles.list;
      if isempty(list)
         fprintf('\n Event list is empty, please create one first! \n')
         ERR=1;
         beep;
      end
end
if ERR==0
   UpdateDisp(hObject, eventdata, handles,list);
end
guidata(hObject,handles);

function PlaySpeed_Callback(hObject, ~, handles)
speed=get(hObject,'Value');
switch speed
   case 1
      handles.PAUSE_T=0;
   case 2
      handles.PAUSE_T=.1;
   case 3
      handles.PAUSE_T=.5;
end
guidata(hObject,handles);

function PlayRange_Callback(hObject,~, handles)
handles.range=get(hObject,'Value');
guidata(hObject,handles);

function SYrange_Callback(hObject, eventdata, handles)
handles.SYrange=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
UpdateDisp(hObject, eventdata, handles,eventNum);
guidata(hObject,handles);

function EventNum_Callback(hObject, eventdata, handles)
eventNum=int16(str2double(get(hObject,'String')));
loadParams=handles.loadParams;
if loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
   source='handles.';
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
   source='handles.Data.';
end

if handles.correct==1
   eval(['var1=',source,'MaxEshift(eventNum);']);
   eval(['var_er1=',source,'MaxEshift_er(eventNum);']);  
   handles = guidata(hObject);
   eval(['var2=',source,'MaxEshift(eventNum);']);
   eval(['var_er2=',source,'MaxEshift_er(eventNum);']);
   
   fprintf('MaxEshift corrected: \n')
   fprintf('%4.2f +- %4.2f changed to %4.2f +- %4.2f! \n',...
      var1,var_er1,var2,var_er2);
   UpdateCor(hObject, eventdata, handles);
end
guidata(hObject,handles);

function WidthEh_Callback(hObject, eventdata, handles)
handles.EmaxSet=str2double(get(hObject,'String'));
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
   source='handles.';
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
   source='handles.Data.';
end
%handles.EmaxSet=0;
if handles.correct==1
   eval(['var1=',source,'MaxEshift(eventNum);']);
   eval(['var_er1=',source,'MaxEshift_er(eventNum);']);
   handles = guidata(hObject);
   eval(['var2=',source,'MaxEshift(eventNum);']);
   eval(['var_er2=',source,'MaxEshift_er(eventNum);']);
   
   fprintf('MaxEshift corrected: \n')
   fprintf('%4.2f +- %4.2f changed to %4.2f +- %4.2f! \n',...
      var1,var_er1,var2,var_er2);
   UpdateCor(hObject, eventdata, handles);
end
handles.EmaxSet=0;
guidata(hObject,handles);

function Slider_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderNum=get(hObject,'Value');
sliderNum=round(sliderNum);

eventNum=int16(str2double(get(handles.EventNum,'String')));
if eventNum==0 && sliderNum==2
   sliderNum=1;
   set(handles.Slider,'Value',1);
end

if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
set(handles.EventNum,'String',num2str(sliderNum,'%g'));
guidata(hObject,handles);

function Ymin_Callback(hObject, eventdata, handles)
handles.Ymin=int16(str2double(get(hObject,'String')));
eventNum=int16(str2double(get(handles.EventNum,'String')));
ROI=handles.ROI;
ROI(3)=handles.Ymin;
handles.ROI=ROI;
UpdateDisp(hObject, eventdata, handles,eventNum);
guidata(hObject,handles);

function Ymax_Callback(hObject, eventdata, handles)
handles.Ymax=int16(str2double(get(hObject,'String')));
eventNum=int16(str2double(get(handles.EventNum,'String')));
ROI=handles.ROI;
ROI(4)=handles.Ymax;
handles.ROI=ROI;
UpdateDisp(hObject, eventdata, handles,eventNum);
guidata(hObject,handles);

function Xmin_Callback(hObject, eventdata, handles)
handles.Xmin=int16(str2double(get(hObject,'String')));
eventNum=int16(str2double(get(handles.EventNum,'String')));
ROI=handles.ROI;
ROI(1)=handles.Xmin;
handles.ROI=ROI;
UpdateDisp(hObject, eventdata, handles,eventNum);
guidata(hObject,handles);

function Xmax_Callback(hObject, eventdata, handles)
handles.Xmax=int16(str2double(get(hObject,'String')));
eventNum=int16(str2double(get(handles.EventNum,'String')));
ROI=handles.ROI;
ROI(2)=handles.Xmax;
handles.ROI=ROI;
UpdateDisp(hObject, eventdata, handles,eventNum);
guidata(hObject,handles);

function Smin_Callback(hObject, eventdata, handles)
handles.sMin=str2double(get(hObject,'String'));
eventNum=int16(str2double(get(handles.EventNum,'String')));
UpdateDisp(hObject, eventdata, handles,eventNum);
guidata(hObject,handles);

function Smax_Callback(hObject, eventdata, handles)
handles.sMax=str2double(get(hObject,'String'));
eventNum=int16(str2double(get(handles.EventNum,'String')));
UpdateDisp(hObject, eventdata, handles,eventNum);
guidata(hObject,handles);

function MakeFit1_Callback(hObject, eventdata, handles)
handles.makeFit1=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function Fit1ModeList_Callback(hObject, eventdata, handles)
handles.fit1Mode=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function MakeFit2_Callback(hObject, eventdata, handles)
handles.makeFit2=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function Fit2ModeList_Callback(hObject, eventdata, handles)
handles.fit2Mode=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function FOMChoice_Callback(hObject, eventdata, handles)
prevChoice=handles.fomChoice;
handles.fomChoice=get(hObject,'Value');
if handles.fomChoice~=prevChoice
   handles.minFOM=-1;
   handles.maxFOM=-1;
   set(handles.MaxFOM,'String','def');
   set(handles.MinFOM,'String','def');
end
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function VarChoice_Callback(hObject, eventdata, handles)
prevChoice=handles.varChoice;
handles.varChoice=get(hObject,'Value');
if handles.varChoice~=prevChoice
   handles.delayMin=-1;
   handles.delayMax=-1;
   set(handles.DelayMin,'String','def');
   set(handles.DelayMax,'String','def');
end
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function HoldPlot_Callback(hObject, eventdata, handles)
handles.holdPlot=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function PlotFit_Callback(hObject, eventdata, handles)
handles.plotFit=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function PlotSec_Callback(hObject, eventdata, handles)
handles.plotSec=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function PlotMain_Callback(hObject, eventdata, handles)
handles.plotMain=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function PlotData_Callback(hObject, eventdata, handles)
handles.plotData=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function PlotWidths_Callback(hObject, eventdata, handles)
handles.plotWidths=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function PeakSep_Callback(hObject, eventdata, handles)
handles.peakSep=str2double(get(hObject,'String'));
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);


function MaxAmp_Callback(hObject, eventdata, handles)
handles.maxAmp=str2double(get(hObject,'String'));
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function MinAmp_Callback(hObject, eventdata, handles)
handles.minAmp=str2double(get(hObject,'String'));
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function LPStrength_Callback(hObject, eventdata, handles)
handles.lpStrength=str2double(get(hObject,'String'));
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function LPFilter_Callback(hObject, eventdata, handles)
handles.lpFilter=get(hObject,'Value');
handles.showLPF=get(hObject,'Value');
set(handles.ShowLPF,'Value',handles.showLPF);

eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function ImgFilter_Callback(hObject, eventdata, handles)
handles.imgFilter=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function XrayFilter_Callback(hObject, eventdata, handles)
handles.xrayFilter=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
UpdateDisp(hObject, eventdata, handles,eventNum);
guidata(hObject,handles);

function AutoStripe_Callback(hObject, eventdata, handles)
handles.autoStripe=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function Debug_Callback(hObject, eventdata, handles)
handles.debug=get(hObject,'Value');
guidata(hObject,handles);

function ShowMean_Callback(hObject, eventdata, handles)
handles.showMean=get(hObject,'Value');
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function ShowStanDev_Callback(hObject, eventdata, handles)
handles.showStanDev=get(hObject,'Value');
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function MaxFOM_Callback(hObject, eventdata, handles)
handles.maxFOM=str2double(get(hObject,'String'));
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function MinFOM_Callback(hObject, eventdata, handles)
handles.minFOM=str2double(get(hObject,'String'));
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function DelayMax_Callback(hObject, eventdata, handles)
handles.delayMax=str2double(get(hObject,'String'));
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function DelayMin_Callback(hObject, eventdata, handles)
handles.delayMin=str2double(get(hObject,'String'));
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function MinSpec_Callback(hObject, eventdata, handles)
handles.minSpec=str2double(get(hObject,'String'));
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function MaxSpec_Callback(hObject, eventdata, handles)
handles.maxSpec=str2double(get(hObject,'String'));
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function MinPix_Callback(hObject, eventdata, handles)
handles.minPix=str2double(get(hObject,'String'));
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function MaxPix_Callback(hObject, eventdata, handles)
handles.maxPix=str2double(get(hObject,'String'));
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
else
   UpdateDisp2(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function Record_Callback(hObject, eventdata, handles)
handles.record=get(hObject,'Value');
guidata(hObject,handles);

function SaveAll_Callback(hObject, eventdata, handles)
handles.saveAll=get(hObject,'Value');
guidata(hObject,handles);

function UseOffRef_Callback(hObject, eventdata, handles)
handles.useOffRef=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
UpdateDisp(hObject, eventdata, handles,eventNum);
guidata(hObject,handles);

function CorFit_Callback(hObject, eventdata, handles)
handles.corFit=get(hObject,'Value');
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function ErrorBar_Callback(hObject, eventdata, handles)
handles.errBar=get(hObject,'Value');
if handles.errBar==1
   handles.errOnly=0;
   set(handles.ErrOnly,'Value',0);
end
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function ErrOnly_Callback(hObject, eventdata, handles)
handles.errOnly=get(hObject,'Value');
if handles.errOnly==1
   handles.errBar=0;
   set(handles.ErrorBar,'Value',0);
end
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function NormMain_Callback(hObject, eventdata, handles)
handles.normMain=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
UpdateDisp(hObject, eventdata, handles,eventNum);
guidata(hObject,handles);

function CenterMain_Callback(hObject, eventdata, handles)
handles.centMain=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
UpdateDisp(hObject, eventdata, handles,eventNum);
guidata(hObject,handles);

function ShowLPF_Callback(hObject, eventdata, handles)
handles.showLPF=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function CheckROI_Callback(hObject, eventdata, handles)
handles.checkROI=get(hObject,'Value');
eventNum=int16(str2double(get(handles.EventNum,'String')));
if handles.loadParams==0
   UpdateDisp(hObject, eventdata, handles,eventNum);
end
guidata(hObject,handles);

function MakeAvgOnOff_Callback(hObject, eventdata, handles)
handles.makeAvgOnOff=get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filepath=['V:\ARDB\E163\DATA\',datestr(now,'yymmdd'),'\'];
if exist(filepath)==0
   mkdir('V:\ARDB\E163\DATA',datestr(now,'yymmdd'));
end

%currentEvent=str2num(get(handles.EventNum,'String'))
currentEvent=int16(str2double(get(handles.EventNum,'String')))
%currentEvent=str2num(handles.EventNum)
%RUNNO=handles.RUNNO;
RUNNO=get(handles.RUNNO,'String');
filename=['spectrum_run',RUNNO,'_event',num2str(currentEvent),'.dat'];
imag=handles.imag;
ROI=handles.ROI;
%singleImag=imag{currentEvent}(ROI(3):ROI(4),ROI(1):ROI(2));
singleImag=imag{currentEvent};
display(['Saving output file ',filepath,filename]);
save([filepath,filename],'singleImag','-ascii');
guidata(hObject,handles);

function SaveParams_Callback(hObject, eventdata, handles)
savefolder=handles.savefolder;
savefolder
% if handles.loadParams==1
%    Data=handles.Data;
%    Data
% end
varList=handles.varList;
modeF1=handles.fit1Mode;
modeF2=handles.fit2Mode;
f1=num2str(modeF1-1,'%g');
f2=num2str(modeF2-1,'%g');
filename=handles.filename(1:end-4);
filename=[savefolder,filename,'_f',f1,'-',f2,'_1'];
fileNum=2;
while exist([filename,'.mat'],'file')~=0
   filename=[filename(1:end-1),num2str(fileNum,'%g')];
   fileNum=fileNum+1;
end

loadParams=handles.loadParams;
if loadParams==0
   
   Data.fit1Mode=modeF1;
   Data.fit2Mode=modeF2;
   Data.peakSep=handles.peakSep;
   Data.ROI=handles.ROI;
   Data.SpreadW=handles.SpreadW;
   Data.lpFilter=handles.lpFilter;
   Data.spectra=handles.spectra;
   Data.spread=handles.spread;
   Data.Ybkgd=handles.Ybkgd;
   Data.Cout1=handles.Cout1;
   Data.Cout2=handles.Cout2;
   for i=1:length(varList)
      vName=char(varList{i});
      eval(['Data.',vName,'=handles.',vName,';']);
      if i<17
         if strcmp(vName,'ChiSqdE')==0 & strcmp(vName,'ChiSqdS')==0
            eval(['Data.',vName,'_er=handles.',vName,'_er;']);
         end
      end
   end
   Data.regen=handles.regen;
else
   Data=handles.Data;
end
Data.omitList=handles.omitList;
Data.reasonList=handles.reasonList;


save(filename,'Data');
guidata(hObject,handles);

function LoadParams_Callback(hObject, eventdata, handles)
handles.loadParams=1;
folder=handles.savefolder;
%folder='V:\ARDB\E163\Data\130416\Fluence1\';
%folder='V:\ARDB\E163\Data\130513\Polarization1\';
%folder='V:\ARDB\E163\Data\130423\Fluence1\';
folder='V:\ARDB\E163\Data\130416\Fluence2\';
[filename,folder] = uigetfile('*.mat;','Choose previously saved parameters',folder);
%    filename = uigetfile('*.mat;','Choose Optics Setting file','V:\ARDB\E163\matlab\OpticsModel\');
eval(['load ',folder,filename]);

if Data.clockT(1)~=0
   clockT0=Data.clockT(1);
   Data.clockT=Data.clockT-clockT0*ones(size(Data.clockT));
end

handles.Data=Data;
omitList=Data.omitList;
handles.omitList=omitList;
handles.reasonList=Data.reasonList;
handles.SpreadW=Data.SpreadW;
handles.filename=filename(1:11);
handles.fit1Mode=Data.fit1Mode;
handles.fit2Mode=Data.fit2Mode;
set(handles.Fit1ModeList,'Value',handles.fit1Mode);
set(handles.Fit2ModeList,'Value',handles.fit2Mode);
set(handles.RUNNO,'String',filename(1:end-4));

UpdateCor(hObject, eventdata, handles);

slist=[];
for i=1:length(omitList)
   slist=[slist,', ',num2str(omitList(i))];
end
set(handles.OmitList,'String',slist);
set(handles.EventList,'String',[]);
axes(handles.Screen_fig)
cla
axes(handles.Spectra_fig)
cla
axes(handles.Spread_fig)
cla
ROI=Data.ROI;
set(handles.Xmin,'String',num2str(ROI(1)),'Value',1);
set(handles.Xmax,'String',num2str(ROI(2)),'Value',1);
set(handles.Ymin,'String',num2str(ROI(3)),'Value',1);
set(handles.Ymax,'String',num2str(ROI(4)),'Value',1);
set(handles.PeakSep,'String',num2str(Data.peakSep));
handles.peakSep=handles.Data.peakSep;
guidata(hObject,handles);

% --- Executes on button press in SaveCorImg.
function SaveCorImg_Callback(hObject, eventdata, handles)
savefolder=handles.savefolder
saveAll=handles.saveAll;
specOnly=handles.specOnly;

f1=num2str(handles.fit1Mode-1,'%g');
f2=num2str(handles.fit2Mode-1,'%g');

nSaves=1;
savePlot=1;
if saveAll==1
   nSaves=13;
   savePlot=[1:5 8:9 13:15];
   if specOnly==1
      nSaves=4;
      savePlot=[1:5 8:9 13:15];
   end
end
%for i=1:nSaves
for i=1:length(savePlot)
   if saveAll==1
      %handles.fomChoice=i;
      handles.fomChoice=savePlot(i);
   end
   figure (99)
   figTitle=UpdateCor(hObject, eventdata, handles,1);
   %set(gcf,'PaperPositionMode','auto')
   %filename=handles.filename;
   filename=[handles.filename(1:end-4),'_f',f1,'-',f2,figTitle,'_1'];
   ImgNum=2;
   
   while exist([savefolder,filename,'.png'],'file')~=0
      filename=[filename(1:end-1),num2str(ImgNum,'%g')];
      ImgNum=ImgNum+1;
   end
   
   enhance_plot;
   print(gcf, '-dpng',[savefolder,filename]);
   %close (99)
end
sstring=['Correlation plot "',filename,'" saved'];
set(handles.Status,'String',sstring);
guidata(hObject,handles);

function SaveSpecImg_Callback(hObject, eventdata, handles)
savefolder=handles.savefolder;
list=handles.list;
varChoice=handles.varChoice;
varList=handles.varList;

%eval(['xVar=handles.',char(varList{varChoice}),';']);

loadParams=handles.loadParams;
if loadParams==1
   source='handles.Data.';
else
   source='handles.';
end
eval(['xVar=',source,char(varList{varChoice}),';']);

makeFit1=handles.makeFit1;
f1=num2str(handles.fit1Mode-1,'%g');

filename=handles.filename(1:end-4);
if makeFit1==1
   titleStr=[filename,', Fit mode ',f1];
   filename=[filename,'_f',f1];
end

if isempty(list)
   EventStr=get(handles.EventNum,'String');
   EventNum=int16(str2double(EventStr));
   list=EventNum;
   filename=[filename,'_Event',EventStr];
   titleStr=[titleStr,', Event ',EventStr];
else
   minD=min(xVar(list));
   maxD=max(xVar(list));
   filename=[filename,'_EventsWith',char(varList{varChoice}),'=',num2str(round(minD),'%g'),'-',num2str(round(maxD),'%g'),'ps'];
   titleStr=[titleStr,', Events between ',char(varList{varChoice}),' = ',num2str(minD,'%3.2f'),' and ',num2str(maxD,'%3.2f'),' ps'];
end

filename=[filename,'_1'];
ImgNum=2;
while exist([savefolder,filename,'.png'],'file')~=0
   filename=[filename(1:end-1),num2str(ImgNum,'%g')];
   ImgNum=ImgNum+1;
end

set(0,'defaultfigureposition',[10 10 560 420]);
figure (99)
if loadParams==0
   UpdateDisp(hObject, eventdata, handles,list,1);
else
   UpdateDisp2(hObject, eventdata, handles,list,1);
end
%titleStr=[handles.filename(1:end-4),titleStr]
title(titleStr);
%enhance_plot;
enhance_plot(0,0,-1);
print(gcf, '-dpng',[savefolder,filename]);
%close (99)
guidata(hObject,handles);

function OnEvents_Callback(hObject, eventdata, handles)
handles.onEvents=get(hObject,'Value');
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function OffEvents_Callback(hObject, eventdata, handles)
handles.offEvents=get(hObject,'Value');
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function Correct_Callback(hObject, eventdata, handles)
handles.correct=get(hObject,'Value');
guidata(hObject,handles);

function AddToList_Callback(hObject, eventdata, handles)
eventNum=int16(str2double(get(handles.EventNum,'String')));

list=handles.list;
if sum(list==eventNum)==0
   list=[list eventNum];
   list=sort(list);
   handles.list=list;
else
   fprintf('\n Event %g already in list \n', eventNum)
   beep;
end
slist=[];
for i=1:length(list)
   if length(list)==1
      slist=num2str(list(i));
   else
      slist=[slist,', ',num2str(list(i))];
   end
end

set(handles.EventList,'String',slist);
guidata(hObject,handles);

function RemFromList_Callback(hObject, eventdata, handles)
eventNum=int16(str2double(get(handles.EventNum,'String')));

list=handles.list;
ind=find(list==eventNum);
list(ind)=[];
handles.list=list;

fprintf('\n Event %g removed from list \n', eventNum)

slist=[];
for i=1:length(list)
   slist=[slist,', ',num2str(list(i))];
end

set(handles.EventList,'String',slist);

guidata(hObject,handles);

function ClearList_Callback(hObject, eventdata, handles)
list=[];
handles.list=list;
set(handles.EventList,'String',num2str(list,'%g'));
guidata(hObject,handles);

function ClearOmit_Callback(hObject, eventdata, handles)
omitList=[];
handles.omitList=[];
handles.reasonList=[];
set(handles.OmitList,'String',num2str(omitList,'%g'));
handles.maxFOM=-1;
handles.minFOM=-1;
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function AddRangeToList_Callback(hObject, eventdata, handles)
onEvents=handles.onEvents;
offEvents=handles.offEvents;

loadParams=handles.loadParams;
if loadParams==1
   source='handles.Data.';
else
   source='handles.';
end

varChoice=handles.varChoice;
fomChoice=handles.fomChoice;
varList=handles.varList;

eval(['yVar=',source,char(varList{fomChoice}),';']);
eval(['xVar=',source,char(varList{varChoice}),';']);
eval(['regen=',source,'regen;']);
% eval(['xVar=handles.',char(varList{varChoice}),';']);
% eval(['yVar=handles.',char(varList{fomChoice}),';']);

xMin=handles.delayMin;
xMax=handles.delayMax;
yMax=handles.maxFOM;
yMin=handles.minFOM;

list=handles.list;
omitList=handles.omitList;

eventNum=1:length(xVar);
eventNum(omitList)=[];
xVar(omitList)=[];
yVar(omitList)=[];
regen(omitList)=[];

if xMin==-1
   xMin=min(xVar);
end

if xMax==-1
   xMax=max(xVar);
end
xMin
xMax
min(xVar)
max(xVar)

if yMin==-1
   yMin=min(yVar);
end
if yMax==-1
   yMax=max(yVar);
end

indExcX=xVar<xMin | xVar>xMax;
eventNum(indExcX)=[];
xVar(indExcX)=[];
yVar(indExcX)=[];
regen(indExcX)=[];

indExcY=yVar<yMin | yVar>yMax;
eventNum(indExcY)=[];
%xVar(indExcY)=[];
%yVar(indExcY)=[];
regen(indExcY)=[];

if onEvents-offEvents==1
   indExcR=(regen==0);
elseif onEvents-offEvents==-1
   indExcR=(regen==1);
end

if onEvents-offEvents~=0
   eventNum(indExcR)=[];
end

list=unique([list, eventNum]);
list=sort(list);
handles.list=list;

slist=[];
for i=1:length(list)
   slist=[slist,', ',num2str(list(i))];
end

set(handles.EventList,'String',slist);
guidata(hObject,handles);

function FixEvent_Callback(hObject, eventdata, handles)
list=handles.list;

fixList=handles.fixList;
delay=handles.delay;
regen=handles.regen;
for i=1:length(list)
   if sum(fixList==list(i))==0
      fprintf('\n Event %g delay changed: from %g to %g \n', list(i),delay(list(i)),delay(list(i)-1))
      delay(list(i))=delay(list(i)-1);
      regen(list(i))=regen(list(i)-1);
      fixList=[fixList list(i)];
   else
      fprintf('\n Event %g already corrected \n', list(i))
      beep;
   end
   
   
end
handles.delay=delay;
handles.regen=regen;

fixList=sort(fixList);
slist=[];
for i=1:length(fixList)
   slist=[slist,', ',num2str(fixList(i))];
end
handles.fixList=fixList;
set(handles.FixList,'String',slist);
UpdateCor(hObject, eventdata, handles)
guidata(hObject,handles);

function OmitEvent_Callback(hObject, eventdata, handles)
list=handles.list;
omitList=handles.omitList;
reasonList=handles.reasonList;
fomChoice=handles.fomChoice;
varList=handles.varList;

reason{1}=char(varList{fomChoice});

for i=1:length(list)
   if sum(omitList==list(i))==0
      %fprintf('\n Event %g has been omitted \n', list(i))
      omitList=[omitList list(i)];
      reasonList=[reasonList reason];
   else
      fprintf('\n Event %g already omitted \n', list(i))
      beep;
   end
end

[omitList ind]=sort(omitList);
reasonList=reasonList(ind);
slist=[];
for i=1:length(omitList)
   slist=[slist,', ',num2str(omitList(i))];
end
handles.omitList=omitList;
set(handles.OmitList,'String',slist);
handles.reasonList=reasonList;
handles.maxFOM=-1;
handles.minFOM=-1;
handles.delayMin=-1;
handles.delayMax=-1;
UpdateCor(hObject, eventdata, handles);
guidata(hObject,handles);

function EventList_Callback(hObject, eventdata, handles)
eventNum=int16(str2double(get(hObject,'String')));
list=handles.list;
if eventNum~=0
   list=sort(unique([list, eventNum]));
   handles.list=list;
end
slist=[];
for i=1:length(list)
   slist=[slist,', ',num2str(list(i))];
end
set(handles.EventList,'String',slist);
guidata(hObject,handles);

function Slider_CreateFcn(hObject, eventdata, handles)
MaxVal=1000;
set(hObject,'Max',MaxVal,'Min',1);
set(hObject,'SliderStep',[1/(MaxVal-1) .1]); %[arrow-step bar-step] in % change
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Slider.
function Slider_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on key press with focus on Slider and none of its controls.
function Slider_KeyPressFcn(hObject, eventdata, handles)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed

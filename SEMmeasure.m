%function [lengths] = image_measure;
% Estimate lengths on an image by point-click.
% [lengths] = image_measure;
%
% This function prompts the user to select an image file to open, then
% opens and plots the image.  The program waits for the user to resize the
% image if desired.  The user then:
% 1.) left-clicks on the endpoints of a known dimension in the image to
% establish a reference length and then
% 2.) types the length of the reference in a dialog box,
% 3.) clicks on as many pairs of points as desired in the image to measure
% their separation distance relative to the reference.
%
CONT=1;

close
%FOLDER='~/Dropbox/Research/AccelFab/SEM/101020b_FEMonQ4(955wBARC)';
%FOLDER='~/Dropbox/Research/AccelFab/SEM/101011_Q3-FEM';
FOLDER='~/Dropbox/Research/AccelFab/SEM/120403_Bonded_V2_1_CrossSec';
[filename,pathname] = uigetfile('*80kx*.JPG; *80kx*.TIF','Pick an Image File',FOLDER); % *.tif;*.jpg;*.png;*.bmp;
im_data = imread([pathname,filename]);

scrsz = get(0,'ScreenSize');
figure('Position',[5 1*scrsz(4)/6 3*scrsz(3)/4 3*scrsz(4)/4])
%figure; 
imagesc(im_data)
%colorbar hsv
colormap gray

%input('Press any key when ready to proceed');

axis image % make sure image is not distorted before starting
meas_num=0;
clear lengths
if ~CONT
%meas_num=0;
[xr,yr] = ginput(2);
line(xr,yr,'Marker','+','Color','c');
l_ref_img = norm([xr(1),yr(1)]-[xr(2),yr(2)]);


% Prompt for length of this line:
    prompt = {'Enter the length of the reference line:'};
    dlg_title = 'Length of Ref. Line';
    num_lines = 1;
    def = {'1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    lscale = str2num(answer{1})/l_ref_img;
end    
    % Displace labels slightly
    text(xr(1)+5,yr(1)+10,num2str(l_ref_img*lscale),'Color','c');

% Zoom out if Right mouse button is depressed.  If the middle button is
% depressed, zoom in again.
loop_break = 0; %meas_num = 1;
while loop_break == 0;
    
  
    [x1,y1,mbutton] = ginput(1);
    line(x1,y1,'Marker','+','Color','y');
    if mbutton == 3;
        loop_break = 1;
        % This signifies that the user is done
    elseif mbutton == 1;
        % Get second point and measure
        [x2,y2,mbutton] = ginput(1);
        if mbutton == 3; loop_break = 1; break; end
        meas_num = meas_num +1;
        line([x1,x2],[y1,y2],'Marker','+','Color','y');
        lengths(meas_num) = norm([x1,y1]-[x2,y2])*lscale;
%         dy=abs(y2-y1);
%         dx=abs(x2-x1);
%         if dx<dy
%         
        %text(x1,y1+2,num2str(lengths(meas_num)),'Color','y');
        text(x1+5,y1+10,num2str(meas_num),'Color','y');
        sprintf('\r %d: %5.2f',meas_num,lengths(meas_num))
        %meas_num = meas_num +1;
    end
end


lengths
%mean(lengths)
%std(lengths)
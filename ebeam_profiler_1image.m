

%cd V:\ARDB\E163\matlab\ControlAndDAQSoftware
VERBOSE=1;
PLOTEACHIMG=0;
SIG_OPTIONS=optimset('Display','off',...
    'TolFun',1e-05,...
    'TolX',1e-05);
showplots=0;
CalX=3.92;%6.45*4;

CalY=3.92;%6.45*4;

cXi=0;
cYi=0;

% filepath='V:\ARDB\E163\Data\140115\';
% img=imread([folder,'2.tif']);

%filepath=['V:\ARDB\E163\DATA\',datestr(now,'yymmdd'),'\'];

%[filename,filepath] = uigetfile({'*.tif;*.jpg;*.png'},'Choose file to process',filepath);

%img=imread([filepath filename]);


img=beam;

ImageAnalyze_E; %gives c1,c2,hor,ver,sigx,sigy

            % add ROI start point to centroid
            cxni = cXi+c1(2);
            cyni = cYi+c2(2);
            
%
        if(VERBOSE)
            figure (1)
            subplot(4,4,[1 10]);
            %figure
            imagesc(img)
            %xlim([100 400])
            %ylim([250 400])
            %title({'Background Subtracted & Rotated Image'})
            
            %title({[filename]})

            %colorbar
            [RN CN]=size(img);
            V=reshape(img,RN*CN,1);
%             if(PLOTEACHIMG)
%                 drawnow;
%             end
            figure (1)
            subplot(4,4,11:12);
            [nb,xb]=hist(double(V),256);
            semilogy(xb,nb,'k.')
%             if bd==1 % 16-bit
%                 xlim([0,4096]);
%             else % 8-bit
                 xlim([0,256]);
%            end
            title('Pixel Intensity Histogram')
            grid on;
%             if(PLOTEACHIMG)
%                 drawnow;
%             end
            
            
            figure (1)
            subplot(4,4,3:4);
            plot(hor);
            hold on;
            sumx=sum(f(c1,1:length(hor)));
            plot(f(c1,1:length(hor)),'r');
            hold on
            line([c1(2)-c1(3)*sqrt(2*log(2)) c1(2)+c1(5)*sqrt(2*log(2))],[c1(4)+c1(1)/2 c1(4)+c1(1)/2], 'LineWidth',2,'Color','k');
            
            
            grid on;
            %title(sprintf('Line Out X. sx = %f',sigx))
            %title({sprintf('FWHM_x = %f [pix] %f [um]',sigx/CalX ,sigx),sprintf('Centroid at %f [pix] %f [um]',cxni,cxni*CalX)});
            ax1 = gca;
            set(ax1,'YAxisLocation','right');
            ylabel('X line out')
            axis tight

            if(PLOTEACHIMG)
                drawnow;
            end
            hold off;
            
            figure (1)
            subplot(4,4,7:8);
            plot(ver);
            hold on;
            sumy=sum(f(c2,1:length(ver)));
            plot(f(c2,1:length(ver)),'g');
            hold on
            line([c2(2)-c2(3)*sqrt(2*log(2)) c2(2)+c2(5)*sqrt(2*log(2))],[c2(4)+c2(1)/2 c2(4)+c2(1)/2],'LineWidth',2,'Color','k');
            
            grid on;
            
            %title(sprintf('Line Out Y. sy = %f',sigy))
                        
            %title({sprintf('FWHM_y = %f [pix] %f [um]',sigy/CalY,sigy),sprintf('Centroid at %f [pix] %f [um]',cyni,cyni*CalY)});
            ax1 = gca;
            set(ax1,'YAxisLocation','right');
            ylabel('Y line out')
            axis tight
            
            if(PLOTEACHIMG)
                drawnow;
            end
            hold off;
          
        end
%                     Xpos=cxni
% Ypos=cyni
% Xwidth=sigx
% Ywidth=sigy


%fprintf('\n%s:\n',filename);
fprintf('X [um] = %f\n',cxni*CalX);
fprintf('Y [um] = %f\n',cyni*CalY);
fprintf('FWHM_x [um] = %f\n',sigx);
fprintf('FWHM_y [um] = %f\n',sigy);

% %%
% figure
% plot(Int)
% figure
% plotyy(1:length(Xpos),Xpos,1:length(Xpos),Ypos)
% 
% figure
% plotyy(1:length(Xpos),Xwidth,1:length(Xpos),Ywidth)

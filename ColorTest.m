%% INDIVIDUAL COLORS

clear all
%close all

%From colorbrewer
col1=[230 97 1]/255;    %dark orange
col2=[253 184 99]/255;  %light orange
col3=[178 171 210]/255; %light purple
col4=[94 60 153]/255;   %dark purple

%Modified for better B&W contrast
col0(1,:)=[1 1 1]; 
col0(2,:)=[255 205 135]/255;  %light orange
col0(3,:)=[166 160 220]/255; %light purple
col0(4,:)=[205 66 0]/255;    %dark orange
col0(5,:)=[47 08 117]/255;   %dark purple
col0(6,:)=[0 0 0];  

grey0=(col0(:,1).^1.8+col0(:,2).^1.8+col0(:,3).^1.8).^(1/1.8);
figure (1)
hold on
plot(grey0,'Color',col0(6,:))

%based on morgenstemning
col(1,:) = [1 1 1];        %white
col(2,:) = [252 229 1]./255;        %yellow
col(3,:)=[253 160 80]/255;          %light orange
col(4,:)=[220 87 1]/255;            %dark orange
col(5,:) = [20 50 95]./255;         %navy blue
col(6,:)=[0 0 0];                   %black

grey=(col(:,1).^1.8+col(:,2).^1.8+col(:,3).^1.8).^(1/1.8);
hold on
plot(grey,'Color',col(2,:))

%% COLORMAPS
N=128;
col1=colormap(morgenstemning(N,'invert',1));
col3=colormap(morgenstemning(3*N,'invert',1));
col=[col1(1:N/8,:);col3(3*N/8+1:3*N*7/8,:);col1(7*N/8+1:end,:)];

figure


subplot(1,3,1)
imagesc(peaks(256))
freezeColors;
cbfreeze(jet)

subplot(1,3,2)
imagesc(peaks(256))
col=colormap(morgenstemning(256,'invert',1));
colormap(col)
cbfreeze(col)
freezeColors;

subplot(1,3,3)
imagesc(peaks(256))
colormap(col)
colorbar
size(col)
%%

colP = [25 53 95]./255;
figure
plot(1:100,exp(1:100),'Color',colP)

%% More plots for testing

map=morgenstemning(7,'invert',1);
map2=map;%[map(1:4,:);.5*map(5,:)+.5*map(6,:);map(7:8,:)];

grey=(col(:,1).^1.8+col(:,2).^1.8+col(:,3).^1.8).^(1/1.8);
grey2=(map2(:,1).^1.8+map2(:,2).^1.8+map2(:,3).^1.8).^(1/1.8);

%%
figure
% 
x=linspace(0,2*pi,100);
% for i=1:8
%     plot(x,cos(x+(i-1)*2*pi/9),'Color',col(i,:),'LineWidth',2*(10-i));
%     hold on
% end
figure
for i=1:6
    plot(x,cos(x+(i-1)*.32),'Color',col0(i,:),'LineWidth',20);
    %plot(x,cos(x+(i-1)*.32),'Color',map2(i,:),'LineWidth',20);
    hold on
end
for i=1:6
    plot(x,cos(x+(i-1)*.4-2.7),'Color',col0(i,:),'LineWidth',10);
    %plot(x,cos(x+(i-1)*.4-2.9),'Color',map2(i,:),'LineWidth',10);
    hold on
end


% for i=1:8
%     plot(x,cos(x+1.6*pi),'Color',col(9-i,:),'LineWidth',2*(9-i));
%     hold on
% end

axis tight

export_fig Test2.eps -cmyk -r300 -painters%-painters
%export_fig Test2.eps -rgb -r300 -painters%-painters

figure
plot((1:length(grey))/length(grey),grey,'b')
hold on
plot((1:length(grey2))/length(grey2),grey2,'r')
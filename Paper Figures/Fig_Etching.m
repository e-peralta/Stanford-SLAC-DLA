clear all
close all

cols(1,:)=[255 205 135]/255;  %light orange
cols(2,:)=[166 160 220]/255; %light purple
cols(3,:)=[205 66 0]/255;    %dark orange
cols(4,:)=[47 08 117]/255;   %dark purple

N=102;
col1=colormap(morgenstemning(N,'invert',1));
col3=colormap(morgenstemning(3*N,'invert',1));
col0=[col1(1:round(N/8),:);col3(round(3*N/8)+1:round(3*N*7/8),:);col1(round(7*N/8)+1:end,:)];
%col0=[col3(1:round(3*N*3/4),:);col1(round(3*N/4)+1:end,:)];
%col0=[col1(1:round(N/2),:);col3(round(3*N/2)+1:round(3*N*7/8),:);col1(round(7*N/8)+1:end,:)];
col0size=size(col0)

%    Q27=struct([]);  
% save([pathname(1:60),'meas.mat']) %laptop
FOLDER='~/Dropbox/Research/AccelFab/SEM/120410_Q20_EtchInTrench/';
load([FOLDER,'meas.mat']) 

FOLDER='~/Dropbox/Research/AccelFab/SEM/120411_Q27_EtchInTrench/';
load([FOLDER,'meas.mat']) 

%

marktype={'^-','o-','*-','v-','x-','<-','+-','s-','>-','p-','h-','d-','.'};
colortype={'r','y','m','g','c','k','b',cols(4,:)};

% figure
% for i=1:size(Q20,2)
% errorbar(Q20(i).pos, Q20(i).depth,Q20(i).d_err,marktype{i},'color',colortype{i});
% hold on
% end
% 
% hold on
% %figure
% for i=1:size(Q27,2)
% errorbar(Q27(i).pos, Q27(i).depth,Q27(i).d_err,marktype{i},'color',colortype{i});
% hold on
% end

% Surface map

Q27m=zeros(54,34);
dmax=0;
k=0;
k2=0;
k3=0;
for i=1:2:length(Q27)
    for j=1:length(Q27(i).pos)
        y=1+8*(i-1);
        %x=-(1+3*(Q27(i).pos(j)-1));
        x=2+3*(Q27(i).pos(j)-1);
        Q27m(y,x)=Q27(i).depth(j);
        k=k+1;
        %Q20_800r(k,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
        Q27_r(k,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j) Q27(i).d_err(j)];
        %Want to use this to exclude the top and bottom most rows        
        if i>1 && i<length(Q27)-2
        k2=k2+1;
        %Q20_800rR(k2,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
        Q27_rR(k2,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j) Q27(i).d_err(j)];
        else
        k3=k3+1;
        %Q20_800rE(k3,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];            
        Q27_rE(k3,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j) Q27(i).d_err(j)];
        end
               
        if Q27(i).depth(j)>dmax
            dmax=Q27(i).depth(j);
        end
    end
end

for i=2:2:length(Q27)
    for j=1:length(Q27(i).pos)
        y=7+8*(i-2);
        %x=-(1+3*(Q27(i).pos(j)-1));
        x=2+3*(Q27(i).pos(j)-1);
        
        Q27m(y,x)=Q27(i).depth(j);
        k=k+1;
        %Q20_1200r(k,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
        Q27_r(k,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j) Q27(i).d_err(j)];
        if i>2 && i<length(Q27)-1
        k2=k2+1;
        %Q20_1200rR(k2,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
        Q27_rR(k2,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j) Q27(i).d_err(j)];
        else
        k3=k3+1;
        %Q20_1200rE(k3,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];            
        Q27_rE(k3,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j) Q27(i).d_err(j)];
        end
        
        if Q27(i).depth(j)>dmax
            dmax=Q27(i).depth(j);
        end
    end
end


Q20m=zeros(54,34);
%dmax=0;
k=0;
k2=0;
k3=0;
for i=1:2:length(Q20)
    for j=1:length(Q20(i).pos)
        y=1+8*(i-1);
        x=1+3*(Q20(i).pos(j)-1);
        Q20m(y,x)=Q20(i).depth(j);
        k=k+1;
        Q20_800r(k,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
        %Q27_r(k,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j) Q27(i).d_err(j)];
        %Want to use this to exclude the top and bottom most rows        
        if i>1 && i<length(Q20)-2
        k2=k2+1;
        Q20_800rR(k2,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
        %Q27_rR(k2,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j) Q27(i).d_err(j)];
        else
        k3=k3+1;
        Q20_800rE(k3,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];            
        %Q27_rE(k3,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j) Q27(i).d_err(j)];
        end
        
        if Q20(i).depth(j)>dmax
            dmax=Q20(i).depth(j);
        end
              
    end
end

Q20_r=Q20_800r;
k0=k;

k=0; 
k2=0;
k3=0;
for i=2:2:length(Q20)
    for j=1:length(Q20(i).pos)
        y=7+8*(i-2);
        x=1+3*(Q20(i).pos(j)-1);
        Q20m(y,x)=Q20(i).depth(j);
        k=k+1;
        Q20_1200r(k,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
        Q20_r(k0+k,:)=Q20_1200r(k,:);
        %Q27_r(k,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j) Q27(i).d_err(j)];
        if i>2 && i<length(Q20)-1
        k2=k2+1;
        Q20_1200rR(k2,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
        %Q27_r(k2,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j) Q27(i).d_err(j)];
        else
        k3=k3+1;
        Q20_1200rE(k3,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];            
        %Q27_r(k3,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j) Q27(i).d_err(j)];
        end
        
        if Q20(i).depth(j)>dmax
            dmax=Q20(i).depth(j);
        end
    end
end



width=4;
height=3.5;
hFig = figure(94);

set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')


imagesc(-34:0,1:55,fliplr(Q27m)./12)
colormap(col0)
%imagesc(Q27m./12)
axis xy
axis image
%title('.8 um and 1.2 um trench')
%title('No trench')
xlabel('X [mm]')
ylabel('Y [mm]')

%colorbar
%caxis([425/12 dmax/12])    
%hc=colorbar();
%hyc=ylabel(hc,'Etch Rate [A/s]')

%insert value next to data point
[Ni Nj]=size(Q27m);
for i=1:Ni
    for j=1:Nj
        if Q27m(i,j)~=0
            if i>28
            text(-j-1,i-2,num2str(int8(Q27m(i,j)/12)),'FontSize',8)
            
            %text(j,i-2,['\color{white}' num2str(int8(Q27m(i,j)/12))])
            %text(j,i-2,['\color{white}' num2str(int8(sqrt(j^2+(i-28)^2)))])
            else
            text(-j-1,i+2,num2str(int8(Q27m(i,j)/12)),'FontSize',8)
            %text(j,i+2,['\color{white}' num2str(int8(Q27m(i,j)/12))])
            %text(j,i+2,['\color{white}' num2str(int8(sqrt(j^2+(i-28)^2)))])
            end
        end
    end
end


%figure
hold on
imagesc(Q20m./12)
colormap(col0)
axis xy
axis image
%title('.8 um and 1.2 um trench')
%title('No trench')
hx=xlabel('X [mm]')
hy=ylabel('Y [mm]')

%colorbar
%caxis([30 51]);%425/12 dmax/12]) 
caxis([425/12 dmax/12]) 

%h=colorbar();
%hyc=ylabel(h,'Etch Rate [A/s]')

%insert value next to data point
[Ni Nj]=size(Q20m);
for i=1:Ni
    for j=1:Nj
        if Q20m(i,j)~=0
            if i>28
            text(j,i-2,num2str(int8(Q20m(i,j)/12)),'FontSize',8)
            %text(j,i-2,['\color{white}' num2str(int8(sqrt(j^2+(i-28)^2)))])
            else
            text(j,i+2,num2str(int8(Q20m(i,j)/12)),'FontSize',8)
            %text(j,i+2,['\color{white}' num2str(int8(sqrt(j^2+(i-28)^2)))])
            end
        end
    end
end

hold on

%insert circle contours
for r=5:5:35
rectangle('Position',[-r,28-r,2*r,2*r],'Curvature',[1,1],'LineWidth',1,'LineStyle',':');
%rectangle('Position',[-r,28-r,2*r,2*r],'Curvature',[1,1],'LineWidth',2,'LineStyle',':','EdgeColor','w')
hold on
end
line([0 0],[0 55],'color','k','LineWidth',1,'LineStyle',:)

set(gca,'FontSize', 10);
set([hx,hy],'FontSize', 10);
tightfig;

%export_fig EtchMap.eps -r300

%

Q27_r=sortrows(Q27_r);
Q20_r=sortrows(Q20_r);

Q20_800r=sortrows(Q20_800r);
Q20_1200r=sortrows(Q20_1200r);


width=3;
height=2.5;
hFig = figure(95);

%set(gcf,'PaperPositionMode','auto')
set(hFig,'ActivePositionProperty','position')
set(hFig,'Units','inches')
set(hFig, 'Position', [0 0 width height])

set(hFig, 'Color', 'w');
set(hFig,'Units','points')


%h1=errorbar(Q27_r(:,1),Q27_r(:,2)/12,Q27_r(:,3)/12,'xr');
h1=errorbar(Q27_r(:,1),Q27_r(:,2)/12,Q27_r(:,3)/12,'o','Color',cols(4,:),...
    'Linewidth',1,'MarkerSize',4,'MarkerFaceColor',cols(2,:));
hold on

% %h2=errorbar(Q20_800r(:,1),Q20_800r(:,2)/12,Q20_800r(:,3)/12,'ob');
h2=errorbar(Q20_r(:,1),Q20_r(:,2)/12,Q20_r(:,3)/12,'s','Color',cols(3,:),...
    'Linewidth',1,'MarkerSize',4,'MarkerFaceColor',cols(1,:));
hold on

% %h3=errorbar(Q20_1200r(:,1),Q20_1200r(:,2)/12,Q20_1200r(:,3)/12,'sg');
% h3=errorbar(Q20_1200r(:,1),Q20_1200r(:,2)/12,Q20_1200r(:,3)/12,':s','Color',cols(3,:),...
%     'Linewidth',1,'MarkerSize',4,'MarkerFaceColor',cols(1,:));

axis tight;

hx=xlabel('R [mm]')
hy=ylabel('Etch rate [A/s]')

hL=legend([h1,h2],{'On surface', 'In trench'});
set(hL,'Location','northwest','Box','Off')
set(gca,'FontSize', 10);
set([hx,hy,hL],'FontSize', 10);
tightfig;

export_fig EtchRates.eps -cmyk -r300 -painters%

% Find mean and std in uniform region

Q_r=[Q20_r;Q27_r];
Q_r=sortrows(Q_r);

figure

h1=errorbar(Q_r(:,1),Q_r(:,2)/12,Q_r(:,3)/12,'o','Color',cols(4,:),...
    'Linewidth',1,'MarkerSize',4,'MarkerFaceColor',cols(2,:));
hold on

mean(Q_r(5:26,2))/12
std(Q_r(5:26,2))/12


%%
% Q20_800rR=sortrows(Q20_800rR);
% Q20_800rE=sortrows(Q20_800rE);
% Q20_1200rR=sortrows(Q20_1200rR);
% Q20_1200rE=sortrows(Q20_1200rE);
% Q27_rR=sortrows(Q27_rR);
% Q27_rE=sortrows(Q27_rE);
% 
% figure
% errorbar(Q20_800rR(:,1),Q20_800rR(:,2)/12,Q20_800rR(:,3)/12,'ob-');
% hold on
% errorbar(Q20_800rE(:,1),Q20_800rE(:,2)/12,Q20_800rE(:,3)/12,'xc:');
% hold on
% errorbar(Q20_1200rR(:,1),Q20_1200rR(:,2)/12,Q20_1200rR(:,3)/12,'or-');
% hold on
% errorbar(Q20_1200rE(:,1),Q20_1200rE(:,2)/12,Q20_1200rE(:,3)/12,'xm:');
% %figure
% hold on
% errorbar(Q27_rR(:,1),Q27_rR(:,2)/12,Q27_rR(:,3)/12,'o-','color',cols(3,:));
% hold on
% errorbar(Q27_rE(:,1),Q27_rE(:,2)/12,Q27_rE(:,3)/12,'xg:');
% xlabel('R [mm]')
% ylabel('Etch rate [A/s]')
% legend('1200nm trench inner','1200nm trench outer',...
%     '800nm trench inner','800nm trench outer',...
%     'no trench inner','no trench outer')


%% Doing same thing but for trenches
% %close all
% clear Q27_r
% 
% Q10=struct();
% scale=27.5;
% Q10(1).depth=[1079 1082]/scale;
% Q10(2).depth=[1068 1081]/scale;
% Q10(1).pos=[1 7];
% Q10(2).pos=[1 7];
% 
% Q10(3).depth=[1122 1151]/scale;
% Q10(4).depth=[1147 1155]/scale;
% Q10(3).pos=[1 5];
% Q10(4).pos=[1 5];
% 
% Q10(5).depth=[1219 1245 1270 1309]/scale;
% Q10(6).depth=[1211 1239 1264 1316]/scale;
% Q10(5).pos=[1 5 11 11];
% Q10(6).pos=[1 5 11 11];
% 
% Q10(7).depth=[1346 1344 1344]/scale;
% Q10(8).depth=[1338 1340 1349]/scale;
% Q10(7).pos=[1 7 7];
% Q10(8).pos=[1 7 7];
% 
% Q27=Q10;
% Q27m=zeros(54,34);
% dmax=0;
% k=0;
% k2=0;
% k3=0;
% %skipped i=1 values something off
% for i=3:2:length(Q27)
%     for j=1:length(Q27(i).pos)
%         y=1+8*(i-1);
%         %x=-(1+3*(Q27(i).pos(j)-1));
%         x=2+3*(Q27(i).pos(j)-1);
%         Q27m(y,x)=Q27(i).depth(j);
%         k=k+1;
%         %Q20_800r(k,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
%         Q27_r(k,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j)];
%         %Want to use this to exclude the top and bottom most rows        
%         if Q27(i).depth(j)>dmax
%             dmax=Q27(i).depth(j);
%         end
%     end
% end
% 
% %skipped i=2 values something off
% for i=4:2:length(Q27)
%     for j=1:length(Q27(i).pos)
%         y=3+8*(i-2);
%         %x=-(1+3*(Q27(i).pos(j)-1));
%         x=2+3*(Q27(i).pos(j)-1);
%         
%         Q27m(y,x)=Q27(i).depth(j);
%         k=k+1;
%         %Q20_1200r(k,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
%         Q27_r(k,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j)];
%         
%         if Q27(i).depth(j)>dmax
%             dmax=Q27(i).depth(j);
%         end
%     end
% end
% 
% figure
% imagesc(Q27m./27.5*scale)
% %caxis([500/27.5*scale dmax/27.5*scale])
% axis xy
% %caxis([425/12 dmax/12]) 
% %imagesc(-34:0,1:55,fliplr(Q27m)./12)
% %insert value next to data point
% [Ni Nj]=size(Q27m);
% for i=1:Ni
%     for j=1:Nj
%         if Q27m(i,j)~=0
%             if i>28
%             %text(-j-1,i-2,num2str(int8(Q27m(i,j)/12)),'FontSize',8)
%             
%             text(j,i-2,['\color{white}' num2str(int8(Q27m(i,j)*scale/27.5))])
%             %text(j,i-2,['\color{white}' num2str(int8(sqrt(j^2+(i-28)^2)))])
%             else
%             %text(-j-1,i+2,num2str(int8(Q27m(i,j)/12)),'FontSize',8)
%             text(j,i+2,['\color{white}' num2str(int8(Q27m(i,j)*scale/27.5))])
%             %text(j,i+2,['\color{white}' num2str(int8(sqrt(j^2+(i-28)^2)))])
%             end
%         end
%     end
% end
% 
% 
% Q27_r=sortrows(Q27_r);
% size(Q27_r)
% 
% figure (12)
% plot(Q27_r(:,1),Q27_r(:,2)*scale/27.5,'x','Color',cols(4,:),...
%     'Linewidth',1,'MarkerSize',4,'MarkerFaceColor',cols(2,:));
% hold on
% 
% Qall=Q27_r;
% clear Q27_r
% Q20a=struct();
% 
% Q20a(1).depth=[1260 1350]/scale;
% Q20a(2).depth=[1282 1305]/scale;
% Q20a(1).pos=[1 7];
% Q20a(2).pos=[1 7];
% Q20a(3).depth=[1112 1195]/scale;
% Q20a(4).depth=[1209 1313]/scale;
% Q20a(3).pos=[1 11];
% Q20a(4).pos=[1 11];
% Q20a(5).depth=[1236 1370]/scale;
% Q20a(6).depth=[1239 1377]/scale;
% Q20a(5).pos=[1 11];
% Q20a(6).pos=[1 11];
% Q20a(7).depth=[1329 1376]/scale;
% Q20a(8).depth=[1354 1413]/scale;
% Q20a(7).pos=[1 7];
% Q20a(8).pos=[1 7];
% 
% Q27=Q20a;
% Q27m=zeros(54,34);
% dmax=0;
% k=0;
% k2=0;
% k3=0;
% for i=1:2:length(Q27)
%     for j=1:length(Q27(i).pos)
%         y=1+8*(i-1);
%         %x=-(1+3*(Q27(i).pos(j)-1));
%         x=2+3*(Q27(i).pos(j)-1);
%         Q27m(y,x)=Q27(i).depth(j);
%         k=k+1;
%         %Q20_800r(k,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
%         Q27_r(k,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j)];
%         %Want to use this to exclude the top and bottom most rows               
%         if Q27(i).depth(j)>dmax
%             dmax=Q27(i).depth(j);
%         end
%     end
% end
% 
% for i=2:2:length(Q27)
%     for j=1:length(Q27(i).pos)
%         y=3+8*(i-2);
%         %x=-(1+3*(Q27(i).pos(j)-1));
%         x=2+3*(Q27(i).pos(j)-1);
%         
%         Q27m(y,x)=Q27(i).depth(j);
%         k=k+1;
%         %Q20_1200r(k,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
%         Q27_r(k,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j)];
%         if Q27(i).depth(j)>dmax
%             dmax=Q27(i).depth(j);
%         end
%     end
% end
% 
% 
% figure
% imagesc(Q27m./27.5*scale)
% %caxis([500/27.5*scale dmax/27.5*scale])
% axis xy
% %caxis([425/12 dmax/12]) 
% %imagesc(-34:0,1:55,fliplr(Q27m)./12)
% %insert value next to data point
% [Ni Nj]=size(Q27m);
% for i=1:Ni
%     for j=1:Nj
%         if Q27m(i,j)~=0
%             if i>28
%             %text(-j-1,i-2,num2str(int8(Q27m(i,j)/12)),'FontSize',8)
%             
%             text(j,i-2,['\color{white}' num2str(int8(Q27m(i,j)*scale/27.5))])
%             %text(j,i-2,['\color{white}' num2str(int8(sqrt(j^2+(i-28)^2)))])
%             else
%             %text(-j-1,i+2,num2str(int8(Q27m(i,j)/12)),'FontSize',8)
%             text(j,i+2,['\color{white}' num2str(int8(Q27m(i,j)*scale/27.5))])
%             %text(j,i+2,['\color{white}' num2str(int8(sqrt(j^2+(i-28)^2)))])
%             end
%         end
%     end
% end
% 
% 
% Q27_r=sortrows(Q27_r);
% size(Q27_r)
% 
% figure (12)
% plot(Q27_r(:,1),Q27_r(:,2)*scale/27.5,'o','Color',cols(3,:),...
%     'Linewidth',1,'MarkerSize',4,'MarkerFaceColor',cols(1,:));
% hold on
% 
% Qall=[Qall;Q27_r];
% clear Q27_r
% scale=18.3;
% Q20b=struct();
% Q20b(1).depth=[795 830]/scale;
% Q20b(2).depth=[806 816]/scale;
% Q20b(1).pos=[1 7];
% Q20b(2).pos=[1 7];
% Q20b(3).depth=[766 871]/scale;
% Q20b(4).depth=[759 871]/scale;
% Q20b(3).pos=[1 11];
% Q20b(4).pos=[1 11];
% Q20b(5).depth=[778 764]/scale;
% Q20b(6).depth=[793 787]/scale;
% Q20b(5).pos=[1 11];
% Q20b(6).pos=[1 11];
% Q20b(7).depth=[873 931]/scale;
% Q20b(8).depth=[899 940]/scale;
% Q20b(7).pos=[1 7];
% Q20b(8).pos=[1 7];
% 
% Q27=Q20b;
% Q27m=zeros(54,34);
% dmax=0;
% k=0;
% k2=0;
% k3=0;
% for i=1:2:length(Q27)
%     for j=1:length(Q27(i).pos)
%         y=5+8*(i-1);
%         %x=-(1+3*(Q27(i).pos(j)-1));
%         x=2+3*(Q27(i).pos(j)-1);
%         Q27m(y,x)=Q27(i).depth(j);
%         k=k+1;
%         %Q20_800r(k,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
%         Q27_r(k,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j)];
%         %Want to use this to exclude the top and bottom most rows        
%         if Q27(i).depth(j)>dmax
%             dmax=Q27(i).depth(j);
%         end
%     end
% end
% 
% for i=2:2:length(Q27)
%     for j=1:length(Q27(i).pos)
%         y=7+8*(i-2);
%         %x=-(1+3*(Q27(i).pos(j)-1));
%         x=2+3*(Q27(i).pos(j)-1);
%         
%         Q27m(y,x)=Q27(i).depth(j);
%         k=k+1;
%         %Q20_1200r(k,:)=[sqrt(x^2+(y-28)^2) Q20(i).depth(j) Q20(i).d_err(j)];
%         Q27_r(k,:)=[sqrt(x^2+(y-28)^2) Q27(i).depth(j)];
%         if Q27(i).depth(j)>dmax
%             dmax=Q27(i).depth(j);
%         end
%     end
% end
% 
% 
% figure
% imagesc(Q27m./18.3*scale)
% %caxis([500/18.3 dmax/18.3])
% axis xy
% %caxis([425/12 dmax/12]) 
% %imagesc(-34:0,1:55,fliplr(Q27m)./12)
% %insert value next to data point
% [Ni Nj]=size(Q27m);
% for i=1:Ni
%     for j=1:Nj
%         if Q27m(i,j)~=0
%             if i>28
%             %text(-j-1,i-2,num2str(int8(Q27m(i,j)/12)),'FontSize',8)
%             
%             text(j,i-2,['\color{white}' num2str(int8(Q27m(i,j)*scale/18.3))])
%             %text(j,i-2,['\color{white}' num2str(int8(sqrt(j^2+(i-28)^2)))])
%             else
%             %text(-j-1,i+2,num2str(int8(Q27m(i,j)/12)),'FontSize',8)
%             text(j,i+2,['\color{white}' num2str(int8(Q27m(i,j)*scale/18.3))])
%             %text(j,i+2,['\color{white}'
%             %num2str(int8(sqrt(j^2+(i-28)^2)))])
%             end
%         end
%     end
% end
% 
% 
% Q27_r=sortrows(Q27_r);
% size(Q27_r)
% figure (12)
% plot(Q27_r(:,1),Q27_r(:,2)*scale/18.3,'s','Color',cols(4,:),...
%     'Linewidth',1,'MarkerSize',4,'MarkerFaceColor',cols(2,:));
% hold on
% 
% Qall=[Qall;Q27_r];
% Qall=sortrows(Qall);
% %%
% figure (95)
% hold on
% h3=plot(Qall(:,1),Qall(:,2),'+k','Linewidth',1,'MarkerSize',4);
% hold on
% h4=plot(Qall(:,1),Qall(:,2),'dk','Linewidth',1,'MarkerSize',4);
% %,'Color',cols(4,:),...
%   % 'Linewidth',1,'MarkerSize',4,'MarkerFaceColor',cols(2,:));
%   
% hL=legend([h1,h2,h3,h4],{'Grat on surface', 'Grat in trench','Trench'});
% set(hL,'Location','northwest','Box','Off')
% set(gca,'FontSize', 10);
% set([hx,hy,hL],'FontSize', 10);
% tightfig;
% 
% export_fig EtchRates2.eps -cmyk -r300 -painters%
% 
%   mean(Qall(1:16,2))
%   std(Qall(1:16,2))
%function f20160720_1_anim_with_contactPoints_1(simY)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
close all; %clear;
nft=3;
nfr=3;
ns=3;
nf=nft+nfr;
gD=.01;

y=simY.signals.values;

tm=length(y(:,1));
T=15;
t=0:T/(length(y(:,1))-1):T; t=t';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Joint variables
figure; hold on;
plot(t, y(:,7:15),'LineWidth',2);
grid on; 
%axis equal; 
set(gca,'FontSize',14); set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');
xlim([0 16]);xticks([0:2:16]);
ylim([-0.0075 0.015]);yticks([-0.0075:0.0025:0.015]);
lgd=legend({'l_{11}','l_{12}','l_{13}','l_{21}','l_{22}','l_{23}','l_{31}','l_{32}','l_{33}'},'NumColumns',3);
%lgd.FontSize = 12;
xlabel('Time [s]');
hold off;
     
% Plotting X,Y,Z of the base
figure;
plot(t, y(:,1:3),'LineWidth',2);
grid on; 
%axis equal; 
set(gca,'FontSize',20); set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');
xlim([0 16]);xticks([0:2:16]);xlabel('Time [s]');ylabel('Base translation [m]');
%ylim([-0.01 0.05]);yticks([-0.01:0.01:0.05]);
lgd=legend('x_{b}','y_{b}','z_{b}');
lgd.FontSize = 20;

% Plotting alpha, beta, gamma, of the base
figure;
plot(t, y(:,4:6),'LineWidth',2);
grid on; 
%axis equal; 
set(gca,'FontSize',20); set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');
xlim([0 16]);xticks([0:2:16]);xlabel('Time [s]'); ylabel('Base angular offsets [rad]');
%ylim([-0.01 0.05]);yticks([-0.01:0.01:0.05]);
lgd=legend('\alpha','\beta','\gamma');
lgd.FontSize = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp=length(y(1,:));
N=(tmp/2-nf)/ns;
dy=y(:,nf+N*ns+1:end);
y=y(:,1:nf+N*ns);
n=length(y(:,1));

size(y)
size(dy)


%writerObj=VideoWriter('anim1.avi');
% open(writerObj);

% for i=1:n
%     %     clc;
%     [i n]
%     clf;
%     q=y(i,:)';
%     dq=dy(i,:)';
%     [pc,vc,Jc] = Jx5(q,dq,nft,nfr,ns,gD);
%     
%     %     cof=.5*[max(pc(1,i))+min(x(1,:));max(x(2,:))+min(x(2,:));max(x(3,:))+min(x(3,:))];
%     
%     xrange=[min(pc(1,:))-.05 max(pc(1,:))+.05];
%     yrange=[min(pc(2,:))-.05 max(pc(2,:))+.05];
%     zrange=[min(pc(3,:))-.05 max(pc(3,:))+.05];
%     
%     
%     
%     % find which points are in contact
%     I=pc(3,:)<0;
%     PC=pc(:,I);
%     pc=pc(:,~I);
%     
%     
%     if i==70
%         fh=figure('Position',[1 1 800 600],'Visible','on');
%         fig1 = subplot(3,3,1:6);
%         hold on
%         
%         plot3(pc(1,:),pc(2,:),pc(3,:),'.','LineWidth',10); set(gca,'FontSize',12); set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');
%         plot3(PC(1,:),PC(2,:),PC(3,:),'.r','LineWidth',10);
%         
%         axis([xrange, yrange, zrange])
%         %xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
%         grid on
%         axis equal
% %         view(30,15);
%         view(-60,15);
%         
%         hold off
%         
%         fig2 = subplot(3,3,7); copyobj(allchild(fig1), fig2); view([-90 90]); grid on; axis equal;     axis([xrange, yrange, zrange]); set(gca,'FontSize',12); set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');%xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
%         fig3 = subplot(3,3,8); copyobj(allchild(fig1), fig3); view([90 0]); grid on; axis equal;    axis([xrange, yrange, zrange]);set(gca,'FontSize',12); set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');%xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
%         fig4 = subplot(3,3,9); copyobj(allchild(fig1), fig4); view([0 0]); grid on; axis equal;    axis([xrange, yrange, zrange]);set(gca,'FontSize',12); set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');%xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
%         
%         set(fh, 'PaperPositionMode', 'auto');
%         cdata = print('-RGBImage');
%         pause()
%     end
%     
% end
% 
% plot()
% plot all points


%     fh=figure('Position',[1 1 800 600],'Visible','on');
%     fig1 = subplot(3,3,1:6);
%     hold on
%
%     plot3(pc(1,:),pc(2,:),pc(3,:),'.','LineWidth',1);
%     plot3(PC(1,:),PC(2,:),PC(3,:),'.r','LineWidth',1);
%
%     axis([xrange, yrange, zrange])
%     xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
%     grid on
%     axis equal
%     view(30,15);
%     hold off
%
%     fig2 = subplot(3,3,7); copyobj(allchild(fig1), fig2); view([90 90]); grid on; axis equal;     axis([xrange, yrange, zrange]);xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
%     fig3 = subplot(3,3,8); copyobj(allchild(fig1), fig3); view([90 0]); grid on; axis equal;    axis([xrange, yrange, zrange]);xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
%     fig4 = subplot(3,3,9); copyobj(allchild(fig1), fig4); view([0 0]); grid on; axis equal;    axis([xrange, yrange, zrange]);xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
%
%     set(fh, 'PaperPositionMode', 'auto');
%     cdata = print('-RGBImage');

% writ the figure to the video object
%writeVideo(writerObj,im2frame(cdata));

%end

%close(writerObj);
%close all;

%end

function f20150825_3_animSimout_1(simY)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


fh=figure('Position',[1 1 400 300],'Visible','off');

y=simY.signals.values;
tmp=length(y(1,:));
N=(tmp/2-3)/2;
y=y(:,1:3+N*2);
n=length(y(:,1));
np=30;
x=zeros(2,np);
xi=linspace(0,N,np);

writerObj=VideoWriter('anim1.avi');
open(writerObj);

for i=1:n
    clf;
    for j=1:np
        % Jx(y,xi,nft,nfr,ns,N)
        [tmp,~,~] = Jx2(y(i,:),xi(j),2,1,2,N);
        x(:,j)=tmp;
%         x(:,j)=p(y(i,:),(j-1)/np);
    end
    plot(x(1,:),x(2,:),'LineWidth',2)
    grid on
    axis equal
    axis([-.5 .5 -.5 .5])
    
    set(fh, 'PaperPositionMode', 'auto');
    cdata = print('-RGBImage');
    
    % writ the figure to the video object
    writeVideo(writerObj,im2frame(cdata));
    
end

close(writerObj);
close all;

end


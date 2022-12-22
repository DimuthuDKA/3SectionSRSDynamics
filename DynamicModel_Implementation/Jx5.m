function [pc,vc,Jc] = Jx5(y,dy,nft,nfr,ns,gD)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

lp=3;
r=.0125;
L=.15;

nf=nft+nfr; % total float dof
N=round((length(y)-nf)/ns); % number of sections
dof=nf+ns*N; % total dof

% determine number of divisions on a plane
Nrz=ceil(2*pi/(gD/r));
nTv=linspace(0,2*pi,Nrz+1);nTv(end)=[];
cos_nTv=cos(nTv);
sin_nTv=sin(nTv);

% determine number of segments
Nz=ceil(L/gD);
nZv=linspace(0,1,Nz+1);
nEls=Nrz*((N-1)*Nz+(Nz+1)); % number of elements
pc=zeros(3,nEls);
vc=pc;
Jc=zeros(3*nEls,dof);

Jv=zeros(lp,dof);
Jw=zeros(lp,lp*dof);
tJ=zeros(lp,dof);
R=eye(3);
p=zeros(3,1);

% floating coordinate frame
% =========================
tmp=Tb(y(1:6));
tR=tmp(1:lp,1:lp);
tp=tmp(1:lp,lp+1);
[d1_Rb,d2_Rb,d3_Rb,~,~,~,~,~,~]=dd_d_Rb(y(4:6));

% Jacobian
Jv(:,1:nft)=eye(lp);
Jw(:,lp*nft+1:lp*nf)=[d1_Rb d2_Rb d3_Rb];
% position and orientation
p=p+R*tp;
R=R*tR;

% main loop
for i=1:N % traverse through all the sections
    nfj=nf+ns*(i-1);
    nfi=nfj+ns;
    q=y(nfi-2:nfi);
    
    if i~=N
        nZT=Nz;
    else
        nZT=Nz+1;
    end
    
    for j=1:nZT
        tmp=T(q,nZv(j));
        tR=tmp(1:lp,1:lp);
        tp=tmp(1:lp,lp+1);
        [d1_p,d2_p,d3_p]=d_px(q,nZv(j));
        [d1_R,d2_R,d3_R]=d_Rx(q,nZv(j));
        
        for k=1:Nrz
            ind=(i-1)*Nz*Nrz+Nrz*(j-1)+k;
            tTzr=[cos_nTv(k) -sin_nTv(k) 0;sin_nTv(k) cos_nTv(k) 0; 0 0 1]*[r;0;0];
            % contact jacobian
            tJ(:,1:nfi)=[Jv(:,1:nfj)+Jw(:,1:lp*nfj)*kron(eye(nfj),tp+tR*tTzr) R*[d1_p+d1_R*tTzr d2_p+d2_R*tTzr d3_p+d3_R*tTzr]];
            Jc(lp*ind-2:lp*ind,:)=tJ;
            vc(:,ind)=tJ*dy;
            % pos update
            tmp=p+R*(tp+tR*tTzr);
            pc(:,ind)=tmp;
        end
    end
    
    if i~=N
        % section updates
        % ===============
        % jacobian
        tmp=T(q,1);
        tR=tmp(1:lp,1:lp);
        tp=tmp(1:lp,lp+1);
        [d1_p,d2_p,d3_p]=d_px(q,1);
        [d1_R,d2_R,d3_R]=d_Rx(q,1);
        
        Jv(:,1:nfi)=[Jv(:,1:nfj)+Jw(:,1:lp*nfj)*kron(eye(nfj),tp) R*[d1_p d2_p d3_p]];
        Jw(:,1:lp*nfi)=[Jw(:,1:lp*nfj)*kron(eye(nfj),tR) R*[d1_R d2_R d3_R]];
        
        % position + orientation
        p=p+R*tp;
        R=R*tR;
    end
end
function Dy=f20160722_2_3secNewCoGwcntcts_1(u)
% compute forward dynamics with contac forces for a single section
% continuum arm with a floating base

N=3; % number of sections
nfr=3;
nft=3;
ns=3;
nf=nfr+nft;
n=nf+ns*N;

F=[zeros(nf,1);u(1:ns*N)];
y=u(ns*N+1:2*ns*N+nf);
dy=u(2*ns*N+nf+1:3*ns*N+2*nf);
Fc=u(3*ns*N+2*nf+1:4*ns*N+3*nf);
tmp=4*ns*N+3*nf;
% K_g=u(tmp+1);
% C_g=u(tmp+2);
% R_n=u(tmp+3);
D_coeff=u(tmp+1); % damping coefficient length variables
D_coeff_float=u(tmp+2);

Dy=zeros(2*n,1);

lp=3;
R=eye(lp);
p=zeros(lp,1);
Jv=zeros(lp,n);
Jw=zeros(lp,lp*n);
Hv=zeros(lp*n,n);
Hw=zeros(lp*n,lp*n);

M=zeros(n);
c=zeros(n,n,n);
C=zeros(n);
V=zeros(n,1);

m=.2; % sectio nmass
mb=.01; % base connecting plate
% C_g=u(1); % Ground damping coefficient
% R_n=u(2); % normal friction coefficient
% mu=.2;

g=[0 0 9.81]'; % gravitational acceleration

D=diag([D_coeff_float*ones(1,nf), D_coeff*ones(1,ns*N)]);

% stiffness regulation parameters
% lmin=0;
% lmax=.06;
Kmin=1750;
% Kmax=100000;
% k=2000;
% Ke=Kmin+.5*Kmax*((1-tanh(k*(y(nf+1:end)-lmin+.006)))+(1+tanh(k*(y(nf+1:end)-lmax-.006))));
Ke=Kmin*ones(ns*N,1);

cScale=[0.6658 - 0.9436i;
    0.9927 - 0.9037i;
    1.2648 - 0.8148i;
    1.5083 - 0.7040i;
    1.7057 - 0.6168i;
    1.8883 - 0.5107i;
    2.0277 - 0.4209i;
    2.1539 - 0.3827i;
    2.2901 - 0.2509i;
    2.3589 - 0.2598i];

for i=0:N % to account for floating base
    
    %%
    % ////////////////////////
    % ////compute MCV loop////
    % ////////////////////////
    if i==0
        tmp=Tb(y(1:6));
        tR=tmp(1:lp,1:lp);
        tp=tmp(1:lp,lp+1);
        [d1_Rb,d2_Rb,d3_Rb,d1d1_Rb,d1d2_Rb,d1d3_Rb,d2d2_Rb,d2d3_Rb,d3d3_Rb]=dd_d_Rb(y(4:6));
        Ix=.25*.0125^2+.01^2/12; Iy=Ix; Iz=.5*.0125^2;
        Ib=mb*[Ix 0 0;0 Iy 0;0 0 Iz];
        msk=[6;7;2];
        t=tR'*d1_Rb;wx=t(msk);
        t=tR'*d2_Rb;wy=t(msk);
        t=tR'*d3_Rb;wz=t(msk);
        % **** M ****
        M(1:3,1:3)=mb*eye(3);
        M(4:6,4:6)=[    wx'*Ib*wx   wx'*Ib*wy   wx'*Ib*wz;...
                        wy'*Ib*wx   wy'*Ib*wy   wy'*Ib*wz;...
                        wz'*Ib*wx   wz'*Ib*wy   wz'*Ib*wz];
        
        % **** dM ****
        t=d1_Rb'*d1_Rb+tR'*d1d1_Rb;xwx=t(msk);
        t=d2_Rb'*d1_Rb+tR'*d1d2_Rb;ywx=t(msk);
        t=d3_Rb'*d1_Rb+tR'*d1d3_Rb;zwx=t(msk);
        t=d1_Rb'*d2_Rb+tR'*d1d2_Rb;xwy=t(msk);
        t=d2_Rb'*d2_Rb+tR'*d2d2_Rb;ywy=t(msk);
        t=d3_Rb'*d2_Rb+tR'*d2d3_Rb;zwy=t(msk);
        t=d1_Rb'*d3_Rb+tR'*d1d3_Rb;xwz=t(msk);
        t=d2_Rb'*d3_Rb+tR'*d2d3_Rb;ywz=t(msk);
        t=d3_Rb'*d3_Rb+tR'*d3d3_Rb;zwz=t(msk);
        
        c(4:6,4:6,4)=c(4:6,4:6,4)+...
            [xwx'*Ib*wx+wx'*Ib*xwx xwx'*Ib*wy+wx'*Ib*xwy xwx'*Ib*wz+wx'*Ib*xwz;...
            xwy'*Ib*wx+wy'*Ib*xwx xwy'*Ib*wy+wy'*Ib*xwy xwy'*Ib*wz+wy'*Ib*xwz;...
            xwz'*Ib*wx+wz'*Ib*xwx xwz'*Ib*wy+wz'*Ib*xwy xwz'*Ib*wz+wz'*Ib*xwz];
        c(4:6,4:6,5)=c(4:6,4:6,5)+...
            [ywx'*Ib*wx+wx'*Ib*ywx ywx'*Ib*wy+wx'*Ib*ywy ywx'*Ib*wz+wx'*Ib*ywz;...
            ywy'*Ib*wx+wy'*Ib*ywx ywy'*Ib*wy+wy'*Ib*ywy ywy'*Ib*wz+wy'*Ib*ywz;...
            ywz'*Ib*wx+wz'*Ib*ywx ywz'*Ib*wy+wz'*Ib*ywy ywz'*Ib*wz+wz'*Ib*ywz];
        c(4:6,4:6,6)=c(4:6,4:6,6)+...
            [zwx'*Ib*wx+wx'*Ib*zwx zwx'*Ib*wy+wx'*Ib*zwy zwx'*Ib*wz+wx'*Ib*zwz;...
            zwy'*Ib*wx+wy'*Ib*zwx zwy'*Ib*wy+wy'*Ib*zwy zwy'*Ib*wz+wy'*Ib*zwz;...
            zwz'*Ib*wx+wz'*Ib*zwx zwz'*Ib*wy+wz'*Ib*zwy zwz'*Ib*wz+wz'*Ib*zwz];
        
        %  **** fill V *****
        V(1:3,1)=mb*eye(3)'*g;
        
        %/////// compute Jacobian/Hessian loop ///////
        %/////////////////////////////////////////////
        % floating base n
        
        % Hessian
        Hw(lp*nft+1:lp*nf,lp*nft+1:lp*nf)=[d1d1_Rb d1d2_Rb d1d3_Rb;d1d2_Rb d2d2_Rb d2d3_Rb;d1d3_Rb d2d3_Rb d3d3_Rb];
        % Jacobian
        Jv(:,1:nft)=eye(3);
        Jw(:,lp*nft+1:lp*nf)=[d1_Rb d2_Rb d3_Rb];
        
        % update p and R
        p=p+R*tp;
        R=R*tR;
        
    else
        nfj=nf+ns*(i-1);
        nfi=nfj+ns;
        q=y(nfi-2:nfi);
        % section n COG, derivatives
        [pc,d1_pc,d2_pc,d3_pc,~,~,~,~,~,~]=dd_d_pc(q,ones(10,1));
        [pBc,d1_pBc,d2_pBc,d3_pBc,d1d1_pBc,d1d2_pBc,d1d3_pBc,d2d2_pBc,d2d3_pBc,d3d3_pBc]=dd_d_pc(q,cScale);
        
        % M
        tJ1=Jv(:,1:nfj);
        tJ2=Jw(:,1:lp*nfj)*kron(eye(nfj),pc);
        tJ3=R*[d1_pc d2_pc d3_pc];
        tJB2=Jw(:,1:lp*nfj)*kron(eye(nfj),pBc);
        tJB3=R*[d1_pBc d2_pBc d3_pBc];
        M(1:nfi,1:nfi)=M(1:nfi,1:nfi)+real([tJ1'*(tJ1+tJ2)+tJ2'*tJ1+tJB2'*tJB2, tJ1'*tJ3+tJB2'*tJB3;(tJ1'*tJ3+tJB2'*tJB3)',tJB3'*tJB3]);
        
        hJ1=[Hv(1:lp*nfj,1:nfj);...
            Jw(:,1:lp*nfj)*kron(zeros(nfj),d1_pc);...
            Jw(:,1:lp*nfj)*kron(zeros(nfj),d2_pc);...
            Jw(:,1:lp*nfj)*kron(zeros(nfj),d3_pc)];
        hJ2=[Hw(1:lp*nfj,1:lp*nfj)*kron(eye(nfj),pc);...
            Jw(:,1:lp*nfj)*kron(eye(nfj),d1_pc);...
            Jw(:,1:lp*nfj)*kron(eye(nfj),d2_pc);...
            Jw(:,1:lp*nfj)*kron(eye(nfj),d3_pc)];
        hJB2=[Hw(1:lp*nfj,1:lp*nfj)*kron(eye(nfj),pBc);...
            Jw(:,1:lp*nfj)*kron(eye(nfj),d1_pBc);...
            Jw(:,1:lp*nfj)*kron(eye(nfj),d2_pBc);...
            Jw(:,1:lp*nfj)*kron(eye(nfj),d3_pBc)];
        hJB3=[reshape(Jw(:,[1:lp:lp*nfj,2:lp:lp*nfj,3:lp:lp*nfj]),[],lp)*[d1_pBc d2_pBc d3_pBc];...
            R*[d1d1_pBc d1d2_pBc d1d3_pBc];...
            R*[d1d2_pBc d2d2_pBc d2d3_pBc];...
            R*[d1d3_pBc d2d3_pBc d3d3_pBc]];
        
        for j=0:i
            if j==0
                for k=1:nf
                    th1=hJ1(lp*k-2:lp*k,:);
                    th2=hJ2(lp*k-2:lp*k,:);
                    thB2=hJB2(lp*k-2:lp*k,:);
                    thB3=hJB3(lp*k-2:lp*k,:);
                end
            else
                for k=1:ns
                    th1=hJ1(lp*(nf+ns*(j-1)+k)-2:lp*(nf+ns*(j-1)+k),:);
                    th2=hJ2(lp*(nf+ns*(j-1)+k)-2:lp*(nf+ns*(j-1)+k),:);
                    thB2=hJB2(lp*(nf+ns*(j-1)+k)-2:lp*(nf+ns*(j-1)+k),:);
                    thB3=hJB3(lp*(nf+ns*(j-1)+k)-2:lp*(nf+ns*(j-1)+k),:);
                end
            end
            tmp=[th1'*(tJ1+tJ2)+th2'*tJ1+thB2'*tJB2, th1'*tJ3+thB2'*tJB3;(th1'*tJ3+thB2'*tJB3)',thB3'*tJB3];
            c(1:nfi,1:nfi,nf+ns*(j-1)+k)=c(1:nfi,1:nfi,nf+ns*(j-1)+k)+real(tmp+tmp');
        end
        
        %  **** fill V *****
        V(1:nfi,1)=V(1:nfi,1)+[tJ1+tJ2 tJ3]'*g;
        
        %/////// compute Jacobian/Hessian loop ///////
        %/////////////////////////////////////////////
        tmp=T(q,1);
        tR=tmp(1:lp,1:lp);
        tp=tmp(1:lp,lp+1);
        
        [d1_p,d2_p,d3_p,d1d1_p,d1d2_p,d1d3_p,d2d2_p,d2d3_p,d3d3_p]=dd_d_p(q);
        [d1_R,d2_R,d3_R,d1d1_R,d1d2_R,d1d3_R,d2d2_R,d2d3_R,d3d3_R]=dd_d_R(q);
        
        % Hessians
        Hv(1:lp*nfi,1:nfi)=[...
            Hv(1:lp*nfj,1:nfj)+Hw(1:lp*nfj,1:lp*nfj)*kron(eye(nfj),tp),...
            reshape(Jw(:,[1:lp:lp*nfj,2:lp:lp*nfj,3:lp:lp*nfj]),[],lp)*[d1_p d2_p d3_p];...
            Jw(:,1:lp*nfj)*kron(eye(nfj),d1_p),...
            R*[d1d1_p d1d2_p d1d3_p];...
            Jw(:,1:lp*nfj)*kron(eye(nfj),d2_p),...
            R*[d1d2_p d2d2_p d2d3_p];...
            Jw(:,1:lp*nfj)*kron(eye(nfj),d3_p),...
            R*[d1d3_p d2d3_p d3d3_p]];
        
        Hw(1:lp*nfi,1:lp*nfi)=[Hw(1:lp*nfj,1:lp*nfj)*kron(eye(nfj),tR),...
            reshape(Jw(:,[1:lp:lp*nfj,2:lp:lp*nfj,3:lp:lp*nfj]),[],lp)*[d1_R d2_R d3_R];...
            Jw(:,1:lp*nfj)*kron(eye(nfj),d1_R),...
            R*[d1d1_R d1d2_R d1d3_R];...
            Jw(:,1:lp*nfj)*kron(eye(nfj),d2_R),...
            R*[d1d2_R d2d2_R d2d3_R];...
            Jw(:,1:lp*nfj)*kron(eye(nfj),d3_R),...
            R*[d1d3_R d2d3_R d3d3_R]];
        
        % Jacobians
        Jv(:,1:nfi)=[Jv(:,1:nfj)+Jw(:,1:lp*nfj)*kron(eye(nfj),tp) R*[d1_p d2_p d3_p]];
        Jw(:,1:lp*nfi)=[Jw(:,1:lp*nfj)*kron(eye(nfj),tR) R*[d1_R d2_R d3_R]];
        
        % update p and R
        p=p+R*tp;
        R=R*tR;
    end
    
    
end

% extract C
for k=1:n
    for j=1:n
        for i=1:n
            C(k,j)=C(k,j)+(c(k,j,i)+c(k,i,j)-c(i,j,k))*dy(i);
        end
    end
end

M=real(m*M); % generalized inertia matrix
C=real(.5*m*C); % Coriolis and centrifugal force matrix
V=real(m*V);
V(nf+1:n)=V(nf+1:n)+Ke.*y(nf+1:n); % conservative force matrix

% solving state equation
Dy(1:n,1)=dy;
Dy(n+1:2*n,1)=M\(F+Fc-(C+D)*dy-V);
% Dy(n+1:2*n,1)=pinv(M)*(F-(C+D)*dy-V);

end
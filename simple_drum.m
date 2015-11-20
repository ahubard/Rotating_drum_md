%% Simple rotating drum simulator based in molecular dynamics.
% Linear aproximation of the foR_drume.
%% Input parameters;
D = 1;                                        %Particle diameter
D_drum = 30*D;                                %Drum diameter
R_drum = D_drum/2;                              %Drum radius
N = ceil((D_drum/D)^2*0.8/2);                 %Number of particles
g = 0.001;                                        %Gravity
K = 1500;                                       %Elastic constant
rot_step = .2;                              %Rotation step
cos_th = cos(rot_step*pi/180);
sin_th = sin(rot_step*pi/180);
Qn = 0.8;                                        %Normal dissipation coefficient
mu = 1;                                        %static friction.
M = 1;


%% Initial conditions
% [x, y] = ndgrid(-R_drum:D:R_drum,-R_drum:D:R_drum);
% r = sqrt(x.^2+y.^2);
% ii = find(r<(R_drum-D/2));
% x = x(ii(1:N))';
% y = y(ii(1:N))';
% vx = randn(1,N)/50;
% vy = randn(1, N)/50;
load('Initial_Conditions.mat','xi','yi','vxi','vyi');
x = xi*cos(.2)-yi*sin(.2);
y = xi*sin(.2)+yi*cos(.2);

%x = xi;
%y = yi;
vx = vxi;
vy = vyi;
ax_old = zeros(1,N);
ay_old = zeros(1,N);
nt = 0;
cool_t = 0;
rotation_marker = [R_drum 0];
%% Simulation Parameters
dt = 0.1;
Nt = 1000000;
NR = ceil(360/rot_step);
Min_cool_down_time = 400;
Min_Energy = 0.0005;
%% Display Parameters
plotit = false;            % do we plot?
Nplotskip = 100;          % number of timesteps to skip before plotting

%% Save parameters
Nsave_skip = 200;
xs = zeros(N,Nt);
ys = zeros(N,Nt);
ts = zeros(1,Nt);
rots = zeros(1,Nt);
% vxs = zeros(N,Nt);
% vys = zeros(N,Nt);
Eks = zeros(1,Nt);
t_rot = zeros(1,NR);
%% Setup Plotting
if(plotit)
    figure(1);
    clf;
    h=zeros(1,N);
    hold on
    for np=1:N
        h(np)=rectangle('Position',[x(np)-.5.*D y(np)-.5.*D D D],'Curvature',[1 1],'edgecolor','b');
    end
    axis('equal');
    rectangle('Position',[-R_drum -R_drum D_drum D_drum],'Curvature',[1 1],'edgecolor','k');
    hc = plot(rotation_marker(1),rotation_marker(2),'xr');
    h_cm = plot(sum(x)/N,sum(y)/N,'mo');
end

%% Main Loop

for nr = 1:NR
    Ek = sum(vx.^2+vy.^2);
    %Rotation
    auxx = rotation_marker(1);
    auxy = rotation_marker(2);
    rotation_marker(1) = auxx*cos_th-auxy*sin_th;
    rotation_marker(2) = auxy*cos_th+auxx*sin_th;
    x_old = x;
    y_old = y;
    x = x_old*cos_th-y_old*sin_th;
    y = x_old*sin_th+y_old*cos_th;
    
    cool_t = 0;
    while(Ek > Min_Energy | cool_t < Min_cool_down_time)
        Ek = sum(vx.^2+vy.^2);
        nt = nt + 1;
        
        if (rem(nt-1,Nsave_skip)==0)
            Eks((nt-1)/Nsave_skip +1) = Ek;
        end
        cool_t = cool_t+1;
        if(plotit)
            if(rem(nt-1,Nplotskip)==0)
                
                figure(1);
                set(hc,'xdata',rotation_marker(1),'ydata',rotation_marker(2));
                set(h_cm,'xdata',sum(x)/N,'ydata',sum(y)/N);
                for np=1:N
                    
                    set(h(np),'Position',[x(np)-.5*D y(np)-.5*D D D]);
                end;
               
                drawnow;
                
                figure(2);
                plot(Eks(1:(nt-1)/Nsave_skip +1));
                axis([1 inf 0 inf])
                
            end;
        else
            if(rem(nt,1000)==0)
                fprintf('%2d.',nr);
                if(rem(nt,30000)==0)
                    fprintf('\n');
                end
            end
        end
        
        
        x = x+vx*dt+ax_old.*dt.^2/2;  % first step in Verlet integration
        y = y+vy*dt+ay_old.*dt.^2/2;
        
        if (rem(nt-1,Nsave_skip)==0)
            xs(:,(nt-1)/Nsave_skip +1) = x;
            ys(:,(nt-1)/Nsave_skip +1) = y;
            ts((nt-1)/Nsave_skip +1) = nt;
            rots((nt-1)/Nsave_skip +1) = nr;
        end
        % Interaction detector
        
        Fx = zeros(1,N);
        Fy = zeros(1,N);
        
        %particle-particle forces
        for nn = 1:N
            for mm = nn+1:N
                dx = x(mm)-x(nn);
                if(abs(dx) < D)
                    dy = y(mm)-y(nn);
                    dnm = dx^2+dy^2;
                    if(dnm < D^2)                   %Check if particles are in contact
                        dnm = sqrt(dnm);              %distancebtw particle n-m
                        
                        dvx = vx(mm)-vx(nn);          %relative velocities
                        dvy = vy(mm)-vy(nn);
                        
                        nx = dx/dnm;                  %Normal to contact line
                        ny = dy/dnm;
                        vn = dvx*nx+dvy*ny;             %Normal velocity
                        
                        delta = abs(D-dnm)/2;
                        
                        Fn = (-K*delta+Qn*vn).*sqrt(delta);  %Normal force of interactions.
                        
                        Fx(nn)=Fx(nn)+Fn.*nx;  % particle-particle Force Law
                        Fx(mm)=Fx(mm)-Fn.*nx;
                        Fy(nn)=Fy(nn)+Fn.*ny;  % particle-particle Force Law
                        Fy(mm)=Fy(mm)-Fn.*ny;
                    end
                end
            end
        end
        
        % Wall interaction;
        
        r2 = x.^2+y.^2;
        ii = find(r2 > (R_drum-D/2)^2); %Particles in contact with walls
        r = sqrt(r2(ii));
        nx = x(ii)./r;
        ny = y(ii)./r;
        vb_n = vx(ii).*nx+vy(ii).*ny;  %normal component of velocity
        delta = r-(R_drum-D/2);
        Fn = -(K*delta+Qn*vb_n).*sqrt(delta); %Normal force
        vb_t = vx(ii).*ny-vy(ii).*nx;   %Tangential velocity
        Ft = mu*sign(vb_t).*(Fn-g*ny);
        %Ft = -mu*sign(vb_t);
        
        
        Fx(ii) = Fx(ii)+Fn.*nx+Ft.*ny;
        Fy(ii) = Fy(ii)+Fn.*ny-Ft.*nx;
        
        ax = Fx./M;
        ay = Fy./M-g;
        
        
        vx = vx+(ax_old+ax).*dt/2;  % second step in Verlet integration
        vy = vy+(ay_old+ay).*dt/2;
        
        
        ax_old = ax;
        ay_old = ay;
        Ek = sum(vx.^2+vy.^2);
       
    end
    t_rot(nr) = nt;
    
end



%           tx = -ny;                       %Tangential
%           ty = nx;
%
%           Vs=-dvx*tx-dvy*ty+w(nn)*R+w(mm)*R;
%           as_old=-vax_old*tx-vay_old*ty+aw_old(nn)*R+aw_old(mm)*R;
%           cmm=find(cc(:,nn)==mm, 1);
%           if(isempty(cmm))
%             cmm=find(cc(:,nn)==0,1);
%             cc(cmm,nn)=mm;
%           end
%           cc_new(cmm,nn)=1;
%           delta=(1-dnm/D);
%           Fn=(-K*delta+Qn*Vn).*sqrt(delta);
%           tmp=xs(cmm,nn)+Vs*dt+as_old*dt^2/2;
%           xs(cmm,nn)=(tmp/xx-fix(tmp/xx))*xx;
%           Fsf=-abs(Fn)*mu*xs(cmm,nn)/xx;
%           Fs=Fsf-Qs*Vs;  %*(Ft>abs(Qs*Vs));
%           Fx(nn)=Fx(nn)+Fn.*nx+Fs*tx;
%           Fy(nn)=Fy(nn)+Fn.*ny+Fs*ty;
%           Tw(nn)=Tw(nn)+Fs*R;
%           Fx(mm)=Fx(mm)-Fn.*nx-Fs*tx;
%           Fy(mm)=Fy(mm)-Fn.*ny-Fs*ty;
%           Tw(mm)=Tw(mm)+Fs*R;
%         end;
%       end;
%     end;
%   end;
%
%   cc=cc.*cc_new;
%   xs=xs.*cc_new;
%
%   %% walls
%   r=sqrt(x.^2+y.^2);
%   nn=find(r>R_drum-R);
%   if(~isempty(nn))
%     nx=x(nn)./r(nn);
%     ny=y(nn)./r(nn);
%     tx=-ny;
%     ty=nx;
%     dvx=R_drum*wc*tx-vx(nn);
%     dvy=R_drum*wc*ty-vy(nn);
%     vax_old=-wc^2/R_drum*nx-ax_old(nn);
%     vay_old=-wc^2/R_drum*ny-ay_old(nn);
%     Vn=dvx.*nx+dvy.*ny;
%     Vs=-dvx.*tx-dvy.*ty+w(nn)*R-wc*R_drum;
%     as_old=-vax_old.*tx-vay_old.*ty+aw_old(nn)*R-wc^2/R_drum*R_drum;
%     delta=(R+r(nn)-R_drum)/D;
%     Fn=(-K*delta+Qn*Vn).*sqrt(delta);
%     xswl(nn)=xswl(nn)+Vs*dt+as_old*dt^2/2;
%     xswl(~nn)=0;
%     xswl(nn)=(xswl(nn)/xx-fix(xswl(nn)/xx))*xx;
%     Fsf=-abs(Fn).*mu.*xswl(nn)/xx;
%     Fs=Fsf-Vs*Qs;
%     Fx(nn)=Fx(nn)+Fn.*nx+Fs.*tx;
%     Fy(nn)=Fy(nn)+Fn.*ny+Fs.*ty;
%     Tw(nn)=Tw(nn)+Fs*R;
%   end
%   %%
%   ax=Fx./M-B*vx;
%   ay=Fy./M-g-B*vy;
%   aw=Tw./I;
%
%   vx=vx+(ax_old+ax).*dt/2;  % second step in Verlet integration
%   vy=vy+(ay_old+ay).*dt/2;
%   w=w+(aw_old+aw).*dt/2;
%
%   px(:,nt) = x;
%   py(:,nt) = y;
%   RR(nt)=-sum(y)-1i*sum(x);
%   Ek(nt)=sum(M.*(vx.^2+vy.^2)/2);
%   Er(nt)=sum(M.*(D*w).^2)/8*I;
% %   Em(nt)=Ek(nt)+Er(nt);
%
%   ax_old=ax;
%   ay_old=ay;
%   aw_old=aw;
%
% end



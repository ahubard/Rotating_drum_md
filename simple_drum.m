%% Simple rotating drum simulator based in molecular dynamics.
% Linear aproximation of the foR_drume.
%% Input parameters;
% D = 1;                                        %Particle diameter
% D_drum = 30*D;                                %Drum diameter
% R_drum = D_drum/2;                              %Drum radius
% N = ceil((D_drum/D)^2*0.8/2);                 %Number of particles
g = 0.001;                                        %Gravity
K = 1500;                                       %Elastic constant
rot_step = .2;                              %Rotation step
cos_th = cos(rot_step*pi/180);
sin_th = sin(rot_step*pi/180);
Qn = 30;                                        %Normal dissipation coefficient
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
load('Initial_Conditions.mat','xi','yi','vxi','vyi','DD','D','N','D_drum','R_drum');
initial_rot = .2;
x = xi*cos(initial_rot)-yi*sin(initial_rot);
y = xi*sin(initial_rot)+yi*cos(initial_rot);
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
Min_Energy = 0.004;
%% Display Parameters
plotit = false;
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
stop_save = 1000;
save_count = 1;
file_save = sprintf('Simple_rotdrum_');
%% Setup Plotting
if(plotit)
    figure(1);
    clf;
    h=zeros(1,N);
    hold on
    for np=1:N
        h(np)=rectangle('Position',[x(np)-.5.*D(np) y(np)-.5.*D(np) D(np) D(np)],'Curvature',[1 1],'edgecolor','b');
    end
    axis('equal');
    rectangle('Position',[-R_drum -R_drum D_drum D_drum],'Curvature',[1 1],'edgecolor','k');
    hc = plot(rotation_marker(1),rotation_marker(2),'xr');
    h_cm = plot(sum(x)/N,sum(y)/N,'mo');
    axis([-R_drum R_drum -R_drum R_drum]);
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
                    
                    set(h(np),'Position',[x(np)-.5*D(np) y(np)-.5*D(np) D(np) D(np)]);
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
                    figure(2);
                plot(Eks(1:(nt-1)/Nsave_skip +1));
                axis([1 inf 0 inf])
                    fprintf('\n');
                end
            end
        end
        
        
        x = x+vx*dt+ax_old.*dt.^2/2;  % first step in Verlet integration
        y = y+vy*dt+ay_old.*dt.^2/2;
        
        if (rem(nt-1,Nsave_skip)==0)
            i_ts = (nt-1)/Nsave_skip +1;
            xs(:,i_ts) = x;
            ys(:,i_ts) = y;
            ts(i_ts) = nt;
            rots(i_ts) = nr;
            if(rem(i_ts,stop_save) == 0)
                save(sprintf('%s%i.mat',file_save,i_ts/stop_save),'xs','ys'...
                    ,'Eks','ts','rots');
            end
            
        end
        % Interaction detector
        
        Fx = zeros(1,N);
        Fy = zeros(1,N);
        
        %particle-particle forces
        for nn = 1:N
            for mm = nn+1:N
                dx = x(mm)-x(nn);
                if(abs(dx) < DD(nn,mm))
                    dy = y(mm)-y(nn);
                    dnm = dx^2+dy^2;
                    if(dnm < DD(nn,mm)^2)                   %Check if particles are in contact
                        dnm = sqrt(dnm);              %distancebtw particle n-m
                        
                        dvx = vx(mm)-vx(nn);          %relative velocities
                        dvy = vy(mm)-vy(nn);
                        
                        nx = dx/dnm;                  %Normal to contact line
                        ny = dy/dnm;
                        vn = dvx*nx+dvy*ny;             %Normal velocity
                        
                        delta = abs(DD(nn,mm)-dnm)/2;
                        
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
        ii = find(r2 > (R_drum-D/2).^2); %Particles in contact with walls
        r = sqrt(r2(ii));
        nx = x(ii)./r;
        ny = y(ii)./r;
        vb_n = vx(ii).*nx+vy(ii).*ny;  %normal component of velocity
        delta = r-(R_drum-D(ii)/2);
        Fn = -(K*delta+Qn*vb_n).*sqrt(delta); %Normal force
        vb_t = vx(ii).*ny-vy(ii).*nx;   %Tangential velocity
        Ft = mu*sign(vb_t).*(Fn-g*ny);
        
        
        
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




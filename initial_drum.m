%%  Particles in a circle under gravity

%% Input parameters;
D = 1;                                         %Particle diameter
D_drum = 50*D;                                 %Drum diameter
R_drum = D_drum/2;                             %Drum radius
N = ceil((D_drum/D)^2/2)+60;                   %Number of particles
g = 0.001;                                     %Gravity
K = 1500;                                      %Elastic constant
Qn = 30;                                      %Normal dissipation coefficient                                       %static friction.
M = 1;
radius_ratio = 1/1.4;
Nsmall = floor(N/2);

%% Initial conditions
[x, y] = ndgrid(-R_drum:D:R_drum,-R_drum:D:R_drum);
r = sqrt(x.^2+y.^2);
ii = find(r < (R_drum-D/2));
x = x(ii(1:N))';
y = y(ii(1:N))';
ii_small = randi([1 N],1,Nsmall);
D = D*ones(1,N);
D(ii_small) = radius_ratio*D(ii_small);
DD = repmat(D,length(D),1)/2;
DD = DD + DD';
vx = randn(1,N)/5000;
vy = randn(1, N)/5000;
ax_old = zeros(1,N);
ay_old = zeros(1,N);
nt = 0;
cool_t = 0;
rotation_marker = [R_drum 0];
%% Simulation Parameters
dt = 0.05;
Nt = 10000;

%% Display Parameters
plotit = true;            % do we plot?
Nplotskip = 300;          % number of timesteps to skip before plotting


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

for nt = 1:Nt
    
    Eks(nt) = sum(vx.^2+vy.^2);
    
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
            plot(Eks(1:nt));
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
    
    
    Fx(ii) = Fx(ii)+Fn.*nx;
    Fy(ii) = Fy(ii)+Fn.*ny;
    
    ax = Fx./M;
    ay = Fy./M-g;
    
    
    vx = vx+(ax_old+ax).*dt/2;  % second step in Verlet integration
    vy = vy+(ay_old+ay).*dt/2;
    
    
    ax_old = ax;
    ay_old = ay;
    
    
end

    




close all; clear all; clc;


%% Physical parameters with dimensions

kB = 1.38e-23;              % Boltzmann constant [J/K]
k = 1e-6;                   % spring potential stiffness [N/m]
rad = 6.6e-6;               % radius of cell [m] (from Leesa)
alpha = 1e-9;               % attractive force coefficient [N.m^2]
T = 300;                    % bath temperature [K]
eta = 1e-3;                 % coefficient of viscosity for the bath [Pa.s]
gamma = 6*pi*rad*eta;       % friction coefficient of the bath [kg/s]
DT = kB*T/gamma;            % translational diffusion coefficient [m^2/s]
gamma_r = 8*pi*eta*rad^3;   % rotational friction coefficient of a sphere [kg.m^3]
DR = kB*T/gamma_r;          % rotational diffusion coefficient [rad^2/s] 
rho = 1000;                 % density of water [kg/m^3]
rho_fe = 5300;              % density of hematite (the material of the colloids) [kg/m^3]
volume = 4/3*pi*rad^3;      % approx volume of particle [m^3]
m = 27e-15;                 % approx mass of cell [kg] (taken from red blood cell number---not used here)
I = 2/5*m*rad^2;            % moment of inertia of a sphere about its central axis [kg.m^2]
kappa = 5e-10;%5e-15;              % flocking torque stiffness constant [N.m]
u = 3.75e-8;%3.75e-8; %2.9e-6*60;    % self-propulsion velocity [m/s] (The avg vel from Leesa was given as 2.25 um/min.)

Nbox_yi = 1.5*rad;           % Range in the y direction (positive and negative) [m] relative to particle i in which possible particle collisions are checked. If particle j, which may collide with i, is not within range Nboxyi of i collision time is not checked.
Nbox_xi = 1.5*rad;           % Range in the x direction (positive and negative) [m] relative to particle i in which possible particle collisions are checked. If particle j, which may collide with i, is not within range Nboxxi of i collision time is not checked.

error = 1e-25;              % [UNUSED] negligible number used in the denominator of various equations to prevent equations from blowing up
smallerr = 1e-5;            % [UNUSED] number used for preventing particles from going through eachother (in case of identical collision times for various particles)
maxfract_p = 0.75;          % [UNUSED] maximum area available for particles other than particle 1 (the large particle)
fracfact = 0.7;             % [UNUSED] fraction factor (it is multiplied with maxfractionp to vary amount of particles placed in system)

N_time = 1e3;               % number of time steps
dt = 60;                    % time step [s]
t = 0:dt:(N_time-1)*dt;     % time array [s]

AoE_rad = 20*rad; % radius within which particles begin orientation "flocking"
                  % (this is an arbitrarily chosen value right now)


%% Initialization

box_dim = 2e-4;             % box dimension [m]
Area_box = box_dim^2;       % [UNUSED] area of box [m^2]
frac = maxfract_p*fracfact; % [UNUSED] fraction of available area for particles
%N_part = floor(Area_box*frac/(pi*rad^2));  % [UNUSED] number of particles placed in the system
N_part = 10;

%--- Translation and rotational position ---%
X = zeros(N_part,N_time);
Y = zeros(N_part,N_time);
Theta = unifrnd(0,2*pi,[N_part 1]); % select angle from uniform distribution

%--- Noise ---%
noise_X = randn(N_part,N_time);
noise_Y = zeros(N_part,N_time);
noise_theta = randn(N_part,N_time);


%================ NON-OVERLAPPING INITIAL POSITIONS ==============%

for i = 1:N_part

    X(i,1) = (2*rand - 1)*(box_dim - rad*1.1);
%     Y(i,1) = -0.5e-3;
    Y(i,1) = 0.2e-4;

    there_is_overlap = false;

    while there_is_overlap == true

        there_is_overlap = false;

        for j = 1:i-1

            dXvec = sqrt((X(i,1) - X(j,1))^2);
            dYvec = sqrt((Y(i,1) - Y(j,1))^2);

            if dXvec^2 + dYvec^2 < 1.01*(2*rad)^2
                there_is_overlap = false;
            end

        end

        if there_is_overlap == true
            X(i,1) = (2*rand - 1)*(box_dim - rad*1.1);
%             Y(i,1) = -0.5e-3;
            Y(i,1) = 0.2e-4;
        end

    end

end


%% Simulation


for tt = 2:N_time


    %================== RELATIVE POSITION ARRAYS =================%

    R = [X(:,tt-1)'; Y(:,tt-1)']; % position vector
    
    %--- 3D matrix approach to displacement vector (easier for python translation)
    %    Dim 1: particle of interest (element num is # particles)
    %    Dim 2: x and y components (element num is # spatial dims; x and y)
    %    Dim 3: particle applying the force to the particle of interest (element num is # particles)
    for ri = 1:N_part
        for rj = 1:N_part 
            script_R(ri,:,rj) = R(:,ri) - R(:,rj); % direction vector from the particle of interest to the source of the force (a different particle)
            mag_script_R(ri,rj) = sqrt(dot(script_R(ri,:,rj),script_R(ri,:,rj))); % magnitude of the direction vector (a 2D matrix rather than a 3D matrix)
            script_R_hat(ri,:,rj) = script_R(ri,:,rj)/mag_script_R(ri,rj); % direction unit vector
        end
    end

    % Replace NaNs with zeros
    mag_script_R(isnan(mag_script_R)) = 0;
    script_R_hat(isnan(script_R_hat)) = 0;


    %================ ATTRACTIVE (CLUSTERING) FORCE ==============%
    
    [row_idx_attract, col_idx_attract] = find(mag_script_R ~= 0 & mag_script_R > 2*rad); % to prevent infinite force
    F_attractive = zeros(N_part, 2, N_part);
    invSq_script_R = script_R.^-2;
    invSq_script_R(isinf(invSq_script_R) | isnan(invSq_script_R)) = 0; % changes elements of inf to zero
    F_attractive(row_idx_attract,:,col_idx_attract) = F_attractive(row_idx_attract,:,col_idx_attract) - alpha./invSq_script_R(row_idx_attract,:,col_idx_attract).*script_R_hat(row_idx_attract,:,col_idx_attract);
    F_attractive(isnan(F_attractive)) = 0; % replace NaNs with zeros
    F_attractive_X = sum(F_attractive(:,1,:),3);
    F_attractive_Y = sum(F_attractive(:,2,:),3);


    %======================= FLOCKING TORQUE =====================%

    %--- find row/col of flocking interacting particle
    [row_idx_flock, col_idx_flock] = find(mag_script_R ~= 0 & mag_script_R <= AoE_rad);
    [row_idx_flock, I] = sort(row_idx_flock);
    col_idx_flock = col_idx_flock(I);
    row_col_flock_mat = [row_idx_flock, col_idx_flock];
    prev_row_idx_flock = 0; % set to -1 in python

    for ii = 1:length(row_idx_flock)

        % The row of indices corresponding to particles exhibiting flocking interactions has repeated values. 
        % To avoid redundant calculations, if a repeated index value is "hit," we skip to the next iteration 
        % of the loop.
        if row_idx_flock(ii) == prev_row_idx_flock
            continue % skip to the next iteration in the loop
        end

        temp_mat = row_col_flock_mat(row_col_flock_mat(:,1) == row_idx_flock(ii),:); % isolate the section of the matrix corresponding to rows with values equaling row_idx_flock(ii) 
        
        %--- find average theta of particles in column
        theta_avg = zeros(N_part, 1);
        theta_avg(row_idx_flock(ii)) = mean(Theta(temp_mat(:,2),tt-1)); % theta average corresponding to torque applied to i-th particle

        prev_row_idx_flock = row_idx_flock(ii);
    end

    T_flocking = -kappa*(Theta(:,tt-1) - theta_avg);


    %====================== ROTATIONAL STEPS ======================% 

    Theta(:,tt) = Theta(:,tt-1) + sqrt(2*DR*dt)*noise_theta(:,tt-1) + T_flocking/gamma_r*dt;
    Theta(:,tt) = (Theta(:,tt)/(2*pi) - floor(Theta(:,tt)/(2*pi))) * 2*pi; % ensures 0 <= Theta <= 2*pi


    %===================== TRANSLATIONAL STEPS ===================%

    % Step size in x direction from different forces
    delX_noise = sqrt(2*DT*dt)*noise_X(:,tt-1);
    delX_active = u*cos(Theta(:,tt-1))*dt;
    delX_attractive = 1/gamma*F_attractive_X*dt;

    % Step size in y direction from different forces    
    delY_noise = sqrt(2*DT*dt)*noise_Y(:,tt-1);
    delY_active = u*sin(Theta(:,tt-1))*dt;
    delY_attractive = 1/gamma*F_attractive_Y*dt;

    X(:,tt) = X(:,tt-1) + delX_noise + delX_active + delX_attractive;
    Y(:,tt) = Y(:,tt-1) + delY_noise + delY_active + delY_attractive;


    %==================== HARD SPHERE REPULSION ==================%
    
    there_is_overlap = false; % THIS WAS CHANGED

    R = [X(:,tt)'; Y(:,tt)']; % position vector
    while there_is_overlap == true

        disp('In While loop.')

        there_is_overlap = false;
    
        
        row_overlap = [];
        col_overlap = [];
    
        %--- 3D matrix approach to displacement vector (easier for python translation)
        %    Dim 1: particle of interest (element num is # particles)
        %    Dim 2: x and y components (element num is # spatial dims; x and y)
        %    Dim 3: particle applying the force to the particle of interest (element num is # particles)
        for ri = 1:N_part
            for rj = 1:N_part 
                script_R(ri,:,rj) = R(:,ri) - R(:,rj); % direction vector from the particle of interest to the source of the force (a different particle)
                mag_script_R(ri,rj) = sqrt(dot(script_R(ri,:,rj),script_R(ri,:,rj))); % magnitude of the direction vector (a 2D matrix rather than a 3D matrix)
                script_R_hat(ri,:,rj) = script_R(ri,:,rj)/mag_script_R(ri,rj); % direction unit vector
                
                % find overlapping particles
                if mag_script_R(ri,rj) > 0 && mag_script_R(ri,rj) < 2*rad
                    there_is_overlap = false;
                    row_overlap = [row_overlap, ri];
                    col_overlap = [col_overlap, rj];
                end
    
            end
        end
    
        for ii = 1:length(row_overlap)
    
            R(:,row_overlap(ii)) = R(:,row_overlap(ii)) + (2*rad - mag_script_R(row_overlap(ii),col_overlap(ii)))*script_R_hat(row_overlap(ii),:,col_overlap(ii))';
    
        end

    end

    disp('Exited While loop.')

    X(:,tt) = R(1,:);
    Y(:,tt) = R(2,:);


    %================ PERIODIC BOUNDARY CONDITIONS ==============%
    
    % Added code
    % Using Pythargorean's Theorem to determine if particle is out of
    % bounds by marking the ones with a radius greater than 1 out of bounds
    
    part_radius = sqrt(power(X(:,tt), 2) + power(Y(:,tt), 2));
    beyond_radius_bound = find(part_radius > box_dim);
    location_x = X(beyond_radius_bound, tt);
    location_y = Y(beyond_radius_bound, tt);
    X(beyond_radius_bound, tt:end) = location_x - (X(beyond_radius_bound, tt:end) - location_x);
    Y(beyond_radius_bound, tt:end) = location_y - (Y(beyond_radius_bound, tt:end) - location_y);

    beyond_lower_bound = find(Y(:,tt) < 0);
    Y(beyond_lower_bound, tt:end) = -Y(beyond_lower_bound, tt:end);
    
% 
%     % Move particles beyond positive x bound
%     idx_beyond_PosX_bound = find(X(:,tt) > box_dim);
%     X(idx_beyond_PosX_bound,tt) = X(idx_beyond_PosX_bound,tt) - 2 * box_dim;
% 
%     % Move particles beyond negative x bound
%     idx_beyond_NegX_bound = find(X(:,tt) < -box_dim);
%     X(idx_beyond_NegX_bound,tt) = X(idx_beyond_NegX_bound,tt) + 2 * box_dim;
% 
%     % Move particles beyond positive x bound
%     idx_beyond_PosY_bound = find(Y(:,tt) > box_dim);
%     Y(idx_beyond_PosY_bound,tt) = Y(idx_beyond_PosY_bound,tt) - 2 * box_dim;
% 
%     % Move particles beyond negative x bound
%     idx_beyond_NegY_bound = find(Y(:,tt) < -box_dim);
%     Y(idx_beyond_NegY_bound,tt) = Y(idx_beyond_NegY_bound,tt) + 2 * box_dim;
%     
%     % NOTE: Currently, distance between particles is calculated only by
%     % considering particles within the box, as if the box were finite and
%     % not subject to the periodic conditions that are meant to imitate an
%     % infinite landscape. I'm not sure if my approach is correct.

end


%% Colors for plotting

% custom colors
red = 1/255*[255 102 102];
blue = 1/255*[0 128 255]; 
teal = 1/255*[0 153 153]; 
magenta = 1/255*[153 0 153]; 
orange = 1/255*[255 128 0];
lightgreen = 1/255*[0 153 76];
lightgray = 1/255*[211 211 211];
lightblue = 1/255*[0 211 255];
violet = 1/255*[127 0 255];

% Matlab colors
ML_blue = [0, 0.4470, 0.7410];
ML_orange = [0.8500, 0.3250, 0.0980];	  
ML_yellow = [0.9290, 0.6940, 0.1250];
ML_purple = [0.4940, 0.1840, 0.5560];
ML_green = 	[0.4660, 0.6740, 0.1880];
ML_cyan = [0.3010, 0.7450, 0.9330];
ML_bergundy = [0.6350, 0.0780, 0.1840];


%% Video

figure, set(gcf, 'Color','black')
set(gca, 'nextplot','replacechildren', 'Visible','off');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
vidObj = VideoWriter('Brownian_sim_test_05.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);



for i = 1:N_time
    hold on  
%     txt = "X: " + X(1,i) + " Y: " + Y(1,i) + " R: " + sqrt(power(X(1,i), 2) + power(Y(1,i), 2));
%     text1 = text(X(1,i),Y(1,i),txt,'FontSize',14,'Color','white');
%     txt = "X: " + X(2,i) + " Y: " + Y(2,i) + " R: " + sqrt(power(X(2,i), 2) + power(Y(2,i), 2));
%     text2 = text(X(2,i),Y(2,i),txt,'FontSize',14,'Color','white');
    h1 = circles(X(:,i),Y(:,i),rad,'edgecolor',1/255*[0 128 0],'facecolor',1/255*[152,251,152]);
    h2 = rectangle('Position',[-box_dim, -box_dim, 2*box_dim, 2*box_dim],'EdgeColor','white','Curvature',[1 1]);
    hold off
    axis equal;
    xlim([-box_dim box_dim]) 
    ylim([0 box_dim])
    writeVideo(vidObj, getframe(gcf));
    delete(h1);
%     delete(text1);
%     delete(text2);
end

close(gcf)
close(vidObj);







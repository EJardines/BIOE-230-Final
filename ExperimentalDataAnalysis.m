
% This script is to plot the MSD of the cell position data given to me by
% Leesa.

clc; clear all; close all;

%% Import and organize data, and calculate MSD

% Import data from Leesa's Excel file
% %                                                                         file name                          sheet name               beginning-to-end columns
% [TrackN, SliceN, X1, Y1, Distance, Velocity, PixelValue] = ImportExcelFile("ExperimentalData.xlsx", "Results from MAX_20210920_Sox17", [1, 3432]);
[TrackN, SliceN, X1, Y1, Distance, Velocity, PixelValue] = ImportExcelFile("ExperimentalData_02.xlsx", "RESULTS", [1, 1818]);


% Remove missing (NaN) entries in the columns
TrackN = rmmissing(TrackN);
SliceN = rmmissing(SliceN);
X1 = rmmissing(X1);
Y1 = rmmissing(Y1);
Distance = rmmissing(Distance);
dt = 1; % experimental time step is one minute

% Separate particle quantities into distinct cell array elements
uniquecells = unique(TrackN);
for i = 1:length(uniquecells)

    idx = find(TrackN == uniquecells(i));
    X{i} = X1(idx);
    Y{i} = Y1(idx);
    del_X{i} = diff(X{i}); % find the delta in X (used to later find theta)
    del_Y{i} = diff(Y{i}); % find the delta in X (used to later find theta)
    t{i} = 0:dt:(length(X{i})-1)*dt;  % time array [min]

    % Set initial distance an displacement values
    distance{i}(1) = 0;
    displacement{i}(1) = 0;

    for j = 1:length(del_X{i})

        % If the displacements are zero, assume the orientation is the same 
        % this frame as in the previous frame.
        if del_Y{i}(j) == 0 && del_X{i}(j) == 0
            Theta{i}(j) = prevTheta;
        else
            Theta{i}(j) = atan2(del_Y{i}(j),del_X{i}(j));
        end

        prevTheta = Theta{i}(j);
    
        distance{i}(j+1) = distance{i}(j) + sqrt(del_X{i}(j)^2 + del_X{i}(j)^2);
        displacement{i}(j+1) = sqrt((X{i}(j+1) - X{i}(1))^2 + (Y{i}(j+1) - Y{i}(1))^2);

    end
    
    ConfinementRatio{i} = displacement{i}./distance{i};

end

% Translational and rotational MSD
[MSD_trans, Tau_trans] = CalcTransMSD(X, Y, dt);
[MSD_rot, Tau_rot] = CalcRotMSD(Theta, dt);

% Average the confinement ratio
NumTraj = length(ConfinementRatio);                                              
for i = 1:NumTraj
    ConfinementRatio_matrix(i,1:length(ConfinementRatio{i})) = ConfinementRatio{i};
    TrajLength_CR(i) = length(ConfinementRatio{i});
end
ConfinementRatio_mean = mean(ConfinementRatio_matrix,1,'omitnan');
maxTrajLength_CR = max(TrajLength_CR);


%% Average translational trajectory MSD

for i = 1:length(X)
    TrajLength_trans(i) = length(MSD_trans{i});
    MaxMSD_trans_vec(i) = max(MSD_trans{i});
end

[TrajMax_trans, idx_trans] = max(TrajLength_trans); % idx is the index corresponding to the longest trajectory
                                                    % In this case, multiple trajectories share the max length.
[MaxMSD_trans, idx_maxMSD_trans] = max(MaxMSD_trans_vec); % idx_maxMSD corresponds to the index of the maximum MSD value   
                                                          % Fortunately, the max MSD is also within an MSD array with the greatest number of elements

NumTraj = length(MSD_trans);                                              
for i = 1:NumTraj
    MSD_trans_matrix(i,1:length(MSD_trans{i})) = MSD_trans{i};
end
MSD_trans_mean = mean(MSD_trans_matrix,1,'omitnan');


%% Average rotational trajectory MSD

for i = 1:length(Theta)
    TrajLength_rot(i) = length(MSD_rot{i});
    MaxMSD_rot_vec(i) = max(MSD_rot{i});
end

[TrajMax_rot, idx_rot] = max(TrajLength_rot); % idx is the index corresponding to the longest trajectory
                                              % In this case, multiple trajectories share the max length.
[MaxMSD_rot, idx_maxMSD_rot] = max(MaxMSD_rot_vec); % idx_maxMSD corresponds to the index of the maximum MSD value   
                                                    % Fortunately, the max MSD is also within an MSD array with the greatest number of elements

NumTraj = length(MSD_rot);                                              
for i = 1:NumTraj
    MSD_rot_matrix(i,1:length(MSD_rot{i})) = MSD_rot{i};
end
MSD_rot_mean = mean(MSD_rot_matrix,1,'omitnan');


%% Linearly fit to MSDs (get diffusion coefficients)

d = 2; % number of dimensions in the simulation

log_MSD_trans_mean = log(MSD_trans_mean);
log_MSD_trans_mean = log_MSD_trans_mean(1:end-1); % removes the NaN at the final element
log_Tau_trans = log(Tau_trans{idx_trans}(1:end-1)); % make size the same as MSD vector
[log_tau_trans_fit, log_MSD_trans_fit, slope_trans, log_y_int_trans] = LinearFit(log_Tau_trans, log_MSD_trans_mean, length(log_MSD_trans_mean));
tau_trans_fit = exp(log_tau_trans_fit);
MSD_trans_fit = exp(log_MSD_trans_fit);

D_trans = exp(log_y_int_trans)/(2*d); % translational "diffusion" coefficient [um^2/min]
print_D_trans = sprintf('D_trans = %d um^2/s', D_trans);
disp(print_D_trans);

log_MSD_rot_mean = log(MSD_rot_mean);
log_MSD_rot_mean = log_MSD_rot_mean(1:end-1); % removes the NaN at the final element
log_Tau_rot = log(Tau_rot{idx_rot}(1:end-1)); % make size the same as MSD vector
[log_tau_rot_fit, log_MSD_rot_fit, slope_rot, log_y_int_rot] = LinearFit(log_Tau_rot, log_MSD_rot_mean, length(log_MSD_rot_mean));
tau_rot_fit = exp(log_tau_rot_fit);
MSD_rot_fit = exp(log_MSD_rot_fit);

D_rot = exp(log_y_int_rot)/(2*d); % rotational "diffusion" coefficient [rad^2/min]
print_D_rot = sprintf('D_rot = %d rad^2/s', D_rot);
disp(print_D_rot);


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


%% Plot

fontsize = 15;
markersize = 3; 
linewidth = 2.5;

% Translational MSD vs time lag for all trajectories
fh1 = figure(1); 
hold on
for i = 1:length(X)
    loglog(Tau_trans{i},MSD_trans{i},'o','MarkerSize', markersize, 'LineWidth', linewidth);
end
hold off
xlabel('$\tau$ [min]','Interpreter','latex')
ylabel('$\langle (\Delta\mathbf{r})^2\rangle$  [$\mu$m$^2$]','Interpreter','latex')
PlotSettingsLogLog(fontsize,fh1)

% Rotational MSD vs time lag for all trajectories
fh2 = figure(2); 
hold on
for i = 1:length(X)
    loglog(Tau_rot{i},MSD_rot{i},'o','MarkerSize', markersize, 'LineWidth', linewidth);
end
hold off
xlabel('$\tau$ [min]','Interpreter','latex')
ylabel('$\langle (\Delta\theta)^2\rangle$  [rad$^2$]','Interpreter','latex')
% PlotSettings(fontsize,fh2)
PlotSettingsLogLog(fontsize,fh2)

% Average translational MSD vs time lag
fh3 = figure(3);
hold on
h1 = loglog(tau_trans_fit,MSD_trans_fit,'-','Color',orange,'MarkerSize', markersize, 'LineWidth', linewidth);
h2 = loglog(Tau_trans{idx_trans},MSD_trans_mean,'o','Color',blue,'MarkerSize', markersize, 'LineWidth', linewidth);
hold off
xlabel('$\tau$ [min]','Interpreter','latex')
ylabel('$\langle (\Delta\mathbf{r})^2\rangle_{\mathrm{avg}}$  [$\mu$m$^2$]','Interpreter','latex')
leg1 = legend([h2(1), h1(1)],...
           'data', 'fit',...
           'Location','northwest');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',fontsize);
legend boxoff
PlotSettingsLogLog(fontsize,fh3)

% Average rotational MSD vs time lag
fh4 = figure(4);
hold on
h1 = loglog(tau_rot_fit,MSD_rot_fit,'-','Color',orange,'MarkerSize', markersize, 'LineWidth', linewidth);
h2 = loglog(Tau_rot{idx_rot},MSD_rot_mean,'o','Color',blue,'MarkerSize', markersize, 'LineWidth', linewidth);
hold off
xlabel('$\tau$ [min]','Interpreter','latex')
ylabel('$\langle (\Delta\theta)^2\rangle_{\mathrm{avg}}$  [rad$^2$]','Interpreter','latex')
leg1 = legend([h2(1), h1(1)],...
           'data', 'fit',...
           'Location','southwest');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',fontsize);
legend boxoff
PlotSettingsLogLog(fontsize,fh4)
% PlotSettings(fontsize,fh4)

% Confinement ratio vs time
fh5 = figure(5);
hold on
for i = 1:length(t)
    h1 = plot(t{i},ConfinementRatio{i},'o','MarkerSize', markersize, 'LineWidth', linewidth);
end
hold off
xlabel('$t$ [min]','Interpreter','latex')
ylabel('C.R.$(t)$','Interpreter','latex')
% leg1 = legend([h2(1), h1(1)],...
%            'data', 'fit',...
%            'Location','northwest');
% set(leg1,'Interpreter','latex');
% set(leg1,'FontSize',fontsize);
% legend boxoff
% PlotSettingsLogLog(fontsize,fh4)
PlotSettings(fontsize,fh5)

fh6 = figure(6);
plot(1:maxTrajLength_CR,ConfinementRatio_mean,'o','Color',blue,'MarkerSize', markersize, 'LineWidth', linewidth)
xlabel('$t$ [min]','Interpreter','latex')
ylabel('C.R.$_{\mathrm{avg}}(t)$','Interpreter','latex')
PlotSettings(fontsize,fh6)


%% Local functions

% Calculate the translation Mean-Squared Displacement (MSD)
function [MSD, Tau] = CalcTransMSD(X, Y, dt)

% X: Cell array for the positions along the x dimension.
% Y: Cell array for the positions along the y dimension.
% dt: Time step size.

    NumParticles = length(X); % number of particles

    for i = 1:NumParticles
        NumTimeSteps = length(X{i}); % number of time steps
        for tau = 1:NumTimeSteps
            DispX = X{i}(1+tau:end) - X{i}(1:end-tau); % displacement along the x direction
            DispY = Y{i}(1+tau:end) - Y{i}(1:end-tau); % displacement along the y direction
            SqDispR = DispX.^2 + DispY.^2; % R^2
            MSD{i}(tau) = mean(SqDispR);
            Tau{i}(tau) = tau*dt;

        end
    end
end

% Calculate the rotational Mean-Squared Displacement (MSD)
function [MSD, Tau] = CalcRotMSD(Theta, dt)

% Theta: Cell array for the angular positions in radians.
% dt:    Time step size.

    NumParticles = length(Theta); % number of particles

    for i = 1:NumParticles
        NumTimeSteps = length(Theta{i}); % number of time steps
        for tau = 1:NumTimeSteps
            Disp = Theta{i}(1+tau:end) - Theta{i}(1:end-tau); % displacement of theta
            SqDisp = Disp.^2;
            MSD{i}(tau) = mean(SqDisp);
            Tau{i}(tau) = tau*dt;

        end
    end
end

% Calculate the linear fit
function [x_out, y_out, slope, y_int] = LinearFit(x_in, y_in, NumFitPts)

    % The below linear fit calculation uses the method of least-squares.

    % x_in:      the input vector of x values
    % y_in:      the input vector of y values
    % NumFitPts: the number of points the fit is desired to have.
    % x_out:     the output vector of x values
    % y_out:     the output linear-fitted vector of y-values
    % slope:     the slope of the output
    % y_int:     the y-intercept of the output
    
    n = numel(x_in);
    xmax = max(x_in);
    xmin = min(x_in);
    
    sum_x = sum(x_in);
    sum_y = sum(y_in);
    sum_xx = sum(x_in.*x_in);
    sum_yy = sum(y_in.^2);
    sum_xy = sum(x_in.*y_in);
    sum_x2 = sum_x * sum_x;
    sum_y2 = sum_y * sum_y;
    
    slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x2); % slope
    y_int = (1 / n) * (sum_y - slope * sum_x); % y-intecept
    
    s = sqrt((sum_yy - sum_y2 / n - slope * (sum_xy - sum_x * sum_y / n)) / (n-2));
    e_m = s * sqrt(n / (n * sum_xx - sum_x2)); % error in the slope
    e_b = s * sqrt(sum_xx / (n * sum_xx - sum_x2)); % error in y-intercept
    r = (sum_xy - sum_x * sum_y / n) *  sqrt((sum_xx - sum_x2 / n) * (sum_yy - sum_y2 / n));
    
    x_out = linspace(xmin, xmax, NumFitPts);
    y_out = slope .* x_out + y_int;

end

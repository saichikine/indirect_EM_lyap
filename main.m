% Sai Chikine
% ASEN 6020 Optimal Trajectories
% Final Project

%% let there be light
clear; close all; clc;
set(groot,'defaultLineLineWidth', 1.5)
set(groot,'defaultAxesFontSize', 16)
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
p = gcp;

%% Load (setup)
load('environment_EM_spacecraft_smallsat.mat');

%%

% ODE tolerances
ode_opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

%% Load ref traj
load('ref_connection.mat');
fprintf('Reference trajectory loaded.\n');

%% Get ref traj with lots of points
tic
[t_ref,X_ref_hist] = ode113(@(t,X) CR3BP(t,X,mu_EM), linspace(0,connection_result.traj{1}(end),10000), connection_result.X0, ode_opts);
toc

ref_traj = {t_ref, X_ref_hist};

% %% Plot reference trajectory (2D)
% 
% figure
% addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
% plot(1-mu_EM, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Moon'); hold on % Moon
% plot(x_L1, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
% plot(x_L2, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
% plot(connection_result.traj{2}(1,1), connection_result.traj{2}(2,1), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
% plot(connection_result.traj{2}(1,:), connection_result.traj{2}(2,:), 'r-','DisplayName', 'Reference Trajectory (Heteroclinic Connection)'); hold on
% title('Reference Trajectory')
% xlabel('X')
% ylabel('Y')
% axis('equal')
% grid on;
% legend();
% hold off

%% 

% From plot, pick start and end states
index_start = 660;
%index_start = 853;
index_end = 1160;

% index_start = 1;
% index_end = length(connection_result.traj{1});

rng(50);

X0 = connection_result.traj{2}(:,index_start) + [normrnd(0,1.3e-6,[2,1]);0;normrnd(0,10e-4,[2,1]);0];
%X0 = connection_result.traj{2}(:,index_start) + [normrnd(0,1e-5,[2,1]);0;normrnd(0,1e-5,[2,1]);0];
t0 = connection_result.traj{1}(index_start);
Xf = connection_result.traj{2}(:,index_end);
tf = connection_result.traj{1}(index_end);

tfmin0 = tf-t0;
cf_guess = 1;

%% Plot reference trajectory (from index_start to index_end)

figure
addToolbarExplorationButtons(gcf)
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot(1-mu_EM, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Moon'); hold on % Moon
plot(x_L1, 0, 'dk', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
plot(x_L2, 0, 'dk' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
plot(connection_result.traj{2}(1,index_start), connection_result.traj{2}(2,index_start), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
plot(connection_result.traj{2}(1,index_start:index_end), connection_result.traj{2}(2,index_start:index_end), 'r-','DisplayName', 'Reference Trajectory (Heteroclinic Connection)'); hold on
plot(connection_result.traj{2}(1,index_end), connection_result.traj{2}(2,index_end), 'ok', 'markerfacecolor', 'r', 'DisplayName', 'Target Point'); hold on
title('Reference Trajectory','FontSize',14)
xlabel('X')
ylabel('Y')
axis('equal')
grid on;
legend('FontSize',12);
%%

% Get initial guess for costates from linear problem (fixed-endpoint LQR)

l = 6; % number of states
R = 2*eye(3);
A_func = @(X) Jacobian(X,mu_EM); % A matrix depends on state (time varying)
B = [zeros(3,3); 47*max_thrust_mag/m_sc*eye(3)]; %dfdu
augSTM_0 = eye(12,12);

X_ref = X0;
Xlarge0 = [X_ref; reshape(augSTM_0,[],1)];
deltaX0 = X_ref - connection_result.traj{2}(:,index_start);

%ode_opts = odeset('RelTol',3e-14,'AbsTol',1e-22);
%[~, Xlarge_loop_hist] = ode113(@(t,Xlarge) CR3BP_costate_STM_dynamics(t,Xlarge,B,R,mu_EM), [t0, t0+cf_guess*tfmin0], Xlarge0, ode_opts);
[~, Xlarge_loop_hist] = ode113(@(t,Xlarge) CR3BP_costate_STM_dynamics(t,Xlarge,B,R,mu_EM), ...
    [t0, t0+cf_guess*tfmin0], [X_ref; reshape(eye(12),[],1)], ode_opts);

deltaXf = Xf - Xlarge_loop_hist(end,1:6)';
STM_t = reshape(Xlarge_loop_hist(end,7:end),12,12);

% Solve for lambda0 (at tnow)
lambda0 = STM_t(1:6,7:12)\(deltaXf - STM_t(1:6,1:6)*deltaX0);

%% Test initial guess for lambda0

m0 = m_sc;
Tmax = max_thrust_mag;
amax = max_thrust_mag/m_sc;

params_accel = struct('m_norm',MU, 'c',exh_vel, 'mu',mu_EM, 'amax',amax, 'Tmax',max_thrust_mag,'L_EM',L_EM, 'T_EM',T_EM, 'force_norm',FU, 'vel_norm',VU, 'accel_norm', AU);

chi0 = [X0; lambda0];
fprintf("Controlled trajectory with guess for lambda0\n")
tic
[t_opt, X_opt] = ode113(@(t,X) state_costate_reduced_dynamics(t,X,m0,@opt_control_min_energy_reduced,params_accel), [t0, t0+cf_guess*tfmin0], chi0, ode_opts);
[t_no, X_no] = ode113(@(t,X) CR3BP(t,X,mu_EM), [t0, t0+cf_guess*tfmin0], X0, ode_opts);
toc

%% Plot controlled traj guess w/ uncontrolled

figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot(1-mu_EM, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Moon'); hold on % Moon
plot(x_L1, 0, 'dk', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
plot(x_L2, 0, 'dk' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
plot(X_opt(1,1), X_opt(1,2), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
plot(Xf(1), Xf(2), 'ok', 'markerfacecolor', 'g', 'DisplayName', 'Target'); hold on
plot(X_opt(:,1), X_opt(:,2), 'r-','DisplayName', 'Initial Guess Controlled Trajectory'); hold on
plot(X_no(:,1), X_no(:,2), 'b-','DisplayName','Uncontrolled Trajectory'); hold off
title('Controlled Trajectory using Linearized Initial Guess','FontSize',14)
xlabel('X')
ylabel('Y')
axis('equal')
grid on;
legend('FontSize',12);
hold off
%% Correct into nonlinear model (without mass)

cfr = 1;
%cfr = 0.9305;
tfmin = tf-t0;
tf0 = tf-t0;
lambda_i = lambda0;
free_vars = lambda_i;

shooting_func_min_energy_reduced_simple = @(free_vars) shooting_func_min_energy_reduced(0,free_vars,X0,m0,Xf,cfr*tfmin,params_accel);

min_energy_reduced_result = shooter(free_vars,shooting_func_min_energy_reduced_simple,'tol',1e-10,'max_iter',5000);


%% Plot results from above

free_vars_converged_reduced = min_energy_reduced_result.free_vars;
lambda0_masspre = free_vars_converged_reduced;
err_mag_hist = min_energy_reduced_result.err_hist;
chi0 = [X0; free_vars_converged_reduced];
fprintf("Controlled trajectory with corrected lambda0\n")
tic
[t_opt_pass1, X_opt_pass1] = ode113(@(t,X) state_costate_reduced_dynamics(t,X,m0,@opt_control_min_energy_reduced,params_accel), [t0, t0+cfr*tfmin0], chi0, ode_opts);

figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot(1-mu_EM, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Moon'); hold on % Moon
plot(x_L1, 0, 'dk', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
plot(x_L2, 0, 'dk' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
plot(X_opt_pass1(1,1), X_opt_pass1(1,2), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
plot(Xf(1), Xf(2), 'ok', 'markerfacecolor', 'g', 'DisplayName', 'Target'); hold on
plot(X_opt_pass1(:,1), X_opt_pass1(:,2), 'r-','DisplayName', 'Converged Controlled Trajectory (Reduced Model)'); hold on
plot(X_no(:,1), X_no(:,2), 'b-','DisplayName','Uncontrolled Trajectory'); hold off
title('Controlled Trajectory using Corrected Initial Costate (Reduced Model)','FontSize',14)
xlabel('X')
ylabel('Y')
axis('equal')
grid on;
legend('FontSize',12);
hold off

figure
addToolbarExplorationButtons(gcf)
semilogy(1:1:length(err_mag_hist), err_mag_hist, '-x','DisplayName', 'Residual (Error) Magnitude');
xlim([1, inf])
xlabel('Iteration','FontSize',14)
ylabel('Residual Magnitude','FontSize',14)
title('Residual Magnitude vs Iteration Number','FontSize',14)
grid on
legend('FontSize',12)

%%%
% chi0 = [X0; free_vars];
% tic
% [t_opt, X_opt] = ode113(@(t,X) state_costate_reduced_dynamics(t,X,m0,@opt_control_min_energy_reduced,params_accel), [t0, t0+cf*tfmin], chi0, ode_opts);
% [t_no, X_no] = ode113(@(t,X) CR3BP(t,X,mu_EM), [t0, tf], X0, ode_opts);
% toc
% 
% figure
% addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
% plot(1-mu_EM, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Moon'); hold on % Moon
% plot(x_L1, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
% plot(x_L2, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
% plot(X_opt(1,1), X_opt(1,2), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
% plot(Xf(1), Xf(2), 'ok', 'markerfacecolor', 'g', 'DisplayName', 'Target'); hold on
% plot(X_opt(:,1), X_opt(:,2), 'r-','DisplayName', 'Initial Guess Controlled Trajectory'); hold on
% plot(X_no(:,1), X_no(:,2), 'b-','DisplayName','Uncontrolled Trajectory'); hold off
% title('Reference Trajectory')
% xlabel('X')
% ylabel('Y')
% axis('equal')
% grid on;
% legend();
% hold off

%% with mass

cff = 1;
%cff = 0.91; 
tfmin = tf-t0;

%Tmax = 0.32; % newtons
%exh_vel = 3000; % m/s
m0 = m_sc;

tf0 = tf-t0;
lambda_i = [lambda0_masspre; 1e-6];
free_vars = lambda_i;

shooting_func_min_energy_simple = @(free_vars) shooting_func_min_energy(0,free_vars,X0,m0,Xf,cff*tfmin,params_accel);

min_energy_full_result = shooter(free_vars,shooting_func_min_energy_simple,'tol',1e-10,'max_iter',5000);


%%

free_vars_converged_full = min_energy_full_result.free_vars;
%% Plot Results

chi0 = [X0; m0; free_vars_converged_full];
tic
[t_opt_pass2, X_opt_pass2] = ode113(@(t,X) state_costate_full_dynamics(t,X,@opt_control_min_energy,params_accel), [t0, t0+cff*tfmin], chi0, ode_opts);
toc

% Compute control history
final_control_results = get_control_hist_min_energy_full(X_opt_pass2,params_accel.c);

% figure
% addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
% plot3(1-mu_EM, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Moon'); hold on % Moon
% plot3(x_L1, 0, 0, 'dk', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
% plot3(x_L2, 0, 0, 'dk' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
% plot3(X0(1), X0(2), X0(3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
% plot3(Xf(1), Xf(2), Xf(3), 'ok', 'markerfacecolor', 'g', 'DisplayName', 'Target'); hold on
% plot3(X_opt_pass2(:,1), X_opt_pass2(:,2), X_opt_pass2(:,3), 'r-','DisplayName', 'Controlled Trajectory w Variable Mass'); hold off
% title('Controlled Trajectory Using Corrected Initial Costate (Full Model with Mass)')
% xlabel('X')
% ylabel('Y')
% grid on;
% % 3D axis equal
% h = get(gca,'DataAspectRatio');
% if h(3)==1
%       set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
% else
%       set(gca,'DataAspectRatio',[1 1 h(3)])
% end
% legend();

% Plot test (2D)

figure
addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
plot(1-mu_EM, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Moon'); hold on % Moon
plot(x_L1, 0, 'dk', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
plot(x_L2, 0, 'dk' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
plot(X0(1), X0(2), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
plot(Xf(1), Xf(2), 'ok', 'markerfacecolor', 'g', 'DisplayName', 'Target'); hold on
plot(X_opt_pass2(:,1), X_opt_pass2(:,2), 'r-','DisplayName', 'Controlled Trajectory w Variable Mass'); hold on
plot(X_no(:,1), X_no(:,2), 'b-','DisplayName','Uncontrolled Trajectory'); hold off
title('Controlled Trajectory Using Corrected Initial Costate (Full Model with Mass)','FontSize',14)
xlabel('X')
ylabel('Y')
grid on;
axis('equal')
legend('FontSize',12);

% Plot control vector and thrust factor history
figure
sgtitle('Control and Mass History','FontSize',14)
subplot(1,2,1)
plot(t_opt_pass2-t0, FU*1e6*Tmax.*final_control_results.u_vec_hist(1,:), 'DisplayName', 'Thrust Vector (x component) [mN]'); hold on
plot(t_opt_pass2-t0, FU*1e6*Tmax.*final_control_results.u_vec_hist(2,:), 'DisplayName', 'Thrust Vector (y component) [mN]'); hold off
title('Control Vector History','FontSize',12)
xlabel('Time [non-dimensional units]','FontSize',12)
ylabel('Control Vector Magnitude [milli-Newtons]','FontSize',12)
grid on
legend('FontSize',12)

subplot(1,2,2)
yyaxis left
plot(t_opt_pass2-t0, final_control_results.u_hist, 'DisplayName', 'Thrust Factor');
ylabel('Thrust Factor','FontSize',12)
yyaxis right
plot(t_opt_pass2-t0, MU*X_opt_pass2(:,7), 'DisplayName', 'Spacecraft Mass');
ylabel('Mass [kg]','FontSize',12)
xlabel('Time [non-dimensional units]','FontSize',12)
title('Thrust Factor and Spacecract Mass vs Time','FontSize',12)
legend('FontSize',12)
grid on

% Plot residual history
figure
addToolbarExplorationButtons(gcf)
semilogy(1:1:length(err_mag_hist), err_mag_hist, '-x','DisplayName', 'Residual (Error) Magnitude');
xlim([1, inf])
xlabel('Iteration','FontSize',14)
ylabel('Residual Magnitude','FontSize',14)
title('Residual Magnitude vs Iteration Number','FontSize',14)
grid on
legend('FontSize',12)

max_thrust_used = FU*1e6*max(Tmax.*final_control_results.u_hist); % [mN]
total_deltav = trapz(TU.*t_opt_pass2,AU*Tmax/m_sc.*final_control_results.u_hist); % [m/s]
fuel_used = MU*(m_sc - X_opt_pass2(end,7)); % [g]
final_mass = MU*X_opt_pass2(end,7); % [kg]
fprintf("Max Thrust Used: %d mN\n",max_thrust_used);
fprintf("Total Delta V: %d m/s\n",total_deltav);
fprintf("Fuel Mass Used: %d g\n",fuel_used);
fprintf("Spacecraft Final mass: %d kg\n",final_mass);
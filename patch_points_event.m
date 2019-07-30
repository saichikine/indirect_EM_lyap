function [zero_val,isterminal,direction] = patch_points_event(t,X, mu, gamma, L_point_num, x0_L_center, s, beta_3, h, eps)

X_approx = linearized_planar_collinear_dynamics(t, mu, gamma, L_point_num, x0_L_center, s, beta_3);

%zero_val = norm(X) - h*norm(X_approx); % The value that we want to be zero
%zero_val = norm(X(1:6)) - norm(X_approx) - eps;
zero_val = norm(X(1:6)-X_approx) - eps;
isterminal = 1;  % Halt integration 
direction = 1;   % Increasing
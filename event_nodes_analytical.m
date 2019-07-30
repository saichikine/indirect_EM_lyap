function [zero_val,isterminal,direction] = event_nodes_analytical(t,X, analytical_sln_handle, eps)

X_approx = analytical_sln_handle(t);

zero_val = norm(X(1:6)-X_approx) - eps;
isterminal = 1;  % Halt integration 
direction = 1;   % Increasing
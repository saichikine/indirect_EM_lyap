function Xlargedot = CR3BP_costate_STM_dynamics(t,Xlarge,B,R,mu)
    % Integrates augmented STM with phidot = (A, -B*R^-1B^T; 0, -A^T)*phi
    % Need to integrate CR3BP simultaneously since A depends on state
    
    %% CR3BP dynamics
    % Unpack CR3BP state
    X = Xlarge(1:6);
    
    % Positions
    x = X(1);
    y = X(2);
    z = X(3);
    % Velocities
    xdot = X(4);
    ydot = X(5);
    zdot = X(6);
    
    % distances from particle to primaries
    r1 = sqrt((x+mu)^2 + y^2 + z^2);
    r2 = sqrt((x-1+mu)^2 + y^2 + z^2);
    
    % Accelerations
    xddot = 2*ydot + x -(1-mu)*((x+mu)/(r1^3)) - mu*(x-1+mu)/(r2^3);
    yddot = -2*xdot + y - (1-mu)*(y/(r1^3)) - mu*(y)/(r2^3);
    zddot = -(1-mu)*(z)/(r1^3) - mu*(z)/(r2^3);
    
    Xdot = [xdot; ydot; zdot; xddot; yddot; zddot];
    
    %% Augmented STM Dynamics
    
    STM = reshape(Xlarge(7:end),12,12);
    
    A = Jacobian(Xlarge,mu);
    dfdX = [A, -B/(R)*B'; zeros(6,6), -A'];
    %dfdX = [A, [zeros(3), zeros(3); zeros(3), eye(3)]; zeros(6), -A'];
    
    STMdot = dfdX*STM;
    
    %% Repack large state
    
    Xlargedot = [Xdot; reshape(STMdot,[],1)];
    
    if any(isnan(Xlargedot))
        warning("NaNs in dynamics.\n")
    end
    
end
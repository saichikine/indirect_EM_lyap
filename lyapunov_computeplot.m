function lyapunov_orbit = lyapunov_computeplot(mu, Ax_km, L_km, L_point_string, varargin)

    %apdgsh
    
    %% Input handling
    
    if nargin>7
        error('Maximum of 7 arguments.')
    end
    
    if ~isstring(L_point_string)
        error('Fourth argument must be a string.')
    end
    
    if (~strcmpi(L_point_string,"l1") && ~strcmpi(L_point_string,"l2"))
        error('Fourth argument must be either "L1" or "L2"')
    end
    
    plot_flag = 1;
    if ~isempty(varargin)>0
        plot_flag = varargin{1};
        if length(varargin)>1
            ydot0 = varargin{2};
            orbit_period = varargin{3};
        end
    end
    
    %% L1 or L2?
    
    L_points = lagrangePoints(mu);
    gammas = CR3BP_L1_L2_gammas(mu);
    xL1 = L_points(1,1);
    xL2 = L_points(1,2);
    
    if strcmpi(L_point_string,'l1')
        L_point = L_points(:,1);
        gamma = gammas.one;
    elseif strcmpi(L_point_string,'l2')
        L_point = L_points(:,2);
        gamma = gammas.two;
    else
        error('You should not be here.')
    end
    
    x0_L_center = Ax_km/L_km;
    %% Get things
    
    if nargin<=5
        lyap_approx = lyapunov_approx_ICs(Ax_km, L_km, mu, L_point_string, L_point, gamma);
        X0 = lyap_approx.ICs;
        orbit_period = lyap_approx.orbit_period;
    else
        if strcmpi(L_point_string,"L1")
            x0 = (1-mu) + gamma*(-x0_L_center - 1);
        elseif strcmpi(L_point_string,"L2")
            x0 = (1-mu) + gamma*(x0_L_center + 1);
        else
            error('This is bad.')
        end
        X0 = [x0; 0; 0; 0; ydot0; 0];
    end
 
    %% Simulate initial guess
    
    ode_opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
    
    %fprintf('Simulating...')
    %tic
    [t_guess, X_guess_hist] = ode113(@(t,X) CR3BP(t,X,mu), [0 orbit_period], X0, ode_opts);
    %fprintf('done.\n')
    %toc
    
    if plot_flag
        
        figure
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar

        subplot(2,1,1)
        plot3(X_guess_hist(:,1), X_guess_hist(:,2), X_guess_hist(:,3), 'DisplayName', 'Approximate Lyapunov Orbit'); hold on
        plot3(X_guess_hist(1,1), X_guess_hist(1,2), X_guess_hist(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
        plot3(-mu, 0, 0, 'ok', 'markerfacecolor', 'k', 'markersize', 10, 'DisplayName', 'Larger Primary'); hold on % Larger primary
        plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 5, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
        plot3(xL1, 0, 0, 'ok', 'markerfacecolor', 'r', 'markersize', 2, 'DisplayName', 'L1 Point'); hold on % L1 location
        plot3(xL2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'markersize', 2, 'DisplayName', 'L2 Point'); hold on % L2 location
        title('Initial (Approximate) Lyapunov Orbit');
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        grid on;
        legend();
        hold off

        % Close up around Lagrange point
        subplot(2,1,2)
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
        plot3(X_guess_hist(1,1), X_guess_hist(1,2), X_guess_hist(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
        plot3(xL1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
        plot3(xL2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
        plot3(X_guess_hist(:,1), X_guess_hist(:,2), X_guess_hist(:,3), 'DisplayName', 'Approximate Lyapunov Orbit'); hold on
        title('Initial (Approximate) Lyapunov Orbit (Close Up)')
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        grid on;
        legend();
        hold off
        
        h = get(gca,'DataAspectRatio');
        if h(3)==1
              set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))]);
        else
              set(gca,'DataAspectRatio',[1 1 h(3)]);
        end

    end
    
    %% Differential Correction
    
    ode_opts_XZ_crossing = odeset('Events',@Findyzero,'RelTol',1e-13,'AbsTol',1e-13);
    
    STM_0 = reshape(eye(6), 36, []);
    X0 = [X0; STM_0];
    
    delta_vec = [100; 100];
    delta = 100;
    counter = 0;
    verbose = 0;
    
    %fprintf('Targeting periodic Lyapunov orbit...\n')
    %tic
    while (norm(delta) > 1e-13)
        counter = counter+1;
        [t_half,X_half] = ode113(@(t,X) CR3BP(t,X,mu), [0 inf], X0, ode_opts_XZ_crossing);

        % State when y=0 (half orbit period)
        dX = X_half(end,1:6);

        % STM at half orbit period
        dPhi = reshape(X_half(end,7:end),6,[]);

        %% dX/dt (X is state)
        x_dot = dX(4);
        y_dot = dX(5);
        z_dot = dX(6);

        r1 = sqrt((dX(1)+mu)^2 + dX(2)^2 + dX(3)^2); %S/C distance to Sun
        r2 = sqrt((dX(1)-1+mu)^2 + dX(2)^2 + dX(3)^2); %S/C distance to Earth

        % Accelerations
        x_ddot = 2*y_dot + dX(1) - (1 - mu)*((dX(1) + mu)/(r1^3)) - mu*(dX(1) - 1 + mu)/(r2^3);
        y_ddot = -2*x_dot + dX(2) - (1 - mu)*(dX(2)/(r1^3)) - mu*(dX(2))/(r2^3);
        z_ddot = -(1 - mu)*(dX(3))/(r1^3) - mu*(dX(3))/(r2^3); 

        % Derivative of new state
        dXdt = [x_dot y_dot z_dot x_ddot y_ddot z_ddot];

        %Update matrix to correct deltaVec
        %Only changes z0 and doty0, holds x0 constant
%         update_mat = [dPhi(4,3)-dPhi(2,3)*(x_ddot/y_dot), dPhi(4,5)-dPhi(2,5)*(x_ddot/y_dot);...
%             dPhi(6,3)-dPhi(2,3)*(z_ddot/y_dot), dPhi(6,5)-dPhi(2,5)*(z_ddot/y_dot)];
% 
%         delta_vec = inv(update_mat)*[-X_half(end,4); -X_half(end,6)];
% 
%         if verbose
%             norm(delta_vec)
%             fprintf('Iteration counter: %d\n', counter)
%         end
        
        delta = (-X_half(end,4))/(dPhi(4,5)-dPhi(2,5)*y_dot/x_ddot);

        %delta_X = [delta_vec(1), 0, 0, 0, delta_vec(2), 0]';
        delta_X = [0, 0, 0, 0, delta, 0]';

        X0(1:6) = X0(1:6) + delta_X;
    end
    %fprintf('done.\n')
    %toc
    
    %% New Corrected Orbit
    
    %Integrate orbit with solved initial conditions
    %fprintf('Corrected periodic orbit.\n')
    %fprintf('Simulating...')
    %tic
    [t_corrected, X_corrected_hist] = ode113(@(t,X) CR3BP(t,X,mu), [0 t_half(end)*2], X0, ode_opts);
    %fprintf('done.\n')
    %toc
    
    lyapunov_orbit{1} = X0(1:6);
    lyapunov_orbit{2} = t_half(end)*2;
    
    %% Plot Corrected Orbit
    
    if plot_flag
        figure
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar

        subplot(2,1,1)
        plot3(X_corrected_hist(:,1), X_corrected_hist(:,2), X_corrected_hist(:,3), 'DisplayName', 'Corrected Lyapunov Orbit'); hold on
        plot3(X_corrected_hist(1,1), X_corrected_hist(1,2), X_corrected_hist(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
        plot3(-mu, 0, 0, 'ok', 'markerfacecolor', 'k', 'markersize', 10, 'DisplayName', 'Larger Primary'); hold on % Larger primary
        plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 5, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
        plot3(xL1, 0, 0, 'ok', 'markerfacecolor', 'r', 'markersize', 2, 'DisplayName', 'L1 Point'); hold on % L1 location
        plot3(xL2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'markersize', 2, 'DisplayName', 'L2 Point'); hold on % L2 location
        title('Corrected Periodic Lyapunov Orbit');
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        grid on;
        legend();
        hold off
        
        h = get(gca,'DataAspectRatio');
        if h(3)==1
              set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))]);
        else
              set(gca,'DataAspectRatio',[1 1 h(3)]);
        end

        % Close up around Lagrange point
        subplot(2,1,2)
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
        plot3(X_corrected_hist(1,1), X_corrected_hist(1,2), X_corrected_hist(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
        plot3(xL1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
        plot3(xL2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
        plot3(X_corrected_hist(:,1), X_corrected_hist(:,2), X_corrected_hist(:,3), 'DisplayName', 'Corrected Lyapunov Orbit'); hold on
        title('Corrected Periodic Lyapunov Orbit (Close Up)')
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        grid on;
        legend();
        hold off
        
        h = get(gca,'DataAspectRatio');
        if h(3)==1
              set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))]);
        else
              set(gca,'DataAspectRatio',[1 1 h(3)]);
        end

    end
end
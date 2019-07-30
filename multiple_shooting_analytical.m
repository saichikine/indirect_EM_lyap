function results = multiple_shooting_analytical(initial_guess_ICs, initial_guess_period, mu, analytical_sln_handle, node_event_handle, varargin)
    
    %% Input handling
    
    default_tol = 1e-14;
    default_max_iter = 30;
    default_plot_flag = 0;
    default_verbose_flag = 1;
    
    p = inputParser;
    valid_scalar_pos_num = @(x) isnumeric(x) && isscalar(x) && (x>0);
    valid_ICs = @(x) all(size(initial_guess_ICs==[6 1]));
    valid_f_handle = @(x) isa(x, 'function_handle');
    valid_bool = @(x) x==1 || x==0;
    valid_int = @(x) mod(x,1)==0;
    
    addRequired(p,'initial_guess_ICs',valid_ICs);
    addRequired(p,'initial_guess_period',valid_scalar_pos_num);
    addRequired(p,'mu',valid_scalar_pos_num);
    addRequired(p,'analytical_sln_handle',valid_f_handle);
    addRequired(p,'node_event_handle',valid_f_handle);
    
    addParameter(p,'plot_flag',default_plot_flag,valid_bool);
    addParameter(p,'verbose_flag',default_verbose_flag,valid_bool);
    addParameter(p,'tol',default_tol,valid_scalar_pos_num);
    addParameter(p,'max_iter',default_max_iter,valid_int);
    
    parse(p,initial_guess_ICs,initial_guess_period,mu,analytical_sln_handle,node_event_handle,varargin{:});
    
    initial_guess_ICs = p.Results.initial_guess_ICs;
    initial_guess_period = p.Results.initial_guess_period;
    mu = p.Results.mu;
    analytical_sln_handle = p.Results.analytical_sln_handle;
    node_event_handle = p.Results.node_event_handle;
    plot_flag = p.Results.plot_flag;
    verbose_flag = p.Results.verbose_flag;
    tol = p.Results.tol;
    max_iter = p.Results.max_iter;
    
    orbit_period = initial_guess_period;
    
    %% Find nodes
    
    t_now = 0;
    X_hist_tot = [];
    t_hist_tot = [];
    patch_point_index = 0;
    arc_initial_states = [];
    arc_final_states = [];
    integration_times = [];
    time_derivs = [];

    ode_opts_patch_points = odeset('RelTol',1e-13,'AbsTol',1e-20,'Events',node_event_handle);
    
    if verbose_flag
        tic
    end
    while t_now < orbit_period

        X0 = [analytical_sln_handle(t_now); reshape(eye(6),[],1)];
        arc_initial_states = [arc_initial_states, X0(1:6)];

        [t_hist,X_hist,te,Xe,ie] = ode113(@(t,X) CR3BP(t,X,mu), [t_now, orbit_period], X0, ode_opts_patch_points);

        patch_point_index = patch_point_index+1;

        arc_final_states = [arc_final_states, reshape(X_hist(end,1:6),6,[])];
        integration_times = [integration_times; t_hist(end)-t_now]; % store duration of arc
        STMs{patch_point_index} = reshape(X_hist(end,7:end),6,6);
        time_derivs = [time_derivs, CR3BP(t_hist(end),X_hist(end,1:6),mu)];

        X_hist_tot = [X_hist_tot; X_hist];
        t_hist_tot = [t_hist_tot; t_hist];
        
        % Keep track of how many integration points each arc is (for last
        % arc check)
        arc_size = length(t_hist);
        
        t_now_prev = t_now;
        t_now = t_hist(end);
    end
    if verbose_flag
        toc
    end
    
    % This makes sure last arc ends at y=0 plane
%     if X_hist(end,2) ~= 0
%         X0 = [analytical_sln_handle(t_now); reshape(eye(6),[],1)];
%         ode_opts_yzero = odeset('RelTol',1e-13,'AbsTol',1e-20,'Events',@Findyzero);
%         [t_hist,X_hist,te,Xe,ie] = ode113(@(t,X) CR3BP(t,X,mu), [t_now_prev, orbit_period], X0, ode_opts_yzero);
%         % Delete last arc time history
%         t_hist_tot(end-(arc_size-1):end) = [];
%         X_hist_tot(end-(arc_size-1):end,:) = [];
%         arc_final_states(:,end) = [];
%         integration_times(end) = [];
%         time_derivs(:,end) = [];
%         % Rewrite last arc time history
%         t_hist_tot = [t_hist_tot; t_hist];
%         X_hist_tot = [X_hist_tot; X_hist];
%         arc_final_states = [arc_final_states, reshape(X_hist(end,:),6,[])];
%         integration_times = [integration_times; t_hist(end) - t_now_prev];
%         time_derivs = [time_derivs, CR3BP(t_hist(end),X_hist(end,1:6),mu)];
%     end
    %% Total linearized dynamics

    X_lin_hist = zeros(6,length(X_hist_tot));

    for i = 1:length(t_hist_tot)
        X_lin_hist(:,i) = analytical_sln_handle(t_hist_tot(i));
    end

    X_lin_hist = X_lin_hist';
    %% Build dg/dX and g(X)

    % Build initial big X
    chi = [reshape(arc_initial_states,[],1); integration_times];

    N = patch_point_index;
    l = 6; % number of states (6)

    %% Break out and integrate forward again

    constraint_norm = 100;
    tol = 1e-14;
    counter = 0;
    max_iter = 30;

    while constraint_norm > tol

        times = chi(end-(N-1):end);
        arc_initial_states = reshape(chi(1:end-N),6,[]);
        % Manually enforce z,zdot=0 constraints for Lyapunov case
        arc_initial_states(3,:) = 0;
        arc_initial_states(6,:) = 0;
        arc_final_states = NaN(size(arc_initial_states));

        time_derivs = NaN(6,N);

        X_hists = {};
        t_hist_total = [];

        parfor i = 1:N
            X0 = [arc_initial_states(:,i); reshape(eye(6),[],1)];

            [t_hist,X_hist,~,~,~] = ode113(@(t,X) CR3BP(t,X,mu), [0 times(i)], X0, ode_opts_patch_points);

            X_hists{i} = X_hist;

            arc_final_states(:,i) = X_hist(end,1:6)';
            STMs{i} = reshape(X_hist(end,7:end),6,6);

            time_derivs(:,i) = CR3BP(t_hist(end),X_hist(end,1:6),mu);

        end
        
        % Patch together X_hist_total
        X_hist_total = [];
        for i = 1:N
            X_hist_total = [X_hist_total; X_hists{i}];
        end
        
        % y=0 and xdot=0 constraints
        dgdX = zeros(N*l+2,N*(l+1));
        % second to last row of jacobian
        dgdX(end-1,2) = 1;
        % last row of jacobian
        dgdX(end,4) = 1;

        for i = 1:N
            dgdX((i-1)*l+1:i*l, (i-1)*l+1:i*l) = STMs{i}; % STMs along main diagonal
            if i<N
                dgdX((i-1)*l+1:i*l, i*l+1:(i+1)*l) = -eye(l);
            elseif i==N
                dgdX((i-1)*l+1:i*l, 1:l) = -eye(l);
            end
            dgdX((i-1)*l+1:i*l, N*l+i) = time_derivs(:,i);
        end

        gX = [];

        for i = 1:N-1
            gX = [gX; arc_final_states(:,i) - arc_initial_states(:,i+1)];
        end

        gX = [gX; arc_final_states(:,N) - arc_initial_states(:,1)];
        gX = [gX; arc_initial_states(2,1)]; % initial y = 0
        gX = [gX; arc_initial_states(4,1)]; % initial xdot = 0

        % Compute and correct

        delta_chi = -dgdX'/(dgdX*dgdX')*1*gX;
        %delta_chi = -pinv(dgdX)*gX;
        constraint_norm_prev = constraint_norm;
        constraint_norm = norm(gX);
        
        norm_ratio = abs(constraint_norm)/constraint_norm_prev^2;

        chi = chi + delta_chi;
        
        counter = counter+1;
        
        if verbose_flag
            fprintf('Iteration count: %i\n',counter);
            fprintf('Ratio: %d\n',norm_ratio);
        end
        if counter > max_iter
            break
        end
    end
    
    if verbose_flag
        fprintf('Final constraint norm: %i\n',constraint_norm);
    end
    total_time = sum(times);
    orbit_IC = reshape(X_hist_total(1,1:6),[],1);
    C = jacobi_constant(orbit_IC,mu);
    %% Plot comparison (starting off)
    
    if plot_flag
        
        L_points = lagrangePoints(mu);
        x_L1 = L_points(1,1);
        x_L2 = L_points(1,2);

        figure
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
        plot3(X_hist_tot(1,1), X_hist_tot(1,2), X_hist_tot(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
        plot3(x_L1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
        plot3(x_L2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
        plot3(X_hist_tot(:,1), X_hist_tot(:,2), X_hist_tot(:,3), 'r-','DisplayName', 'Approximate Lyapunov Orbit'); hold on
        plot3(X_lin_hist(:,1), X_lin_hist(:,2), X_lin_hist(:,3),'b-','DisplayName','Linearized Dynamics'); hold on
        title('Multiple Shooting Initial Guess')
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        grid on;
        legend();
        hold off

        %% Plot comparison (after done)

        figure
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
        plot3(X_hist_total(1,1), X_hist_total(1,2), X_hist_total(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
        plot3(x_L1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
        plot3(x_L2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
        plot3(X_hist_total(:,1), X_hist_total(:,2), X_hist_total(:,3), 'r-','DisplayName', 'Converged Lyapunov Orbit'); hold on
        plot3(X_lin_hist(:,1), X_lin_hist(:,2), X_lin_hist(:,3),'b-','DisplayName','Linearized Dynamics'); hold on
        title('Multiple Shooting Converged')
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        grid on;
        legend();
        hold off
    end
    
    %% Save results
    
    results = struct('X0',orbit_IC,'period',total_time,'jacobi_constant',C,'free_vars',chi);
end
function results = multiple_shooting_hetero(initial_nodes, initial_times, mu, varargin)
    
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
    
    addRequired(p,'initial_nodes',@isnumeric);
    addRequired(p,'initial_times',@isnumeric);
    addRequired(p,'mu',valid_scalar_pos_num);
    
    addParameter(p,'plot_flag',default_plot_flag,valid_bool);
    addParameter(p,'verbose_flag',default_verbose_flag,valid_bool);
    addParameter(p,'tol',default_tol,valid_scalar_pos_num);
    addParameter(p,'max_iter',default_max_iter,valid_int);
    
    parse(p,initial_nodes,initial_times,mu,varargin{:});
    
    initial_nodes = p.Results.initial_nodes;
    initial_times = p.Results.initial_times;
    mu = p.Results.mu;
    plot_flag = p.Results.plot_flag;
    verbose_flag = p.Results.verbose_flag;
    tol = p.Results.tol;
    max_iter = p.Results.max_iter;
 
    %% Setup

    % Build initial big X
    chi = [reshape(initial_nodes,[],1); initial_times];

    N = size(initial_nodes,2); % number of nodes
    l = 6; % number of states (6)

    %% Break out and integrate forward again

    constraint_norm = 100;
    counter = 0;
    
    ode_opts = odeset('RelTol',5e-14,'AbsTol',1e-20);

    while constraint_norm > tol

        times = chi(end-(N-1)+1:end);
        %arc_initial_states = initial_nodes;
        arc_initial_states = reshape(chi(1:N*l),6,[]);
        arc_final_states = NaN(size(arc_initial_states));

        time_derivs = NaN(6,N);

        X_hists = {};
        t_hist_total = [];

        parfor i = 1:N-1
            X0 = [arc_initial_states(:,i); reshape(eye(6),[],1)];
            
            if times(i)==0
                X_hist = X0';
                t_hist = [0];
            else
                [t_hist,X_hist] = ode113(@(t,X) CR3BP(t,X,mu), [0 times(i)], X0, ode_opts);
            end

            X_hists{i} = X_hist;

            arc_final_states(:,i) = X_hist(end,1:6)';
            STMs{i} = reshape(X_hist(end,7:end),6,6);

            time_derivs(:,i) = CR3BP(t_hist(end),X_hist(end,1:6),mu);

        end
        
        % Patch together X_hist_total
        X_hist_total = [];
        for i = 1:N-1
            X_hist_total = [X_hist_total; X_hists{i}];
        end
        
        % Build dgdX
        dgdX = zeros((N-1)*l,N*(l+1)-1);

        for i = 1:(N-1)
            dgdX((i-1)*l+1:i*l, (i-1)*l+1:i*l) = STMs{i}; % STMs along main diagonal
            if i<=(N-1)
                dgdX((i-1)*l+1:i*l, i*l+1:(i+1)*l) = -eye(l);
            end
            if i<=(N-1)
                dgdX((i-1)*l+1:i*l, N*l+i) = time_derivs(:,i);
            end
        end

        gX = [];

        for i = 1:N-1
            gX = [gX; arc_final_states(:,i) - arc_initial_states(:,i+1)];
        end

        % Compute and correct

        delta_chi = -dgdX'/(dgdX*dgdX')*1*gX;
        %delta_chi = -pinv(dgdX)*gX;
        constraint_norm = norm(gX);

        chi = chi + delta_chi;
        
        counter = counter+1;
        
        if verbose_flag
            fprintf('Iteration count: %i\n',counter);
            fprintf('Constraint norm: %d\n',constraint_norm);
        end
        if counter > max_iter
            break
        end
    end

    total_time = sum(times);
    orbit_IC = reshape(X_hist_total(1,1:6),[],1);
    C = jacobi_constant(orbit_IC,mu);
    %% Plot comparison (starting off)
    
    if plot_flag
        
        L_points = lagrangePoints(mu);
        x_L1 = L_points(1,1);
        x_L2 = L_points(1,2);

        %% Plot comparison (after done)
        
        % 3D Plot
%         figure
%         addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
%         plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
%         plot3(X_hist_total(1,1), X_hist_total(1,2), X_hist_total(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
%         plot3(x_L1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
%         plot3(x_L2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
%         arcs = {};
%         for i = 1:N-1
%             arcs{i} = scatter3(X_hists{i}(:,1), X_hists{i}(:,2), X_hists{i}(:,3),'.'); hold on
%         end
%         title('Multiple Shooting Converged')
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
%         grid on;
%         legend();
%         hold off
        
        % 2D Plot
        figure
        addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
        plot(1-mu, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
        plot(X_hist_total(1,1), X_hist_total(1,2), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
        plot(x_L1, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
        plot(x_L2, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
        arcs = {};
        for i = 1:N-1
            arcs{i} = scatter(X_hists{i}(:,1), X_hists{i}(:,2), '.'); hold on
        end
        title('Multiple Shooting Converged')
        xlabel('X')
        ylabel('Y')
        axis('equal')
        grid on;
        legend();
        hold off

%         figure
%         addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
%         plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'm', 'markersize', 8, 'DisplayName', 'Smaller Primary'); hold on % Smaller primary
%         plot3(X_hist_total(1,1), X_hist_total(1,2), X_hist_total(1,3), 'ok', 'markerfacecolor', 'y', 'DisplayName', 'Initial Point'); hold on
%         plot3(x_L1, 0, 0, 'ok', 'markerfacecolor', 'r', 'DisplayName', 'L1 Point'); hold on % L1 location
%         plot3(x_L2, 0, 0, 'ok' , 'markerfacecolor', 'b', 'DisplayName', 'L2 Point'); hold on % L2 location
%         plot3(X_hist_total(:,1), X_hist_total(:,2), X_hist_total(:,3), 'r-','DisplayName', 'Converged Connection'); hold on
%         title('Multiple Shooting Converged')
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
%         grid on;
%         legend();
%         hold off
    end
    
    %% Save results
    
    results = struct('X0',orbit_IC,'total_time',total_time,'jacobi_constant',C,'free_vars',chi);
end
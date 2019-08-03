function jacobian = numerical_jacobian(func,X,varargin)

    % Computes numerical jacobian of 'func' w.r.t. X
    
    %% Input Handling
    default_rel_eps = 1e-7;
    default_abs_eps = 1e-8;
    
    valid_f_handle = @(x) isa(x, 'function_handle');
    valid_scalar_pos_num = @(x) isnumeric(x) && isscalar(x) && (x>0);
    
    p = inputParser;
    
    addRequired(p,'func',valid_f_handle);
    addRequired(p,'X',@isvector);
    
    addParameter(p,'RelEps',default_rel_eps,valid_scalar_pos_num);
    addParameter(p,'AbsEps',default_abs_eps,valid_scalar_pos_num);
    
    parse(p,func,X,varargin{:});
    
    func = p.Results.func;
    X = p.Results.X;
    
    rel_eps = p.Results.RelEps;
    abs_eps = p.Results.AbsEps;
    
    X = reshape(X,[],1);
    if ~isvector(func(X))
        error("Function must be a scalar or vector.\n")
    end
    
    N = length(X); % Length of domain vector
    M = length(func(X)); % Length of codomain vector
    
    %% Compute Numerical Jacobian
    jacobian = NaN(M,N);
    func_output_nominal = func(X); % precompute for speed
    
    parfor i = 1:N
        
        Xi_perturbation = X(i)*rel_eps + abs_eps;
        X_perturbed = X;
        X_perturbed(i) = X(i) + Xi_perturbation;
        
        jacobian(:,i) = (func(X_perturbed) - func_output_nominal)./Xi_perturbation;
    end
end
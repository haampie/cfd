function solution = central(params, method)
    dx = params.L / method.N;
    
    % eta = params.u * method.dt / dx;
    % d = 2 * params.k * method.dt / dx ^ 2;
    % assert(d > 0);
    % assert(d < 1);
    % assert(eta^2 < d);

    ssub_diag = zeros(method.N - 1, 1);
    sub_diag = zeros(method.N - 1, 1);
    the_diag = zeros(method.N - 1, 1);
    sup_diag = zeros(method.N - 1, 1);
    ones_vec = ones(method.N - 1, 1);

    % Switch spatial discretizations
    switch method.method
        case 'central'
            if method.uniform
                sub_diag = (params.k / dx ^ 2 + params.u / (2 * dx)) * ones_vec;
                the_diag = -2 * params.k / dx ^ 2 * ones_vec;
                sup_diag = (params.k / dx ^ 2 - params.u / (2 * dx)) * ones_vec;
            else
                h = method.grid(2 : method.N + 1) - method.grid(1 : method.N);
                h_prev = h(1 : method.N - 1);
                h_next = h(2 : method.N);
                h_sum = h_prev + h_next;
                
                sub_diag = 2 * params.k ./ (h_prev .* h_sum) + params.u ./ h_sum;
                the_diag = -2 * params.k ./ (h_prev .* h_next);
                sup_diag = 2 * params.k ./ (h_next .* h_sum) - params.u ./ h_sum;
            end
        case 'upwind'
            ssub_diag = zeros(method.N - 1, 1);
            sub_diag = (params.k / dx ^ 2 + params.u / dx) * ones_vec;
            the_diag = (-2 * params.k / dx ^ 2 - params.u / dx) * ones_vec;
            sup_diag = params.k / dx ^ 2 * ones_vec;
        case 'b3'
            ssub_diag = -params.u / (2 * dx) * ones_vec;
            sub_diag = (params.k / dx ^ 2 + 2 * params.u / dx) * ones_vec;
            the_diag = (-2 * params.k / dx ^ 2 - 3 * params.u / (2 * dx)) * ones_vec;
            sup_diag = params.k / dx ^ 2 * ones_vec;

            ssub_diag(1) = 0;
            sub_diag(1) = params.k / dx ^ 2 + params.u / dx;
            the_diag(1) = -2 * params.k / dx ^ 2 - params.u / dx;
            sup_diag(1) = params.k / dx ^ 2;
        otherwise
            error('Not supported method');
    end


    steps = floor(method.T / method.dt);
    solution = zeros(steps, method.N + 1);
    solution(1, :) = params.initial_condition;

    % Switch between explicit and implicit
    if method.explicit
        for step = 2 : steps
            solution(step, 1) = params.phi_left;
            solution(step, 2) = ...
                    solution(step - 1, 1) * sub_diag(1) * method.dt + ...
                    solution(step - 1, 2) * (1 + the_diag(1) * method.dt) + ...
                    solution(step - 1, 3) * sup_diag(1) * method.dt;
            
            for idx = 3 : method.N
                solution(step, idx) = ...
                    solution(step - 1, idx - 2) * ssub_diag(idx - 1) * method.dt + ...
                    solution(step - 1, idx - 1) * sub_diag(idx - 1) * method.dt + ...
                    solution(step - 1, idx) * (1 + the_diag(idx - 1) * method.dt) + ...
                    solution(step - 1, idx + 1) * sup_diag(idx - 1) * method.dt;
            end
            
            solution(step, method.N + 1) = params.phi_right;
        end
    else
        for step = 2 : steps
            rhs = solution(step - 1, 2 : method.N);
            rhs(1) = rhs(1) + method.dt * sub_diag(1) * params.phi_left;
            rhs(2) = rhs(2) + method.dt * ssub_diag(2) * params.phi_left;
            rhs(method.N - 1) = rhs(method.N - 1) + method.dt * sup_diag(method.N - 1) * params.phi_right;

            solution(step, 1) = params.phi_left;
            solution(step, 2 : method.N) = pdma( ...
                -method.dt * ssub_diag, ...
                -method.dt * sub_diag, ...
                1 - method.dt * the_diag, ...
                -method.dt * sup_diag, ...
                zeros(method.N - 1, 1), ...
                rhs, ...
                method.N - 1 ...
            );
            solution(step, method.N + 1) = params.phi_right;
        end
    end
end
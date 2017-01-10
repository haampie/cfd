function ex2
    % Set some default values.
    N = 32;

    params = struct( ...
        'u', 10, ...
        'k', 1, ...
        'L', 1, ...
        'initial_condition', linspace(0, 1, N + 1), ...
        'phi_left', 0, ...
        'phi_right', 1 ...
    );

    method = struct( ...
        'method', 'central', ...
        'uniform', false, ...
        'grid', linspace(0, 1, N + 1), ...
        'N', N, ...
        'T', 1, ...
        'dt', 0.01, ...
        'explicit', false ...
    );

    %% ORDER OF TEMPORAl ACCURACY 

    % x_values holds the numerical value of the solution at (x = 0.75, t = 0.25)
    y_values = [];
    time_steps_values = 2 .^ (3 : 16);

    for time_steps = time_steps_values
        method.dt = method.T / time_steps;
        solution = unsteady_conv_diff(params, method);

        % x = 0.75 and t = 0.25.
        t_idx = time_steps / 4 + 1;
        x_idx = 3 * method.N / 4 + 1;
        y_values(end + 1) = solution(t_idx, x_idx);
    end
    
    figure;
    orders = order_of_accuracy(y_values, 2);
    semilogx(time_steps_values(3 : end), orders, '-*')
    title('Order of temporal accuracy');
    xlabel('1 / dt');
    ylabel('Order of accuracy');
    grid on;

    %% ORDER OF SPATIAL ACCURACY
    y_values = [];
    N_values = 2 .^ (3 : 12);

    for N = N_values
        method.N = N;
        method.dt = 0.01;
        method.grid = linspace(0, 1, N + 1);
        params.initial_condition = linspace(0, 1, N + 1);
        solution = unsteady_conv_diff(params, method);
                
        % x = 0.75 and t = 10.
        x_idx = 3 * N / 4 + 1; 
        y_values(end + 1) = solution(end, x_idx);
    end
    
    figure;
    orders = order_of_accuracy(y_values, 2);
    semilogx(N_values(3 : end), orders, '-*')
    title('Order of spatial accuracy');
    xlabel('1 / dx');
    ylabel('Order of accuracy');
    grid on;
end

function order = order_of_accuracy(values, factor)
    diff_next = abs(values(3 : end) - values(2 : end - 1));
    diff_prev = abs(values(2 : end - 1) - values(1 : end - 2));
    order = abs(log(diff_next ./ diff_prev)) ./ log(factor);
end
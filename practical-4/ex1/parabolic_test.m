function parabolic_test
    N = 20;
    initial_condition = linspace(0, 1, N + 1);
    [grid, solution] = parabolic_dirichlet(2, 1, 0.05, N, 4, 0.04, initial_condition, 0, 1, 'upwind');
        
    plot(grid, solution)
end
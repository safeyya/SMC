function single_link_manipulator_simulation()
    % System parameters 
    
    I = 1.0;     
    J = 1.0;     
    m0 = 1.0;     
    g = 9.81;  
    L = 1.0;     
    k = 1.0;     
    a = m0 * g * L / I; 
    b = k / I;   
    c = k / J;    
    d = 1 / J;    

    % Control parameters from Theorem 1
    K = 18;       % Control gain 
    delta = 1.1;  % Uncertainty parameter from Theorem 1

  
    
    % Define lambda bound from Equation (7)
    q = 3; 
    n = 4;  % Number of states
    lambda_bound = 1 / (2 * (q + n - 1)); 
    lambda = 0.9 * lambda_bound;  

  
    
    % Values for r_i and s_l from Equation (9)
    r_1 = 4.4;  
    r_2 = 3.3;  
    r_3 = 2.2;  
    r_i = [r_1, r_2, r_3];

    s_1 = 2.2;
    s_2 = 6.6;
    s_3 = 6.6;
    s_l = [s_1, s_2, s_3];  

  
    
    % Compute epsilon from Equation (8)
    epsilon = (lambda * K) / (2 * max([r_i, s_l]));

    % Uncertainty bounds for omega 
    underline_omega_0 = .5; overline_omega_0 = 3;
    underline_omega_1 = 0.5; overline_omega_1 = 2;
    underline_omega_2 = 0.2; overline_omega_2 = 1.5;
    underline_omega_3 = 0.1; overline_omega_3 = 1;

    % Define nu and w 

    nu = 0.001;  % event-triggering threshold

    
    w = 0.0005;  % Chosen value for w, smaller than nu

    %  M for sliding surface 
    M = [1, 3, 3];  

    % Initial conditions (z(0))
    z0 = [pi/6; 0; pi/6; 0];  
    t_span = [0, 50];  

    % Solve 
    [t, x] = ode45(@(t, x) manipulator_dynamics(t, x, M, K, lambda, epsilon, delta, nu, w, a, b, c, d, ...
                                                 underline_omega_0, overline_omega_0, ...
                                                 underline_omega_1, overline_omega_1, ...
                                                 underline_omega_2, overline_omega_2, ...
                                                 underline_omega_3, overline_omega_3), ...
                   t_span, z0);

 
    % Plot the results
    figure('Position', [100, 100, 1600, 1200], 'Color', 'w');
    plot(t, x(:,1), 'r', 'LineWidth', 4);
    hold on;
    plot(t, x(:,3), 'b--', 'LineWidth', 4);
    xlabel('Time (s)', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel('Angle', 'FontSize', 24, 'FontWeight', 'bold');
    title('Angular Positions \theta_1 and \theta_2 Over Time', 'FontSize', 30, 'FontWeight', 'bold');
    legend('\theta_1', '\theta_2', 'Location', 'Best', 'FontSize', 24, 'FontWeight', 'bold');
    set(gca, 'FontSize', 18, 'FontWeight', 'bold');
    grid on;
    ylim([-1, 1.2]);
    yticks(-1:0.5:1); 
    ax = gca;
    ax.YAxis.FontSize = 24; 
    ax.YAxis.FontWeight = 'bold'; 
    set(gca, 'YTickLabel', get(gca,'YTickLabel'), 'FontSize', 24, 'FontWeight', 'bold');
    print(gcf, 'high_res_figure.png', '-dpng', '-r300');
end

function dx = manipulator_dynamics(t, x, M, K, lambda, epsilon, delta, nu, w, a, b, c, d, ...
                                   underline_omega_0, overline_omega_0, ...
                                   underline_omega_1, overline_omega_1, ...
                                   underline_omega_2, overline_omega_2, ...
                                   underline_omega_3, overline_omega_3)
    % Sliding surface
    sigma = compute_sigma(x, M, lambda);

    % State-dependent uncertainties (72)
    [omega_0, omega_1, omega_2, omega_3] = compute_omega(x, a, b, c, d);

    % Event-triggered control mechanism
    u = event_triggered_control(t, x, sigma, K, epsilon, delta, nu, w, omega_0, omega_1, omega_2, omega_3);

    % System dynamics update
    dx = [x(2);
          -a * sin(x(1)) - b * (x(1) - x(3));
          x(4);
          c * (x(1) - x(3)) + d * u];
end

function sigma = compute_sigma(x, M, lambda)
    % Sliding surface calculation 
    sigma = lambda * (M(1)*x(1) + M(2)*x(2) + M(3)*x(3) + x(4));  
end

function u = event_triggered_control(t, x, sigma, K, epsilon, delta, nu, w, omega_0, omega_1, omega_2, omega_3)
    % Control law with event-triggering from Theorem 1
    sign_sigma = sign(sigma);

    % Compute \mathcal{A}_i(x) for each index i
    A1 = compute_A1(x, epsilon);
    A2 = compute_A2(x, epsilon);
    A3 = compute_A3(x, epsilon);
    
    % Compute \mathcal{B}_\ell(x) for each index l
    x_sharp = compute_x_sharp(x);  %  x_sharp in \mathcal{B}
    B1 = compute_B1(x_sharp, epsilon);
    B2 = compute_B2(x_sharp, epsilon);
    B3 = compute_B3(x_sharp, epsilon);
    
    % Control law from Equation 21
    u = -K * sign_sigma + abs(sign_sigma) * (A1 + B1 + A2 + B2 + A3 + B3);

    % Apply event-triggering mechanism (22)-(24)
    if abs(sigma) >= nu  % Condition from Equation (22)
        u = u;  % Apply control when sigma exceeds nu
    elseif abs(sigma) == w  
        u = u;  
    elseif abs(sigma) < w  
        u = 0;  % keep the control input 
    end
end

function x_sharp = compute_x_sharp(x)
    % Define x_sharp 
    x_sharp = [x(2), x(3), x(4)];  
end

% Functions for \mathcal{A}_i(x) for i = 1, 2, 3
function A1 = compute_A1(x, epsilon)
    g_i1 = sin(x(1));  
    if abs(g_i1) >= epsilon
        A1 = -4.4 * sign(g_i1);  % Using r_1 = 4.4
    else
        A1 = (-g_i1 / epsilon);  
    end
end

function A2 = compute_A2(x, epsilon)
    g_i2 = x(2)^2;  
    if abs(g_i2) >= epsilon
        A2 = -3.3 * sign(g_i2);  % Using r_2 = 3.3
    else
        A2 = (-g_i2 / epsilon);  
    end
end

function A3 = compute_A3(x, epsilon)
    g_i3 = x(3);  
    if abs(g_i3) >= epsilon
        A3 = -2.2 * sign(g_i3);  % Using r_3 = 2.2
    else
        A3 = (-g_i3 / epsilon);  
    end
end

% Functions for \mathcal{B}_\ell(x) for l = 1, 2, 3
function B1 = compute_B1(x_sharp, epsilon)
    g_l1 = abs(x_sharp(1));  
    if abs(g_l1) >= epsilon
        B1 = 2.2 * sign(g_l1);  % Using s_1 = 2.2
    else
        B1 = (g_l1 / epsilon);  
    end
end

function B2 = compute_B2(x_sharp, epsilon)
    g_l2 = abs(x_sharp(2));  
    if abs(g_l2) >= epsilon
        B2 = 6.6 * sign(g_l2);  % Using s_2 = 6.6
    else
        B2 = (g_l2 / epsilon);  
    end
end

function B3 = compute_B3(x_sharp, epsilon)
    g_l3 = abs(x_sharp(3));  
    if abs(g_l3) >= epsilon
        B3 = 6.6 * sign(g_l3);  % Using s_3 = 6.6
    else
        B3 = (g_l3 / epsilon);  
    end
end

% Function for computing the omega values (Equation 72)
function [omega_0, omega_1, omega_2, omega_3] = compute_omega(x, a, b, c, d)
    theta1 = x(1);
    omega_0 = b * d;  
    omega_1 = -a * c;  
    omega_2 = a * sin(theta1);  
    omega_3 = -(a * cos(theta1) + b + c);  
end

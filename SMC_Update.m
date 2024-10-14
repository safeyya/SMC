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
    K = 8;       % Control gain for sliding mode control
    delta = 1.1;  % Uncertainty parameter from Theorem 1

    % Define lambda bound from Equation (7)
    q = 3; 
    n = 4;  % Number of states
    lambda_bound = 1 / (2 * (q + n - 1)); 
    lambda = 0.9 * lambda_bound;  

    % Values for r_i and s_i from Equation (9)
    r_i = [1.5, 2.0, 1.8]; 
    s_i = [0.9, 1.0, 1.2];  

    % Compute epsilon from Equation (8)
    epsilon = (lambda * K) / (2 * max([r_i, s_i]));

    % Uncertainty bounds for omega 
    underline_omega_0 = .5; overline_omega_0 = 3;
    underline_omega_1 = 0.5; overline_omega_1 = 2;
    underline_omega_2 = 0.2; overline_omega_2 = 1.5;
    underline_omega_3 = 0.1; overline_omega_3 = 1;

    % Define nu and w from the draft (based on the event-triggering condition)
    nu = 0.001;  % event-triggering threshold
    w = 0.0005;  % Chosen value for w, smaller than nu

    % Matrix M for sliding surface (from draft)
    M = [1, 3, 3];  % Example values for the sliding surface weights

    % Initial conditions (z(0))
    z0 = [pi/6; 0; pi/6; 0];  
    t_span = [0, 30];  

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

    %  sigma_dot
    sigma_dot = compute_sigma_dot(x, M, a, b, c);  % Time derivative of sliding surface

    %  state-dependent uncertainties (72)
    [omega_0, omega_1, omega_2, omega_3] = compute_omega(x, a, b, c, d);

    %  event-triggered control mechanism
    u = event_triggered_control(t, x, sigma, sigma_dot, K, epsilon, delta, nu, w, omega_0, omega_1, omega_2, omega_3);

    % System dynamics update
    dx = [x(2);
          -a * sin(x(1)) - b * (x(1) - x(3));
          x(4);
          c * (x(1) - x(3)) + d * u];
end

function sigma = compute_sigma(x, M, lambda)
    % Sliding surface calculation (sigma) with M
    sigma = lambda * (M(1)*x(1) + M(2)*x(2) + M(3)*x(3) + x(4));  
end

function sigma_dot = compute_sigma_dot(x, M, a, b, c)
    % Time derivative of sliding surface calculation (sigma_dot)
 
    
   sigma_dot = M(1)*x(2) + M(2)*(-a * sin(x(1)) - b * (x(1) - x(3))) + M(3)*x(4);



end 

function u = event_triggered_control(t, x, sigma, sigma_dot, K, epsilon, delta, nu, w, omega_0, omega_1, omega_2, omega_3)
    % Control law with event-triggering from Theorem 1
    sign_sigma = sign(sigma);

    % Compute \mathcal{A}_i(x) 
    A = compute_A(x, epsilon);  
    
    % Compute \mathcal{B}
    x_sharp = compute_x_sharp(x);  
    B = compute_B(x_sharp, epsilon);  
    
    % Control law from Equation 22
    u = -K * sign_sigma + abs(sign_sigma) * (A + B);

    % Apply event-triggering mechanism 
    if abs(sigma) > nu && abs(sigma_dot) > w * delta  % Case 1: sigma > nu 
        u = u + 8;  % Control self adjustment during event-trigger
    
    elseif abs(sigma) == nu  % Case 2: sigma equals nu
        % Keep the control input until sigma hits w 
   
    elseif abs(sigma) < w  % Case 3: sigma is less than w
        % Wait for sigma to hit nu again
    end


end

function x_sharp = compute_x_sharp(x)
   
    x_sharp = [x(2), x(3), x(4)];  
end

function A_val = compute_A(x, epsilon)

    g_i = [x(1), x(2), x(3)];  

    % Approximate using Equation (12)


    if abs(g_i(1)) >= epsilon
        A_val = 2 * abs(g_i(1));  
    else
        A_val = (-1) * (g_i(1) / epsilon);  
    end
end

function B_val = compute_B(x_sharp, epsilon)
    % Piecewise function \mathcal B based on x_sharp
    if abs(x_sharp(1)) >= epsilon
        B_val = 0.8 * abs(x_sharp(1)) + 0.5 * abs(x_sharp(2));  % Apply original B function
    else
        B_val = (-1) * (x_sharp(1) / epsilon);  % Approximate using Equation (12)
    end
end

function [omega_0, omega_1, omega_2, omega_3] = compute_omega(x, a, b, c, d)
    theta1 = x(1);
    theta1_dot = x(2);
    
    % Compute omega values as per Equation (72) in the draft
  
    
    omega_0 = b * d;  
    omega_1 = -a * c;  
    omega_2 = a * sin(theta1);  
    omega_3 = -(a * cos(theta1) + b + c);  


end 

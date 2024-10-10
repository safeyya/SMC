function single_link_manipulator_simulation()
    % Parameters from draft 
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

    % Control parameters from Theorem 1 and Lemma 1
    K = 10;       % Control gain for sliding mode control
   
   
    nu = 0.001;    % Defined as per draft assumptions
   
    epsilon = 0.01; 
   
   
    
    gamma = 1.01;  

    % Ensure that w is smaller than nu based on draft
    w = nu / 2;  

    % Assumptions for lambda, q, and n based on the draft
    q = 3; 
    n = 4; % n is the number of states of the system

  
    lambda_bound = w / (2 * w * (q + n - 1)); 
    lambda = 0.9 * lambda_bound;  % lambda within the appropriate bound

     M = [1, 3, 3]; 


    delta = 1.1;   % Uncertainty parameter from Theorem 1
    underline_omega_0 = 1; overline_omega_0 = 3;
    underline_omega_1 = 0.5; overline_omega_1 = 2;
    underline_omega_2 = 0.2; overline_omega_2 = 1.5;
    underline_omega_3 = 0.1; overline_omega_3 = 1;

   

% Compute state-dependent uncertainties 
function [omega_0, omega_1, omega_2, omega_3] = compute_omega(x)
    theta1 = x(1);
    theta1_dot = x(2);
    omega_0 = 1 + 0.5 * sin(theta1);   
    omega_1 = 0.5 * cos(theta1);       
    omega_2 = 1 + theta1_dot^2;        
    omega_3 = 0.1 * abs(theta1_dot);  
end





    z0 = [pi/6; 0; pi/6; 0]; %  initial condition z(0)
   
    t_span = [0, 30];

       % A and B terms for the control law
    A = [1, 1, 1]; 
    B = [1, 1];    
    g = [1, 1, 1];
    x_sharp = [1, 1]; 





    % Solve
    [t, x] = ode45(@(t, x) manipulator_dynamics(t, x, M, K, nu, lambda, w, delta, ...
            A, B, g, x_sharp, a, b, c, d), t_span, z0);

    % Plot the angular positions 
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
    yticks([-1:0.5:1]); 
    
 
    ax = gca; 
    ax.YAxis.FontSize = 24; 
    ax.YAxis.FontWeight = 'bold'; 
    set(gca, 'YTickLabel', get(gca,'YTickLabel'), 'FontSize', 24, 'FontWeight', 'bold');  

   
    print(gcf, 'high_res_figure.png', '-dpng', '-r300');
end


  



function dx = manipulator_dynamics(t, x, M, K, nu, lambda, w, delta, ...
                                   A, B, g, x_sharp, ...
                                   a, b, c, d)

    sigma = compute_sigma(x, M, lambda);  
    sigma_dot = compute_sigma_dot(x, M, a, b, c);

    % Implement event-triggered control based on the sliding surface σ
    u = event_triggered_control(x, sigma, sigma_dot, nu, w, K, delta, A, B, g, x_sharp, t);

 

    
    
    
    % Dynamics update 
    dx = [x(2);
          -a * sin(x(1)) - b * (x(1) - x(3));
          x(4);
          c * (x(1) - x(3)) + d * u]; 
end


 


function u = event_triggered_control(x, sigma, sigma_dot, nu, w, K, delta, ...
                                                A, B, g, x_sharp, t)
    
  
    sign_sigma = sign(sigma);

    A_term = sum(A .* g); % Sum for q components
    B_term = sum(B .* x_sharp); % Sum for n-1 components

    % Implement control law from Equation 22
    u = -K * sign_sigma + abs(sign_sigma) * (A_term + B_term);

    % Check triggering conditions based on equations (23) - (27)
    if abs(sigma) > nu  % Case 1: σ(x(tr)) > ν
        if abs(sigma_dot) > w * delta
            % Triggering occurs
        end
    elseif abs(sigma) == w % Case 2: σ(x(tr)) = w
        
    else  % Case 3: σ(x(tr)) < w
        % No triggering until σ hits ν again
    end
end



function sigma = compute_sigma(x, M, lambda)
    sigma = lambda * (M(1)*x(1) + M(2)*x(2) + M(3)*x(3) + x(4));  
end

%  compute time derivative of sliding surface
function sigma_dot = compute_sigma_dot(x, M, a, b, c)
    sigma_dot = M(1)*x(2) + M(2)*(-a * sin(x(1)) - b * (x(1) - x(3))) + M(3)*x(4);
end 


% Interations to equal outlet temperature to 4 decimal places
% Initialize variables
max_iterations = 1000; % Maximum number of iterations to prevent infinite loop
tolerance = 1e-4;      % Tolerance for the error difference (4 decimal places)

% Initial Assumptions
previous_value = 304.423;           % Inlet temperature of coolant, K                
current_value = previous_value -1;  % T_Outlet_Coolant = T_inlet_coolant - 1;
T_inlet_air = 303.403;          % Inlet temperature of air, K
T_Outlet_Air = T_inlet_air + 1; % Outlet T of Air, K
Epsilon = 0.4;                  % Heat Changer Effectiveness

% Loop until the error difference is within the specified tolerance
iteration = 0;
while abs(current_value - previous_value) >= tolerance
    % Increment the iteration counter
    iteration = iteration + 1;
    
    % Check if the maximum number of iterations has been reached
    if iteration > max_iterations
        warning('Maximum number of iterations reached. The solution may not have converged.');
        break;
    end
    
    % Update previous_value with current_value
    previous_value = current_value;
    
    % Calculate the new current_value using the example function or some process
    [current_value, T_Outlet_Air, Epsilon] = evrad_calc(previous_value,T_Outlet_Air,Epsilon);
    
    % Display the iteration number and current error difference
    % fprintf('Iteration %d: Current Value = %.6f, Error Difference = %.6f\n', iteration, current_value, abs(current_value - previous_value));
end

% Display the final result
fprintf('Final Value of Coolant Outlet Temp: %.4f K and Air Outlet Temp: %.4f K after %d iterations with Error Difference = %.6f\n', current_value, T_Outlet_Air, iteration, abs(current_value - previous_value));
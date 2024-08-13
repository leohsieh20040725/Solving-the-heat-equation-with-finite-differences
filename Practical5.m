%Question 4

% parameters
alpha = 1e-6; % thermal diffusivity in m^2/s
L = 20; % length in meters
TO = 0; % temperature in degrees Celsius
Td = 1000; % temperature in degrees Celsius
w = 20e-2; % width in meters
tfinal = 100*24*60*60; % time in seconds

% discretization
I = 201; % number of spatial points
x = linspace(0, L, I); % spatial grid
dx = x(2) - x(1); % step size
dtmax = 0.5*((dx^2)/alpha); % maximum time step size based on stability criterion
dt = 0.5*dtmax; % chosen time step size
t = 0:dt:tfinal; % time vector

% Initialize temperature matrix
N = length(t); % number of time steps
T = zeros(N, I); % temperature matrix
Tinit = TO + (Td - TO).*exp(-(x.^2)/(2*w.^2)); % initial temperature distribution
T(1, :) = Tinit; % set initial temperature distribution
T(1, I) = TO;  % enforce boundary condition at x = L

% loop
for n = 1:N-1
    % Boundary condition at x = L
    T(n+1, I) = TO; % enforce boundary condition at x = L
    
    % Boundary condition at x = 0: Zero flux (Neumann boundary condition)
    T(n+1, 1) = T(n, 1) + (alpha*dt/dx^2)*(2*T(n, 2) - 2*T(n, 1)); % No flux at x = 0
    
    % Interior points
    for i = 2:I-1
        % Finite difference scheme for heat equation
        T(n+1, i) = T(n, i) + (alpha*dt/dx^2)*(T(n, i+1) - 2*T(n, i) + T(n, i-1)); 
    end
end

%Question 5
% Plotting contour of temperature with default axis range
figure;
[C, h] = contour(t, x, T');
colorbar();
xlabel('Time (seconds)');
ylabel('Distance (m)');
title('Temperature in Cooling Dyke (Contour)');
clabel(C, h); % Add temperature labels

% Plotting contour of temperature with adjusted x and y axis range
figure;
[C, h] = contour(t, x, T');
colorbar();
xlabel('Time (seconds)');
ylabel('Distance (m)');
title('Temperature in Cooling Dyke (Adjusted x and y axis range) (Contour)');
xlim([0, 20e5]); % Adjust x-axis range
ylim([0, 2]);   % Adjust y-axis range
clabel(C, h); % Add temperature labels

%Question 6

% Initialize an array to store maximum temperatures at each depth
max_temperatures = zeros(1, I);

% Loop through each depth
for i = 1:I

% Extract temperatures at depth i
depth_temperatures = T(:, i);

% Find the maximum temperature reached at depth i
max_temperatures(i) = max(depth_temperatures);
end

% Plot maximum temperature reached at each depth
figure;
plot(x, max_temperatures, 'LineWidth', 2);
xlabel('Distance (m)');
ylabel('Maximum Temperature (°C)');
title('Maximum Temperature Reached at Each Depth');
grid on;



%Question 8
% Define ranges for w and alpha
w_values = [0.3, 0.4]; % example values for w
alpha_values = [1e-6, 2e-6]; % example values for alpha
% Loop over each combination of w and alpha
for w_index = 1:length(w_values)
for alpha_index = 1:length(alpha_values)
% Get current values of w and alpha
ww = w_values(w_index);
alphas = alpha_values(alpha_index);
% Copy the code from Question 4 for the calculation
% Discretization
II= 201;
xx = linspace(0, L, II);
dxx = xx(2) - xx(1);
dtmaxI = 0.5*((dxx^2)/alphas);
dtt = 0.5*dtmaxI;
tt = 0:dtt:tfinal;
% Initialize temperature matrix
NN = length(tt);
TT = zeros(NN, II);
Tinitt = TO + (Td - TO).*exp(-(xx.^2)/(2*ww.^2));
TT(1, :) = Tinitt;
TT(1, II) = TO;
% Loop
for n = 1:NN-1
TT(n+1, II) = TO;
TT(n+1, 1) = TT(n, 1) + (alphas*dtt/dxx^2)*(2*TT(n, 2) - 2*TT(n, 1));
for i = 2:II-1
TT(n+1, i) = TT(n, i) + (alphas*dtt/dxx^2)*(TT(n, i+1) - 2*TT(n, i) + TT(n, i-1)); 
end
end
% Plotting
% Plotting contour of temperature with default axis range
figure;
[c, h] = contour(tt, xx, TT');
colorbar();
xlabel('Time (seconds)');
ylabel('Distance (m)');
title(['Temperature in Cooling Dyke (Contour) (w = ', num2str(ww), ', alpha = ', num2str(alphas), ')']);
xlim([0, 90e5]); % Adjust x-axis range
ylim([0, 5]);
clabel(c, h)
end
end










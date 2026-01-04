%% PowerCast-Motion Lite: (x,y) Only with Forecasting
% Simplified version for garment worker analysis using only x,y coordinates
% Features: Idle/slow detection + Position forecasting

clear; clc; close all;

%% ==================== STEP 1: LOAD OR GENERATE DATA ====================
fprintf('Loading/generating hand position data...\n');

% OPTION 1: Load your actual data (uncomment and modify)
data1 = load('Garmentstimexy.mat'); % Should contain: time, x, y
data = data1.data;
time = data(:,1);
x = data(:,2);
y = data(:,3);

% OPTION 2: Generate synthetic sewing motion data
% [time, x, y] = generate_sewing_motion_data();

% % % Plot raw data
figure('Position', [100, 100, 1200, 400]);
subplot(1,3,1);
plot(time, x, 'b-', 'LineWidth', 1.5); hold on;
plot(time, y, 'r-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Position');
title('Raw Hand Position Data');
legend({'X position', 'Y position'});
grid on;

%% ==================== STEP 2: EXTRACT MOTION FEATURES ====================
fprintf('Extracting motion features from (x,y)...\n');

% Calculate velocity (2D)
dt = mean(diff(time));
vx = gradient(x, dt);
vy = gradient(y, dt);
velocity = sqrt(vx.^2 + vy.^2);

% Calculate motion "energy" signature
window_size = round(2/dt); % 2-second windows
energy = movmean(velocity.^2, window_size); % Kinetic energy proxy

% Calculate direction consistency
dx = diff(x); dy = diff(y);
angles = atan2(dy, dx);
angle_changes = abs(diff(unwrap(angles))); % Unwrap for continuity
angle_changes = [angle_changes(1);  angle_changes; angle_changes(end)]; % Pad to match length

% Simple periodicity detection using autocorrelation
[acf_x, lags] = autocorr(x, 100);
[acf_y, ~] = autocorr(y, 100);
periodicity = (max(acf_x(10:end)) + max(acf_y(10:end))) / 2;

% % % velocity = velocity(1:end-1);
% % % velocity = velocity(1:end-1);

% Compile features matrix (time × features)
features = [velocity, energy(1:length(velocity)), ...
            angle_changes(1:length(velocity))];
feature_names = {'Velocity', 'Energy', 'DirectionChange'};

% Plot features
% subplot(1,3,2);
% plot(time, velocity, 'b-', 'LineWidth', 1.5); hold on;
% plot(time, energy(1:length(time)), 'r-', 'LineWidth', 1.5);
% plot(time, angle_changes(1:length(time)), 'g-', 'LineWidth', 1.5);
% xlabel('Time (s)'); ylabel('Feature Value');
% title('Extracted Motion Features');
% legend({'Velocity', 'Energy', 'Direction Change'});
% grid on;

%% ==================== STEP 3: IDLE/SLOW DETECTION ====================
fprintf('Detecting idle and slow work periods...\n');

% Adaptive thresholding (percentile-based)
velocity_th_idle = prctile(velocity, 30); % Bottom 10% = idle
velocity_th_slow = prctile(velocity, 50); % Bottom 30% = potentially slow

energy_th_idle = prctile(energy, 30);
direction_th_high = prctile(angle_changes, 90); % Top 10% = erratic

% Classify each time point
work_states = cell(length(time), 1);
confidence = zeros(length(time), 1);

for i = 1:length(time)
    v = velocity(i);
    e = energy(i);
    d = angle_changes(i);
    
    % Rule-based classification
    if v < velocity_th_idle && e < energy_th_idle
        work_states{i} = 'IDLE';
        confidence(i) = 0.95;
    elseif v < velocity_th_slow || d > direction_th_high
        work_states{i} = 'SLOW';
        confidence(i) = 0.8;
    else
        work_states{i} = 'NORMAL';
        confidence(i) = 0.9;
    end
end

% Calculate efficiency metrics
idle_percent = 100 * sum(strcmp(work_states, 'IDLE')) / length(work_states);
slow_percent = 100 * sum(strcmp(work_states, 'SLOW')) / length(work_states);

fprintf('Efficiency Analysis:\n');
fprintf('  Idle time: %.1f%%\n', idle_percent);
fprintf('  Slow time: %.1f%%\n', slow_percent);
fprintf('  Normal work: %.1f%%\n', 100 - idle_percent - slow_percent);

% Plot work states
subplot(1,3,3);
hold on;
colors = containers.Map;
colors('IDLE') = [1, 0.7, 0.7];
colors('SLOW') = [1, 1, 0.7];
colors('NORMAL') = [0.7, 1, 0.7];

% Create colored regions
prev_state = work_states{1};
start_idx = 1;
for i = 1:length(work_states)
    if ~strcmp(work_states{i}, prev_state) || i == length(work_states)
        % Fill region with state color
        area_x = [time(start_idx), time(i), time(i), time(start_idx)];
        area_y = [min(velocity), min(velocity), max(velocity), max(velocity)];
        fill(area_x, area_y, colors(prev_state), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        start_idx = i;
        prev_state = work_states{i};
    end
end

plot(time, velocity, 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Velocity');
title(sprintf('Work States (Idle: %.1f%%, Slow: %.1f%%)', idle_percent, slow_percent));
legend({'Velocity'}, 'Location', 'best');
grid on;

%% ==================== STEP 4: BUILD FORECASTING MODEL ====================
fprintf('\nBuilding forecasting model for (x,y) positions...\n');

% Prepare data for forecasting
forecast_horizon = 50; % Predict 50 steps ahead
lookback_window = 100; % Use 100 past points

% Split data: 80% training, 20% testing
split_idx = floor(0.8 * length(time));
train_time = time(1:split_idx);
train_x = x(1:split_idx);
train_y = y(1:split_idx);

test_time = time(split_idx+1:end);
test_x = x(split_idx+1:end);
test_y = y(split_idx+1:end);

% Train simple forecasting models (AR-like for x and y)
fprintf('Training position forecasters...\n');

% Function to forecast using linear regression on recent history
forecast_x = forecast_position(train_x, test_x, lookback_window, forecast_horizon);
forecast_y = forecast_position(train_y, test_y, lookback_window, forecast_horizon);

% Calculate forecast errors
test_start = split_idx - lookback_window + 1;
true_x_test = x(test_start+lookback_window:test_start+lookback_window+forecast_horizon-1);
true_y_test = y(test_start+lookback_window:test_start+lookback_window+forecast_horizon-1);

mse_x = mean((forecast_x - true_x_test).^2);
mse_y = mean((forecast_y - true_y_test).^2);

fprintf('Forecast Performance:\n');
fprintf('  X position MSE: %.6f\n', mse_x);
fprintf('  Y position MSE: %.6f\n', mse_y);
fprintf('  Combined RMSE: %.4f units\n', sqrt(mse_x + mse_y));

%% ==================== STEP 5: VISUALIZE FORECASTING ====================
figure('Position', [100, 100, 1200, 500]);

% Plot 1: X position forecast
subplot(2,2,1);
hold on;
plot(train_time(end-200:end), train_x(end-200:end), 'b-', 'LineWidth', 1.5);
plot(test_time(1:forecast_horizon), forecast_x, 'r--', 'LineWidth', 2);
plot(test_time(1:forecast_horizon), true_x_test, 'k-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('X Position');
title('X Position Forecasting');
legend({'Past', 'Forecast', 'Actual'}, 'Location', 'best');
grid on;

% Plot 2: Y position forecast
subplot(2,2,2);
hold on;
plot(train_time(end-200:end), train_y(end-200:end), 'b-', 'LineWidth', 1.5);
plot(test_time(1:forecast_horizon), forecast_y, 'r--', 'LineWidth', 2);
plot(test_time(1:forecast_horizon), true_y_test, 'k-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Y Position');
title('Y Position Forecasting');
legend({'Past', 'Forecast', 'Actual'}, 'Location', 'best');
grid on;

% Plot 3: 2D trajectory forecast
subplot(2,2,3);
hold on;
plot(train_x(end-100:end), train_y(end-100:end), 'b-', 'LineWidth', 1.5);
plot(forecast_x, forecast_y, 'r--', 'LineWidth', 2);
plot(true_x_test, true_y_test, 'k-', 'LineWidth', 1);
scatter(forecast_x(1), forecast_y(1), 100, 'go', 'filled');
scatter(forecast_x(end), forecast_y(end), 100, 'ro', 'filled');
xlabel('X Position'); ylabel('Y Position');
title('2D Trajectory Forecast');
legend({'Past Path', 'Forecast', 'Actual', 'Start', 'End'});
grid on; axis equal;

% Plot 4: Forecasting error over horizon
subplot(2,2,4);
forecast_steps = 1:forecast_horizon;
error_x = abs(forecast_x - true_x_test);
error_y = abs(forecast_y - true_y_test);

plot(forecast_steps, error_x, 'b-', 'LineWidth', 1.5); hold on;
plot(forecast_steps, error_y, 'r-', 'LineWidth', 1.5);
plot(forecast_steps, sqrt(error_x.^2 + error_y.^2), 'k-', 'LineWidth', 2);
xlabel('Forecast Steps Ahead'); ylabel('Absolute Error');
title('Forecast Error vs Horizon');
legend({'X error', 'Y error', 'Total error'});
grid on;

%% ==================== STEP 6: REAL-TIME DETECTION & FORECASTING ====================
fprintf('\nSimulating real-time detection and forecasting...\n');

% Select a test point
test_point = split_idx + 50;

% Extract recent history
recent_x = x(test_point-lookback_window+1:test_point);
recent_y = y(test_point-lookback_window+1:test_point);
recent_time = time(test_point-lookback_window+1:test_point);

% Current features for detection
current_v = velocity(test_point);
current_e = energy(test_point);
current_d = angle_changes(test_point);

% Detect current state
if current_v < velocity_th_idle && current_e < energy_th_idle
    current_state = 'IDLE';
elseif current_v < velocity_th_slow || current_d > direction_th_high
    current_state = 'SLOW';
else
    current_state = 'NORMAL';
end

% Forecast from current point
future_x = forecast_position(x(1:test_point), [], lookback_window, forecast_horizon);
future_y = forecast_position(y(1:test_point), [], lookback_window, forecast_horizon);

fprintf('Real-time Analysis at t=%.1fs:\n', time(test_point));
fprintf('  Current state: %s\n', current_state);
fprintf('  Velocity: %.4f\n', current_v);
fprintf('  Forecasted position in %d steps: (%.3f, %.3f)\n', ...
        forecast_horizon, future_x(end), future_y(end));

%% ==================== SUPPORTING FUNCTIONS ====================

function [time, x, y] = generate_sewing_motion_data()
    % Generate synthetic (x,y) sewing motion data
    duration = 300; % 5 minutes
    fs = 30; % 30 Hz
    N = duration * fs;
    time = (0:N-1)' / fs;
    
    % Base sewing pattern (rhythmic motion)
    t = time;
    
    % Main sewing motion (elliptical pattern)
    freq1 = 2.0; % 2 Hz - fast stitch
    freq2 = 0.5; % 0.5 Hz - fabric movement
    
    x = 0.2 * sin(2*pi*freq1*t) + 0.05 * sin(2*pi*freq2*t) + 0.01 * randn(N,1);
    y = 0.15 * cos(2*pi*freq1*t) + 0.03 * cos(2*pi*freq2*t + pi/4) + 0.01 * randn(N,1);
    
    % Add idle periods (stationary)
    idle_starts = [50, 150, 250] * fs; % At 50s, 150s, 250s
    idle_durations = [5, 8, 3] * fs; % 5s, 8s, 3s idle
    
    for i = 1:length(idle_starts)
        start_idx = idle_starts(i);
        end_idx = min(start_idx + idle_durations(i), N);
        x(start_idx:end_idx) = x(start_idx);
        y(start_idx:end_idx) = y(start_idx);
    end
    
    % Add slow work periods (reduced amplitude/frequency)
    slow_starts = [80, 200] * fs;
    slow_durations = [10, 15] * fs;
    
    for i = 1:length(slow_starts)
        start_idx = slow_starts(i);
        end_idx = min(start_idx + slow_durations(i), N);
        segment_len = end_idx - start_idx + 1;
        t_slow = (0:segment_len-1)' / fs;
        
        % Slower, smaller motions
        x_slow = 0.05 * sin(2*pi*0.8*freq1*t_slow) + 0.01 * randn(segment_len,1);
        y_slow = 0.04 * cos(2*pi*0.8*freq1*t_slow) + 0.01 * randn(segment_len,1);
        
        x(start_idx:end_idx) = x(start_idx) + x_slow;
        y(start_idx:end_idx) = y(start_idx) + y_slow;
    end
    
    % Smooth the data
    x = smooth(x, 5);
    y = smooth(y, 5);
end

function forecasts = forecast_position(train_data, test_data, lookback, horizon)
    % Simple position forecasting using linear regression on recent history
    % train_data: historical positions
    % test_data: future positions (for validation, can be empty for prediction)
    % lookback: number of past points to use
    % horizon: number of steps to forecast
    
    if isempty(test_data)
        % Prediction mode: forecast from end of train_data
        all_data = train_data;
        forecast_start = length(all_data);
    else
        % Validation mode: forecast from split point
        all_data = [train_data; test_data];
        forecast_start = length(train_data);
    end
    
    forecasts = zeros(horizon, 1);
    
    % Forecast one step at a time
    for h = 1:horizon
        % Prepare features: recent history
        start_idx = forecast_start + h - lookback;
        end_idx = forecast_start + h - 1;
        
        if start_idx < 1
            % Pad with earliest value if not enough history
            padded_data = [repmat(all_data(1), 1-start_idx, 1); all_data(1:end_idx)];
            X = padded_data';
        else
            X = all_data(start_idx:end_idx)';
        end
        
        % Target: next value
        if end_idx + 1 <= length(all_data)
            y_true = all_data(end_idx + 1);
        else
            y_true = NaN; % Unknown future
        end
        
        % Simple linear predictor: weighted average of recent points
        % (More sophisticated: AR model, LSTM, etc.)
        weights = exp(-0.1 * (lookback-1:-1:0))'; % Exponential decay weights
        weights = weights / sum(weights);
        
        forecast = sum(X .* weights');
        
        forecasts(h) = forecast;
        
        % For multi-step forecasting, append forecast to history
        if h < horizon && end_idx + 1 > length(all_data)
            all_data = [all_data; forecast];
        end
    end
end

function [acf, lags] = autocorr(x, maxlag)
    % Simple autocorrelation function
    N = length(x);
    x = x - mean(x);
    
    if nargin < 2
        maxlag = N-1;
    end
    
    acf = zeros(maxlag+1, 1);
    lags = 0:maxlag;
    
    for k = 0:maxlag
        acf(k+1) = sum(x(1:N-k) .* x(1+k:N)) / (N * var(x));
    end
end

%% ==================== SUMMARY AND USAGE ====================
fprintf('\n==================================================\n');
fprintf('POWERCAST-MOTION LITE: (x,y) ONLY WITH FORECASTING\n');
fprintf('==================================================\n\n');
fprintf('CAPABILITIES:\n');
fprintf('1. IDLE/SLOW DETECTION:\n');
fprintf('   - Uses velocity, energy, and direction changes\n');
fprintf('   - Adaptive percentile-based thresholds\n');
fprintf('   - Real-time state classification\n\n');
fprintf('2. POSITION FORECASTING:\n');
fprintf('   - Predicts future (x,y) positions\n');
fprintf('   - Multi-step forecasting capability\n');
fprintf('   - Visualizes trajectory predictions\n\n');
fprintf('3. EFFICIENCY METRICS:\n');
fprintf('   - Calculates %% idle/slow time\n');
fprintf('   - Provides work state timeline\n\n');

fprintf('TO USE WITH YOUR DATA:\n');
fprintf('1. Replace generate_sewing_motion_data() with your data loader\n');
fprintf('2. Ensure your data has columns: [time, x, y]\n');
fprintf('3. Adjust fs (sampling frequency) if different from 30Hz\n');
fprintf('4. Tune lookback_window and forecast_horizon for your needs\n\n');

fprintf('NEXT STEPS FOR IMPROVEMENT:\n');
fprintf('1. Replace simple forecaster with ARIMA/LSTM for better accuracy\n');
fprintf('2. Add more sophisticated features (curvature, smoothness)\n');
fprintf('3. Implement online learning for adaptive thresholds\n');
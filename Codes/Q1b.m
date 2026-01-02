
Fs = 100;         
T = 1;            
N = Fs * T;       
t = (0:N-1) / Fs; 
B_max = 20;       
k_cutoff = round(B_max * N / Fs) + 1;

% x_n = 2*sin(pi*t) + 5*cos(7*pi*t);
x_n = cos(5*pi*t) + sin(10*pi*t);
%x_n = sin(4*pi*t) + 0.5*cos(6*pi*t);
%x_n = cos(3*pi*t) + sin(8*pi*t);

function [x_reconstructed, error_history] = reconstructer(x_n, available, num_iterations, k_cutoff, N)
    % Initialize x_est with available samples and zero-fill missing ones
    x_est = x_n;
    x_est(~available) = 0;
  
    x_known_values = x_n(available);
    error_history = zeros(1, num_iterations);

    for k = 1:num_iterations

        X_est = fft(x_est); 
        
        %Remove Frequencies above the Cutoff
        X_est(k_cutoff+1 : N-k_cutoff+1) = 0;
        
        x_filtered = ifft(X_est, 'symmetric'); 
        error_history(k) = mean((x_n - x_filtered).^2);
        
        %Replace the estimated signal with the known values
        x_est = x_filtered;
        x_est(available) = x_known_values;
    end
    x_reconstructed = x_est;
end

%available = rand(size(x_n)) > p;

%% Single Example
p = 0.05;          
num_iterations = 40; 

available = rand(size(x_n)) > p;



num_missing = sum(~available);


[x_reconstructed, error] = reconstructer(x_n, available, num_iterations, k_cutoff, N);


% PLOTS

figure('Name', 'POCS Signal Reconstruction - Single Case', 'Position', [100, 100, 1000, 800]);

% Time Domain Comparison
subplot(4, 1, 1);
plot(t, x_n, 'LineWidth', 2, 'DisplayName', 'Original Bandlimited Signal', 'Color', [0 0.5 0]);
hold on;
plot(t(available), x_n(available), 'o', 'MarkerSize', 4, 'DisplayName', 'Corrupted Samples (Known)', 'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', [0.5 0 0]);
plot(t, x_reconstructed, '--', 'LineWidth', 1.5, 'DisplayName', sprintf('Reconstructed Signal (Iters: %d)', num_iterations), 'Color', [0 0 1]);
hold off;
title(sprintf('Time Domain Comparison (p = %.2f, Missing: %.1f%%)', p, (num_missing/N)*100));
xlabel('Time (s)'); ylabel('Amplitude');
legend('show', 'Location', 'NorthEast');
grid on;

% Original vs. Corrupted Spectrum (NEW PLOT)
subplot(4, 1, 2);
X_corr_fft = abs(fftshift(fft(x_corr))); 

plot(f_axis, X_n_fft, 'LineWidth', 2, 'DisplayName', 'Original Spectrum', 'Color', [0 0.5 0]);
hold on;
plot(f_axis, X_corr_fft, 'LineWidth', 1.5, 'DisplayName', 'Corrupted Spectrum (Zero-filled)', 'Color', [0.8 0 0]);
line([-B_max, -B_max], ylim, 'Color', 'r', 'LineStyle', ':', 'DisplayName', 'Bandwidth Limit');
line([B_max, B_max], ylim, 'Color', 'r', 'LineStyle', ':', 'HandleVisibility', 'off');
hold off;
title('2. Original vs. Corrupted Frequency Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
legend('show', 'Location', 'NorthEast');
xlim([-Fs/2, Fs/2]);
grid on;

% Frequency Domain Comparison
subplot(4, 1, 3);
X_n_fft = abs(fftshift(fft(x_n)));
X_reconstructed_fft = abs(fftshift(fft(x_reconstructed)));
f_axis = (-N/2 : N/2-1) * Fs / N; % Frequency axis

plot(f_axis, X_n_fft, 'LineWidth', 2, 'DisplayName', 'Original Spectrum', 'Color', [0 0.5 0]);
hold on;
plot(f_axis, X_reconstructed_fft, '--', 'LineWidth', 1.5, 'DisplayName', 'Reconstructed Spectrum', 'Color', [0 0 1]);
line([-B_max, -B_max], ylim, 'Color', 'r', 'Linewidth',3, 'LineStyle', ':', 'DisplayName', 'Bandwidth Limit');
line([B_max, B_max], ylim, 'Color', 'r','Linewidth',3,  'LineStyle', ':', 'HandleVisibility', 'off');
hold off;
title('Frequency Domain Analysis');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
legend('show', 'Location', 'NorthEast');
xlim([-Fs/2, Fs/2]);
grid on;

% Error Convergence
subplot(4, 1, 4);
plot(1:num_iterations, 20*log(error), 'LineWidth', 2, 'Color', [0.8 0 0]);
title('Reconstruction Error Convergence (MSE)');
xlabel('Iteration Number'); ylabel('Mean Squared Error (MSE) (in dB)');
grid on;
set(gca, 'YScale', 'linear'); 


%% Accurqacy as a Function of P

num_runs = 10;    
num_iterations_sweep = 50; 
P_values = 0.01:0.01:0.1; 
mae_results = zeros(length(P_values), 1);



for idx = 1:length(P_values)
    p_current = P_values(idx);
    current_mae_sum = 0;
    
    for run = 1:num_runs 

        available = rand(size(x_n)) > p_current;
        [x_reconstructed, ~] = reconstructer(x_n, available, num_iterations_sweep, k_cutoff, N);

        current_mae_sum = current_mae_sum + mean(abs(x_n - x_reconstructed));
    end
    
    % Store average MAE for this probability p
    mae_results(idx) = current_mae_sum / num_runs;
end

% PLOT: MAE v P
figure('Name', 'POCS Reconstruction Performance Sweep', 'Position', [100, 100, 800, 500]);
plot(P_values, mae_results, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0 0.4 0.7], 'Color', [0.1 0.1 0.1]);
title(sprintf('Mean Absolute Error (MAE) vs. P (Avg. %d Runs)', num_runs));
xlabel('Probability of Missing Sample (p)');
ylabel('Mean Absolute Error (MAE)');
grid on;
set(gca, 'FontSize', 12);
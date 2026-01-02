%% Reading the File

[x, fs] = audioread("File4.wav"); % Change File to see for different audio files
x = 0.5 * (x(:,1) + x(:,2)); % Average out L and R channels         

%% Plotting x vs time
N = length(x);
t = (0:N-1)/fs;        

figure;
plot(t, x);
xlabel("Time (seconds)");
ylabel("Amplitude");
title("Audio Signal x(t)");
grid on;

%% Finding peaks
windowLen = round(0.02*fs);  

energy = movmean(abs(x), windowLen); 
energy = energy / max(energy); % Normalized

[pks, locs] = findpeaks(energy, ...
    'MinPeakHeight', 0.05, ...
    'MinPeakDistance', round(0.1*fs), ...
    'MinPeakProminence', 0.04);

main_thresh = 0.25 * median(pks);
valid = pks > main_thresh; % Eliminates noise

pks  = pks(valid);
locs = locs(valid);

beat_times = locs / fs;

figure;
plot(t, x); hold on;
plot(t, energy * max(abs(x)), 'g', 'LineWidth', 1.5);

stem(beat_times, ones(size(beat_times))*max(abs(x)), 'r', 'filled');

legend("Original", "Smoothened Signal", "Detected Beats");
xlabel("Time(s)");
ylabel("Amplitude");
title("Drum Beat Detection");
grid on;


%% Finding duration of beats

thresh_ratio = 0.05;  % If lesser than this, we'll consider beat 'off'
numBeats = length(locs);

hit_durations = zeros(numBeats,1);

for k = 1:numBeats
    peak_idx = locs(k);

    % Local peak amplitude from raw signal
    peak_val = abs(x(peak_idx));
    thresh = thresh_ratio * peak_val;

    % find either thresh_ratio or next beat as rightmost part of current beat 
    right = peak_idx;
    while right < length(x) && abs(x(right)) > thresh
        right = right + 1;
    end

    hit_durations(k) = (right - peak_idx) / fs;
end

%% Detecting Instruments

% Take small range for just the beat
pre_ms  = 10;
post_ms = 200;

preN  = round(pre_ms  * fs / 1000);
postN = round(post_ms * fs / 1000);

numHits = length(locs);
hit_snips = cell(numHits,1);

for k = 1:numHits
    L = max(1, locs(k) - preN);
    R = min(length(x), locs(k) + postN);
    hit_snips{k} = x(L:R);
end

% Extract duration of beat and mean of frequency components of beat, so make a features array
% to keep track
features = zeros(numHits, 2);

for k = 1:numHits
    h = hit_snips{k};

    %  A weighted average (by freq. amplitude) of frequency components
    H = abs(fft(h));
    f = (0:length(H)-1)*(fs/length(H));
    centroid = sum(f'.*H) / sum(H);

    % Duration
    dur = hit_durations(k);

    features(k,:) = [centroid, dur];
end

features = normalize(features);

% Estimating no. of instruments
Kmax = 8; % (Assuming that we have a max of 8 instruments)
sil = zeros(Kmax,1);

for K = 2:Kmax
    labels_tmp = kmeans(features, K, 'Replicates', 10);
    sil(K) = mean(silhouette(features, labels_tmp)); % Find how well a 'hit' fits into its own cluster compared to other clusters.
end

[~, num_instruments] = max(sil);

disp("Estimated number of instruments: " + num_instruments);

% Final Grouping
instrument_labels = kmeans(features, num_instruments, 'Replicates', 20);
instrument_groups = cell(num_instruments,1);

for k = 1:num_instruments
    instrument_groups{k} = find(instrument_labels == k);
end

figure;
gscatter(features(:,1), features(:,2), instrument_labels);
xlabel("Weighted Average of Frequency Components");
ylabel("Hit Duration");
title("Drum Types Grouping");
grid on;

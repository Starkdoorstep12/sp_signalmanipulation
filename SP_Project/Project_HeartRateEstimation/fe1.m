
sampling_rate = 128; 
threshold = 0.17;     % Threshold for R-peak detection
min_distance = 0.3 * sampling_rate; % Minimum distance between peaks (in samples)

data = load('E1.mat'); 
E1 = data.E1;      

time = (0:length(E1)-1) / sampling_rate;

subplot(1,2,1);
plot(time, E1);
xlabel('Time (s)');
xlim([0 10]);
ylabel('E1');
title('Complete E1 Signal with Respect to Time');

n = length(E1);
f = (-n/2:n/2-1) * (sampling_rate / n);
subplot(1,2,2);
plot(f, abs(fftshift(abs(fft(E1)))));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum of E1 Signal (Shifted)');
grid on;


[pks, locs] = findpeaks(E1, 'MinPeakHeight', threshold, ...
                             'MinPeakDistance', min_distance);

RR_intervals = diff(locs) / sampling_rate;

HR = 60 ./ RR_intervals;

figure;
plot((1:length(E1)) / sampling_rate, E1, 'b'); % Plot entire ECG signal
hold on;
plot(locs / sampling_rate, pks, 'ro'); 
xlabel('Time (s)');
xlim([0 10]);
ylabel('ECG Amplitude');
title('Detected R-Peaks in Entire ECG Signal');
grid on;
legend('ECG Signal', 'R-Peaks');


avg_HR = mean(HR);
fprintf('Average Heart Rate: %.2f bpm\n', avg_HR);

num_minutes=10;
HRmin=[];
for k=1:num_minutes
    start_idx = (k-1)*60 + 1;
    end_idx = min(k*60, length(HR));
    temp = mean(HR(start_idx:end_idx));
    HRmin = [HRmin temp];
end

minutes = 1:num_minutes;

% Plot bar graph of mean HR vs. minutes
figure;
bar(minutes, HRmin);
xlabel('Minutes');
ylabel('Mean Heart Rate (bpm)');
title('Mean Heart Rate per Minute');
grid on;
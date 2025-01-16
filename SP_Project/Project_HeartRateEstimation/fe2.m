
sampling_rate = 128;
threshold = 0.38;     % Threshold for R-peak detection
min_distance = 0.3 * sampling_rate; % Minimum distance between peaks (in samples)

data = load('E2.mat'); 
E2 = data.E2;         

filter_order = 4; 

[b, a] = butter(filter_order, [5 45] / (sampling_rate / 2), 'bandpass');

E2_filtered = filtfilt(b, a, E2);

time = (0:length(E2)-1) / sampling_rate;


subplot(2,2,1);
plot(time, E2);
xlabel('Time (s)');
xlim([0 10]);
ylabel('E2');
title('Complete E2 Signal with Respect to Time');

n = length(E2_filtered); % Number of samples

subplot(2,2,2);
plot(time, E2_filtered);
xlabel('Time (s)');
xlim([0 10]);
ylabel('E2 (Filtered)');
title('Complete E2 Filtered Signal with Respect to Time');


f = (-n/2:n/2-1) * (sampling_rate / n);

subplot(2,2,3);

plot(f, abs(fftshift(abs(fft(E2)))));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ylim([0 10000]);
title('Magnitude Spectrum of  ECG Signal');

subplot(2,2,4);
plot(f, abs(fftshift(abs(fft(E2_filtered)))));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ylim([0 10000]);
title('Magnitude Spectrum of Filtered ECG Signal');



[pks, locs] = findpeaks(E2_filtered, 'MinPeakHeight', threshold, ...
                             'MinPeakDistance', min_distance);

RR_intervals = diff(locs) / sampling_rate;

HR = 60 ./ RR_intervals;


figure;
plot((1:length(E2_filtered)) / sampling_rate, E2_filtered, 'b'); 
hold on;
plot(locs / sampling_rate, pks, 'ro'); % Mark detected R-peaks
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



time = (0:length(E2)-1) / sampling_rate;


figure;
plot(time, E2);
xlabel('Time (s)');
xlim([0 10]);
ylabel('E2');
title('Complete E2 Signal with Respect to Time');


f = (-n/2:n/2-1) * (sampling_rate / n);

figure;
plot(f, abs(fftshift(abs(fft(E2)))));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ylim([0 10000]);
title('Magnitude Spectrum of  ECG Signal');




n = length(E2_filtered); % Number of samples

figure;
plot(time, E2_filtered);
xlabel('Time (s)');
xlim([0 10]);
ylabel('E2 (Filtered)');
title('Complete E2 Filtered Signal with Respect to Time');


f = (-n/2:n/2-1) * (sampling_rate / n);

figure;
plot(f, abs(fftshift(abs(fft(E2_filtered)))));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ylim([0 10000]);
title('Magnitude Spectrum of Filtered ECG Signal');
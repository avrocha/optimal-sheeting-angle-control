clear; clc; close all; 
fig_cnt = 1;

% Select data source
data_source = 'tacking';

switch data_source
    case 'tacking'
        filename      = 'awa_pm_45\awa_data_45.txt';
    case 'awa_100'
        filename      = 'awa_100\awa_data_100.txt';
    otherwise
        disp('Error: Select valid data source.\n')
end

fid = fopen(filename, 'r');
fgets(fid); % Skip title
out = fscanf(fid, ['%f', ',', '%f', ',', '%f', ',', '%f'], [4, inf]);
out = out';

time    = out(:, 1);
awa     = out(:, 2);
awa_hat = out(:, 3);
heading = out(:, 4);

% Small deviation from 5Hz
fs_data = 1/(time(2) - time(1));
n       = size(out, 1);

% Edit time to obtain constant sampling frequency from the sensor
fs_data = 5; % Both datasets are approximately 5 Hz
time    = (0:1/fs_data:(1/fs_data)*(n-1))';

% Frequency analysis from data
fx        = (0:n-1)*(fs_data/n);
y         = fft(awa);
y_hat     = fft(awa_hat);
y_abs     = abs(y).^2/n;
y_hat_abs = abs(y_hat).^2/n;

figure(fig_cnt); clf(fig_cnt);
subplot(2, 1, 1); hold on;
plot(fx, y_abs);
title('Frequency spectrum of AWA raw data')
xlabel('f [Hz]')
ylabel('|awa(jw)|')

subplot(2, 1, 2); hold on;
plot(fx, y_hat_abs);
title('Frequency spectrum of AWA filtered data')
xlabel('f [Hz]')
ylabel('|awa(jw)|')
fig_cnt = fig_cnt + 1;

% LPF - IIR design
bworder  = 5;
fc       = 0.01;
dt       = 1 / fs_data;
[b, a]   = butter(bworder, fc*dt, 'low');
M = bworder + 1; % filter length

figure(fig_cnt); clf(fig_cnt);
freqz(b, a)
fig_cnt = fig_cnt + 1;

% Add prefix
x    = [awa(1)*ones(1, bworder+1), awa'];
y_bw = [awa(1)*ones(1, bworder+1), zeros(1, n)];

% Flip coeffs for vector-wise multiplication
a = fliplr(a);
b = fliplr(b);

for i = M:n+M
    y_bw(i) = 1/a(end) * (-a(1:end-1) * y_bw(i-M+1:i-1)' + b * x(i-M+1:i)');
end

% EMA filter
y_ema  = [awa(1), zeros(1, n-1)];
alpha  = 0.005; 

for i = 2:n
    y_ema(i) = (1-alpha)*y_ema(i-1) + alpha*x(i);
end

figure(fig_cnt); clf(fig_cnt);
freqz([alpha, 0], [1, -(1-alpha)])
fig_cnt = fig_cnt + 1;

figure;
hold on;
plot(1:n, awa)
plot(1:n, y_bw(M+1:n+M), '--', 'LineWidth', 1.2)
plot(1:n, y_ema, '--', 'LineWidth', 1.2)
legend('awa', 'Butterworth LPF', 'EMA LPF')
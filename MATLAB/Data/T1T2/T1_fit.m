clc;clear;close all;
tic
load T1data.mat

T1Data = NmrData.SPECTRA;
if size(T1Data,1)>size(T1Data,2)
    T1Data = T1Data.';  
end  
t = NmrData.d2;
cs=NmrData.Specscale;     % chemical shift

T1Data = real(T1Data);
T1Data = T1Data / max(T1Data(:));

[peaks, locs] = findpeaks(T1Data(1, :), 'MinPeakHeight', 0.01);
S = T1Data(:, locs);
S = 1 - S ./ S(10, :);
NormS = S ./ S(1, :);
figure(1);
plot(1 : length(T1Data(1, :)), T1Data(1, :));
hold on;
plot(locs, peaks, '*', 'Color','r');

figure(2);
plot(t, NormS)

figure(3)
plot(t, S)

ppm = cs;
b = t';
S = S';
idx_peaks = locs;
save('T1_net_input.mat', 'S', 'ppm', 'b', 'idx_peaks', '-mat');


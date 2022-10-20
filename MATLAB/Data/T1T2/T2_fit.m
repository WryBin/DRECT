clc;clear;close all;
tic
load T2data.mat

T2Data = NmrData.SPECTRA;
if size(T2Data,1)>size(T2Data,2)
    T2Data = T2Data.';  
end  
t = NmrData.bigtau;
cs=NmrData.Specscale;     % chemical shift

T2Data = real(T2Data);
T2Data = T2Data / max(T2Data(:));

[peaks, locs] = findpeaks(T2Data(1, :), 'MinPeakHeight', 0.01);
S = T2Data(:, locs);
NormS = S ./ S(1, :);
figure(1);
plot(1 : length(T2Data(1, :)), T2Data(1, :));
hold on;
plot(locs, peaks, '*', 'Color','r');

figure(2);
plot(t, NormS)

ppm = cs;
b = t';
S = S';
idx_peaks = locs;
b(end) = 12.8;
save('T2_net_input.mat', 'S', 'ppm', 'b', 'idx_peaks', '-mat');

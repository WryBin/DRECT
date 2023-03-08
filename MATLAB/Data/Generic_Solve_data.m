clc;clear;close all;

% Parameters
type = 'GSP';  % Sample name
In_file_path = 'GSP\GSP.mat';% You can specify your own In_file_path here
result_path = ''; % Specify the directory to save the file
Threshold = 0.03; 
find_peaks = 'False'; % Take either the peak value or all values above a certain threshold.
Denoise = 'False'; % Whether to take the peak value along with the average 
% of the left and right points to achieve a simple noise reduction effect."


% (Currently, only GSP.mat and QGC.mat are at your disposal. 
% We do not possess the authority to publicly disclose any other data. 
% Hence, we solely offer fundamental information extracted from 
% the source data for the purpose of graphical representation. 
% This encompasses b, HNMR, idx_peaks, and ppm.)

% Load file
load(In_file_path)

DOSYData = NmrData.SPECTRA;% DOSYData    
 if size(DOSYData,1)>size(DOSYData,2)
     DOSYData = DOSYData.';  
 end

g=100*NmrData.Gzlvl; % gradient values
BD=NmrData.DELTAOriginal; % diffusion time
LD=NmrData.deltaOriginal; % diffusion encoding time
cs=NmrData.Specscale;     % chemical shift
gamma = 4257.7;
g2 = (2*pi*gamma*g*LD).^2*(BD-LD/3)*1e4;
b = g2*1e-10;

% For the sample VD, only data within a limited range are processed
if strcmp(type, 'VD')
    DOSYData = DOSYData(:, 2991:3277); % VD:fn of Fourier transform is 8192
end

% Get the real part of a signal
DOSYData = real(DOSYData);

% Normalize
DOSYData = DOSYData / max(DOSYData(:));

% Preprocess to remove the irrelevant values.
if strcmp(find_peaks, 'True')
    [peaks, idx_peaks] = findpeaks(DOSYData(1, :), 'MinPeakHeight', Threshold, 'MinPeakProminence',0.001);
else
    idx_peaks = find(DOSYData(1, :) >= Threshold);
end
ppm = cs;

if strcmp(type, 'VD')
    ppm = ppm(2991:3277); % VD
end

if strcmp(Denoise, 'True')
    S = DOSYData(:, idx_peaks) + DOSYData(:, idx_peaks+1) + DOSYData(:, idx_peaks-1);
else 
    S = DOSYData(:, idx_peaks);
end

% Create a plot to assist in parameter selection
whole_spec = zeros([length(b), length(ppm)]);
whole_spec(:, idx_peaks) = S;

figure(1)
subplot(2, 1, 1)
plot(ppm, DOSYData(1, :))
hold on
if strcmp(find_peaks, 'True')
    plot(ppm(idx_peaks), peaks, '*', 'Color','r', MarkerSize=1.5);
end
title("Original ^1H NMR")
subplot(2, 1, 2)
plot(ppm, whole_spec(1, :))
S_new = S ./ S(1, :);
title("The values above the specified threshold")

figure(2)
plot(b, S_new(:, :))
title("The normalized decay signal")
S = S';

% Save the results in .mat format, which can be used as input for the model
save([result_path, type, '_net_input.mat'], 'S', 'ppm', 'b', 'idx_peaks', '-mat');
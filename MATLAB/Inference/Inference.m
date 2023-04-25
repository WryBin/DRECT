clc;clear;close all;
clear classes

%% Parameters
% Parameters of preprocessing
type = 'GSP';
In_file_path = '..\Data\GSP\GSP.mat';% You can specify your own In_file_path here
result_path = ''; % Specify the directory to save the file
Threshold = 0.03; 
find_peaks = 'False'; % Take either the peak value or all values above a certain threshold.
Denoise = 'False'; % Whether to take the peak value along with the average 
% of the left and right points to achieve a simple noise reduction effect."

% Parameter for DRILT
model_path = "..\..\Result\DOSY\last.ckpt";

% (Currently, only GSP.mat and QGC.mat are at your disposal. 
% We do not possess the authority to publicly disclose any other data. 
% Hence, we solely offer fundamental information extracted from 
% the source data for the purpose of graphical representation. 
% This encompasses b, HNMR, idx_peaks, and ppm.)

% Parameters of plotting
scale = 1.0;  % Parameters used to adjust the calculation range of D. When the required value of the diffusion coefficient exceeds 
% the computational boundary, increasing this parameter enables obtaining a complete output result. However, it is important to 
% note that this may affect the processing effect of data with diffusion coefficient values that are close to each other.

dc1 = 1.5;  % dc1, dc2: Display range of diffusion coefficient dimension              
dc2 = 5;


%% Preprocessing
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

% interpolate
NmrDatai = zeros(size(S, 1), 30);
for i = 1:size(S, 1)
    NmrDatai(i,:) = interp1(b(1,:), S(i,:), linspace(0, max(b(1,:)), 30), 'linear', 'extrap');
end

S = NmrDatai;
b = linspace(0, max(b(1,:)), 30);

SPECTRA = zeros(length(ppm), length(S(1, :)));
SPECTRA(idx_peaks, :) = S;
NmrData.SPECTRA = S;
NmrData.arraydim = length(S(1, :));
NmrData.ngrad = length(S(1, :));
NmrData.Ppmscale = ppm;
NmrData.np = length(ppm);
NmrData.b = b; % b was changed because intepolated
NmrData.Gzlvl=sqrt(b/(NmrData.dosyconstant*1e-10));% for interpolated

%% DRILT
% Configuration for running python files in matlab
if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end

mod = py.importlib.import_module('DRILT');
py.importlib.reload(mod);



max_D_DRILT = 14*(0.8*scale/max(b));

% Get the result
DRILT_result = double(py.DRILT.test(scale, py.numpy.array(S), py.numpy.array(b), model_path))';

%% Plot

% ploting parameters
cs1 = ppm(idx_peaks(1)-fix(length(ppm)/100));  % cs1, cs2: Display range of chemical shift dimension
cs2 = ppm(idx_peaks(end)+fix(length(ppm)/100));

figure()
t = tiledlayout(4, 1);
t.TileSpacing = 'tight';
% t.Padding = 'tight';

% Sparsity
DRILT_result(DRILT_result < repmat((max(DRILT_result, [], 1) * 0.7), [140, 1])) = 0;

nexttile()  %1H NMR
cs_spec = zeros([length(ppm), 1]);
cs_spec(idx_peaks, :) = S(:, 1);
plot(ppm,cs_spec, "Color",'k');set(gca,'Xdir','reverse');axis off;
xlim([cs1,cs2]);

nexttile(2, [3, 1])
result_plot(DRILT_result, ppm, idx_peaks, dc1, dc2, cs1, cs2,max_D_DRILT);

function result_plot(result, ppm, idx_peaks,dc1, dc2, cs1, cs2, max_D)

    contour_level = 40; % Determines the number and positions of the contour lines / regions.

    spec_whole = zeros([length(result(:, 1)), length(ppm)]);
    spec_whole(:, idx_peaks) = result;
    decay_range = linspace(0, max_D, length(result(:, 1)));
    contour(ppm,decay_range, real(spec_whole), contour_level);
    
    set(gca,'Ydir','reverse','Xdir','reverse'); 
    
    xlim([cs1,cs2]);
    ylim([dc1,dc2]);
end

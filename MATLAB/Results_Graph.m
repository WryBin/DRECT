% An example showing how to plot the result
close all;clear;clc
load ./Result/QGC_Result.mat  % Change to your own filename
load ./Data/QGC/QGC_net_input.mat

% parameters
dc1 = 0;  % dc1, dc2: Display range of diffusion coefficient dimension              
dc2 = calculated_max_D;
cs1 = ppm(idx_peaks(1)-fix(length(ppm)/100));  % cs1, cs2: Display range of chemical shift dimension
cs2 = ppm(idx_peaks(end)+fix(length(ppm)/100));
contour_level = 40; % Determines the number and positions of the contour lines / regions.

% plot
figure()
t = tiledlayout(3, 1);
t.TileSpacing = 'tight';
t.Padding = 'tight';

nexttile()
cs_spec = zeros([length(ppm), 1]);
cs_spec(idx_peaks, :) = S(:, 1);
plot(ppm,cs_spec, "Color",'k');set(gca,'Xdir','reverse');axis off;
xlim([cs1,cs2]);

nexttile(2, [2, 1])

spec_whole = zeros([length(Z(1, :)), length(ppm)]);
spec_whole(:, idx_peaks) = Z.';
decay_range = linspace(0, calculated_max_D, length(Z(1, :)));
contour(ppm,decay_range,spec_whole,contour_level);

set(gca,'Ydir','reverse','Xdir','reverse'); 

xlim([cs1,cs2]);
ylim([dc1,dc2]);
xlabel('Chemical Shift(ppm)');
ylabel('Diffusion Coefficien“t(10^{-10}m^2/s)');
% An example showing how to plot the result
close all;clear;clc
load ./Result/QGC_Result.mat  % Change to your own filename
load ./Data/QGC/QGC_net_input.mat

% parameters
cs1 = 4;  % cs1, cs2: Display range of chemical shift dimension
cs2 = 12.5;
dc1 = 1;  % dc1, dc2: Display range of diffusion coefficient dimension              
dc2 = 12;

contour_level = 40; % Determines the number and positions of the contour lines / regions.
max_b_train = 0.8; % Consistent with the max_b set in the file config.py


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
decay_range = linspace(0, (length(Z(1, :))-1)/10, length(Z(1, :)));
contour(ppm,decay_range*(max_b_train/b(end)),spec_whole,contour_level);
set(gca,'Ydir','reverse','Xdir','reverse'); 

xlim([cs1,cs2]);
ylim([dc1,dc2]);
xlabel('Chemical Shift(ppm)');
ylabel('Diffusion Coefficien“t(10^{-10}m^2/s)');
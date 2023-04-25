clc;clear;close all;
clear classes

%% Preprocessing
In_file_path = '..\Data\GSP\GSP.mat';% You can specify your own In_file_path here
Threshold = 0.03; 
find_peaks = 'False'; % Take either the peak value or all values above a certain threshold.
Denoise = 'False'; % Whether to take the peak value along with the average 
% of the left and right points to achieve a simple noise reduction effect."

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

% Get the real part of a signal
DOSYData = real(DOSYData);

% Normalize
DOSYData = DOSYData / max(DOSYData(:));

% Threshold
if strcmp(find_peaks, 'True')
    [peaks, idx_peaks] = findpeaks(DOSYData(1, :), 'MinPeakHeight', Threshold, 'MinPeakProminence',0.001);
else
    idx_peaks = find(DOSYData(1, :) >= Threshold);
end
ppm = cs;

if strcmp(Denoise, 'True')
    S = DOSYData(:, idx_peaks) + DOSYData(:, idx_peaks+1) + DOSYData(:, idx_peaks-1);
else 
    S = DOSYData(:, idx_peaks);
end

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

%% SCORE
% default parameter in GNAT
ncomp = 3; 
speclim = [-2.5, 12.5, 20];
NmrData.SCOREopts=[0 0 0 0 0 0 0 0 0];
NmrData.nug = [1, 0, 0, 0];
fixed=[];

scoredata=score_mn(NmrData,ncomp,speclim,NmrData.SCOREopts,NmrData.nug,fixed);

component1 = (normpdf(1:140, scoredata.Dval(1)*10, 1)/max(normpdf(1:140, scoredata.Dval(1)*10, 1)))'*scoredata.COMPONENTS(1,:);
component2 = (normpdf(1:140, scoredata.Dval(2)*10, 1)/max(normpdf(1:140, scoredata.Dval(2)*10, 1)))'*scoredata.COMPONENTS(2,:);
component3 = (normpdf(1:140, scoredata.Dval(3)*10, 1)/max(normpdf(1:140, scoredata.Dval(3)*10, 1)))'*scoredata.COMPONENTS(3,:);
SCORE_result_GSP = component1+component2+component3;
%% DECRA

ncomp = 4;  % if set to 3, only two components can be shown.
decradata=decra_mn(NmrData,ncomp,speclim);

component1 = (normpdf(1:140, decradata.Dval(1)*10^11, 1)/max(normpdf(1:140, decradata.Dval(1)*10^11, 1)))'*decradata.COMPONENTS(1,:);
component2 = (normpdf(1:140, decradata.Dval(2)*10^11, 1)/max(normpdf(1:140, decradata.Dval(2)*10^11, 1)))'*decradata.COMPONENTS(2,:);
component3 = (normpdf(1:140, decradata.Dval(3)*10^11, 1)/max(normpdf(1:140, decradata.Dval(3)*10^11, 1)))'*decradata.COMPONENTS(3,:);
DECRA_result_GSP = component1+component2+component3;

%% CONTIN
% rilt.m is an emulation of S. Provencher CONTIN program
% We made some necessary modifications to accommodate the data to be processed

G = zeros([length(idx_peaks), 20]);
f = waitbar(0, 'Sovling:  0%', 'Name', 'Please Wait','CreateCancelBtn','setappdata(gcbf,''canceling'',1)'); % revised 20220328 by chen
setappdata(f,'canceling',0); 

for i = 1:length(idx_peaks)
    if getappdata(f, 'canceling')
        break
    end
    str=['CONTIN: Sovling:  ',num2str(round(i/(length(idx_peaks))*100)),'%']; 
    waitbar(i / (length(idx_peaks)), f, str); 

    data = S(i, :)';

    alpha = 1e-2;

    % Guess s space and g function
    s = linspace(0, 14, 20)';
    g = ones(size(s));
    [g,yfit,cfg] = rilt(b',data,s,g,alpha,'linear',[],[],[],{'g>0'},[],[]);

    G(i, :) = g';
end
CONTIN_result_GSP = G';

delete (f)
%% DRILT
% Configuration for running python files in matlab
if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end

mod = py.importlib.import_module('DRILT');
py.importlib.reload(mod);

% Parameter for DRILT
model_path = "..\..\Result\DOSY\last.ckpt";

scale = 1.0;  % Parameters used to adjust the calculation range of D. When the required value of the diffusion coefficient exceeds 
% the computational boundary, increasing this parameter enables obtaining a complete output result. However, it is important to 
% note that this may affect the processing effect of data with diffusion coefficient values that are close to each other.
max_D_DRILT = 14*(0.8*scale/max(b));

% Get the result
DRILT_result_GSP = double(py.DRILT.test(scale, py.numpy.array(S), py.numpy.array(b), model_path))';

%% Plot

% ploting parameters
dc1 = 1;  % dc1, dc2: Display range of diffusion coefficient dimension              
dc2 = 5.5;
cs1 = ppm(idx_peaks(1)-fix(length(ppm)/100));  % cs1, cs2: Display range of chemical shift dimension
cs2 = ppm(idx_peaks(end)+fix(length(ppm)/100));

figure()
t = tiledlayout(9, 1);
t.TileSpacing = 'tight';
t.Padding = 'tight';

% Sparsity
SCORE_result_GSP(SCORE_result_GSP < repmat((max(SCORE_result_GSP, [], 1) * 0.7), [140, 1])) = 0;
DECRA_result_GSP(DECRA_result_GSP < repmat((max(DECRA_result_GSP, [], 1) * 0.7), [140, 1])) = 0;
CONTIN_result_GSP(CONTIN_result_GSP < repmat((max(CONTIN_result_GSP, [], 1) * 0.7), [20, 1])) = 0;
DRILT_result_GSP(DRILT_result_GSP < repmat((max(DRILT_result_GSP, [], 1) * 0.7), [140, 1])) = 0;

nexttile()  %1H NMR
cs_spec = zeros([length(ppm), 1]);
cs_spec(idx_peaks, :) = S(:, 1);
plot(ppm,cs_spec, "Color",'k');set(gca,'Xdir','reverse');axis off;
xlim([cs1,cs2]);

nexttile(2, [2, 1])
DiffCoef = [2.22, 3.03, 4.17];
result_plot(DRILT_result_GSP, DiffCoef, ppm, idx_peaks, t, dc1, dc2, cs1, cs2,max_D_DRILT);
set(gca,'xtick',[])

nexttile(4, [2, 1])
DiffCoef = [2.06, 2.97, 4.71];
result_plot(SCORE_result_GSP, DiffCoef, ppm, idx_peaks, t, dc1, dc2, cs1, cs2, 14);
set(gca,'xtick',[])


nexttile(6, [2, 1])
DiffCoef = [1.95, 2.76, 4.68];
result_plot(DECRA_result_GSP, DiffCoef, ppm, idx_peaks, t, dc1, dc2, cs1, cs2, 14);
set(gca,'xtick',[])

nexttile(8, [2, 1])
DiffCoef = [2.2, 3.3, 4.0];
result_plot(CONTIN_result_GSP, DiffCoef, ppm, idx_peaks, t,dc1, dc2, cs1, cs2, 14);

xlabel(t, 'Chemical Shift(ppm)');
ylabel(t, 'Diffusion Coefficien"t(10^{-10}m^2/s)');

function result_plot(result_GSP, DiffCoef, ppm, idx_peaks, t, dc1, dc2, cs1, cs2, max_D)

    contour_level = 40; % Determines the number and positions of the contour lines / regions.

    spec_whole = zeros([length(result_GSP(:, 1)), length(ppm)]);
    spec_whole(:, idx_peaks) = result_GSP;
    decay_range = linspace(0, max_D, length(result_GSP(:, 1)));
    contour(ppm,decay_range, real(spec_whole), contour_level);
    
    set(gca,'Ydir','reverse','Xdir','reverse'); 
    
    xlim([cs1,cs2]);
    ylim([dc1,dc2]);
    
    text(5.4, 1.8, 'a', FontSize=12)
    text(5.0, DiffCoef(1), "PEG600", FontSize=8, Color=[1, 0.07, 0.07])
    text(5.0, DiffCoef(2), "sucrose", FontSize=8, Color=[1, 0.07, 0.07])
    text(5.0, DiffCoef(3), "glucose", FontSize=8, Color=[1, 0.07, 0.07])
    
    for i = 1:length(DiffCoef)
        l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
        uistack(l, "bottom");
    end
    set(gca,'YTick',unique(DiffCoef) );
end
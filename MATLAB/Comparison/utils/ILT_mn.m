function [ILTdata]=ILT_mn(pfgnmrdata,PeakPick,thresh,Specrange,diffrange, ILTOpts,nugflag, nugc)
%  [ILTdata]=ILT_mn(pfgnmrdata,PeakPick,thresh,Specrange,diffrange, ILTOpts,nugflag, nugc)
%   ILT DOSY (diffusion-ordered spectroscopy) fitting of
%   PFG-NMR diffusion data (aka DOSY data)
%   
%   -------------------------------INPUT--------------------------------------
%   pfgnmrdata      Data structure containing the PFGNMR experiment.
%                   containing the following members:
%                       filename: the name and path of the original file
%                       np: number of (real) data points in each spectrum
%                       wp: width of spectrum (ppm)
%                       sp: start of spectrum (ppm)
%                       dosyconstant: gamma.^2*delts^2*DELTAprime
%                       Gzlvl: vector of gradient amplitudes (T/m)
%                       ngrad: number of gradient levels
%                       Ppmscale: scale for plotting the spectrum
%                       SPECTRA: spectra (processed)
%   ---------------------------INPUT - OPTIONAL-------------------------------
%   thresh =        Noise threshold. Any data point below won't be used in
%                   the fitting. Expressed as % of highest peak. Default 5
%                   (%).
%   Specrange       The spectral range (in ppm) in which the score fitting
%                   will be performed. if set to [0] defaults will be used.
%                   DEFAULT is [sp wp+sp];
%   Diffrange =     The range in diffusion coeeficients that will be
%                   calculated AND the number of data points in the
%                   diffusion dimension (fn1). DEFAULT is [0 20 128];
%                     Diffrange(1) - Min diffusion coefficient plotted.
%                              (0)   DEFAULT. (*10e-10 m2/s)
%                     Diffrange(2) - Max diffusion coefficient plotted.
%                              (20)  DEFAULT  (*10e-10 m2/s)
%                     Diffrange(3) - Number of datapoints in the diffusion
%                                    dimension
%    
%   Peak picking
%                    (0) DEFAULT does peak picking;
%                    (1) fits each frequency individually
%                         
%
%   ILTOpts =       
%
%                    ILTOpts(1) - Regularisation Method
%                    (0) No regularisation;
%                    (1) Tickonov 
%
%                    ILTOpts(2) - Regularisation Parameter Optimisation
%                    (0) Manual
%                    (1) L-curve
%                    (2) GCV
%
%                    ILTOpts(3) - Smoothing
%                    (0) No smoothing
%                    (1) First derivative smoothing
%                    (2) Second derivative smoothing
%
%                    ILTOpts(4) - Constraint
%                    (0) No constraint
%                    (1) Non negativity
%
%   nugflag          (0) pure exponential
%                    (1) NUG correction using coefficients in 'nugc'
%
%
%   nugc =          Vector contaning coefficients for non-uniform field
%                   gradient correction. If not supplied default values
%                   from a Varian ID 5mm probe (Manchester 2006) will be
%                   used.
%
%   --------------------------------------OUTPUT------------------------------
%   ILTdata        Structure containg the data obtained after dosy with
%                  the follwing elements:
%                     
%                  
%                  fit_time:   The time it took to run the fitting
%                  Gzlvl: vector of gradient amplitudes (T/m)
%                  wp: width of spectrum (ppm)
%                  sp: start of spectrum (ppm)
%                  Ppmscale: scale for plotting the spectrum
%                  filename: the name and path of the original file
%                  Options: options used for the fitting
%                  Nug: coeffients for the non-uniform field gradient
%                       compensation
%                  FITSTATS: amplitude and diffusion coeffcient with 
%                            standard error for the fitted peaks
%                  FITTED: the fitted decys for each peak
%                  freqs: frequencies of the fitted peaks (ppm)
%                  RESIDUAL: difference between the fitted decay and the
%                            original (experimental) for all peaks
%                  ORIGINAL: experimental deacy of all fitted peaks
%                  type: type of data (i.e. dosy or something else)
%                  DOSY: matrix containing the DOSY plot
%                  Spectrum: the least attenuated spectrum      
%                  fn1: number of data points (digitisation) in the
%                       diffusion dimension
%                  dmin: plot limit in the diffusion dimension 10-10 m^2/s)
%                  dmax: plot limit in the diffusion dimension 10-10 m^2/s)
%                  Dscale: Plot scale for the diffusion dimension       
%                  threshold: theshold (%) over which peaks are included in
%                             the fit
%                   
%   Example:  
%
%   See also: dosy_mn, score_mn, decra_mn, mcr_mn, varianimport, 
%             brukerimport, jeolimport, peakpick_mn, dosyplot_mn, 
%             dosyresidual, dosyplot_gui, scoreplot_mn, decraplot_mn,
%             mcrplot_mn
%             
%   References:
%   D. Iain J., Journal of Magnetic Resonance 211 (2011) 178-185.
%   P. Hansen, Numerical Algorithms 46 (2007) 189-194.


%   
%   
%   This is a part of the GNAT        
%   Copyright  2017  <Mathias Nilsson>
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License along
%   with this program; if not, write to the Free Software Foundation, Inc.,
%   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
%   Dr. Mathias Nilsson
%   School of Chemistry, University of Manchester,
%   Oxford Road, Manchester M13 9PL, UK
%   Telephone: +44 (0) 161 306 4465
%   Fax: +44 (0)161 275 4598
%   mathias.nilsson@manchester.ac.uk

% MN      12 July 2017      Fixed correct integrals.

th=thresh;
t_start=cputime;
RegularisationMethod=ILTOpts(1);
OptimLambda=ILTOpts(2);
Smooth=ILTOpts(3);
Constraint=ILTOpts(4);
RegLambda=ILTOpts(5);

dmin=diffrange(1);
dmax=diffrange(2);
dres=diffrange(3);




% if Specrange==0
%     %do nothing
% else
%     if length(Specrange)~=2
%         error('DOSY: Specrange should have excatly 2 elements')
%     end
%     if Specrange(1)<pfgnmrdata.sp
%         disp('DOSY: Specrange(1) is too low. The minumum will be used')
%         Specrange(1)=pfgnmrdata.sp;
%     end
%     if Specrange(2)>(pfgnmrdata.wp+pfgnmrdata.sp)
%         disp('DOSY: Specrange(2) is too high. The maximum will be used')
%         Specrange(2)=pfgnmrdata.wp+pfgnmrdata.sp;
%     end
%     for k=1:length(pfgnmrdata.Ppmscale)
%         if (pfgnmrdata.Ppmscale(k)>Specrange(1))
%             begin=k-1;
%             k1=begin;
%             break;
%         end
%     end
%     
%     for k=begin:length(pfgnmrdata.Ppmscale)
%         if (pfgnmrdata.Ppmscale(k)>=Specrange(2))
%             endrange=k;
%             break;
%         end
%     end
%     %make a new stucture
%     pfgnmrdata.sp=pfgnmrdata.Ppmscale(k1);
%     pfgnmrdata.wp=pfgnmrdata.Ppmscale(endrange)-pfgnmrdata.Ppmscale(k1);
%     pfgnmrdata.Ppmscale=pfgnmrdata.Ppmscale(k1:endrange);
%     pfgnmrdata.SPECTRA=pfgnmrdata.SPECTRA(k1:endrange,:);
%     pfgnmrdata.np=length(pfgnmrdata.Ppmscale);
% end


pfgnmrdata.Gzlvl=pfgnmrdata.Gzlvl.^2;
[nspec,~]=size(pfgnmrdata.SPECTRA);

DOSY=zeros(nspec,dres);
if PeakPick==0 %Use peak-picking
    [peak]=peakpick_mn(pfgnmrdata.SPECTRA(:,1),th);
    hp = waitbar(0,'ILT: Fitting after peak picking');
    [npeaks]=length(peak);
    disp(['ILT: Fitting on: ' num2str(npeaks) ' peaks'])
elseif (PeakPick==1)
    peak(pfgnmrdata.np)=struct('max',0,'start',0,'stop',0);
    %npeaks=nspec;
    npeaks=0;
    for k=1:pfgnmrdata.np
        if pfgnmrdata.SPECTRA(k,1)>max(pfgnmrdata.SPECTRA(:,1))*th/100
            npeaks=npeaks+1;
            peak(npeaks).max=k;
            peak(npeaks).start=k;
            peak(npeaks).stop=k;
        end
    end
    hp = waitbar(0,'ILT: Fitting each frequency');
    elseif (PeakPick==2) %use integrals.

    peak=pfgnmrdata.integral_peaks;
    [npeaks]=numel(peak);
    for k=1:npeaks        
        [~, Index]=min(abs(pfgnmrdata.Ppmscale - peak(k).max));
        peak(k).max=Index;
        
        [~, Index]=min(abs(pfgnmrdata.Ppmscale - peak(k).start));
        peak(k).start=Index;
        
        [~, Index]=min(abs(pfgnmrdata.Ppmscale - peak(k).stop));
        peak(k).stop=Index;
        
        
    end
    disp('Hi')
    disp(['ILT: Fitting on: ' num2str(npeaks) ' integral regions'])
    hp = waitbar(0,'ILT: Fitting integrals');
else
    error('ILT: Illegal PeakPick')
end


dline=linspace(dmin,dmax,dres)*1e-10;

grad2=pfgnmrdata.Gzlvl*pfgnmrdata.dosyconstant;

if nugflag==0
    %grad2 is fine
else
    %for NUG
    grad2=nugc(1)*(1e-10*grad2).^1+...
        nugc(2)*(1e-10*grad2).^2+...
        nugc(3)*(1e-10*grad2).^3+...
        nugc(4)*(1e-10*grad2).^4;
    grad2=grad2/1e-10;
end



[dgrid,ggrid]=meshgrid(dline,grad2);



stmat=ggrid.*dgrid;
%Add NUG
exp(-stmat);

STmatrix=exp(-ggrid.*dgrid);

ILTdata.RESIDUAL=zeros(npeaks,length(pfgnmrdata.Gzlvl));
ILTdata.FITTED=zeros(npeaks,length(pfgnmrdata.Gzlvl));
ILTdata.ORIGINAL=zeros(npeaks,length(pfgnmrdata.Gzlvl));
[~,n]=size(STmatrix);
switch Smooth
    case 0 %none
        L=[];
        b0=[];
    case 1 %1st derivative
        L=get_l(n,1);
        b0=zeros(n-1,1);
    case 2 %2nd derivative
        L=get_l(n,2);
        b0=zeros(n-2,1);
    otherwise
end

if issparse(L)
    L=full(L);
    disp('sparse')
end




for i=1:npeaks
    disp(['DOSY: Peak: ' num2str(i) ])
    waitbar(i/npeaks);
    [U,sigma,~] = csvd(STmatrix);
    decay= pfgnmrdata.SPECTRA(peak(i).max,:)';
   
    switch OptimLambda
        case 0 %manual
            RegLambda=ILTOpts(5);
        case 1 %l-curve
            [RegLambda,~,~,~] = l_curve(U,sigma,decay,'Tikh');
        case 2 %gcv
            [RegLambda,~,~] = gcv(U,sigma,decay,'Tikh');
        otherwise
    end

    DiffusionMatrix=cat(1,STmatrix,RegLambda*L);
    DecayArray=cat(1,decay,b0);
    
    switch RegularisationMethod
        case 0  %none
            switch Constraint
                case 0 %none
                    DiffusionMatrix=cat(1,STmatrix,RegLambda*L);
                    DiffusionDistribution=pinv(DiffusionMatrix)*DecayArray;
                case 1 %non-negativity
                    DiffusionMatrix=cat(1,STmatrix,RegLambda*L);
                    DiffusionDistribution=lsqnonneg(DiffusionMatrix,DecayArray);
                otherwise
            end
        case 1 %Tikhonov
            switch Constraint
                case 0 %none
                    DiffusionMatrix=cat(1,STmatrix,L);
                    [U,sigma,V] = csvd(DiffusionMatrix);
                    [DiffusionDistribution,~,~] = tikhonov(U,sigma,V,DecayArray,RegLambda);
                    
                case 1 %non-negativity
                    DiffusionDistribution=lsqnonneg(DiffusionMatrix,DecayArray);
                otherwise
            end
        otherwise
    end
    
    DOSY(peak(i).start:peak(i).stop,:)=pfgnmrdata.SPECTRA(peak(i).start:peak(i).stop,1)*DiffusionDistribution';
    ILTdata.freqs(i)=pfgnmrdata.Ppmscale(peak(i).max);
    ILTdata.ORIGINAL(i,:)=pfgnmrdata.SPECTRA(peak(i).max,:);
    for k=1:dres
        ILTdata.FITTED(i,:)=ILTdata.FITTED(i,:) + DiffusionDistribution(k)*exp(-dline(k) *grad2);
    end
    ILTdata.RESIDUAL(i,:)=ILTdata.ORIGINAL(i,:) - ILTdata.FITTED(i,:);
    ILTdata.FITSTATS(i,:)=[0 0 ];
    
end

endt=cputime-t_start;
h=fix(endt/3600);
m=fix((endt-3600*h)/60);
s=endt-3600*h - 60*m;
fit_time=[h m s];
disp(['ILT: Fitting time was: ' num2str(fit_time(1)) ' h  ' num2str(fit_time(2)) ' m and ' num2str(fit_time(3),2) ' s']);

ILTdata.type='ILT';
ILTdata.fit_time=fit_time;
ILTdata.DOSY=DOSY;
ILTdata.Gzlvl=pfgnmrdata.Gzlvl;
ILTdata.Spectrum=pfgnmrdata.SPECTRA(:,1);
ILTdata.Ppmscale=pfgnmrdata.Ppmscale';
% ILTdata.sp=pfgnmrdata.sp;
% ILTdata.wp=pfgnmrdata.wp;
ILTdata.np=pfgnmrdata.np;
ILTdata.fn1=dres;
ILTdata.dmin=dmin;
ILTdata.dmax=dmax;
ILTdata.Dscale=linspace(dmin,dmax,dres);
ILTdata.Options=[];
ILTdata.threshold=th;
ILTdata.Nug=[];
%ILTdata
close(hp);
% dosyplot_gui(ILTdata,5,th);



end
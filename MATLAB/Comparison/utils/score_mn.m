function [scoredata]=score_mn(pfgnmrdata,ncomp,Specrange,Options,Nug,Fixed,Dinit)
%   [scoredata]=score_mn(pfgnmrdata,ncomp,Specrange,Options,Nug,Fixed,Dinit)
%   SCORE (Speedy Component Resolution) fitting of
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
%
%   ncomp           Number of componets (spectra) to fit.
%
%
%   ---------------------------INPUT - OPTIONAL-------------------------------
%
%   Specrange       The spectral range (in ppm) in which the score fitting
%                   will be performed. if set to [0] defaults will be used.
%                   DEFAULT is [sp wp+sp];
%   Options         Optional parameters%
%                   Options(1): Method to obtain the initial guesses
%                       0 (DEFAULT) obtained by monoexponential fitting of
%                         the summed amplitudes. Centered on the fitted
%                         value and differing by a factor of 1.5
%                       1 random values between 0 and 25 (e-10 m^2/s).
%                       2 User supplied in the vector Dinit.%
%                   Options(2): Constraints on the spectra
%                       0 No constraint (DEFAULT)
%                       1 Non negativity
%                       2 Non negativity - fast algorithm. Requires the
%                         N-way Toolbox by R. Bro (http://www.models.kvl.dk)
%                   Options(3): Shape of the decay
%                       0 (DEFAULT) pure exponetial
%                       1 correcting for non-uniform field gradients.
%                         Supply your own coefficients or use the default
%                         (Varian ID 5mm probe, Manchester 2006). Presently
%                         a maximum of 4 coefficients is supported.
%                   Options(4): plotting
%                        0 (DEFAULT) plot the result using scoreplot
%                        1 Plot the results + diagnostics
%                        2 no plot 
%                   Options(5): Show iteration info
%                        0 (DEFAULT) show info for each iteration
%                        1  show info only after completed fit
%                        2  no info%
%   Nug           array contaning coefficients for non-uniform field
%                 gradient. Set nug=[] for default values.%
%   Fixed         vector containg diffusion coefficients kept fixed
%                 (e.g [1.22 4.51]). Not included in the number of
%                 components (ncomp). An empty vector [] fixes no values.%
%   Dinit        Starting guesses for diffusion coefficients in 1e-10 m2/s
%                 Will only be used in Options(1)=2.
%
%   --------------------------------------OUTPUT------------------------------
%   scoredata      Structure containg the data obtained after score with
%                  the follwing elements:
%
%                  COMPONENTS: The fitted spectra. In the optimal case
%                              representing the true spectra of the
%                              indiviual componets in the mixture
%                  DECAYS:     The decys of the fitted spectra
%                  fit_time:   The time it took to run the fitting
%                  Gzlvl: vector of gradient amplitudes (T/m)
%                  wp: width of spectrum (ppm)
%                  sp: start of spectrum (ppm)
%                  Ppmscale: scale for plotting the spectrum
%                  filename: the name and path of the original file
%                  Dval: fitted diffusion coeffcients (10^-10 m2s-1)
%                  ncomp: number of fitted components
%                  relp: relavtive proportion of the fitted components
%                  Options: options used for the fitting
%                  nug: coeffients for the non-uniform field gradient
%                       compensation
%                  resid: residual sum of squares
%                  nfix: number of componets with a fixed diffusion
%                        coeffcient
%                  Dfix: fixed diffusion coefficients
%
%   Example:
%
%   See also: dosy_mn, score_mn, decra_mn, mcr_mn, varianimport,
%             brukerimport, jeolimport, peakpick_mn, dosyplot_mn,
%             dosyresidual, dosyplot_gui, scoreplot_mn, decraplot_mn,
%             mcrplot_mn
%
%   References:
%   Speedy Component Resolution: An Improved Tool for Processing
%   Diffusion-Ordered SPectroscopy Data. Anal. Chem., 2008 (In Press)
%
%   This is a part of the GNAT
%   Copyright  2017  <Mathias Nilsson>%
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
%
% SCORE processing of
% PFG-NMR diffusion data (aka DOSY data)
%
% see also 'core_mn,decra_mn,dosy_mn and mcr_mn'
%
%     score_mn.m SCORE processing of NMR data
%     Copyright (C) <2007>  <Mathias Nilsson>
%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
%     Dr. Mathias Nilsson
%     School of Chemistry, University of Manchester,
%     Oxford Road, Manchester M13 9PL, UK
%     Telephone: +44 (0) 161 275 4668
%     Fax: +44 (0)161 275 4598
%     mathias.nilsson@manchester.ac.uk


tic;
%% Initial check
if nargin==0
    disp(' ')
    disp(' SCORE ')
    disp(' ')
    disp(' Type <<help score_mn>> for more info')
    return
elseif nargin<2
    error(' The inputs pfgnmrdata and ncomp must be given')
elseif nargin >7
    error('score_mn takes a maximum of 7 arguments')
end

%% Defaults and initialisations
Opts=[0 0 0 0 0 0 0 0 0];
diffcoef=0;
resid=[];
dscale=1e-10;
expfactor=pfgnmrdata.dosyconstant*dscale;
nugc=[9.280636e-1 -9.789118e-3 -3.834212e-4 2.51367e-5]; % Default coefficients (Varian ID 5mm probe, Manchester 2006)
fixed=[];

% if nargin>2
%     if Specrange==0
%         %do nothing
%     else
% %         if length(Specrange)~=3 && length(Specrange)~=2
% %             error('SCORE: Specrange should have 2 or 3 elements')
% %         end
% %         if Specrange(1)<pfgnmrdata.sp
% %             disp('SCORE: Specrange(1) is too low. The minumum will be used')
% %             Specrange(1)=pfgnmrdata.sp;
% %         end
% %         if Specrange(2)>(pfgnmrdata.wp+pfgnmrdata.sp)
% %             disp('SCORE: Specrange(2) is too high. The maximum will be used')
% %             Specrange(2)=pfgnmrdata.wp+pfgnmrdata.sp;
% %         end
% %         for k=1:length(pfgnmrdata.Ppmscale)
% %             if (pfgnmrdata.Ppmscale(k)>Specrange(1))
% %                 begin=k-1;
% %                 break;
% %             end
% %         end
% %         
% %         for k=begin:length(pfgnmrdata.Ppmscale)
% %             if (pfgnmrdata.Ppmscale(k)>=Specrange(2))
% %                 endrange=k;
% %                 break;
% %             end
% %         end
%         %         begin
%         %         endrange
%         %make a new stucture
% %         pfgnmrdata.sp=pfgnmrdata.Ppmscale(begin);
% %         pfgnmrdata.wp=pfgnmrdata.Ppmscale(endrange)-pfgnmrdata.Ppmscale(begin);
% %         load('QGC_net_input.mat', 'ppm')
%         load('QGC_net_input_interpolated.mat', 'ppm')
% 
%         pfgnmrdata.Ppmscale = ppm;
% %         pfgnmrdata.Ppmscale=pfgnmrdata.Ppmscale(begin:endrange);
%         pfgnmrdata.SPECTRA=pfgnmrdata.SPECTRA(1:16384,:);
%         pfgnmrdata.np=length(pfgnmrdata.Ppmscale);
%     end
% end
%% Options
tic

Fixed=Fixed';                                                           % To be able to use 'comma' between fixed D-values (value1,value2) in the GNAT GUI

if nargin>3
    %user Options
    if length(Options)<6
        error('SCORE: Options is a vector of length 6')
    end
    Opts=Options;
end

% For running SCORE in 3D datasets
% gzlvltemp=pfgnmrdata.Gzlvl;
% spectratemp=pfgnmrdata.SPECTRA;
ngrad=pfgnmrdata.ngrad;
dosyrange=[fix((pfgnmrdata.flipnr-1)/ngrad)*ngrad+1:fix((pfgnmrdata.flipnr-1)/ngrad)*ngrad+ngrad];
pfgnmrdata.SPECTRA=pfgnmrdata.SPECTRA(:,dosyrange);
if length(dosyrange)~=length(pfgnmrdata.Gzlvl) %Added to avoid problems with DRONE experiments
    pfgnmrdata.Gzlvl=pfgnmrdata.Gzlvl(dosyrange);
else
    pfgnmrdata.Gzlvl=pfgnmrdata.Gzlvl;
end
                                                                                   
switch Opts(1)

    case 0
        disp('SCORE: Using monoexponential fitting for initialisation of difussion coefficients')              
        opts_fmin=optimset('fminsearch');
        opts_fmin=optimset(opts_fmin,'Display','Final');
        opts_fmin=optimset(opts_fmin,'MaxFunEvals',10000);
%         
%         size(pfgnmrdata.SPECTRA)
%         size(pfgnmrdata.Gzlvl)
%         
        
        tmp=sum(pfgnmrdata.SPECTRA);
        fitparm(1)=tmp(1);
        fitparm(2)= 4;
        fitparm=real(fitparm);
        [fitd]=...
            fminsearch(@pureinit2,fitparm,opts_fmin,pfgnmrdata.Gzlvl.^2,expfactor,...
            sum(pfgnmrdata.SPECTRA));
        
        
        dcenter=fitd(2);
        nFix=length(Fixed);
        if mod(ncomp-nFix,2)==0 %even number
            dval=linspace(dcenter/((ncomp-nFix)/2*1.5),dcenter*((ncomp-nFix)/2*1.5),(ncomp-nFix));                            % To have flexible code, i.e. to be able to run with or without fixed D-value(s), therefore are the (-nFix) parts in the code
        elseif (ncomp-nFix) > 2
            dval=linspace(dcenter/(((ncomp-1)-nFix)/2*1.5),dcenter*(((ncomp-1)-nFix)/2*1.5),(ncomp-nFix));
        else
            dval=dcenter;
        end
        
    case 1
        disp('Random values centred on mono-exponential fitting')              
        opts_fmin=optimset('fminsearch');
        opts_fmin=optimset(opts_fmin,'Display','Final');
        opts_fmin=optimset(opts_fmin,'MaxFunEvals',10000);
%         
%         size(pfgnmrdata.SPECTRA)
%         size(pfgnmrdata.Gzlvl)
%         
        
        tmp=sum(pfgnmrdata.SPECTRA);
        fitparm(1)=tmp(1);
        fitparm(2)= 4;
        fitparm=real(fitparm);
        [fitd]=...
            fminsearch(@pureinit2,fitparm,opts_fmin,pfgnmrdata.Gzlvl.^2,expfactor,...
            sum(pfgnmrdata.SPECTRA));
        
        
        dcenter=fitd(2);

        nguesses=Opts(7);                                                                 % Implementation for the GNAT script
        nFix=length(Fixed);
        sigma=Opts(9);                                                                    % Implementation for the GNAT script
        dval=abs(normrnd(dcenter,sigma,[nguesses ncomp-nFix]));                           % Make normal distribution of random numbers centred around D-value from mono-exponential fitting, with only accepting positive values (abs)
        if nFix==ncomp
            
        else
            figure
            histogram(dval);
            xlabel('diffusion coefficients IN');
            ylabel('occurence');
            title('Normally distributed random numbers')
        end

    case 2
        disp('SCORE: Using user supplied values for initialisation of difussion coefficients')
        if nargin<7
            error('SCORE: No user supplied values supplied in Dinit')
        end
        
        dval=Dinit;
    otherwise
        error('SCORE: Illegal Options(1)')
end
switch Opts(2)                                                              % The non-negativity constraint, having (0) without and (1) with this constraint
    
    case 0
        disp('SCORE: No constraints')
        nnegflag=0;
    case 1
        disp('SCORE: Non negativity constraint')
        nnegflag=1;
    case 2
        if (exist('fastnnls')==2) %#ok
            disp('SCORE: Non negativity constraint (fast algorithm)')
            nnegflag=2;
        else
            disp('SCORE: The N-way Toolbox by R. Bro (http://www.models.kvl.dk) is needed for the fast algorithm')
            disp('SCORE: Using the MATLAB standard instead')
            nnegflag=1;
        end
    case 3
        disp(['SCORE: Constraining D from 0 to ' num2str(Specrange(3))])
        nnegflag=0;
    otherwise
        error('SCORE: Illegal Options(2)')
end
switch Opts(3)                                                              % Third option regarding having an exponential fit or an NUG
    
    case 0
        disp('SCORE: Decays as pure exponentials')
        model=@purexp;
        
    case 1
        disp('SCORE: Decays as power series of exponetials (NUG)')
        if nargin<5
            disp('SCORE: Using default NUG coefficients')
            nugc=[9.280636e-1 -9.789118e-3 -3.834212e-4 2.51367e-5]; % Default coefficients (Varian ID 5mm probe, Manchester 2006)
        end
        model=@nugexp;
        
    otherwise
        error('SCORE: Illegal Options(3)')
end
switch Opts(4)                                                              % Plotting option (diagnostics)
    
    case 0
        disp('SCORE: Plotting data after fit')
        
    case 1
        disp('SCORE: Plotting results and diagnostics')
        
    case 2
        disp('SCORE: No plotting')
        
        
    otherwise
        error('SCORE: Illegal Options(4)')
end
switch Opts(5)                                                              % This option, whether to show the iteration information (regarding the fminsearch algorihm?)
    
    case 0
        disp('SCORE: Information after each iteration of the simplex algorithm')
        
    case 1
        disp('SCORE: Information at the end of fit (simplex)')
        
    case 2
        disp('SCORE: No information of the simplex fit')
        
    otherwise
        error('SCORE: Illegal Options(5)')
end
switch Opts(6)                                                              % This option regarding chosing between SCORE (residual min, as inner loop) and OUTSCORE (cross-talk min, as inner loop)
    
    case 0
        disp('SCORE: Using SCORE (residuals minimisation) inner loop')
        
    case 1
        disp('SCORE: Using OUTSCORE (cross-talk minimisation) inner loop')
                
    case 2
        disp('SCORE: Using HYSCORE (Hybrid SCORE-OUTSCORE minimisation) inner loop')
        alfa=Opts(8);                                                                % Implementation for the GNAT script
        
    otherwise
        error('SCORE: Illegal Options(6)')
end

if nargin>4                                                                 % What is nargin?
    if Opts(3)==1
        if ( (exist('Nug','var')==0) || (max(size(Nug)==[0 0])) )
            if isempty(Nug)
                disp('SCORE: Using default NUG coefficients')
            end
            nugc=[9.280636e-1 -9.789118e-3 -3.834212e-4 2.51367e-5]; % Default coefficients (Varian ID 5mm probe, Manchester 2006)
        else
            nugc=zeros(1,4);
            nugc(1:length(nugc))=Nug;
            disp('SCORE: Using user supplied NUG coefficients')
        end
    end
end
if nargin>5
    disp(['SCORE: Using: ' num2str(length(Fixed)) ' fixed value(s)']);
    fixed=Fixed;
end

if nargin>6
    if Opts(2)==2
        if length(Dinit)~=ncomp
            error('SCORE: The number of starting values must be equal to ncomp')
        elseif Opts(1)==2
            dval=Dinit;
        end
    else
        disp('SCORE: Set Options(1)==2 to use user supplied values for starting guesses')
    end
end

Opts(5)=2;
% dval
% dvsl=[4.28 7.29 11.5 15.8]
opts1=optimset('fminsearch');                                               % Option (1) utilizes fminsearch, which is the Nelder-Mead simplex method 
if Opts(5)==0
    opts1=optimset(opts1,'Display','Iter');
elseif Opts(5)==1
    opts1=optimset(opts1,'Display','Final');
elseif Opts(5)==2
    opts1=optimset(opts1,'Display','Notify');
else
    disp('SCORE: Warning! Unknown Options(5) - using default (0)')
    opts1=optimset(opts1,'Display','Iter');
end
opts1=optimset(opts1,'TolX',1e-6);
opts1=optimset(opts1,'TolFun',1e-6);
opts1=optimset(opts1,'MaxIter',10000);
if ncomp==0 && nnegflag~=0
    opts1=optimset(opts1,'MaxIter',1);
end
opts1=optimset(opts1,'MaxFunEvals',10000);

% Options(2)=3;

if Opts(1)==0
    nguesses=1;

end


if ncomp==0
    % no point in minimising
else
    if Options(2)==3                        %AC using fmincon??
        lb=zeros(1,length(fixed)+ncomp);    %AC Set lower and upper bounds??
        ub=ones(1,length(fixed)+ncomp)*Specrange(3);
        for j=1:length(fixed)               %AC do fixed first.??
            lb(j)=fixed(j)-0.5;
            ub(j)=fixed(j)+0.5;
        end
        lb(length(fixed)+1:end)=0;
        
        [diffcoef, resid]=fmincon(@score_inner,dval,[],[],[],[],lb,ub,[],opts1,pfgnmrdata,ncomp,model,nugc,nnegflag,fixed);
    else
        nfix=length(fixed);
        
        if Opts(6)==2
             if nfix==ncomp
                diffcoef1=zeros(ncomp,nfix);                                                % With 'alfa', if the number of components = number of fixed D-values, only use fixed D-values
                resid1=zeros(1,nfix);
                dval=fixed;
                for p=1:nfix
                    [diffcoef1(:,p), resid1(:,p)]=fminsearch(@score_inner,dval(p,:),opts1,pfgnmrdata,ncomp,model,nugc,nnegflag,fixed,Opts(6),alfa);
                end
            else    
                diffcoef1=zeros(ncomp-nfix,nguesses);                                       % 'With alfa', if the number of components do not equal to the number of fixed D-values, but the latter are non-zero, the fit only some of the diffcoef(s)
                resid1=zeros(1,nguesses);
                for p=1:nguesses
                    [diffcoef1(:,p), resid1(:,p)]=fminsearch(@score_inner,dval(p,:),opts1,pfgnmrdata,ncomp,model,nugc,nnegflag,fixed,Opts(6),alfa);
                end
            end 
        else

            if nfix==ncomp
                diffcoef1=zeros(ncomp,nfix);
                resid1=zeros(1,nfix);                                                        % Some as previous code, but without 'alfa', to have flexibility in the code
                dval=fixed;
                for p=1:nfix
                    [diffcoef1(:,p), resid1(:,p)]=fminsearch(@score_inner,dval(p,:),opts1,pfgnmrdata,ncomp,model,nugc,nnegflag,fixed,Opts(6));
                end
            else    
                diffcoef1=zeros(ncomp-nfix,nguesses);
                resid1=zeros(1,nguesses);
                for p=1:nguesses
                    [diffcoef1(:,p), resid1(:,p)]=fminsearch(@score_inner,dval(p,:),opts1,pfgnmrdata,ncomp,model,nugc,nnegflag,fixed,Opts(6));
                end
            end
        end
        
        if nguesses>1
            if nfix==ncomp
                
            else
%                 figure
%                 edges=[0.0043895323632206:0.000000000000000001:0.0043895323632207];
%                 histogram(resid1, edges);
%                 xlabel('residuals');
%                 ylabel('occurence');
%                 title('Residuals');
% 
%                 figure
%                 xtest=1:nguesses;
%                 scatter(xtest,resid1(1,:), 'x');
%                 xlabel('number of guesses');
%                 ylabel('Residuals'); 

%                 figure
%                 edges=[3.535782:0.0000001:3.535787];                      % First D-value
%                 histogram(diffcoef1, edges);
%                 xlabel('diffusion coefficient');
%                 ylabel('occurence'); 
%                 title('Diffusion coefficients OUT');
%                 
%                 figure
%                 edges=[2.00014634:0.0000000001:2.000146358];                        % Second D-value
%                 histogram(diffcoef1, edges);
%                 xlabel('diffusion coefficient');
%                 ylabel('occurence');
%                 title('Diffusion coefficients OUT');

%                 figure
%                 xtest=1:nguesses;
%                 scatter(xtest,diffcoef1(1,:), 'x');
%                 xlabel('number of guesses');
%                 ylabel('diffusion coefficients'); 
            end

        end
            [resid1, index]=min(resid1);                                                % Here choosing only the minimum residual value out of all other residuals obtained by the different number of guesses
            diffcoef=diffcoef1(:,index);

    end
end
nfix=length(fixed); 
if nfix==ncomp
        C=zeros(pfgnmrdata.ngrad,nfix);                                                 % If the number of fixed D-values = number of components, then only use fixed D-values
        S=zeros(nfix,pfgnmrdata.np);
        m=1;
        if fixed>0
            for k=1:(nfix)
                inval(1)=1;
                inval(2)=fixed(m);
                C(:,k)=model(inval,pfgnmrdata.Gzlvl.^2,expfactor,nugc);
                m=m+1;
            end
        end
else
        C=zeros(pfgnmrdata.ngrad,ncomp);                                                % Else, if the number of fixed D-values not equal and non-zero, then use some of the diffcoef(s)
        S=zeros(ncomp,pfgnmrdata.np); 

        for k=1:(ncomp-nfix)
            inval(1)=1; 
            inval(2)=diffcoef(k);                                                       % Where here we use some of the diffcoef(s)
            C(:,k)=model(inval,pfgnmrdata.Gzlvl.^2,expfactor,nugc, pfgnmrdata.b);
        end
        m=1;
        if fixed>0
            for k=((ncomp-nfix)+1):(ncomp)
                inval(1)=1;
                inval(2)=fixed(m);                                                      % And here we use some of the fixed D-value(s)
                C(:,k)=model(inval,pfgnmrdata.Gzlvl.^2,expfactor,nugc);
                m=m+1;
            end
        else
            for k=(ncomp+1):(ncomp+nfix)
                inval(1)=1;                                                             % If the number of fixed D-value(s) is zero, use original code
                inval(2)=fixed(m);
                C(:,k)=model(inval,pfgnmrdata.Gzlvl.^2,expfactor,nugc);
                m=m+1;
            end
        end
end

X=pfgnmrdata.SPECTRA;                                                       % Defining X as the DOSY dataset
if nnegflag==1
    for k=1:pfgnmrdata.np
        S(:,k)=lsqnonneg(C,X(k,:)');                                        % lsqnoneg= non-linear square
    end
elseif nnegflag==2
    for k=1:pfgnmrdata.np
        S(:,k)=fastnnls(C,X(k,:)');
    end
else
    
    S=C\X';                                                                 % Defining spectra (S) by the previously fitted decaying shapes by the ratio to the data-set (X) - How good does (C) describes (X) in order to obtain good (S)
end

endt=toc;
h=fix(endt/3600);                                                           % Defining how long did the calculation take, in hours, minutes and seconds( the following lines)
m=fix((endt-3600*h)/60);
s=endt-3600*h - 60*m;
fit_time=[h m s];
disp(['SCORE: Fitting time was: ' num2str(fit_time(1)) ' h  ' num2str(fit_time(2)) ' m and ' num2str(fit_time(3),2) ' s']);     % And displaying those 'hours, minutes, seconds'.

disp('SCORE: Finding the relative contributions');
%tic
%Sort the components according to size

if diffcoef==0
    diffcoef=[];
end

if nfix==ncomp
    [scoredata.Dval, sortindex]=sort([fixed]);                             % To have flexibilty in the code, with or without fixed D-values
else
    [scoredata.Dval, sortindex]=sort([diffcoef;fixed]); 
end
% assignin('base','diffcoef_sortline', scoredata.Dval)
% assignin('base','fixed', fixed)
 scoredata.DECAYS=C(:,sortindex);
 scoredata.COMPONENTS=S(sortindex,:);
scoredata.RESIDUALS=X-S'*C';
scoredata.filename=pfgnmrdata.filename;
scoredata.Gzlvl=pfgnmrdata.Gzlvl;
% scoredata.COMPONENTS=S;
% scoredata.DECAYS=C;
scoredata.fit_time=fit_time;
scoredata.Ppmscale=pfgnmrdata.Ppmscale;
scoredata.ncomp=ncomp;
scoredata.nfix=nfix;
scoredata.Options=Opts;
% scoredata.Dval=diffcoef;
scoredata.Dfix=fixed;
scoredata.nug=nugc;
scoredata.resid=resid;
%find the relative contributions
relp=sum(abs(scoredata.COMPONENTS),2);
scoredata.relp=relp;
scoredata.np=pfgnmrdata.np;
% scoredata.sp=pfgnmrdata.sp;
% scoredata.wp=pfgnmrdata.wp;
scoredata.ngrad=ngrad;
scoredata.dosyconstant=pfgnmrdata.dosyconstant;
scoredata.SPECTRA=pfgnmrdata.SPECTRA;
% if Opts(4)==0
%     %default: plot result
%     scoredata.diags=0;
%     scoreplot(scoredata);
% elseif Opts(4)==1
%     %plot result + diagnostics
%     scoredata.diags=1;
%     scoreplot(scoredata);
% else
%     %do nothing
% end

toc

end
%-------------------------------AUXILLIARY FUNCTIONS-----------------------

function y = pureinit(a,xdata,expfactor)
% for fitting a  pure exponential
y =  a(1)*exp(-a(2)*expfactor.*xdata);
end

function sse = pureinit2(a,xdata,expfactor,ydata)
% for fitting a  pure exponential
y =  a(1)*exp(-a(2)*expfactor.*xdata);
sse=ydata-y;
sse=sum(sse.^2,2);
end

function [resid]=score_inner(Dval,pfgnmrdata,ncomp,model,nugc,nnegflag,fixed,OUTSCORE,alfa)     % Whether to use OUTSCORE?
nfix=length(fixed);                                                         % Number of fixed diffusion coefficients
dscale=1e-10;                                                               % The scale of D-values (in 1x10^-10)
expfactor=pfgnmrdata.dosyconstant*dscale;                                   % dosyconstant = gamma.^2*delts^2*DELTAprime - so expfactor is this dosyconstant times the 1x10^10 factor?
C=zeros(pfgnmrdata.ngrad,ncomp);                                            % Making C-matrix of zeros (to pre-define the matrix, easier for calculations), defined as 'the number of gradient levels (e.g. 16) as rows and 'the number of components to fit + number of fixed D-values' as columns
S=zeros(ncomp+nfix,pfgnmrdata.np);                                          % Making S-matrix of zeros (defined as 'the number of components to fit + number of fixed D-values' as row  and 'the number of (real) data points in each spectrum' as the columns.
X=pfgnmrdata.SPECTRA;                                                       % Making X-matrix as the original DOSY dataset
%nnegflag=3;

nfix=length(fixed);                    
if nfix==ncomp
        C=zeros(pfgnmrdata.ngrad,nfix);               
        S=zeros(nfix,pfgnmrdata.np);                                      
        m=1;
        if fixed>0
            for k=1:(nfix)                                               
                inval(1)=1;
                inval(2)=fixed(m);
                C(:,k)=model(inval,pfgnmrdata.Gzlvl.^2,expfactor,nugc, pfgnmrdata.b);
                m=m+1;
            end
        end
else
        C=zeros(pfgnmrdata.ngrad,ncomp);                                    
        S=zeros(ncomp,pfgnmrdata.np);                                       

        for k=1:(ncomp-nfix)
            inval(1)=1;                                                           
            inval(2)=Dval(k);
            C(:,k)=model(inval,pfgnmrdata.Gzlvl.^2,expfactor,nugc, pfgnmrdata.b);
        end
        m=1;
        if fixed>0
            for k=((ncomp-nfix)+1):(ncomp)                                              
                inval(1)=1;
                inval(2)=fixed(m);
                C(:,k)=model(inval,pfgnmrdata.Gzlvl.^2,expfactor,nugc, pfgnmrdata.b);
                m=m+1;
            end
        else
            for k=(ncomp+1):(ncomp+nfix)                                              
                inval(1)=1;
                inval(2)=fixed(m);
                C(:,k)=model(inval,pfgnmrdata.Gzlvl.^2,expfactor,nugc);
                m=m+1;
            end
        end
end
% C=fliplr(C);

if nnegflag==1
    for k=1:pfgnmrdata.np
        S(:,k)=lsqnonneg(C,X(k,:)');
    end
elseif nnegflag==2
    for k=1:pfgnmrdata.np
        S(:,k)=fastnnls(C,X(k,:)');
    end
elseif nnegflag==3
    %NOTE: CORE STYLE not non negativity. NOT available as an option
    for k=1:pfgnmrdata.np
        S(:,k)=C\X(k,:)';
        %S(:,k)=linsolve(C,X(k,:)');
    end
else
    S=C\X';
end

%% Minimisation Criteria
if ncomp<=1  % Can't OUTSCORE one component, revert to SCORE
    OUTSCORE=0;
end    
Sstart=S;

switch OUTSCORE
    case 1 % area normalisation (OUTSCORE)

        for m=1:ncomp    
            S(m,:)=abs(S(m,:))./sum(abs(S(m,:)));
        end
        CrossTalk=ones(1,pfgnmrdata.np);

        CompPairs=nchoosek(1:ncomp,2);
        [rows, columns]=size(CompPairs);
        for n=1:rows % the sum of pairwise .* of spectra
            CrossTalk=CrossTalk+(S(CompPairs(n,1),:).*S(CompPairs(n,2),:));
        end
        resid = sum(CrossTalk);
        
    case 0 % Residuals Minimisation (SCORE)
        M=C*S;
        resid=sum(sum(((M'-X).^2)));
        
    case 2 % Hybrid SCORE-OUTSCORE minimsation (HYSCORE)
        
        for m=1:ncomp    
            S_end(m,:)=abs(Sstart(m,:))./sum(abs(Sstart(m,:)));
        end
        CrossTalk=ones(1,pfgnmrdata.np);

        CompPairs=nchoosek(1:ncomp,2);
        [rows, columns]=size(CompPairs);
        for n=1:rows % the sum of pairwise .* of spectra
            CrossTalk=CrossTalk+(S_end(CompPairs(n,1),:).*S_end(CompPairs(n,2),:));
        end
        resid_OUTSCORE=sum(CrossTalk);                                                  % Defining residuals only from OUTSCORE

        M=C*Sstart;
        resid_SCORE=sum(sum(((M'-X).^2)));                                              % Defining residuals only from SCORE

        resid=((alfa/100)*resid_SCORE)+(resid_OUTSCORE*((100/100)-(alfa/100)));         % Combining both SCORE/OUTSCORE-residuals with a weighting factor (alfa)
end

end
function y = purexp(a,xdata,expfactor,nugc,b) %#ok
% for fitting a number of pure exponentials
y = a(1)*exp(-a(2)*b);
% y = a(1)*exp(-a(2)*expfactor.*xdata);
end
function y = nugexp(a,xdata,expfactor,nugc)
%fitting non-uniform field gradient decays (NUG)
y = a(1)*exp( - nugc(1)*(a(2)*expfactor.*xdata).^1 -...
    nugc(2)*(a(2)*expfactor.*xdata).^2 -...
    nugc(3)*(a(2)*expfactor.*xdata).^3 -...
    nugc(4)*(a(2)*expfactor.*xdata).^4);
end




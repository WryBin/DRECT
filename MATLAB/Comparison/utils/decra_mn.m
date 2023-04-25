
function [decradata]=decra_mn(pfgnmrdata,ncomp,Specrange,plot)

%   DECRA (Direct Exponential Curve Resolution Algorithm) fitting of
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
%   Specrange       The spectral range (in ppm) in which the decra fitting
%                   will be performed. if set to [0] defaults will be used.
%                   DEFAULT is [sp wp+sp];
%
%   --------------------------------------OUTPUT------------------------------
%   decradata      Structure containg the data obtained after decra with
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
%   Example:  
%
%   See also: dosy_mn, score_mn, decra_mn, mcr_mn, varianimport, 
%             brukerimport, jeolimport, peakpick_mn, dosyplot_mn, 
%             dosyresidual, dosyplot_gui, scoreplot_mn, decraplot_mn,
%             mcrplot_mn
%             
%   Based on the script in Antalek, B. Conc.Magn.Reson. 2002, 14(4),
%   225-258.
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


% Explanation of Variables in the Matlab Script
% pspec                   -resolved spectra normalized to the g = 0 point in the experiment
%                         and representative of composition
% diff                    -diffusion coefficients calculated from the eigenvalues
% a                       -relative amounts of each resolved component scaled to g = 0
% fname                   -name of file (phasefile in Varian VNMR)
% constant                -(DELTA - delta/3)*gamma^2*delta^2 (rad s T-1)^2
% grad2                   -gradient values (T m-1)^2 -NOT gradient squared!
% spec                    -number of spectra in total
% startspec and endspec   -range of data to be used (gradient levels)
% startpoint and endpoint -expansion of data to be used (spectral points)
% ncom                    -number of components chosen
% incr                    -difference in g2 values (T m_1)2

if nargin==0
   disp(' ')
   disp(' DECRA ')
   disp(' ')
   disp(' Type <<help decra_mn>> for more info')
   return
elseif nargin<2
   error(' The inputs pfgnmrdata and ncomp must be given')
end

% if nargin>2
%     if Specrange==0
%         %do nothing
%     else
%         if length(Specrange)~=2
%             error('DECRA: Specrange should have excatly 2 elements')
%         end
%         if Specrange(2)<Specrange(1)
%             error('DECRA: Second element of Specrange must be larger than the first')
%         end
%         if ((Specrange(1)<pfgnmrdata.sp) || Specrange(2)>(pfgnmrdata.sp+pfgnmrdata.wp))
%             error('DECRA: Specrange exceeds the spectral width (sp - sp+wp)i')
%         end
%         for k=1:length(pfgnmrdata.Ppmscale)
%             if (pfgnmrdata.Ppmscale(k)>Specrange(1))
%                 begin=k-1;
%                 break;
%             end
%         end
%         for k=begin:length(pfgnmrdata.Ppmscale)
%             if (pfgnmrdata.Ppmscale(k)>=Specrange(2))
%                 endrange=k;
%                 break;
%             end
%         end
%         %make a new stucture
%         pfgnmrdata.sp=pfgnmrdata.Ppmscale(begin);
%         pfgnmrdata.wp=pfgnmrdata.Ppmscale(endrange)-pfgnmrdata.Ppmscale(begin);
%         pfgnmrdata.Ppmscale=pfgnmrdata.Ppmscale(begin:endrange);
%         pfgnmrdata.SPECTRA=pfgnmrdata.SPECTRA(begin:endrange,:);
%         pfgnmrdata.np=length(pfgnmrdata.Ppmscale);     
%     end
% end

if nargin>3
    %loco plot options
else
    plot=1;
end

tic;
data=pfgnmrdata.SPECTRA';
grad2=pfgnmrdata.Gzlvl;
constant=pfgnmrdata.dosyconstant;
startspec=1;
endspec=pfgnmrdata.ngrad;

% Making sure the data is in the right format
grad2=grad2(:);
grad2=grad2.^2;
% Assuming more spectral points than gradient levels

%a constant difference between g2^2 is necessary
incr=grad2(2)-grad2(1);

%define two ranges (split the data set in two)
range1=(startspec:endspec-1);
range2=(startspec+1:endspec);

%create a common base for the two data sets using SVD
[vc,sc,uc]=svd(data(range1,:)',0);
sc=sc(1:ncomp,1:ncomp);
uc=uc(:,1:ncomp);
vc=vc(:,1:ncomp);

%project the two data sets onto the common base
%auv=sc;
buv=uc'*data(range2,:)*vc;

%solve the generalized eigenvalue problem
[v,s]=eig(buv*inv(sc)); %#ok<MINV>
%ev=diag(s);
[ev,sortindex]=sort(diag(s));
v=v(:,sortindex);

%calculate spectra and concentrations
pspec=pinv(vc*inv(sc)*v); %#ok<MINV>
pint=uc*v;

%scale spectra and concentrations
total=sum(data(range1,:),2);
scalefactor=pint\total;
pint=pint*diag(scalefactor);
pspec=diag(1./scalefactor)*pspec;

%calculate proper composition
pint2=pint*diag(ev);
nrows=size(pint,1);
pintcomb=[pint(1,:);(pint(2:nrows,:)+pint2(1:nrows-1,:))/2; pint2(nrows,:)];
b=log(ev)/incr;
diff=-b/constant;
grad21=grad2(range1);
grad22=grad2(range2);
grad2comb=[grad21;grad22(nrows)];
a=zeros(1,length(ev));
for i=1:length(ev)
    expest0=exp(grad2comb*b(i));
    a(i)=expest0\pintcomb(:,i);
end
pspec=diag(a)*pspec;
pint=pintcomb;
%fit_time=toc;
endt=toc;
h=fix(endt/3600);
m=fix((endt-3600*h)/60);
s=endt-3600*h - 60*m;
fit_time=[h m s];
disp(['DECRA: Fitting time was: ' num2str(fit_time(1)) ' h  ' num2str(fit_time(2)) ' m and ' num2str(fit_time(3),2) ' s']);

[decradata.Dval, sortindex]=sort(diff);
decradata.DECAYS=pint(:,sortindex);
decradata.COMPONENTS=pspec(sortindex,:);
decradata.relp=a(:,sortindex);

if isreal(decradata.COMPONENTS)==0
    disp('COMPLEX')
end



decradata.filename=pfgnmrdata.filename;
decradata.Gzlvl=pfgnmrdata.Gzlvl; 
%decradata.COMPONENTS=pspec;
%decradata.Dval=diff;
%decradata.DECAYS=pint;
decradata.fit_time=fit_time;
% decradata.wp=pfgnmrdata.wp;
% decradata.sp=pfgnmrdata.sp;
decradata.Ppmscale=pfgnmrdata.Ppmscale;
decradata.ncomp=ncomp;
decradata.np=pfgnmrdata.np;

% if plot==1
% decraplot(decradata);
% end
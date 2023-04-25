function [peak]=peakpick_mn(spectrum,threshold)

%   [peak]=peakpick_mn(spectrum,threshold)
%   finds and defines the peaks over a certain threshold.
%   ------------INDATA---------------------------------
%   spectrum     Vector containng the spectrum
%   threshold    Threshold in % below which peaks are ignored
%
%   ------------OUTDATA---------------------------------
%   peak         Structure of peaks containing start, stop and max%
%   Example:
%   See also:
%
%   This is a part of the DOSYToolbox
%   Copyright 2007-2008  <Mathias Nilsson>

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


th=threshold*max(spectrum)/100;


%get all extreme values

[ymax,imax,ymin,imin] = peakfind_mn(spectrum);



% all peak maxima that are positive are positive peaks
% all peak min that are negative are negative peaks
%only use those over the threshold
peakmax=find(ymax>th);
imax=imax(peakmax);
ymax=ymax(peakmax);
up= ymax>0;
uppeaks=imax(up);


peakmin=find(ymin>th);
imin=imin(peakmin);
ymin=ymin(peakmin);
down= ymin<0;
downpeaks=imin(down);

peak(1)=struct('max',[],'start',[],'stop',[]);
% for k=1:length(uppeaks)
%     peak(k).max=uppeaks(k)
% end

%Lets find the boundaries for the positive peaks first

for p=1:length(uppeaks)
    peak(p).max=uppeaks(p);
    peak(p).start=peak(p).max;
    peak(p).stop=peak(p).max;
    for k=(peak(p).stop+1):(length(spectrum)+1)
        if peak(p).start==length(spectrum)
            %we have gone to far so this is end of peak
            %disp('1')
            peak(p).stop=length(spectrum);
            break
        end
        if sum((k+1)==imax)
            %one max after another - a rare even, but lets take the end of
            %peak as halfway between
            peak(p).stop=round((peak(p).max+k-1)/2);
        end
        if sum(k==imin)
            %min point  so this is end of peak          
            peak(p).stop=k-1;
            break
        end
        if  spectrum(k)<=th
            %below threshold  so this is end of peak
            peak(p).stop=k-1;            
            break
        end
        if  spectrum(k)<=0
            %below zero so this is end of peak
            peak(p).stop=k-1;           
            break
        end
        if  spectrum(k)<=0.001*spectrum(peak(p).max)
            %below 0.1% of peak height
            peak(p).stop=k-1;           
            break
        end

    end

    for k=(peak(p).start+1):(length(spectrum)+1)
        peak(p).start=peak(p).start-1;
        %test=peak(p).start
        if peak(p).start==0
            %we have gone to far so this is end of peak           
            peak(p).start=1;
            break
        end
        if sum((peak(p).start)==imax)
            %one max after another - a rare even, but lets take the end of
            %peak as halfway between
            peak(p).start=round((peak(p).max+k-1)/2);          
        end
        if sum(peak(p).start==imin)
            %min point  so this is end of peak
            peak(p).start=peak(p).start+1;                    
            break
        end
        if spectrum(peak(p).start)<=th
            %below threshold  so this is end of peak
            peak(p).start=peak(p).start+1;          
            break
        end
        if spectrum(peak(p).start)<=0
            %below zero so this is end of peak
            peak(p).start=peak(p).start+1;
            break
        end
        if  spectrum(peak(p).start)<=0.001*spectrum(peak(p).max) 
            %below 0.1% of peak height
            peak(p).start=peak(p).start+1; 
            break
        end

    end
end
%Lets find the boundaries for the negative peaks now
q=0;
for p=(length(peak)+1):(length(peak)+length(downpeaks))
    q=q+1;
    peak(p).max=downpeaks(q); 
    peak(p).start=peak(p).max; %#ok<*AGROW>
    peak(p).stop=peak(p).max;
    for k=(peak(p).stop+1):(length(spectrum)+1)
        if peak(p).stop==length(spectrum)
            %we have gone to far so this is end of peak
            peak(p).stop=length(spectrum);
            break
        end
        if k>length(spectrum)
            %we have gone to far so this is end of peak
            peak(p).stop=length(spectrum);
            break
        end
        if sum((k+1)==imin)
            %one max after another - a rare even, but lets take the end of
            %peak as halfway between
            peak(p).stop=round((peak(p).max+k-1)/2);
        end
        if sum(k==imax)
            %max point  so this is end of peak
            peak(p).stop=k-1;
            break
        end
        if  spectrum(k)<=th
            %below threshold  so this is end of peak
            peak(p).stop=k-1;
            break
        end
        if  spectrum(k)>=0
            %above zero so this is end of peak
            peak(p).stop=k-1;
            break
        end
        if  spectrum(k)>=0.001*spectrum(peak(p).max)
            %below 0.1% of peak height
            peak(p).stop=k-1;
            break
        end

    end

    for k=(peak(p).start+1):(length(spectrum)+1)
        peak(p).start=peak(p).start-1;
        %test=peak(p).start
        if peak(p).start==0
            %we have gone to far so this is end of peak
            peak(p).start=1;
            break
        end
        if sum((peak(p).start)==imin)
            %one max after another - a rare even, but lets take the end of
            %peak as halfway between
            peak(p).start=round((peak(p).max+k-1)/2);
        end
        if sum(peak(p).start==imin)
            %max point  so this is end of peak
            peak(p).start=peak(p).start+1;
            break
        end
        if spectrum(peak(p).start)<=th
            %below threshold  so this is end of peak
            peak(p).start=peak(p).start+1;
            break
        end
        if spectrum(peak(p).start)>=0
            %above zero so this is end of peak
            peak(p).start=peak(p).start+1;
            break
        end
        if  spectrum(peak(p).start)>=0.001*spectrum(peak(p).max)
            %below 0.1% of peak height
            peak(p).start=k-1;
            break
        end

    end
   % peak(p)
end



if (isempty(imax))
    peak = [];
end

end

function [maxValues,maxIndex,minValues,minIndex]=peakfind_mn(spectrum)

%  [maxValues,maxIndex,minValues,minIndex]=peakfind_mn(spectrum)
%
%  Finds the max and min points of a spectrum

%NOTE: plateau and end points are not handled at the moment
%

maxIndex=[]; %#ok<NASGU>
minIndex=[]; %#ok<NASGU>

%first derivative
dx=gradient(spectrum);


%min point begins when dx is positive
minP=(dx > 0);
minP= diff(minP);

maxIndex=find(minP == -1) + 1;
minIndex=find(minP == +1) + 1;


%max/min points at the ends
minIndex=[1 ;minIndex ;length(spectrum)];

maxValues=spectrum(maxIndex);
minValues=spectrum(minIndex);

end

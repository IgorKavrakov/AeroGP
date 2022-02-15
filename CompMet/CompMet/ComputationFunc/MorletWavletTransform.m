function [WT, FreqBins, Scales,FreqRes] = MorletWavletTransform(Sig, nLevel, fmax, fmin, fc, sampl,beta,Padding)
% Wavelet Transform based on Morlet Wavelet
% This file is part of CompMet.
% CompMet is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%  CompMet is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with CompMet.  If not, see <https://www.gnu.org/licenses/>.

% Copyright (c) Igor Kavrakov, Ahsan Kareem, Guido Morgenthal 2020

SigLen = length(Sig);                             % obtain the length of the signal
FreqBins = linspace(fmin,fmax,nLevel);            % divide the frequency band into equivalent segments
Scales = fc./ FreqBins;                           % calculate the scale vector           
dt=1/sampl;                                       % calculate the time increment
WT = zeros(nLevel, SigLen);                       % for the storage of wavelet coefficients
FreqRes=1./Scales./2/pi/sqrt(2)./beta;

TimeResolution=Scales./sqrt(2).*beta;
%Padding the signal
if Padding==1
SignalExtension=min(fix(TimeResolution(end)/dt),SigLen);
SignalPadded=[fliplr(Sig(2:SignalExtension)) Sig fliplr(Sig(end-SignalExtension+1:end-1))];
elseif Padding==-1
SignalExtension=min(fix(TimeResolution(end)/dt),SigLen);
SignalPadded=[-fliplr(Sig(2:SignalExtension)) Sig -fliplr(Sig(end-SignalExtension+1:end-1))];
else
SignalPadded=Sig;
SignalExtension=1;
end

for m = 1:nLevel
    a = Scales(m);
    TimeStepsForMorlet=fix(TimeResolution(m)/dt);
    t = -TimeStepsForMorlet*dt:dt:TimeStepsForMorlet*dt;                     
    Morl = pi^(-1/4)*conj(exp(1j*2*pi*fc*(t/a)).*exp(-(t/a).^2/2));               % Morlet wavelet at Scales(m)
    temp = conv(SignalPadded,Morl,'same')/sqrt(a)*dt;%*2*pi;                                     % calculate the wavelet coefficients using convolution  
    WT(m,:) = temp(SignalExtension:SignalExtension+SigLen-1);   % obtain the central elements with a length of SigLen
end
end

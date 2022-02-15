function [f1,WavletAutoBiSpectrum, WavletAutoBiCoherence,WaveletAutoBiNormalziation,NonZeroVals] = WaveletBiCoherence(WT,FreqBins,dt,T,NIntervals,Fact,RandomzFlag)
% Calcuate wavelet bi-coherence
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


%% Input
% WT - Wavlet coefficients;
% f - Frequencies
% t -  Time vector
% T - time for averaging [TStart TEnd]
% NIntervals - Number of intervals to be splitted
% RandomzFlag - If RandomzFlag =0, No phaze randomization will be performed
if nargin==3; T(1)=0; T(2)=dt*size(WT,2);NIntervals=1; Fact=10*pi; RandomzFlag=1; end
if nargin==4; NIntervals=1; Fact=10*pi; RandomzFlag=1; end
if nargin==5; Fact=10*pi; RandomzFlag=1; end
if nargin==6; RandomzFlag=1; end

if size(T)==1; T(1)=0; T(2)=dt*size(WT,2); end


TStart=T(1);TEnd=T(2);
StepStart=max(1,floor(TStart/dt));StepEnd=min(floor(TEnd/dt),size(WT,2));
NSteps=StepEnd-StepStart+1;
WT=WT(:,StepStart:1:StepEnd);
IntervalSteps=floor(NSteps/NIntervals);


BinsConsidered=find(FreqBins<FreqBins(end)/2);
BinsConsidered=BinsConsidered(end);
WavletAutoBiSpectrum=zeros(BinsConsidered,BinsConsidered,NIntervals);
WavletAutoBiCoherence=zeros(BinsConsidered,BinsConsidered,NIntervals);
WavletInt=zeros(BinsConsidered,BinsConsidered);
WavletInt2=zeros(BinsConsidered,BinsConsidered);
WaveletAutoBiNormalziation=zeros(BinsConsidered,BinsConsidered,NIntervals);

Mplyer=1;
f1=FreqBins(1:BinsConsidered); %Equal to f2
indx=hankel([1:BinsConsidered],[BinsConsidered,1:BinsConsidered-1])+1;
indx=reshape(indx,[BinsConsidered.^2,1]);
for j=1:NIntervals

    for k=1:IntervalSteps-1
            CC=reshape(WT(indx,k+(j-1)*IntervalSteps),[BinsConsidered BinsConsidered]);
            CC2=WT(1:BinsConsidered,k+(j-1)*IntervalSteps)*(WT(1:BinsConsidered,k+(j-1)*IntervalSteps).');
            CohOneStep=CC2.*conj(CC);
            if RandomzFlag==1
                Phase=angle(CohOneStep); 
                Random=(rand(size(Phase))-0.5).*Fact*2; 
                CohOneStep=CohOneStep.*exp(1j.*Random.*Phase); 
            end   
            WavletInt=WavletInt+Mplyer.*dt./2.*abs(CC2).^2;
            WavletInt2=WavletInt2+Mplyer.*dt./2.*abs(CC).^2;
            WavletAutoBiSpectrum(:,:,j)=WavletAutoBiSpectrum(:,:,j)+Mplyer.*CohOneStep.*dt./2;
            Mplyer=2;
    end

        Mplyer=1;
        CC=reshape(WT(indx,IntervalSteps+(j-1)*IntervalSteps),[BinsConsidered BinsConsidered]);
        CC2=WT(1:BinsConsidered,IntervalSteps+(j-1)*IntervalSteps)*(WT(1:BinsConsidered,IntervalSteps+(j-1)*IntervalSteps).');
        CohOneStep=CC2.*conj(CC);
        if RandomzFlag==1
            Phase=angle(CohOneStep);
            Random=(rand(size(Phase))-0.5).*Fact*2;
            CohOneStep=CohOneStep.*exp(1j.*Random.*Phase); 
        end        
        WavletInt=WavletInt+Mplyer.*dt./2.*abs(CC2).^2;
        WavletInt2=WavletInt2+Mplyer.*dt./2.*abs(CC).^2;
        
        WaveletAutoBiNormalziation(:,:,j)=WavletInt.*WavletInt2;
        WavletAutoBiSpectrum(:,:,j)=WavletAutoBiSpectrum(:,:,j)+Mplyer.*CohOneStep.*dt./2;
        WavletAutoBiCoherence(:,:,j)=abs(WavletAutoBiSpectrum(:,:,j)).^2./WaveletAutoBiNormalziation(:,:,j);

       %Reset
    WavletInt=zeros(BinsConsidered,BinsConsidered);
    WavletInt2=zeros(BinsConsidered,BinsConsidered);
    % Consider lower triangle only
    % Consider lower triangle only or half of it
    if f1(end)<=(1/dt/4)
    a=ones(BinsConsidered,BinsConsidered);           
    a=triu(a);
    else
    warning('Only half of the inner triangle will be considered');
    a=zeros(BinsConsidered,BinsConsidered);
    LastBin=find(f1>(1/dt/4));LastBin=LastBin(1);
    b=ones(LastBin,LastBin);           
    b=triu(b);a(1:LastBin,1:LastBin)=b;    
    end
    WavletAutoBiSpectrum(:,:,j)=WavletAutoBiSpectrum(:,:,j).*a;
    WavletAutoBiCoherence(:,:,j)=WavletAutoBiCoherence(:,:,j).*a;
    NonZeroVals=find(a~=0);
        
    
end

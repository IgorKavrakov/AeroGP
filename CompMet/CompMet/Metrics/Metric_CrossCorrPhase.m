function [OutStruct1,OutStruct2] = Metric_CrossCorrPhase(X1,X2,dt,TPhase)
% Calcuate phase metric
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


[OutStruct1.Value,lag]=xcorr(X1,X2,'coeff');
OutStruct1.Name='Cross correlation';
[~,I] = max(abs((OutStruct1.Value)));
n_start = lag(I); clear lag;

OutStruct2.Value=exp(-abs(n_start.*dt)./TPhase); %Sprague and Geers
OutStruct2.Name='Phase error (Phi)';


end


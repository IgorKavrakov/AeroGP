function [ OutStruct ] = Metric_RMS(X1,X2)
% Calcuate RMS Metric 
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
OutStruct.Name='RMS Metric (RMS)';
OutStruct.Value=exp(-abs(rms(X1)-rms(X2))./rms(X1));


end


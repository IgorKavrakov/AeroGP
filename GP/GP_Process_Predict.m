function [Mean_Pred,Var,lik] = GP_Process_Predict(x,y,x_predict,hyp,L_Kern)
% Prediction for a GP - similar as GP_Process.m, excluding the computation
% of kernel and likelihood.

% By Igor Kavrakov

%%%%%%%%% COPYRIGHT NOTICE %%%%%%%%% 
%  This file is part of AeroGP.
%  AeroGP is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  AeroGP is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with AeroGP.  If not, see <https://www.gnu.org/licenses/>.

% Copyright (c) Igor Kavrakov, Allan McRobie, Guido Morgenthal 2022

%%
Kern_predict_r=Kern_Exp(x_predict,x,hyp);

alpha=L_Kern'\(L_Kern\y);
Mean_Pred=Kern_predict_r*alpha;

if nargout>1
Kern_predict=Kern_Exp(x_predict,x_predict,hyp); %Get the kernel at prediction points    
v=L_Kern\Kern_predict_r';       %Predictive varianse
Var=Kern_predict-v'*v;          %Predictive variances
lik=-1/2*y'*alpha-sum(log(diag(L_Kern)))-length(y)/2*log(2*pi); %Likelihood   
end

end


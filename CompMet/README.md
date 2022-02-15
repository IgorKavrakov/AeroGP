# CompMet
CompMet: Metrics for Comparison of Time Histories

Cite:
Kavrakov, I., Kareem, A., and Morgenthal, G. 2020. Comparison Metrics for Time-histories: Application to Bridge Aerodynamics. J. Eng. Mech., 146 (9), 040200093. https://doi.org/10.1061/(ASCE)EM.1943-7889.0001811

CompMet is a Matlab-based computer code that computes metrics for comparison of two time-histories.
The files constituting this code are included in CompMet.
Given to time-histories, X1 and X2, the code is simply called as:

[M]=CompMet(X1,X2,Properties);

where, "Properties" is an input structure the Metric Properties.
"M" is the output structure, containing all relevant information on the metrics, including all values that are used for plotting.
To understand the structure of "Properties", see GenericSignals.m file. This includes comparisons of generic-time histories. The details on of the utilized methods are given in the abovementioned article.

%%%%%%%%% COPYRIGHT NOTICE %%%%%%%%% 

CompMet is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CompMet is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CompMet.  If not, see <https://www.gnu.org/licenses/>.
    
Copyright (c) Igor Kavrakov, Ahsan Kareem, Guido Morgenthal 2020

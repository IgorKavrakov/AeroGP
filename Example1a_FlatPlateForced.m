%% by Igor Kavrakov (igor.kavrakov@uni-weimar.de; igor.kavrakov@gmail.com)
clear all;  clc;  restoredefaultpath; matlabrc; close all;
addpath(genpath('GP'),'Example1_FlatPlateAnalytical');

% AeroGP: Data-driven model for aerodynamic analyses of structures using Gaussian Processes (GP)
% Please cite our work when you are you are using our software in your research or publications:

% Kavrakov, I., McRobie, A., and Morgenthal, G. 2022. Data-driven aerodynamic analysis of structures using Gaussian Processes. 
% J. Wind Eng. Ind. Aerodyn., 222, 104911. 
% https://doi.org/10.1016/j.jweia.2022.104911

% Accepted manuscript on arXiv:
% https://arxiv.org/abs/2103.13877

% The script includes a forced-vibration analysis for the flat plate example for using a GP-NFIR model.
% The script is based on the aforementioned article
% This script includes the results from Sec. 3 (Fundamental Application: Flat Plate)

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
fprintf(['AeroGP: Aerodynamic Analyses of Structures using Gaussian Processes \nIgor Kavrakov, Allan McRobie, Guido Morgenthal 2022 (c) \nCite as:\n Kavrakov, I., McRobie, A., and Morgenthal, G. 2022.\n Data-driven aerodynamic analysis of structures using Gaussian Processes.\n J. Wind Eng. Ind. Aerodyn., 222, 104911.\n https://doi.org/10.1016/j.jweia.2022.104911\n\n']);

%% Control
%GP properties
GP.Par.Train=0;         %If Train=0 - no training; Hyperparameters are loaded (default as in Github)
GP.Par.Lag=200;         %Regressors for S_alpha - alpha,. Note if S_alpha=1, it means we take 1 previous i.e. Alpha(i-1).
GP.Par.jitter=eps;      %Jitter term
GP.Par.SNR=20;          %Signal to nose ratio
GP.Par.SubSet=3;        %Subset of regressors (F)
GP.Par.Noise=-1;        %Optimise for noise? -1 >yes
GP.Par.OptiControl1=optimset('GradObj','on',...        %Optimiser properties
                             'TolX',1e-8,...
                             'Display','iter',...
                             'MaxFunEval',200000,...
                             'MaxIter', 500);

%Flat plate input & training properties
FP.Train.Alpha_STD=.1;      %Standard deviation of the effecitve angle - training
FP.Train.Alpha_AmpRel=.05;  %Relative Fourier amplitude - DoE for training (r_l)
FP.Train.Vr=[2 14];         %Input angle magnitude. Vr-Reduced velocity
FP.Train.Cycl=20;           %Cycles - number of cycles for training
FP.Par.ds=0.05;             %Reduced time-step (Warning: use 0.05 i.e. be consistent if using FP_RandInp)
%% Training
rng(1); %Reproducibility
if GP.Par.Train  
    input('!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!\n\nThe training of the GP model is computationally expensive (~3-4h - Linux server).\nFor prediction, use the determined hyperparameters by setting GP.Par.Train=0.\n\nPress any key to continue or cancel (CTRL+C).');    
    FP.Pred.Excitation='Train'; % Train   
    FP = Example1_FlatPlateAnalytical(FP,GP.Par.Lag); %FP model analytical
    GP.Par.Indx=sort(randperm(FP.Train.Samp,floor(FP.Train.Samp/GP.Par.SubSet))); %Subset of regressors for training  

    %Lift
    GP.Train(1).x=[FP.Train.alpha_h  FP.Train.alpha_a  FP.Train.alpha_h_d  FP.Train.alpha_a_d]; %Organise training vector
    GP.Train(1).x_max=max(max(abs(GP.Train(1).x))); %Normalization value 
    GP.Train(1).x=GP.Train(1).x/GP.Train(1).x_max;  %Normalize
 
    GP.Train(1).y_max=max(max(abs(FP.Train.CL))); %Normalize force 
    GP.Train(1).y=FP.Train.CL/GP.Train(1).y_max;  %Normalization parameters
    GP.Train(1).noise=randn([FP.Train.Samp,1])/GP.Par.SNR*std(GP.Train(1).y); %Add noise for training
    GP.Train(1).y=GP.Train(1).y+GP.Train(1).noise; %Add noise

    hyp=log(rand(6+2*GP.Par.Lag,1)); %Prior hyper parameters [Amp_kernel,Noise,alpha_h+Lag,alpha_a+Lag,alpha_h_d,alpha_a_d] in vectors.
    f=@(hyp)GP_Process_Opt(GP.Train(1).x(GP.Par.Indx,:),GP.Train(1).y(GP.Par.Indx),hyp,GP.Par.jitter,GP.Par.Noise); %Set up likelihood function to minimise
    GP.Train(1).hyp=fminunc(f,hyp,GP.Par.OptiControl1); % Minimize > get yperparameters

    %Moment
    GP.Train(2).x=GP.Train(1).x; %Same training input for moment as for lift
    GP.Train(2).x_max=GP.Train(1).x_max; %Normalization value > same as for moment

    GP.Train(2).y_max=max(max(abs(FP.Train.CM))); %Normalize force 
    GP.Train(2).y=FP.Train.CM./GP.Train(2).y_max; %Normalization parameters
    GP.Train(2).noise=randn([FP.Train.Samp,1])/GP.Par.SNR*std(GP.Train(2).y); %Add noise for training.
    GP.Train(2).y=GP.Train(2).y+GP.Train(2).noise;%Add noise

    f=@(hyp)GP_Process_Opt(GP.Train(2).x(GP.Par.Indx,:),GP.Train(2).y(GP.Par.Indx),hyp,GP.Par.jitter,GP.Par.Noise); %Set up likelihood function to minimise
    GP.Train(2).hyp=fminunc(f,hyp,GP.Par.OptiControl1); % Minimize > get yperparameters
    save('Example1_FlatPlateAnalytical/Example1_Train','GP');%Save GP hyperparameters
else
    load('Example1_FlatPlateAnalytical/Example1_Train','GP');%Load GP hyperparameters
end
%%  Prediction at Training Input
%Lift
FP.Pred.Excitation='Train'; %Prediction at training excitation
FP = Example1_FlatPlateAnalytical(FP,GP.Par.Lag); %FP model analytical
GP.Pred(1).x=GP.Train(1).x; 
GP.Pred(1).y_targ=GP.Train(1).y-GP.Train(1).noise;
[GP.Pred(1).y,GP.Pred(1).y_var,~,GP.Train(1).L_Kern]=GP_Process(GP.Train(1).x,GP.Train(1).y,GP.Pred(1).x,GP.Train(1).hyp,GP.Par.jitter,GP.Par.Noise); %Predictions

%Moment
GP.Pred(2).x=GP.Pred(1).x;
GP.Pred(2).y_targ=GP.Train(2).y-GP.Train(2).noise;
[GP.Pred(2).y,GP.Pred(2).y_var,~,GP.Train(2).L_Kern]=GP_Process(GP.Train(2).x,GP.Train(2).y,GP.Pred(2).x,GP.Train(2).hyp,GP.Par.jitter,GP.Par.Noise); %Predictions

%% Prediction at Random Input
FP.Pred.Excitation='Rand'; %Random excitation
FP = Example1_FlatPlateAnalytical(FP,GP.Par.Lag); %FP model analytical

GP.Train(3)=GP.Train(2); %Setup prediction model input & output
GP.Pred(3).x=[FP.Pred.alpha_h  FP.Pred.alpha_a  FP.Pred.alpha_h_d  FP.Pred.alpha_a_d];
GP.Pred(3).x=GP.Pred(3).x./GP.Train(3).x_max; %Normalize input             
GP.Pred(3).y_targ=FP.Pred.CM./GP.Train(3).y_max; %Normalize output
[GP.Pred(3).y,GP.Pred(3).y_var]=GP_Process_Predict(GP.Train(3).x,GP.Train(3).y,GP.Pred(3).x,GP.Train(3).hyp,GP.Train(3).L_Kern); %Predictions

%% Prediction at Sinusoidal Input
FP.Pred.Vr=[1:1:18];      %Predict at sinusoidal output - Vr range
FP.Pred.Alpha_Amp=.1;     %Prediction angle amplitude
FP.Pred.Cycl=ones(1,length(FP.Pred.Vr)).*6; %P Cycles - number of cycles for prediction
FP.Pred.Cycl(1)=12; % 12 cycles for Vr=1 to account for all lags.
StepsVr=cumsum([1 FP.Pred.Cycl.*FP.Pred.Vr./FP.Par.ds+1]); %Multi-step-ahead prediction

% Forced Vertical Displacements
FP.Pred.Excitation='SinH'; %Sin For sinusoidal input
FP = Example1_FlatPlateAnalytical(FP,GP.Par.Lag); %FP model analytical
 %Lift
 GP.Train(4)=GP.Train(1);
 GP.Pred(4).x=[FP.Pred.alpha_h  FP.Pred.alpha_a  FP.Pred.alpha_h_d  FP.Pred.alpha_a_d];
 GP.Pred(4).x=GP.Pred(4).x./GP.Train(4).x_max;              
 GP.Pred(4).y_targ=FP.Pred.CL./GP.Train(4).y_max;
 GP.Pred(4).y=GP.Pred(4).y_targ*0;
 for i=1:length(FP.Pred.Vr)
 [GP.Pred(4).y(StepsVr(i):StepsVr(i+1)-1)]=GP_Process_Predict(GP.Train(4).x,GP.Train(4).y,GP.Pred(4).x((StepsVr(i):StepsVr(i+1)-1),:),GP.Train(4).hyp,GP.Train(4).L_Kern); %Predictions
 end

 %Moment
 GP.Train(5)=GP.Train(2);
 GP.Pred(5).x=GP.Pred(4).x;
 GP.Pred(5).y_targ=FP.Pred.CM./GP.Train(5).y_max;
 GP.Pred(5).y=GP.Pred(5).y_targ*0;
 for i=1:length(FP.Pred.Vr)
 [GP.Pred(5).y(StepsVr(i):StepsVr(i+1)-1)]=GP_Process_Predict(GP.Train(5).x,GP.Train(5).y,GP.Pred(5).x((StepsVr(i):StepsVr(i+1)-1),:),GP.Train(5).hyp,GP.Train(5).L_Kern); %Predictions
 end

% Forced Rotation
FP.Pred.Excitation='SinA'; %Sin For sinusoidal input
FP = Example1_FlatPlateAnalytical(FP,GP.Par.Lag); %FP model analytical

 %Lift
 GP.Train(6)=GP.Train(1);
 GP.Pred(6).x=[FP.Pred.alpha_h  FP.Pred.alpha_a  FP.Pred.alpha_h_d  FP.Pred.alpha_a_d];
 GP.Pred(6).x=GP.Pred(6).x./GP.Train(6).x_max;              
 GP.Pred(6).y_targ=FP.Pred.CL./GP.Train(6).y_max;
 GP.Pred(6).y=GP.Pred(6).y_targ*0;
 for i=1:length(FP.Pred.Vr)
 [GP.Pred(6).y(StepsVr(i):StepsVr(i+1)-1)]=GP_Process_Predict(GP.Train(6).x,GP.Train(6).y,GP.Pred(6).x((StepsVr(i):StepsVr(i+1)-1),:),GP.Train(6).hyp,GP.Train(6).L_Kern); %Predictions
 end
 
%Moment
 GP.Train(7)=GP.Train(2);
 GP.Pred(7).x=GP.Pred(6).x;
 GP.Pred(7).y_targ=FP.Pred.CM./GP.Train(7).y_max;
 GP.Pred(7).y=GP.Pred(7).y_targ*0;
 for i=1:length(FP.Pred.Vr)
 [GP.Pred(7).y(StepsVr(i):StepsVr(i+1)-1)]=GP_Process_Predict(GP.Train(7).x,GP.Train(7).y,GP.Pred(7).x((StepsVr(i):StepsVr(i+1)-1),:),GP.Train(7).hyp,GP.Train(7).L_Kern); %Predictions
 end
 
%% Prediction at Static Input
FP.Pred.Excitation='Stat';          %Statatic wind coefficients
FP.Pred.StatAlpha=[-1 -.5 0 .5 1];  %Prediction angle
FP.Pred.StatTau=20;                 %Nondimensional time for prediction
Inc=(FP.Pred.StatTau./FP.Par.ds+1); %Increment in the prediction vector

FP = Example1_FlatPlateAnalytical(FP,GP.Par.Lag); %FP model analytical

GP.Train(8)=GP.Train(1);
GP.Pred(8).x=[FP.Pred.alpha_h  FP.Pred.alpha_a  FP.Pred.alpha_h_d  FP.Pred.alpha_a_d];
GP.Pred(8).x=GP.Pred(8).x./GP.Train(8).x_max;              
GP.Pred(8).y_targ=FP.Pred.CL./GP.Train(8).y_max;
GP.Pred(8).y=GP.Pred(8).y_targ*0;
for i=1:length(FP.Pred.StatAlpha)
[GP.Pred(8).y(1+(i-1)*Inc:i*Inc)]=GP_Process_Predict(GP.Train(8).x,GP.Train(8).y,GP.Pred(8).x((1+(i-1)*Inc:i*Inc),:),GP.Train(8).hyp,GP.Train(8).L_Kern); %Predictions
end

GP.Train(9)=GP.Train(2);
GP.Pred(9).x=[FP.Pred.alpha_h  FP.Pred.alpha_a  FP.Pred.alpha_h_d  FP.Pred.alpha_a_d];
GP.Pred(9).x=GP.Pred(9).x./GP.Train(9).x_max;              
GP.Pred(9).y_targ=FP.Pred.CM./GP.Train(9).y_max;
GP.Pred(9).y=GP.Pred(9).y_targ*0;
for i=1:length(FP.Pred.StatAlpha)
[GP.Pred(9).y(1+(i-1)*Inc:i*Inc)]=GP_Process_Predict(GP.Train(9).x,GP.Train(9).y,GP.Pred(9).x((1+(i-1)*Inc:i*Inc),:),GP.Train(9).hyp,GP.Train(9).L_Kern); %Predictions
end

% save('Example1_FlatPlateAnalytical/Ex1a_Out','GP','FP','-v7.3');
%% Plots
% load('Example1_FlatPlateAnalytical/Ex1a_Out','GP','FP');
Example1a_Plots(GP,FP); 

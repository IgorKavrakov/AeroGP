 function [M] = CompMet(X1,X2,CalculationProperties)
% CompMet is a Matlab-based computer code that computes metrics for comparison of two time-histories.
% The software is by:
% Igor Kavrakov (Chair of Modeling and Simulation of Structures, Bauhaus University Weimar)
% Ahsan Kareem (NatHaz Modeling Labratory, University of Notre Dame)
% Guido Morgenthal (Chair of Modeling and Simulation of Structures, Bauhaus University Weimar) 
% 
% Please cite our work when you are you are using our software in your research or publications:
% Kavrakov, I., Kareem, A., and Morgenthal, G. 2020. Comparison Metrics for Time-histories: Application to Bridge Aerodynamics. J. Eng. Mech., 146 (9), 040200093. https://doi.org/10.1061/(ASCE)EM.1943-7889.0001811

% The code compares two time histories (X1 and X2), taking one as a reference (X1).
% The CalculationProperties structure is explained in one of the examples.

%%%%%%%%% COPYRIGHT NOTICE %%%%%%%%% 
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

fprintf(['CompMet: Comparison Metrics for Time-Histories \nIgor Kavrakov, Ahsan Kareem, Guido Morgenthal 2020 (c) \nCite as: Kavrakov, I., Kareem, A., and Morgenthal, G. 2020. Comparison Metrics for Time-histories: Application to Bridge Aerodynamics. J. Eng. Mech., 146 (9), 040200093. https://doi.org/10.1061/(ASCE)EM.1943-7889.0001811\n\n']);
%% Assigmment

if nargin==3
if ~isfield(CalculationProperties,'Metrics'); CalculationProperties.Metrics=[]; M=[]; else M=[]; end
if ~isfield(CalculationProperties,'Plots'); CalculationProperties.Plots=[]; end

if ~isfield(CalculationProperties,'WavletProperties') && isfield(CalculationProperties.Metrics,'Wavelet') && CalculationProperties.Metrics.Wavelet; error('Please input correct WaveletProperties structure'); end
if ~isfield(CalculationProperties,'StationartyProperties') && isfield(CalculationProperties.Metrics,'Stationarty') && CalculationProperties.Metrics.Stationarty; error('Please input correct StationartyProperties structure'); end
if ~isfield(CalculationProperties,'WaveletBicoherenceProperties') && isfield(CalculationProperties.Metrics,'WaveletBicoherence') && CalculationProperties.Metrics.WaveletBicoherence; error('Please input correct WaveletBicoherenceProperties structure'); end
if ~isfield(CalculationProperties,'PDFProperties') && isfield(CalculationProperties.Metrics,'PD') && CalculationProperties.Metrics.PDF; error('Please input correct PDFProperties structure'); end

if isfield(CalculationProperties,'dt'); dt=CalculationProperties.dt;else dt=1; end
if isfield(CalculationProperties,'Name'); Name=CalculationProperties.Name;else Name='NameCase'; end
if isfield(CalculationProperties,'TPhase'); TPhase=CalculationProperties.TPhase;else TPhase=1; end

if ~isfield(CalculationProperties.Plots,'FigWidth'); CalculationProperties.Plots.FigWidth=19.5; end
if ~isfield(CalculationProperties.Plots,'FigDepth'); CalculationProperties.Plots.FigDepth=7.5; end

% FIX X for easier everything with fft
if mod(length(X1),2)==1;    X1=X1(1:end-1); X2=X2(1:end-1); end

else
    error('Insufficient Input');
end
%if isdir(Name);rmdir(Name,'s');end; mkdir(Name);
set(0,'DefaultTextInterpreter','LaTeX');
MetricsPlot=[];Names={};ID=0;CalculationProperties.Plots.Name=Name;
%% General properties
N=length(X1);
t=dt.*(0:1:N-1);
X1=X1(:); X2=X2(:);
M.Signals.X1=X1; M.Signals.X2=X2; M.dt=dt;
%% Metrics

% Cross correlation & %% Phase error based on the lag
if isfield(CalculationProperties.Metrics,'Phase')&&CalculationProperties.Metrics.Phase
    [M.Phase.CrossCorr,M.Phase] = Metric_CrossCorrPhase(X1,X2,dt,TPhase);
    ID=ID+1;MetricsPlot(ID)=M.Phase.Value;Names{ID}='$M_\phi$';
    fprintf('Phase Metric (Phi): %.2f \n', M.Phase.Value);    
end

% Peak Metric
if isfield(CalculationProperties.Metrics,'Peak')&&CalculationProperties.Metrics.Peak
    [M.Peak] = Metric_Peak(X1,X2);
    ID=ID+1;MetricsPlot(ID)=M.Peak.Value;Names{ID}='$M_p$';
    fprintf('Peak Metric (p): %.2f \n', M.Peak.Value);       
end

% RMS Metric
if isfield(CalculationProperties.Metrics,'RMS')&&CalculationProperties.Metrics.RMS
    [M.RMS] = Metric_RMS(X1,X2);
    ID=ID+1;MetricsPlot(ID)=M.RMS.Value;Names{ID}='$M_{\mathrm{rms}}$';
    fprintf('Root Mean Square Metric (rms): %.2f\n', M.RMS.Value);           
end
% Magnitude Metric
if isfield(CalculationProperties.Metrics,'MagnitudeWarp')&&CalculationProperties.Metrics.MagnitudeWarp
    [M.WrapedSignals.Dist_D,~,~,M.WrapedSignals.X2_Wrap,M.WrapedSignals.X1_Wrap]=dtw(X2,X1); % Discrete Time Wraped Signals - we use reverse (X2,X1) as X1 is referent signal! 
    [M.MagnitudeMetricWrap]=Metric_MagnitudeWrap(M.WrapedSignals.X1_Wrap,M.WrapedSignals.X2_Wrap,t,CalculationProperties.Plots);
    ID=ID+1;MetricsPlot(ID)=M.MagnitudeMetricWrap.Value;Names{ID}='$M_m$';
    fprintf('Warped Magnitude Metric (m): %.2f\n', M.MagnitudeMetricWrap.Value);               
end

% Magnitude Metric Unwarped
if isfield(CalculationProperties.Metrics,'Magnitude')&&CalculationProperties.Metrics.Magnitude
    [M.MagnitudeMetricUnWrap]=Metric_Magntitude(X1,X2);
    ID=ID+1;MetricsPlot(ID)=M.MagnitudeMetricUnWrap.Value;Names{ID}='$M_{muw}$';
    fprintf('Un-warped Magnitude Metric (muw): %.2f\n', M.MagnitudeMetricUnWrap.Value);                   
end

%Probability Distribution metric
if isfield(CalculationProperties.Metrics,'PDF')&&CalculationProperties.Metrics.PDF
    [M.PDF]=Metric_PDF(X1,X2,dt,CalculationProperties.PDFProperties,CalculationProperties.Plots);
    ID=ID+1;MetricsPlot(ID)=M.PDF.Value;Names{ID}='$M_{\mathrm{pdf}}$';
    fprintf('Probability Distribution Function Metric (pdf): %.2f\n', M.PDF.Value);                       
end

% Wavelet metric
if isfield(CalculationProperties.Metrics,'Wavelet')&&CalculationProperties.Metrics.Wavelet
    [M.Wavelet] = Metric_Wavelet(X1,X2,dt,CalculationProperties.WavletProperties,t,CalculationProperties.Plots);
    ID=ID+1;MetricsPlot(ID)=M.Wavelet.Value(1);Names{ID}='$M_{w}$';
    ID=ID+1;MetricsPlot(ID)=M.Wavelet.Value(2);Names{ID}='$M_{wf}$';
    fprintf('Wavelet Metric (w): %.2f\n', M.Wavelet.Value(1));        
    fprintf('Frequency-normalized Wavelet Metric (wf): %.2f\n', M.Wavelet.Value(2));        

    %Stationarity metric
   if isfield(CalculationProperties.Metrics,'Stationarity')&&CalculationProperties.Metrics.Stationarity
    [M.StationarityAnalysis]=Metric_Stationarity(X1,X2,t,dt,M.Wavelet,CalculationProperties.WavletProperties,CalculationProperties.StationartyProperties,[],CalculationProperties.Plots);
    ID=ID+1;MetricsPlot(ID)=M.StationarityAnalysis.Value;Names{ID}='$M_{s}$';
    fprintf('Stationarity Metric (s): %.2f\n', M.StationarityAnalysis.Value);        
   end
   
   %Bicoherence metric
    if isfield(CalculationProperties.Metrics,'WaveletBicoherence')&&CalculationProperties.Metrics.WaveletBicoherence
    [M.WaveletBicoherenceAnalysis,M.StationarityAnalysis]=Metric_WaveletBicoherence(X1,X2,dt,M.Wavelet,CalculationProperties.WavletProperties,M.StationarityAnalysis,CalculationProperties.WaveletBicoherenceProperties,CalculationProperties.Plots);
    ID=ID+1;MetricsPlot(ID)=M.WaveletBicoherenceAnalysis.Value;Names{ID}='$M_{b}$';
    fprintf('Bicoherence Metric (b): %.2f\n', M.WaveletBicoherenceAnalysis.Value);        
    end
end

%% General plots
if CalculationProperties.Plots.General
    FigWidth=CalculationProperties.Plots.FigWidth;FigDepth=CalculationProperties.Plots.FigDepth;
    Blue=[0 155/255 180/255];Red=[183/255 26/255 73/255];Green=[0.4660 0.6740 0.1880];Purp=[0.4940;0.1840;0.5560];Orange=[0.8500 0.3250  0.0980]; Gray=[0.5 0.5 0.5];

    % Figure 1 - Signal
    [positions] = subplot_pos(FigWidth,FigDepth,1.2,0.5,1,0.6,1,1,1,1);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [FigWidth FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
    ax=axes('position',positions{1,1},'Layer','top');
    hold on
     plot(t,X1,'color',Blue);
     plot(t,X2,'color',Red);

     l=legend({'$X_1$','$X_2$'},'location','south');  set(l, 'Interpreter', 'latex','FontSize',8)
     xlabel('$t$','Interpreter','latex','FontSize',9);ylabel('Amplitude','Interpreter', 'latex','FontSize',9);
     set(gca,'FontSize',8,'FontName','Times');
     title('Original Signals');
     box on; grid on;
     print(gcf,[Name '/'  Name '_Signals'],'-dpdf');%savefig(gcf,[Name '/'  Name '_Signals']);

    % Figure 2 - Metrics (Spider Plot)
    [positions] = subplot_pos(FigWidth,FigWidth,0.5,0.5,1,0.6,1,1,1,1);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [FigWidth FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
    ax=axes('position',positions{1,1},'Layer','top');
    SpiderPlot(MetricsPlot,Names,{Blue})    
    print(gcf,[Name '/'  Name '_Metrics'],'-dpdf');
end
%save([Name '/'  Name],'M','-v7.3');fprintf(['\n\n']);
end



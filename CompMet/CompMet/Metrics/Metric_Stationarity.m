function [ OutStruct ] = Metric_Stationarity(X1,X2,t,dt,WaveletAnalysis,WavletProperties,StationartyProperties,WaveletBicoherenceProperties,Plots)
% Calcuate Stationarity Metric 
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
    OutStruct.Name='Stationarity Metric (ST)';
    PlotHandle=0;PlotHandle1=0;
    df=WaveletAnalysis.FreqBins1(2)-WaveletAnalysis.FreqBins1(1);

    
    N_Average=ones(length(WaveletAnalysis.FreqBins1),1).*length(t);
    Vals=ones(length(WaveletAnalysis.FreqBins1),length(t));
    ColumnsOutInfluence=1:1:length(t);
    NColumnAverage=length(ColumnsOutInfluence);
    
    %Signal 1
    OutStruct.SpectraAverage1=sum(abs(WaveletAnalysis.WT1).^2.*Vals,2)./N_Average;
    OutStruct.SpectraAverage1(isnan(OutStruct.SpectraAverage1))=0;
    D_LSD=(abs(WaveletAnalysis.WT1).^2.*Vals)./repmat(OutStruct.SpectraAverage1,[1 length(t)]);
    D_LSD(isnan(D_LSD))=0;D_LSD=abs(log10(D_LSD)); D_LSD(isinf(D_LSD))=0;
    D_LSD=sum(D_LSD,1).*df; %log spectral deviation  
    Mean_LSD=sum(D_LSD(ColumnsOutInfluence))./NColumnAverage;
    OutStruct.Theta_LSD1=1./NColumnAverage.*sum((D_LSD(ColumnsOutInfluence)-Mean_LSD).^2);
    
    H_Norm=OutStruct.SpectraAverage1./sum(OutStruct.SpectraAverage1);G_Norm=abs(WaveletAnalysis.WT1.*Vals).^2./repmat(sum(abs(WaveletAnalysis.WT1.*Vals).^2,1),[length(WaveletAnalysis.Scales1) 1]);
    D_KL=(G_Norm-repmat(H_Norm,[1 length(t)])).*log10(G_Norm./repmat(H_Norm,[1 length(t)]));
    D_KL(isnan(D_KL))=0;D_KL(isinf(D_KL))=0;
    D_KL=sum(D_KL,1).*df; % Kullback Leibler (Symetric! Bassevele 1989)
    Mean_KL=sum(D_KL(ColumnsOutInfluence))./NColumnAverage;
    OutStruct.Theta_KL1=1./NColumnAverage.*sum((D_KL(ColumnsOutInfluence)-Mean_KL).^2);      

    OutStruct.Theta_KLG1=OutStruct.Theta_KL1.*(1+OutStruct.Theta_LSD1);
        
    %Surogates
    OutStruct.Theta_LSD_Sur1=zeros(StationartyProperties.NSurogates,1);
    OutStruct.Theta_KL_Sur1=zeros(StationartyProperties.NSurogates,1);
    [Sur1,~,~,~] = SurrogateSignal(X1,dt,StationartyProperties.NSurogates,0);
    OutStruct.WT_M1=zeros(length(WaveletAnalysis.FreqBins1),length(t)); %Average Scalogram
    OutStruct.SP_M1=zeros(length(WaveletAnalysis.FreqBins1),1); %Average Scalogram
   

        
    for i=1:StationartyProperties.NSurogates
    [WT,~,~,~] = MorletWavletTransform(Sur1(:,i),  WavletProperties.nLevel, WavletProperties.fmax, WavletProperties.fmin, WavletProperties.f0, 1/t(2),WavletProperties.beta,WavletProperties.Padding);

    OutStruct.WT_M1=OutStruct.WT_M1+abs(WT).^2;
    OutStruct.SP_M1=OutStruct.SP_M1+mean(abs(WT).^2,2); %Average Scalogram

    SpectraAverage=sum(abs(WT).^2.*Vals,2)./N_Average;
    SpectraAverage(isnan(SpectraAverage))=0;
    D_LSD=(abs(WT).^2.*Vals)./repmat(SpectraAverage,[1 length(t)]);D_LSD(isnan(D_LSD))=0;
    D_LSD=abs(log10(D_LSD));D_LSD(isinf(D_LSD))=0;D_LSD=sum(D_LSD,1).*df; 
    Mean_LSD=sum(D_LSD(ColumnsOutInfluence))./NColumnAverage;
    OutStruct.Theta_LSD_Sur1(i)=1./NColumnAverage.*sum((D_LSD(ColumnsOutInfluence)-Mean_LSD).^2);
       
    H_Norm=SpectraAverage./sum(SpectraAverage);G_Norm=abs(WT.*Vals).^2./repmat(sum(abs(WT.*Vals).^2,1),[length(WaveletAnalysis.Scales1) 1]);
    D_KL=(G_Norm-repmat(H_Norm,[1 length(t)])).*log10(G_Norm./repmat(H_Norm,[1 length(t)]));
    D_KL(isnan(D_KL))=0;D_KL(isinf(D_KL))=0;
    D_KL=sum(D_KL,1).*df; % Kullback Leibler (Symetric! Bassevele 1989)
    Mean_KL=sum(D_KL(ColumnsOutInfluence))./NColumnAverage;
    OutStruct.Theta_KL_Sur1(i)=1./NColumnAverage.*sum((D_KL(ColumnsOutInfluence)-Mean_KL).^2);        
    end
    OutStruct.WT_M1=OutStruct.WT_M1./StationartyProperties.NSurogates;
    OutStruct.SP_M1=OutStruct.SP_M1./StationartyProperties.NSurogates;
    
    OutStruct.Theta_KLG_Sur1=OutStruct.Theta_KL_Sur1.*(1+OutStruct.Theta_LSD_Sur1);
            
    if ~isfield(StationartyProperties,'Dist')||strcmpi(StationartyProperties.Dist,'gamma');
        % Just for fitting - divide by std (avoid numerical error)
        OutStruct.Theta_LSD_Sur1_ToFIT=OutStruct.Theta_LSD_Sur1./std(OutStruct.Theta_LSD_Sur1);
        OutStruct.Theta_LSD1_ToFIT=OutStruct.Theta_LSD1./std(OutStruct.Theta_LSD_Sur1);
        OutStruct.Theta_KL_Sur1_ToFIT=OutStruct.Theta_KL_Sur1./std(OutStruct.Theta_KL_Sur1);
        OutStruct.Theta_KL1_ToFIT=OutStruct.Theta_KL1./std(OutStruct.Theta_KL_Sur1);
        OutStruct.Theta_KLG_Sur1_ToFIT=OutStruct.Theta_KLG_Sur1./std(OutStruct.Theta_KLG_Sur1);
        OutStruct.Theta_KLG1_ToFIT=OutStruct.Theta_KLG1./std(OutStruct.Theta_KLG_Sur1);        
        
    OutStruct.DistParm1=fitdist(OutStruct.Theta_LSD_Sur1_ToFIT,'gamma');
    OutStruct.ProbExced1=gamcdf(OutStruct.Theta_LSD1_ToFIT,OutStruct.DistParm1.a,OutStruct.DistParm1.b);
    OutStruct.DistParmKL1=fitdist(OutStruct.Theta_KL_Sur1_ToFIT,'gamma');
    OutStruct.ProbExcedKL1=gamcdf(OutStruct.Theta_KL1_ToFIT,OutStruct.DistParmKL1.a,OutStruct.DistParmKL1.b);
    OutStruct.DistParmKLG1=fitdist(OutStruct.Theta_KLG_Sur1_ToFIT,'gamma');
    OutStruct.ProbExcedKLG1=gamcdf(OutStruct.Theta_KLG1_ToFIT,OutStruct.DistParmKLG1.a,OutStruct.DistParmKLG1.b);   
    end

    %Signal 2
    OutStruct.SpectraAverage2=sum(abs(WaveletAnalysis.WT2).^2.*Vals,2)./N_Average;
    OutStruct.SpectraAverage2(isnan(OutStruct.SpectraAverage2))=0;
    D_LSD=(abs(WaveletAnalysis.WT2.^2).*Vals)./repmat(OutStruct.SpectraAverage2,[1 length(t)]);
    D_LSD(isnan(D_LSD))=0;D_LSD=abs(log10(D_LSD)); D_LSD(isinf(D_LSD))=0;
    D_LSD=trapz(D_LSD,1).*df;
    Mean_LSD=sum(D_LSD(ColumnsOutInfluence))./NColumnAverage;
    OutStruct.Theta_LSD2=1./NColumnAverage.*sum((D_LSD(ColumnsOutInfluence)-Mean_LSD).^2);
    
    H_Norm=OutStruct.SpectraAverage2./sum(OutStruct.SpectraAverage2);G_Norm=abs(WaveletAnalysis.WT2.*Vals).^2./repmat(sum(abs(WaveletAnalysis.WT2.*Vals).^2,1),[length(WaveletAnalysis.Scales2) 1]);
    D_KL=(G_Norm-repmat(H_Norm,[1 length(t)])).*log10(G_Norm./repmat(H_Norm,[1 length(t)]));
    D_KL(isnan(D_KL))=0;D_KL(isinf(D_KL))=0;
    D_KL=sum(D_KL,1).*df; % Kullback Leibler (Symetric! Bassevele 1989)
    Mean_KL=sum(D_KL(ColumnsOutInfluence))./NColumnAverage;
    OutStruct.Theta_KL2=1./NColumnAverage.*sum((D_KL(ColumnsOutInfluence)-Mean_KL).^2);  
    
    OutStruct.Theta_KLG2=OutStruct.Theta_KL2.*(1+OutStruct.Theta_LSD2);   

    %Surogates
    OutStruct.Theta_LSD_Sur2=zeros(StationartyProperties.NSurogates,1);
    OutStruct.Theta_KL_Sur2=zeros(StationartyProperties.NSurogates,1);   
    [Sur2,~,~,~] = SurrogateSignal(X2,dt,StationartyProperties.NSurogates,0); 
    OutStruct.WT_M2=zeros(length(WaveletAnalysis.FreqBins2),length(t)); %Average Scalogram
    OutStruct.SP_M2=zeros(length(WaveletAnalysis.FreqBins2),1); %Average Scalogram


    for i=1:StationartyProperties.NSurogates
    [WT,~,~,~] = MorletWavletTransform(Sur2(:,i),  WavletProperties.nLevel, WavletProperties.fmax, WavletProperties.fmin, WavletProperties.f0, 1/t(2),WavletProperties.beta,WavletProperties.Padding);
    OutStruct.WT_M2=OutStruct.WT_M2+abs(WT).^2;
    OutStruct.SP_M2=OutStruct.SP_M2+mean(abs(WT).^2,2); %Average Scalogram

    SpectraAverage=sum(abs(WT).^2.*Vals,2)./N_Average;
    SpectraAverage(isnan(SpectraAverage))=0;
    D_LSD=(abs(WT).^2.*Vals)./repmat(SpectraAverage,[1 length(t)]);D_LSD(isnan(D_LSD))=0;
    D_LSD=abs(log10(D_LSD));D_LSD(isinf(D_LSD))=0;D_LSD=trapz(D_LSD,1).*df; 
    Mean_LSD=sum(D_LSD(ColumnsOutInfluence))./NColumnAverage;
    OutStruct.Theta_LSD_Sur2(i)=1./NColumnAverage.*sum((D_LSD(ColumnsOutInfluence)-Mean_LSD).^2);

    H_Norm=SpectraAverage./sum(SpectraAverage);G_Norm=abs(WT.*Vals).^2./repmat(sum(abs(WT.*Vals).^2,1),[length(WaveletAnalysis.Scales2) 1]);
    D_KL=(G_Norm-repmat(H_Norm,[1 length(t)])).*log10(G_Norm./repmat(H_Norm,[1 length(t)]));
    D_KL(isnan(D_KL))=0;D_KL(isinf(D_KL))=0;
    D_KL=sum(D_KL,1).*df; % Kullback Leibler (Symetric! Bassevele 1989)
    Mean_KL=sum(D_KL(ColumnsOutInfluence))./NColumnAverage;
    OutStruct.Theta_KL_Sur2(i)=1./NColumnAverage.*sum((D_KL(ColumnsOutInfluence)-Mean_KL).^2);           
    end
    OutStruct.WT_M2=OutStruct.WT_M2./StationartyProperties.NSurogates;
    OutStruct.SP_M2=OutStruct.SP_M2./StationartyProperties.NSurogates;
    OutStruct.Theta_KLG_Sur2=OutStruct.Theta_KL_Sur2.*(1+OutStruct.Theta_LSD_Sur2);
    
    if ~isfield(StationartyProperties,'Dist')||strcmpi(StationartyProperties.Dist,'gamma');
                % Just for fitting - divide by std (avoid numerical error)
        OutStruct.Theta_LSD_Sur2_ToFIT=OutStruct.Theta_LSD_Sur2./std(OutStruct.Theta_LSD_Sur2);
        OutStruct.Theta_LSD2_ToFIT=OutStruct.Theta_LSD2./std(OutStruct.Theta_LSD_Sur2);
        OutStruct.Theta_KL_Sur2_ToFIT=OutStruct.Theta_KL_Sur2./std(OutStruct.Theta_KL_Sur2);
        OutStruct.Theta_KL2_ToFIT=OutStruct.Theta_KL2./std(OutStruct.Theta_KL_Sur2);
        OutStruct.Theta_KLG_Sur2_ToFIT=OutStruct.Theta_KLG_Sur2./std(OutStruct.Theta_KLG_Sur2);
        OutStruct.Theta_KLG2_ToFIT=OutStruct.Theta_KLG2./std(OutStruct.Theta_KLG_Sur2);        

    OutStruct.DistParm2=fitdist(OutStruct.Theta_LSD_Sur2_ToFIT,'gamma');
    OutStruct.ProbExced2=gamcdf(OutStruct.Theta_LSD2_ToFIT,OutStruct.DistParm2.a,OutStruct.DistParm2.b);
    OutStruct.DistParmKL2=fitdist(OutStruct.Theta_KL_Sur2_ToFIT,'gamma');
    OutStruct.ProbExcedKL2=gamcdf(OutStruct.Theta_KL2_ToFIT,OutStruct.DistParmKL2.a,OutStruct.DistParmKL2.b); 
    OutStruct.DistParmKLG2=fitdist(OutStruct.Theta_KLG_Sur2_ToFIT,'gamma');
    OutStruct.ProbExcedKLG2=gamcdf(OutStruct.Theta_KLG2_ToFIT,OutStruct.DistParmKLG2.a,OutStruct.DistParmKLG2.b);   
    end
    %% Finally check stationarity
    % Check LG
    if isfield(StationartyProperties,'ConfidenceLevel')&&StationartyProperties.ConfidenceLevel>=0; ConfidenceLevel=(1-StationartyProperties.ConfidenceLevel); else ConfidenceLevel=0.05; end
    if OutStruct.ProbExced1<1-ConfidenceLevel; OutStruct.StationarityLG1=0; else OutStruct.StationarityLG1=1; end
    if OutStruct.ProbExced2<1-ConfidenceLevel; OutStruct.StationarityLG2=0; else OutStruct.StationarityLG2=1; end
    % Check KL
    if OutStruct.ProbExcedKL1<1-ConfidenceLevel; OutStruct.StationarityKL1=0; else OutStruct.StationarityKL1=1; end
    if OutStruct.ProbExcedKL2<1-ConfidenceLevel; OutStruct.StationarityKL2=0; else OutStruct.StationarityKL2=1; end
    % Check KLG
    if OutStruct.ProbExcedKLG1<1-ConfidenceLevel; OutStruct.StationarityKLG1=0; else OutStruct.StationarityKLG1=1; end
    if OutStruct.ProbExcedKLG2<1-ConfidenceLevel; OutStruct.StationarityKLG2=0; else OutStruct.StationarityKLG2=1; end
    
    if ~isfield(StationartyProperties,'DiscriminatingStatistic')||strcmpi(StationartyProperties.DiscriminatingStatistic,'LG')    
    OutStruct.Stationarity1=OutStruct.StationarityLG1;OutStruct.Stationarity2=OutStruct.StationarityLG2;% One if both are the same;
    elseif strcmpi(StationartyProperties.DiscriminatingStatistic,'KL')
    OutStruct.Stationarity1=OutStruct.StationarityKL1;OutStruct.Stationarity2=OutStruct.StationarityKL2;% One if both are the same;       
    elseif strcmpi(StationartyProperties.DiscriminatingStatistic,'KLG')
    OutStruct.Stationarity1=OutStruct.StationarityKLG1;OutStruct.Stationarity2=OutStruct.StationarityKLG2;% One if both are the same;               
    end
    OutStruct.Value=OutStruct.Stationarity1==OutStruct.Stationarity2;% One if both are the same;
    
    clear Mean LSD Indx D_LSD D_LSD D_LSD WT Vals NColumnAverage ColumnsOutInfluence SpectraAverage N_Average D_KL G_Norm H_Norm 
   % OutStruct.Value==1&&OutStruct.Stationarity1==1&&
 %% Surogates - Local nonstationarity  
        if (OutStruct.Value==1&&OutStruct.Stationarity1==1)||(isfield(StationartyProperties,'LocalAnalysis')&&StationartyProperties.LocalAnalysis==1);
            %Signal 1
            OutStruct.WT_STD1=zeros(length(WaveletAnalysis.FreqBins1),length(t)); %STD Scalogram
            OutStruct.SP_STD1=zeros(length(WaveletAnalysis.FreqBins1),1); %STD SPECTRA
            
            for i=1:StationartyProperties.NSurogates
            [WT, ~, ~,~] = MorletWavletTransform(Sur1(:,i),WavletProperties.nLevel, WavletProperties.fmax, WavletProperties.fmin, WavletProperties.f0, 1/t(2),WavletProperties.beta,WavletProperties.Padding);    
            OutStruct.WT_STD1=OutStruct.WT_STD1+(abs(WT).^2-OutStruct.WT_M1).^2;
            OutStruct.SP_STD1=OutStruct.SP_STD1+(mean(abs(WT).^2,2)-OutStruct.SP_M1).^2; %STD Scalogram            
            end
            OutStruct.WT_STD1=sqrt(OutStruct.WT_STD1./(StationartyProperties.NSurogates-1));
            OutStruct.SP_STD1=sqrt(OutStruct.SP_STD1./(StationartyProperties.NSurogates-1));          

            %g factor - for treshold
            if isfield(StationartyProperties,'g')&&StationartyProperties.g>=0; g=StationartyProperties.g; else g=1.8; end

            WT_Treshold=OutStruct.WT_M1+g.*OutStruct.WT_STD1;
            WT_Treshold=(repmat(max(WT_Treshold,[],2),[1 length(t)]));
            OutStruct.WT1_Filtered=abs(WaveletAnalysis.WT1);
            OutStruct.WT1_Filtered(find(OutStruct.WT1_Filtered.^2<WT_Treshold))=0;

            %Signal 2    
            OutStruct.WT_STD2=zeros(length(WaveletAnalysis.FreqBins2),length(t)); %STD Scalogram
            OutStruct.SP_STD2=zeros(length(WaveletAnalysis.FreqBins2),1); %STD SPECTRA          
            
            for i=1:StationartyProperties.NSurogates
            [WT, ~, ~,~] = MorletWavletTransform(Sur2(:,i),WavletProperties.nLevel, WavletProperties.fmax, WavletProperties.fmin, WavletProperties.f0, 1/t(2),WavletProperties.beta,WavletProperties.Padding);    
            OutStruct.WT_STD2=OutStruct.WT_STD2+(abs(WT).^2-OutStruct.WT_M2).^2;
            OutStruct.SP_STD2=OutStruct.SP_STD2+(mean(abs(WT).^2,2)-OutStruct.SP_M2).^2; %STD Scalogram            
            end
            OutStruct.WT_STD2=sqrt(OutStruct.WT_STD2./(StationartyProperties.NSurogates-1));
            OutStruct.SP_STD2=sqrt(OutStruct.SP_STD2./(StationartyProperties.NSurogates-1));          
            
            WT_Treshold=OutStruct.WT_M2+g.*OutStruct.WT_STD2;
            WT_Treshold=(repmat(max(WT_Treshold,[],2),[1 length(t)]));            
            OutStruct.WT2_Filtered=abs(WaveletAnalysis.WT2);
            OutStruct.WT2_Filtered(find(OutStruct.WT2_Filtered.^2<WT_Treshold))=0;
            
            OutStruct.Value=exp(-sqrt(sum(sum((abs(OutStruct.WT1_Filtered.*(WaveletAnalysis.ConeInfluence))-abs(OutStruct.WT2_Filtered.*(WaveletAnalysis.ConeInfluence))).^2)))/sqrt(sum(sum((OutStruct.WT1_Filtered.*(WaveletAnalysis.ConeInfluence)).^2))));
            PlotHandle=1;
        end %Local analysis

%% Plots
if ~isfield(WavletProperties,'MaxNormPlot'); MaxNormPlot=1.5; else MaxNormPlot=WavletProperties.MaxNormPlot; end
    maximumX1=max(max(abs(WaveletAnalysis.WT1)));

if isfield(Plots,'StatSurrMean')&&Plots.StatSurrMean==1; 

    [positions] = subplot_pos(Plots.FigWidth,Plots.FigDepth, 1.2, 0.5, 1.7, 0.6, 2, 1,1.5,1.5);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);

    ax=axes('position',positions{1,1},'Layer','top');
    hold on
    contourf(t,WaveletAnalysis.FreqBins1,sqrt(OutStruct.WT_M1)./maximumX1,1024/20, 'edgecolor','none');
    shading flat
    line(sqrt(2).*WaveletAnalysis.Scales1,WaveletAnalysis.FreqBins1,max(max(abs(sqrt(OutStruct.WT_M1))./maximumX1))*ones(1,length(WaveletAnalysis.Scales1)),'color','r','linewidth',2)     
    line(t(end)-sqrt(2).*WaveletAnalysis.Scales1,WaveletAnalysis.FreqBins1,max(max(abs(sqrt(OutStruct.WT_M1))./maximumX1))*ones(1,length(WaveletAnalysis.Scales1)),'color','r','linewidth',2)
    caxis([0 MaxNormPlot]);
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$t$ [s]','Interpreter', 'latex','FontSize',9);
    ylabel('$f$ [Hz]','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Mean Scalogram of Surrogates - Xs1 (normalized w.r.t. X1 (original))');
    box on; grid on; hold off

    ax=axes('position',positions{2,1},'Layer','top');
    hold on
    contourf(t,WaveletAnalysis.FreqBins2,sqrt(OutStruct.WT_M2)./maximumX1,1024/20, 'edgecolor','none');
    shading flat
    line(sqrt(2).*WaveletAnalysis.Scales2,WaveletAnalysis.FreqBins2,max(max(abs(sqrt(OutStruct.WT_M2))./maximumX1))*ones(1,length(WaveletAnalysis.Scales2)),'color','r','linewidth',2)     
    line(t(end)-sqrt(2).*WaveletAnalysis.Scales2,WaveletAnalysis.FreqBins2,max(max(abs(sqrt(OutStruct.WT_M2))./maximumX1))*ones(1,length(WaveletAnalysis.Scales2)),'color','r','linewidth',2)
    caxis([0 MaxNormPlot]);
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$t$ [s]','Interpreter', 'latex','FontSize',9);
    ylabel('$f$ [Hz]','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Mean Scalogram of Surrogates - Xs2 (normalized w.r.t. X1 (original))');
    box on; grid on; hold off

    print(gcf,[Plots.Name '\'  Plots.Name '_ScalogramSurrogates'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_ScalogramSurrogates']);
end

if isfield(Plots,'StatScalogramFiltered')&&Plots.StatScalogramFiltered==1&&PlotHandle==1;

    [positions] = subplot_pos(Plots.FigWidth,Plots.FigDepth, 1.2, 0.5, 1.7, 0.6, 2, 1,1.5,1.5);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);

    ax=axes('position',positions{1,1},'Layer','top');
    hold on
    contourf(t,WaveletAnalysis.FreqBins1,abs(OutStruct.WT1_Filtered)./maximumX1,1024/20, 'edgecolor','none');
    shading flat
    line(sqrt(2).*WaveletAnalysis.Scales1,WaveletAnalysis.FreqBins1,max(max(abs(OutStruct.WT1_Filtered)./maximumX1))*ones(1,length(WaveletAnalysis.Scales1)),'color','r','linewidth',2)     
    line(t(end)-sqrt(2).*WaveletAnalysis.Scales1,WaveletAnalysis.FreqBins1,max(max(abs(OutStruct.WT1_Filtered)./maximumX1))*ones(1,length(WaveletAnalysis.Scales1)),'color','r','linewidth',2)
    caxis([0 MaxNormPlot]);
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$t$ [s]','Interpreter', 'latex','FontSize',9);
    ylabel('$f$ [Hz]','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Filtered Scalogram - X1 (normalized w.r.t. X1)');
    box on; grid on; hold off

    ax=axes('position',positions{2,1},'Layer','top');
    hold on
    contourf(t,WaveletAnalysis.FreqBins2,abs(OutStruct.WT2_Filtered)./maximumX1,1024/20, 'edgecolor','none');
    shading flat
    line(sqrt(2).*WaveletAnalysis.Scales2,WaveletAnalysis.FreqBins2,max(max(abs(OutStruct.WT2_Filtered)./maximumX1))*ones(1,length(WaveletAnalysis.Scales2)),'color','r','linewidth',2)     
    line(t(end)-sqrt(2).*WaveletAnalysis.Scales2,WaveletAnalysis.FreqBins2,max(max(abs(OutStruct.WT2_Filtered)./maximumX1))*ones(1,length(WaveletAnalysis.Scales2)),'color','r','linewidth',2)
    caxis([0 MaxNormPlot]);
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$t$ [s]','Interpreter', 'latex','FontSize',9);
    ylabel('$f$ [Hz]','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Filtered Scalogram - X2 (normalized w.r.t. X1)');
    box on; grid on; hold off

    print(gcf,[Plots.Name '\'  Plots.Name '_ScalogramFiltered'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_ScalogramFiltered']);
end


end



function [OutStruct,StationartyAnalysis] = Metric_WaveletBicoherence(X1,X2,dt,WaveletAnalysis,WavletProperties,StationartyAnalysis,WaveletBicoherenceProperties,Plots)
% Calcuate Bicoherence Metric 
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

  OutStruct.Name='Wavelet AutobiSpectra Metric (B)';

if isfield(WaveletBicoherenceProperties,'NDivis')&&WaveletBicoherenceProperties.NDivis>=0; NDivis=WaveletBicoherenceProperties.NDivis; else NDivis=1; end
if isfield(WaveletBicoherenceProperties,'TRange')&&length(WaveletBicoherenceProperties.TRange)>=1; TRange=WaveletBicoherenceProperties.TRange; else TRange=0; end
if isfield(WaveletBicoherenceProperties,'Rand')&&WaveletBicoherenceProperties.Rand~=1; Rand=WaveletBicoherenceProperties.Rand; else Rand=1; end
if isfield(WaveletBicoherenceProperties,'Fact')&&WaveletBicoherenceProperties.Fact>0; Fact=WaveletBicoherenceProperties.Fact; else Fact=10*pi; end
if isfield(WaveletBicoherenceProperties,'Eps')&&WaveletBicoherenceProperties.Eps>0; Eps=WaveletBicoherenceProperties.Eps; else Eps=0.2; end

[OutStruct.FreqBins1,OutStruct.WaveletBiSpectra1,OutStruct.WaveletBiCoherence1,~,NonZeroVals]=WaveletBiCoherence(WaveletAnalysis.WT1,WaveletAnalysis.FreqBins1,dt,TRange,NDivis,Fact,Rand);
[OutStruct.FreqBins2,OutStruct.WaveletBiSpectra2,OutStruct.WaveletBiCoherence2]=WaveletBiCoherence(WaveletAnalysis.WT2,WaveletAnalysis.FreqBins2,dt,TRange,NDivis,Fact,Rand);

%For pratical reasons
OutStruct.WaveletBiSpectra1=abs(OutStruct.WaveletBiSpectra1).^2;
OutStruct.WaveletBiSpectra2=abs(OutStruct.WaveletBiSpectra2).^2;

Length=length(OutStruct.FreqBins1); % The random is only for plotting

if isfield(Plots,'WaveletBicoherence')&&Plots.WaveletBicoherence==1; 

    [positions] = subplot_pos(Plots.FigWidth,Plots.FigDepth, 1.2, 0.5, 1.7, 0.6, 2, 1,1.5,1.5);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);

    ax=axes('position',positions{1,1},'Layer','top');
    hold on
    contourf(OutStruct.FreqBins1,OutStruct.FreqBins1,mean(OutStruct.WaveletBiCoherence1,3)+rand(Length,Length).*0.001,1024/20, 'edgecolor','none');
    shading flat;
    caxis([0 1]); axis equal; xlim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]); ylim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]);
    set(gcf,'rend','painters') 
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$f_1$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f_2$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Wavelet Squared Bicoherence - X1');
    box on; grid on; hold off

    ax=axes('position',positions{2,1},'Layer','top');
    hold on
    contourf(OutStruct.FreqBins2,OutStruct.FreqBins2,mean(OutStruct.WaveletBiCoherence2,3)+rand(Length,Length).*0.001,1024/20, 'edgecolor','none');
    shading flat
    caxis([0 1]); axis equal; xlim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]); ylim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]);
    xlabel('$f_1$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f_2$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Wavelet Squared Bicoherence - X2');
    box on; grid ; hold off

    print(gcf,[Plots.Name '\'  Plots.Name '_WaveletBicoherence'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_WaveletBicoherence']);
end

if ~isfield(WaveletBicoherenceProperties,'MaxNormPlot'); MaxNormPlot=1.5; else MaxNormPlot=WaveletBicoherenceProperties.MaxNormPlot; end

if isfield(Plots,'WaveletBispectrum')&&Plots.WaveletBispectrum==1; 

    [positions] = subplot_pos(Plots.FigWidth,Plots.FigDepth, 1.2, 0.5, 1.7, 0.6, 2, 1,1.5,1.5);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);

    ax=axes('position',positions{1,1},'Layer','top');
    hold on
    contourf(OutStruct.FreqBins1,OutStruct.FreqBins1,mean(OutStruct.WaveletBiSpectra1,3)./max(max(mean(OutStruct.WaveletBiSpectra1,3)))+rand(Length,Length).*0.001,1024/20, 'edgecolor','none');
    shading flat;
    caxis([0 MaxNormPlot]); axis equal; xlim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]); ylim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]);
    set(gcf,'rend','painters') 
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$f_1$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f_2$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Wavelet Squared Bispectrum - X1 (normalized w.r.t. X1)');
    box on; grid on; hold off

    ax=axes('position',positions{2,1},'Layer','top');
    hold on
    contourf(OutStruct.FreqBins2,OutStruct.FreqBins2,mean(OutStruct.WaveletBiSpectra2,3)./max(max(mean(OutStruct.WaveletBiSpectra1,3)))+rand(Length,Length).*0.001,1024/20, 'edgecolor','none');
    shading flat
    caxis([0 MaxNormPlot]); axis equal; xlim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]); ylim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]);
    xlabel('$f_1$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f_2$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Wavelet Squared Bispectrum - X2 (normalized w.r.t. X1)');
    box on; grid ; hold off

    print(gcf,[Plots.Name '\'  Plots.Name '_WaveletBispectrum'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_WaveletBispectrum']);
end

%% Filtering
if isfield(WaveletBicoherenceProperties,'SoftTresholding')&&WaveletBicoherenceProperties.SoftTresholding==1
                %g factor - for treshold
        if isfield(WaveletBicoherenceProperties,'g')&&WaveletBicoherenceProperties.g>=0; g_bi=WaveletBicoherenceProperties.g; else g_bi=1.8; end
        if isfield(WaveletBicoherenceProperties,'NSurogates')&&WaveletBicoherenceProperties.NSurogates>=0; NSurogates_bi=floor(WaveletBicoherenceProperties.NSurogates); else NSurogates_bi=1; end
        
            WaveletBiSpectra_Surrogates=zeros(Length,Length,NSurogates_bi);
            for i=1:NSurogates_bi;
            [~,SurSpc]=WaveletBiSpectrumConstantPhaseRandomisation(WaveletAnalysis.WT1,WaveletAnalysis.FreqBins1,dt,TRange,NDivis,1,Fact,Eps);
            WaveletBiSpectra_Surrogates(:,:,i)=abs(SurSpc);            
            end
            OutStruct.WaveletBiSpectra_M1=mean(WaveletBiSpectra_Surrogates,3);
            OutStruct.WaveletBiSpectra_STD1=std(WaveletBiSpectra_Surrogates,[],3);

            WaveletBiSpectra_Surrogates=zeros(Length,Length,NSurogates_bi);
            for i=1:NSurogates_bi;
            [~,SurSpc]=WaveletBiSpectrumConstantPhaseRandomisation(WaveletAnalysis.WT2,WaveletAnalysis.FreqBins2,dt,TRange,NDivis,1,Fact,Eps);
            WaveletBiSpectra_Surrogates(:,:,i)=abs(SurSpc);
            end
            OutStruct.WaveletBiSpectra_M2=mean(WaveletBiSpectra_Surrogates,3);
            OutStruct.WaveletBiSpectra_STD2=std(WaveletBiSpectra_Surrogates,[],3);
            clear Sur1 WT WaveletBiSpectra_Surrogates

        Treshold_1=OutStruct.WaveletBiSpectra_M1+g_bi.*OutStruct.WaveletBiSpectra_STD1;
        Treshold_2=OutStruct.WaveletBiSpectra_M2+g_bi.*OutStruct.WaveletBiSpectra_STD2;
        
        OutStruct.WaveletBiSpectra1_Unfiltered=OutStruct.WaveletBiSpectra1;
        OutStruct.WaveletBiSpectra2_Unfiltered=OutStruct.WaveletBiSpectra2;
        for i=1:size(OutStruct.WaveletBiSpectra1,3)
        TempSpectra=OutStruct.WaveletBiSpectra1(:,:,i);TempSpectra(TempSpectra<=max(max(Treshold_1.^2)))=0;
        OutStruct.WaveletBiSpectra1(:,:,i)=TempSpectra;
        TempSpectra=OutStruct.WaveletBiSpectra2(:,:,i);TempSpectra(TempSpectra<=max(max(Treshold_2.^2)))=0;        
        OutStruct.WaveletBiSpectra2(:,:,i)=TempSpectra;
        end
        clear TempSpectra
        OutStruct.WaveletBiCoherence1_Unfiltered=OutStruct.WaveletBiCoherence1;
        OutStruct.WaveletBiCoherence2_Unfiltered=OutStruct.WaveletBiCoherence2;
        OutStruct.WaveletBiCoherence1(find(OutStruct.WaveletBiSpectra1==0))=0;
        OutStruct.WaveletBiCoherence2(find(OutStruct.WaveletBiSpectra2==0))=0;
end

OutStruct.WaveletBiCoherence1=mean(OutStruct.WaveletBiCoherence1,3);
OutStruct.WaveletBiCoherence2=mean(OutStruct.WaveletBiCoherence2,3);
OutStruct.WaveletBiSpectra1=mean(OutStruct.WaveletBiSpectra1,3);
OutStruct.WaveletBiSpectra2=mean(OutStruct.WaveletBiSpectra2,3);

if isfield(WaveletBicoherenceProperties,'HardTresholding')&&1>WaveletBicoherenceProperties.HardTresholding>0
OutStruct.WaveletBiSpectra1(find(OutStruct.WaveletBiCoherence1<=WaveletBicoherenceProperties.HardTresholding.^2))=0;      
OutStruct.WaveletBiSpectra2(find(OutStruct.WaveletBiCoherence2<=WaveletBicoherenceProperties.HardTresholding.^2))=0;      
OutStruct.WaveletBiCoherence1(find(OutStruct.WaveletBiCoherence1<=WaveletBicoherenceProperties.HardTresholding.^2))=0;      
OutStruct.WaveletBiCoherence2(find(OutStruct.WaveletBiCoherence2<=WaveletBicoherenceProperties.HardTresholding.^2))=0;      
end

if max(max(OutStruct.WaveletBiSpectra1))>0&&max(max(OutStruct.WaveletBiSpectra2))>0; 
OutStruct.Value=exp(-sqrt(sum(sum((sqrt(OutStruct.WaveletBiSpectra1(NonZeroVals))-sqrt(OutStruct.WaveletBiSpectra2(NonZeroVals))).^2)))/sqrt(sum(sum(OutStruct.WaveletBiSpectra1(NonZeroVals)))));   
elseif max(max(OutStruct.WaveletBiSpectra1))==0&&max(max(OutStruct.WaveletBiSpectra2))==0
OutStruct.Value=1;    
else
OutStruct.Value=0;
end 
%% Plots
% Prop.Plots.SurrogateBicoherence=1;  %Plot Bicoherence for filtering

if isfield(Plots,'WaveletBicoherenceFilterd')&&Plots.WaveletBicoherenceFilterd==1; 
    [positions] = subplot_pos(Plots.FigWidth,Plots.FigDepth, 1.2, 0.5, 1.7, 0.6, 2, 1,1.5,1.5);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);
    ax=axes('position',positions{1,1},'Layer','top');
    hold on
    contourf(OutStruct.FreqBins1,OutStruct.FreqBins1,OutStruct.WaveletBiCoherence1+rand(Length,Length).*0.001,1024/20, 'edgecolor','none');
    shading flat;
    caxis([0 1]); axis equal; xlim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]); ylim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]);
    set(gcf,'rend','painters') 
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$f_1$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f_2$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Filtered Squared Wavelet Bicoherence - X1');
    box on; grid on; hold off

    ax=axes('position',positions{2,1},'Layer','top');
    hold on
    contourf(OutStruct.FreqBins2,OutStruct.FreqBins2,OutStruct.WaveletBiCoherence2+rand(Length,Length).*0.001,1024/20, 'edgecolor','none');
    shading flat
    caxis([0 1]); axis equal; xlim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]); ylim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]);
    xlabel('$f_1$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f_2$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Filtered Squared Wavelet Bicoherence - X2');
    box on; grid ; hold off

    print(gcf,[Plots.Name '\'  Plots.Name '_WaveletBicoherenceFiltered'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_WaveletBicoherenceFiltered']);
end

if isfield(Plots,'WaveletBispectrumFilterd')&&Plots.WaveletBispectrumFilterd==1; 

    [positions] = subplot_pos(Plots.FigWidth,Plots.FigDepth, 1.2, 0.5, 1.7, 0.6, 2, 1,1.5,1.5);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);

    ax=axes('position',positions{1,1},'Layer','top');
    hold on
    contourf(OutStruct.FreqBins1,OutStruct.FreqBins1,OutStruct.WaveletBiSpectra1./max(max(OutStruct.WaveletBiSpectra1)),1024/20);
    shading flat;
    caxis([0 MaxNormPlot]); axis equal; xlim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]); ylim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))], 'edgecolor','none');
    set(gcf,'rend','painters') 
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$f_1$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f_2$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Filterd Wavelet Squared Bispectrum - X1 (normalized w.r.t. X1)');
    box on; grid on; hold off

    ax=axes('position',positions{2,1},'Layer','top');
    hold on
    contourf(OutStruct.FreqBins2,OutStruct.FreqBins2,OutStruct.WaveletBiSpectra2./max(max(OutStruct.WaveletBiSpectra1)),1024/20, 'edgecolor','none');
    shading flat
    caxis([0 MaxNormPlot]); axis equal; xlim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]); ylim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]);
    xlabel('$f_1$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f_2$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Filterd Wavelet Squared Bispectrum - X2 (normalized w.r.t. X1)');
    box on; grid ; hold off

    print(gcf,[Plots.Name '\'  Plots.Name '_WaveletBispectrumFilterd'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_WaveletBispectrumFilterd']);
end


if isfield(Plots,'WaveletBispectrumSurogate')&&Plots.WaveletBispectrumSurogate==1; 
    [positions] = subplot_pos(Plots.FigWidth,Plots.FigDepth, 1.2, 0.5, 1.7, 0.6, 2, 1,1.5,1.5);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);

    ax=axes('position',positions{1,1},'Layer','top');
    hold on
    contourf(OutStruct.FreqBins1,OutStruct.FreqBins1,Treshold_1.^2./max(max(OutStruct.WaveletBiSpectra1)),1024/20, 'edgecolor','none');
    shading flat;
    caxis([0 MaxNormPlot]); axis equal; xlim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]); ylim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]);
    set(gcf,'rend','painters') 
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$f_1$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f_2$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Surrogate Wavelet Squared Bispectrum - X1 (normalized w.r.t. X1)');
    box on; grid on; hold off

    ax=axes('position',positions{2,1},'Layer','top');
    hold on
    contourf(OutStruct.FreqBins2,OutStruct.FreqBins2,Treshold_2.^2./max(max(OutStruct.WaveletBiSpectra1)),1024/20, 'edgecolor','none');
    shading flat
    caxis([0 MaxNormPlot]); axis equal; xlim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]); ylim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]);
    xlabel('$f_1$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f_2$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Surrogate Wavelet Squared Bispectrum - X2 (normalized w.r.t. X1)');
    box on; grid ; hold off

    print(gcf,[Plots.Name '\'  Plots.Name '_WaveletBispectrumSurogate'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_WaveletBispectrumSurogate']);
end

if isfield(Plots,'WaveletBispectrumDiff')&&Plots.WaveletBispectrumDiff==1; 
    [positions] = subplot_pos(Plots.FigWidth/2,Plots.FigDepth,1.2,0.5,1.7,0.6,1,1,1.5,1.5);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);
    Length=length(OutStruct.FreqBins1); % The random is only for plotting
    ax=axes('position',positions{1,1},'Layer','top');
    hold on
    contourf(OutStruct.FreqBins1,OutStruct.FreqBins1,abs(OutStruct.WaveletBiSpectra1-OutStruct.WaveletBiSpectra2),1024/20, 'edgecolor','none');
    shading flat;
    caxis([0 MaxNormPlot]); axis equal; xlim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]); ylim([floor(OutStruct.FreqBins1(1)) ceil(OutStruct.FreqBins1(end))]);
    set(gcf,'rend','painters') 
    h=colorbar('location','southoutside','position',[0.1286    0.0525    0.8214    0.0435]);
    xlabel('$f_1$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f_2$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Difference Bicoherence - abs(B(X1(f,f))-B(X2(f,f))','Interpreter', 'none');
    box on; grid on; hold off

    print(gcf,[Plots.Name '\'  Plots.Name 'WaveletBispectrumDiff'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_WaveletBispectrumDiff']);

end

end


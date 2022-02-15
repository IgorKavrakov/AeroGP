function [ OutStruct ] = Metric_PDF( X1,X2,dt,PDProperties,Plots )
% Calcuate PDF Metric 
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

% This code uses the function KDE by Zdravko Botev (2015). See kde.m and
% akde1d for more information (incl. licence)

  OutStruct.Name='Probability Distribution Metric (PD)';
  
    if isfield(PDProperties,'StandardScore')&&PDProperties.StandardScore~=0 
    X1=(X1-mean(X1))./std(X1);
    X2=(X2-mean(X2))./std(X2);    
    end
    NSamp=length(X1);

 StartB=min(min(X1),min(X2));EndB=max(max(X1),max(X2));
     if isfield(PDProperties,'MinMax1')&&length(PDProperties.MinMax1)==2 
       StartB_1=PDProperties.MinMax1(1);
       EndB_1=PDProperties.MinMax1(2);
     else
       StartB_1=StartB;
       EndB_1=EndB;         
     end
      if isfield(PDProperties,'MinMax2')&&length(PDProperties.MinMax2)==2 
       StartB_2=PDProperties.MinMax2(1);
       EndB_2=PDProperties.MinMax2(2);
     else
       StartB_2=StartB;
       EndB_2=EndB;         
      end
  if  isfield(PDProperties,'KernelType')&&strcmpi(PDProperties.KernelType,'Adaptive')
      OutStruct.XVals_1=(StartB_1:(EndB_1-StartB_1)/PDProperties.Kerneldiscretization:EndB_1)';
      OutStruct.XVals_2=(StartB_2:(EndB_2-StartB_2)/PDProperties.Kerneldiscretization:EndB_2)';      
      if isfield(PDProperties,'Gamma')&&PDProperties.Gamma>0 
      [OutStruct.X1_PDF]=akde1d(X1,OutStruct.XVals_1,PDProperties.Gamma);
      [OutStruct.X2_PDF]=akde1d(X2,OutStruct.XVals_2,PDProperties.Gamma);
      else   
      [OutStruct.X1_PDF,OutStruct.XVals_1]=akde1d(X1,OutStruct.XVals_1);
      [OutStruct.X2_PDF,OutStruct.XVals_2]=akde1d(X2,OutStruct.XVals_2);         
      end
  else
      [OutStruct.KernelBandwith1,OutStruct.X1_PDF,OutStruct.XVals_1]=kde(X1,PDProperties.Kerneldiscretization,StartB_1,EndB_1);    %Kernel approximation
      [OutStruct.KernelBandwith2,OutStruct.X2_PDF,OutStruct.XVals_2]=kde(X2,PDProperties.Kerneldiscretization,StartB_2,EndB_2);   
  end
% Just for plotting
    if isfield(PDProperties,'NBins')&&PDProperties.NBins~=0 
       OutStruct.NBins1=PDProperties.NBins;
       OutStruct.NBins2=PDProperties.NBins;
    else
    RangeMin=min([ceil(2*NSamp^(1/3)),ceil(log2(NSamp)+1),ceil(sqrt(NSamp)),...
                            max(ceil(3.5*std(X1)*NSamp^(-1/3)),ceil(3.5*std(X2)*length(X2)^(-1/3)))]);
    RangeMax=max([ceil(2*NSamp^(1/3)),ceil(log2(NSamp)+1),ceil(sqrt(NSamp)),...
                            max(ceil(3.5*std(X1)*NSamp^(-1/3)),ceil(3.5*std(X2)*length(X2)^(-1/3)))]); 
    OutStruct.Range=RangeMin:1:RangeMax;         
            OutStruct.J_X1=zeros(length(OutStruct.Range),1);
            OutStruct.J_X2=zeros(length(OutStruct.Range),1);
            for i=1:length(OutStruct.Range)
                RangVal=abs(max(X1)-min(X1));
                X1Hist=hist(X1, OutStruct.Range(i));
                h=RangVal/OutStruct.Range(i);
                OutStruct.J_X1(i)=2/((NSamp-1)*h)-(NSamp+1)/(NSamp^2*(NSamp-1)*h).*sum((X1Hist).^2);
                RangVal=abs(max(X2)-min(X2));
                X2Hist=hist(X2, OutStruct.Range(i));
                h=RangVal/OutStruct.Range(i);
                OutStruct.J_X2(i)=2/((NSamp-1)*h)-(NSamp+1)/(NSamp^2*(NSamp-1)*h).*sum((X2Hist).^2);
            end
            [~,idx_1]=min(OutStruct.J_X1);[~,idx_2]=min(OutStruct.J_X2);
            OutStruct.NBins1=OutStruct.Range(idx_1);
            OutStruct.NBins2=OutStruct.Range(idx_2);
    end
            [OutStruct.X1_PDF_Hist,OutStruct.Cent1]=hist(X1,OutStruct.NBins1);OutStruct.X1_PDF_Hist=OutStruct.X1_PDF_Hist./(sum(OutStruct.X1_PDF_Hist)*(OutStruct.Cent1(2)-OutStruct.Cent1(1))); %Normalize!
            [OutStruct.X2_PDF_Hist,OutStruct.Cent2]=hist(X2,OutStruct.NBins2);OutStruct.X2_PDF_Hist=OutStruct.X2_PDF_Hist./(sum(OutStruct.X2_PDF_Hist)*(OutStruct.Cent2(2)-OutStruct.Cent2(1))); %Normalize!

             if isfield(Plots,'PDF')&&Plots.PDF==1
               Blue=[0 155/255 180/255];Red=[183/255 26/255 73/255];  
               [positions] = subplot_pos(Plots.FigWidth/2,Plots.FigDepth,1.2,0.5,1.7,0.6,1,1,1.5,1.5);
               f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);
               ax=axes('position',positions{1,1},'Layer','top');
               hold on
               plot(OutStruct.Cent1,OutStruct.X1_PDF_Hist,'o','color',Blue);
               plot(OutStruct.Cent2,OutStruct.X2_PDF_Hist,'s','color',Red);
               plot(OutStruct.XVals_1,OutStruct.X1_PDF,'-','color',Blue);
               plot(OutStruct.XVals_2,OutStruct.X2_PDF,'-','color',Red);               
               xlabel('$X_1$,$X_2,$','Interpreter', 'latex','FontSize',9);
               ylabel('PDF','Interpreter', 'latex','FontSize',9);
               set(gca,'FontSize',8,'FontName','Times');
               title(['Density of X1 and X2 '],'Interpreter', 'none');
               l=legend({'$X_1$','$X_2$','Approx $X_1$','Approx $X_2$'},'location','NorthEast');  set(l, 'Interpreter', 'latex','FontSize',8)               
               box on; grid on; hold off
               print(gcf,[Plots.Name '\'  Plots.Name '_PDF'],'-dpdf');
             end
% Integration - setting a uniform distance if discretisation of kernels is not the same
     if StartB_1~=StartB_2||EndB_1~=EndB_2
         h1=OutStruct.XVals_1(2)-OutStruct.XVals_1(1);h2=OutStruct.XVals_2(2)-OutStruct.XVals_2(1);
         h=min(h1,h2); %Minimum step
         OutStruct.Range_Intg=max(StartB_1,StartB_2):h:min(EndB_1,EndB_2); %Only where there are values
         OutStruct.X1_PDF_Intg=interp1(OutStruct.XVals_1,OutStruct.X1_PDF,OutStruct.Range_Intg);
         OutStruct.X2_PDF_Intg=interp1(OutStruct.XVals_2,OutStruct.X2_PDF,OutStruct.Range_Intg);
     else
         h=OutStruct.XVals_1(2)-OutStruct.XVals_1(1);
         OutStruct.X1_PDF_Intg=OutStruct.X1_PDF;
         OutStruct.X2_PDF_Intg=OutStruct.X2_PDF;         
     end
             


BCoefficient=sum(sqrt(OutStruct.X1_PDF_Intg.*OutStruct.X2_PDF_Intg)*h); %Bhattacharyya coefficient
OutStruct.Value=exp(-(-log(BCoefficient))); %h since integration over probabiltiy d. function!


end


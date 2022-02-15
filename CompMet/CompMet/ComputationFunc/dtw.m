function [Dist,DistMat,Path,X1_Wrap,X2_Wrap]=dtw(X1,X2)
%UNTITLED9 Summary of this function goes here
% Dynamic Time Warping
% This file is part of CompMetTH.
% CompMetTH is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%  CompMetTH is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with CompMetTH.  If not, see <https://www.gnu.org/licenses/>.

% Copyright (c) Igor Kavrakov, Ahsan Kareem, Guido Morgenthal 2020


L1=length(X1); L2=length(X2);
Norm2=(repmat(X1,1,L2)-repmat(X2',L1,1)).^2; %Euclidian norm squared

DistMat=zeros(size(Norm2)); %Distance matrix
DistMat(:,1)=cumsum(Norm2(:,1)); %First column
DistMat(1,:)=cumsum(Norm2(1,:)); %Second column

for i=2:L1
    for j=2:L2
        DistMat(i,j)=Norm2(i,j)+min(DistMat(i-1,j),min(DistMat(i-1,j-1),DistMat(i,j-1))); %Min c
    end
end
Dist=DistMat(end,end); %Total distance

% now try to get the path sequence
i=L2;j=L1; %Initiate!
Path=[L1 L2];
while (i+j)~=2
    if (i-1)==0
        j=j-1;
    elseif (j-1)==0
        i=i-1;
    else 
      [~,CC]=min([DistMat(j-1,i),DistMat(j,i-1),DistMat(j-1,i-1)]); %Check where should we go
      if CC==1
        j=j-1; %Vertical decent
      elseif CC==2
        i=i-1; %Horizontal decent
      else
        j=j-1;i=i-1; %Diagonal decent
      end
    end
    Path=[j i; Path]; 
end
X1_Wrap=X1(Path(:,1));%Warped signals
X2_Wrap=X2(Path(:,2));%Warped signal


function MapColor=colorMap(Value,ThreeDimention)
% ThreeDimention: default:0, if >0 makes sure the output has 3 dimension
if nargin<2 ThreeDimention=0; end
[mm, nn]=size(Value);
if mm==1 & nn==1
    if      (Value<0)       MapColor=[0 0 1];    
    elseif  (Value<=1/3)    MapColor = [0, 3*Value , 1];
    elseif  (Value<=0.50)   MapColor = [0, 1, 6*(0.5-Value)];
    elseif  (Value<=2/3)    MapColor = [6*(Value-0.5), 1, 0];
    elseif  (Value<=1)      MapColor = [1, 3*(1-Value), 0];
    else                    MapColor = [1 0 0]; 
    end
    if ThreeDimention MapColor=reshape(MapColor,1,1,3); end
else
    %MapColor=ones(mm,nn,3);
    MapColor1 = ones(mm,nn);
    MapColor2 = zeros(mm,nn);
    MapColor3 = zeros(mm,nn);
    
    ind=find(Value<=1);
    %MapColor1(ind)=1;
    MapColor2(ind)=3*(1-Value(ind));
    %MapColor3(ind)=0;
    
    ind=find(Value<=2/3);
    MapColor1(ind)=6*(Value(ind)-0.5);
    MapColor2(ind)=1;
    %MapColor3(ind)=0;
    
    ind=find(Value<=0.50);
    MapColor1(ind)=0;
    %MapColor2(ind)=1;
    MapColor3(ind)=6*(0.5-Value(ind));
    
    ind=find(Value<=1/3);
    %MapColor1(ind)=0;
    MapColor2(ind)=3*Value(ind);
    MapColor3(ind)=1;
    
    ind=find(Value<0);
    %MapColor1(ind)=0;
    MapColor2(ind)=0;
    %MapColor3(ind)=1;    
    
    if nn==1 & ~ThreeDimention 
        MapColor = [MapColor1 MapColor2 MapColor3]; 
    elseif mm==1 & ~ThreeDimention 
        MapColor = [MapColor1; MapColor2; MapColor3];
    else
        MapColor(:,:,1) = MapColor1; 
        MapColor(:,:,2) = MapColor2;
        MapColor(:,:,3) = MapColor3;
    end
end
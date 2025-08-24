function [x,y,uwarp,vwarp,iter]= FindDisplacements (I1,I2,wxy,sxy)

%   Last Update: 21-10-14

%   2-D Particle Image Velocimetry
%
%   wxy= [wx,wy]  interrogation window HALF size
%   sxy= [sx,sy]  window spacing in x and y
%

disp('Finding displacements...')


MAXITER= 20;

if size(I1) ~= size(I2)
    disp ('Images should be same size ...');
    return
end

if wxy(1) < 1 | wxy(2) < 1
    disp ('Wrong window size wxy ...');
    return
end

if sxy(1) < 1 | sxy(2) < 1
    disp ('Wrong window spacing sxy ...');
    return
end

[height,width]=size(I1);

% make sure images are in DOUBLE format
I1=double(I1);
I2=double(I2);

FI1= griddedInterpolant(I1);
FI2= griddedInterpolant(I2);

% loop over all subwindows: start with on offset of wxy2 so as to allow
% for displacements at the image outer edge
iout=0;
for xoff=1+wxy(1):sxy(1):width-wxy(1)+1-wxy(1)
    
    iout= iout + 1;
    jout= 0;
    
    for yoff=1+wxy(2):sxy(2):height-wxy(2)+1-wxy(2)
    
        jout= jout + 1;
         
%         [peak,warp,shift,iter(jout,iout)]= fast_warp (I1,I2,[xoff,yoff],[0,0],wxy,MAXITER,'off');               
        [peak,warp,shift,mse,iter(jout,iout)]= fast_warp_BOS_2way ...
                        (I1,I2,FI1,FI2,[xoff,yoff],[0,0],wxy,MAXITER,false);   
        % [peak,warp,shift,mse,iter(jout,iout)]= fast_warp_BOS_2way ...
        %     (I1,I2,FI1,FI2,[xoff,yoff],[0,0],wxy,MAXITER,true);
		% Copy result into the data set 
		x(jout,iout)= xoff + wxy(1);
		y(jout,iout)= yoff + wxy(2);
        uwarp(jout,iout)= shift(1);
        vwarp(jout,iout)= shift(2);
        
    end
end

return


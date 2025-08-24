function [peak,warp,shift,iter]= fast_warp (map1,map2,F1,F2,xyoff,dxy,wxy,maxiter,display)
%
% last update:  25-Apr-2019
%
% symmetric warp projection, using COMPLETE images to avoid boundary
% clipping
%
% transform is  x'= Ax + b;
%               A= eye(2,2) + warp;     b= shift + dxy
%
% the tranform is applied in the forward direction on one image, in the
% backward direction on the other ...
%
% map1, map2:       complete input images
% F1, F2:           gridded interpolant structures for map1, map2
% xyoff(2):         input vector defining the center point of the match
% dxy(2):           input vector defining the initial HALF separation 
%                   (+/- dxy) of the patterns
% wxy(2):           HALF size of pattern region (full size is 2*wxy+1)
% maxiter:          max. number of iterations allowed
% display:          activates progess display ("on", "off")
%
% peak:             correlation value between matched patterns
% warp(2,2):        deformation coefficients in warp matrix
% shift(2):         displacement vector between patterns (includes offset)
% iter:             number of iterations until convergence
%

warp= zeros(2,2);
shift= zeros(2,1);
peak= 0;

map1= double(map1);
map2= double(map2);

% make sure we're dealing with integer coordinates ...
wxy =round(wxy);
dxy= round(dxy);
xyoff= round(xyoff);

% create centered coordinate system in search window
[x,y]= meshgrid(-wxy(1):wxy(1),-wxy(2):wxy(2));

% compute fixed part of forward and reverse coordinate transforms
xyp0= eye(2,2)*[x(:)';y(:)'] + repmat(xyoff',1,numel(x));

try
    % extract data for subwindows (could fail if pixels are outside image)
    xc1= xyoff(1)-dxy(1)-wxy(1);
    yc1= xyoff(2)-dxy(2)-wxy(2);
    newmap1= map1(yc1:yc1+2*wxy(2),xc1:xc1+2*wxy(1));
    xc2= xyoff(1)+dxy(1)-wxy(1);
    yc2= xyoff(2)+dxy(2)-wxy(2);
    newmap2= map2(yc2:yc2+2*wxy(2),xc2:xc2+2*wxy(1));
catch
    iter= maxiter;
    return
end

% prepare display
if strcmp(display,'on')
    handle= figure (1);
    clf
    hold off;
    set (gcf,'color',[1,1,1]);
    subplot (1,3,1);
    cc(:,:,1)= newmap1;
    cc(:,:,2)= newmap2;
    cc(:,:,3)= 0.;
    cc= cc/ max(cc(:));
    imshow (cc,[]);
    title ('Match');
    axis off
    subplot (1,3,2);
    imshow (zeros(wxy(2),wxy(1)),[0 8]);
    title ('MSE');
    subplot (1,3,3);
    rectx= [-.5,.5,.5,-.5,-.5];
    recty= [-.5,-.5,.5,.5,-.5];
    plot (rectx,recty,'r');
    xlim ([-2,2]);
    ylim ([-2,2]);
    axis square;
    grid on
    title ('Transform (rel. to dxy)');
    hold on;
    drawnow;
end
   
% preallocate some arrays ...
nx= 2*wxy(1)+1;
ny= 2*wxy(2)+1;
A= zeros([nx*ny,6]);
kx= zeros(ny,nx);    
ky= zeros(nx,nx); 

iter= 0;
sol= ones(6,1);
while norm(sol)>1.E-3 && iter<maxiter
	
    iter= iter + 1;
    
    % reshape into 2D arrays for gradient computation
    if iter > 1
        newmap1= reshape (newmap1,ny,nx);
        newmap2= reshape (newmap2,ny,nx);
    end
    
    % compute gradients kx, ky of kxy
    [kx,ky]= gradient(newmap1 + newmap2);
    
    A(:,1)= kx(:).*x(:);
    A(:,2)= kx(:).*y(:);
    A(:,3)= ky(:).*x(:);
    A(:,4)= ky(:).*y(:);
    A(:,5)= kx(:);
    A(:,6)= ky(:);
    
    % solve ...
    sol= A \ (newmap1(:) - newmap2(:));

    % copy into output variables
	warp= warp + [[sol(1),sol(2)] ; [sol(3),sol(4)]];
	shift= shift + [sol(5);sol(6)];
    
    % compute new, interpolated maps
    shiftnxy= repmat(shift+dxy',1,nx*ny);
    xyp1= xyp0 - warp*[x(:)';y(:)'] - shiftnxy;
    newmap1= F1(xyp1(2,:),xyp1(1,:));
    xyp2= xyp0 + warp*[x(:)';y(:)'] + shiftnxy;
    newmap2= F2(xyp2(2,:),xyp2(1,:));
      
    if strcmp(display,'on')
        figure (handle)
        subplot (1,3,1);
        cc(:,:,1)= reshape(newmap2,[ny,nx])/max(newmap2(:));
        cc(:,:,2)= reshape(newmap1,[ny,nx])/max(newmap1(:));
        cc(:,:,3)= 0;
        imshow (cc,[]);
        title ('Match');
        axis off
        subplot (1,3,2);
        if iter == 1 
            hmax= max(newmap2(:)-newmap1(:));
        end
        imshow (abs(reshape(newmap2-newmap1,ny,nx)),[0 hmax]);
        title ('MSE');
        subplot (1,3,3);
        newrectxy= (eye(2,2)+warp)*[rectx(:)';recty(:)'] + repmat(shift,1,5);
        plot (newrectxy(1,:),newrectxy(2,:),'b');
        drawnow;
    end
    
end

% compute new correlation peak value
cmat= cov(newmap1,newmap2);
peak= cmat(1,2) / sqrt(cmat(1,1)*cmat(2,2));

% must re-scale results because of symmetric transform
warp= 2*warp;
shift= 2*shift' + 2*dxy;

return




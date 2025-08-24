function [peak,warp,shift,mse,iter]= fast_warp_BOS_2way (ref,map,Fref,Fmap,xypos,xyoffset,wxy,maxiter,display)
%
% last update:  16-Oct-2020
%
% symmetric warp projection, using COMPLETE images to avoid boundary
% clipping
%
% this version for BOS data assumes that the transform matrox A is symmetric, 
% implying that the transforming field is curl-free (i.e. derived from a potential)
%
% transform is  x'= Ax + b;
%               A= eye(2,2) + warp;     b= shift + dxy
%
% the tranform is applied in the forward direction on one image, in the
% backward direction on the other ...
%
% ref, map:         complete reference and deformed inut images
% Fref, Fmap:       gridded interpolant structures for ref, map
% xypos(2):         input vector defining the center point of the match
% xyoffset(2):      input vector defining the initial HALF separation 
%                   (+/- xyoffset) of the patterns
% wxy(2):           HALF size of pattern region (full size is 2*wxy+1)
% maxiter:          max. number of iterations allowed
% display:          activates progess display (true / false)
%
% peak:             correlation value between matched patterns
% warp(2,2):        deformation coefficients in warp matrix
% shift(2):         displacement vector between patterns (includes offset)
% mse               mean squared error between matched maps
% iter:             number of iterations until convergence
%

warp= zeros(2,2);
shift= zeros(2,1);
peak= 0;
mse= NaN;

ref= double(ref);
map= double(map);

% make sure we're dealing with integer coordinates ...
wxy= round(wxy);
xyoffset= round(xyoffset);
xypos= round(xypos);

% create centered coordinate system in search window
[xwin,ywin]= meshgrid(-wxy(1):wxy(1),-wxy(2):wxy(2));

% compute fixed part of forward and reverse coordinate transforms
xyp0= eye(2,2)*[xwin(:)';ywin(:)'] + repmat(xypos',1,numel(xwin));

try
    % extract data for subwindows (could fail if pixels are outside image)
    xc1= xypos(1) - xyoffset(1) - wxy(1);
    yc1= xypos(2) - xyoffset(2) - wxy(2);
    newref= ref(yc1:yc1+2*wxy(2),xc1:xc1+2*wxy(1));
    xc2= xypos(1) + xyoffset(1) - wxy(1);
    yc2= xypos(2) + xyoffset(2) - wxy(2);
    newmap= map(yc2:yc2+2*wxy(2),xc2:xc2+2*wxy(1));
catch
    iter= NaN;
    return
end

% prepare display
if display
    handle= figure (1);
    clf
    hold off;
    set (gcf,'color',[1,1,1]);
    subplot (1,3,1);
    cc(:,:,1)= newref;
    cc(:,:,2)= newmap;
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
A= zeros([nx*ny,5]);
kx= zeros(ny,nx);    
ky= zeros(nx,nx); 

iter= 0;
sol= ones(5,1);
while norm(sol) > 1.E-3 && iter < maxiter
	
    iter= iter + 1;
    
    % reshape into 2D arrays for gradient computation
    if iter > 1
        newref= reshape (newref,2*wxy(2)+1,2*wxy(1)+1);
        newmap= reshape (newmap,2*wxy(2)+1,2*wxy(1)+1);
    end
    
    % compute gradients kx, ky of kxy
    [kx,ky]= gradient(newref + newmap);

    A(:,1)= kx(:).*xwin(:);
    A(:,2)= kx(:).*ywin(:) + ky(:).*xwin(:);
    A(:,3)= ky(:).*ywin(:);
    A(:,4)= kx(:);
    A(:,5)= ky(:);

    % solve the associated least squares problem ...
    % (note implicit minus-sign for r.h.s.)
    sol= A \ (newref(:) - newmap(:));
 
    % copy into output variables
	warp= warp + [[sol(1),sol(2)] ; [sol(2),sol(3)]];
	shift= shift + [sol(4);sol(5)];
    
    % compute new, interpolated maps
    shiftnxy= repmat(shift+xyoffset',1,nx*ny);
    xyp1= xyp0 - warp*[xwin(:),ywin(:)]' - shiftnxy;
    newref= Fref(xyp1(2,:),xyp1(1,:));
    xyp2= xyp0 + warp*[xwin(:),ywin(:)]' + shiftnxy;
    newmap= Fmap(xyp2(2,:),xyp2(1,:));
          
    if strcmp(display,'on')
        figure (handle)
        subplot (1,3,1);
        cc(:,:,1)= reshape(newref,[2*wxy(2)+1,2*wxy(1)+1])/max(newref(:));
        cc(:,:,2)= reshape(newmap,[2*wxy(2)+1,2*wxy(1)+1])/max(newmap(:));
        cc(:,:,3)= 0;
        imshow (cc,[]);
        title ('Match');
        axis off
        subplot (1,3,2);
        if iter == 1 
            hmax= max(newmap(:)-newref(:));
        end
        imshow (abs(reshape(newmap-newref,ny,nx)),[0 hmax]);
        title ('MSE');
        subplot (1,3,3);
        newrectxy= (eye(2,2)+warp)*[rectx(:)';recty(:)'] + repmat(shift,1,5);
        plot (newrectxy(1,:),newrectxy(2,:),'b');
        drawnow;
    end
    
end

% compute mean squared error between matched maps
mse= mean((newmap(:)-newref(:)).^2);

% compute new correlation peak value
cmat= cov(newmap,newref);
peak= cmat(1,2) / sqrt(cmat(1,1)*cmat(2,2));

% must re-scale results because of symmetric transform
warp= 2*warp;
shift= 2*shift(:) + 2*xyoffset(:);

return

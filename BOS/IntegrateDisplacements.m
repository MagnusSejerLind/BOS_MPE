function [phi, err]= IntegrateDisplacements (u,v,delta_xy,left,right,top,bottom)

% latest update: 3-Nov-2020

% solves BOS gradient / Poisson problem using direct differentiation fomula 
% 5 point polynomial stencil is assumed in gradient computation

% input:    (u,v)           = 2D gradient field with horizontal (u) and vertical (v) vector components
%           delta_xy        = conversion between grid units and world coordinates [(m/pixel)*(pixel/gridunit)]
%           left, right ... = vectors with fixed values on boundaries
%                             (use NaN if not used, Inf to enforce zero normal slope)
% output:   phi             = array of same size as (u,v); good points are finite; masked points are NaN

% data is assumed to be stored in column-first (y) order, 2D array

[ny,nx]= size(u);

if nx==1 | ny==1
    disp ('ERROR: grid data has to be two-dimensional');
    return;
end

if nx<5 | ny<5
    disp ('ERROR: data fields are too narrow (< 5 points)');
    return;
end

% first: generate masks to control horizontal and vertical processing
% (use u velocity component to detect NaN points
mask= u;
[hmask,vmask]= MakeMasks (mask);

% number of unknown points
npts= nx*ny;

% number of equations
nboundary= max(sum(isfinite(left)) + sum(isfinite(right)) + sum(isfinite(top)) + sum(isfinite(bottom)) , 1); 
neqns= nboundary + length(hmask(:)) + length(find(isfinite(vmask(:))));

% number of non-zero matrix elements
nh0= length(find(hmask(:)==0));
nh1= length(find(hmask(:)~=0 & isfinite(hmask(:))));
nh2= length(find(isnan(hmask(:))));
nzh= 4*nh0 + 5*nh1 + nh2;
nv0= length(find(vmask(:)==0));
nv1= length(find(vmask(:)~=0 & isfinite(vmask(:))));
nzv= 4*nv0 + 5*nv1;
nonzeros= nzh+nzv+nboundary;

% allocate sparse arrays
% a= spalloc(neqns,npts,nonzeros);
a= spalloc(npts,neqns,nonzeros);
b= zeros(neqns,1);

k= 1;

% insert Dirichlet boundary conditions with boundary value(s) specified

% left edge: refractive index
for kk=1:ny
    if isfinite(left(kk))
        a(kk,k)= 12;
        b(k)= left(kk);
        k= k + 1;
    end
end

% right edge: refractive index
for kk=1:ny
    if isfinite(right(kk))
        a(kk+(nx-1)*ny,k)= 12;
        b(k)= right(kk);
        k= k + 1;
    end
end

% top edge: refractive index
for kk=1:nx
    if isfinite(top(kk))
        a(1+(kk-1)*ny,k)= 12;
        b(k)= top(kk);
        k= k + 1;
    end
end

% bottom edge: refractive index
for kk=1:nx
    if isfinite(bottom(kk))
        a(ny+(kk-1)*ny,k)= 12;
        b(k)= bottom(kk);
        k= k + 1;
    end
end

% if no boundary values are given, insert condition for zero surface value
% at corner (1,1)
if k == 1
	a(1,k)= 12;
	b(k)= 0;
	k= k + 1;
end

% insert v. Neumann Boundary conditions
uu= u;
vv= v;

% left edge
for kk=1:ny
    if isinf(left(kk))
        uu(kk,1)= 0.;
    end
end

% right edge
for kk=1:ny
    if isinf(right(kk))
        uu(kk,end)= 0.;
    end
end

% top edge
for kk=1:nx
    if isinf(top(kk))
        vv(1,kk)= 0.;
    end
end

% bottom edge
for kk=1:nx
    if isinf(bottom(kk))
        vv(end,kk)= 0.;
    end
end

% fill matrix for horizontal gradient equations
for i=1:nx
    for j=1:ny
        ij= (i-1)*ny+j;
        switch hmask(j,i)
        case 0
            % central difference
            a(ij-2*ny:ny:ij+2*ny,k)= [1,-8,0,8,-1];
            b(k)= uu(j,i) * delta_xy;
            k= k + 1;
        case -1
            a(ij-ny:ny:ij+3*ny,k)= [-3,-10,18,-6,1];
            b(k)= uu(j,i) * delta_xy;
            k= k + 1;
        case 1
            a(ij-3*ny:ny:ij+ny,k)= [-1,6,-18,10,3];
            b(k)= uu(j,i) * delta_xy;
            k= k + 1;
        case -2
            a(ij:ny:ij+4*ny,k)= [-25,48,-36,16,-3];
            b(k)= uu(j,i) * delta_xy;
            k= k + 1;
        case 2
            a(ij-4*ny:ny:ij,k)= [3,-16,36,-48,25];
            b(k)= uu(j,i) * delta_xy; 
            k= k + 1;
        otherwise
            % points inside masked area are irrelevant but kept in matrix to simplify cartesian addressing
            a(ij,k)= 12;
            b(k)= 0;
            k= k + 1;
        end
    end
end

% fill matrix for vertical gradient equations
for i=1:nx
    for j=1:ny
        ij= (i-1)*ny+j;
        switch vmask(j,i)
        case 0
            % central difference
            a(ij-2:ij+2,k)= [1,-8,0,8,-1];
            b(k)= vv(j,i) * delta_xy;
            k= k + 1;
        case -1
            a(ij-1:ij+3,k)= [-3,-10,18,-6,1];
            b(k)= vv(j,i) * delta_xy;
            k= k + 1; 
        case 1
            a(ij-3:ij+1,k)= [-1,6,-18,10,3];
            b(k)= vv(j,i) * delta_xy;
            k= k + 1;
        case -2
            a(ij:ij+4,k)= [-25,48,-36,16,-3];
            b(k)= vv(j,i) * delta_xy;
            k= k + 1;
        case 2
            a(ij-4:ij,k)= [3,-16,36,-48,25];
            b(k)= vv(j,i) * delta_xy;
            k= k + 1;    
        otherwise
            % NaN point ignored
        end
    end
end

% all coefficients in matrix "a" are scaled by a factor 12 ...
% due to filling sequence of matrix a, must use its transpose a'
a= a' / 12;

% overdetermined system: MATLAB solves automatically by least squares ...
phi= full(a \ b);

% this should also work, but svd() is not implemented for sparse matrices :-(
% phi= pinv(a)*b;

% could use SVD to kill off small eigenvalues -  if it weren't so 
% computationally expensive :-(
% [u,s,v]= svds(a,npts);
% phi= v*((u'*b)./diag(s));

% standard, explicit least squares solution :-}
% overdetermined system: solve by least squares ...
% ata= a' * a;
% ab= a' * b;
% phi= full(ata \ ab);

% reshuffle into square grid; fill bad points with NaNs
phi= reshape(phi,ny,nx);
phi(find(isnan(mask)))= NaN;

% compute error measure
[px,py]= gradient(phi);
err= sqrt( (px-u*delta_xy).^2 + (py-v*delta_xy).^2 );

return

%---------------------------------------------------------------------------------------------------------------

function [hmask,vmask]= MakeMasks (maskpattern)

% prepares masks for GradientDirectFilterStep ...
% masks are made based on information in image 'maskpattern',
% where a 'NaN' pixel means a bad data point and anything else is ok ...

[ny,nx]= size(maskpattern);

% by default: all points are "bad"; must find valid neighbours...
hmask= repmat(NaN,ny,nx);
vmask= repmat(NaN,ny,nx);

% horizontal processing

% first: check for good central points
for i=3:nx-2
    for j=1:ny
        l2= maskpattern(j,i-2);
        l1= maskpattern(j,i-1);
        c= maskpattern(j,i);
        r1= maskpattern(j,i+1);
        r2= maskpattern(j,i+2);
        if (isfinite(l2) & isfinite(l1) & isfinite(c) & isfinite(r1) & isfinite(r2))
            hmask(j,i)= 0;
        end
    end
end

% next: check remainder for good points with one NaN neighbour left
for i=2:nx-3
    for j=1:ny
        if isnan(hmask(j,i))
            l1= maskpattern(j,i-1);
            c= maskpattern(j,i);
            r1= maskpattern(j,i+1);
            r2= maskpattern(j,i+2);
            r3= maskpattern(j,i+3);
            if isfinite(l1) & isfinite(c) & isfinite(r1) & isfinite(r2) & isfinite(r3)
                hmask(j,i)= -1;
            end
        end
    end
end

% next: check remainder for good points with one NaN neighbour right
for i=4:nx-1
    for j=1:ny
        if isnan(hmask(j,i))
            l3= maskpattern(j,i-3);
            l2= maskpattern(j,i-2);
            l1= maskpattern(j,i-1);
            c= maskpattern(j,i);
            r1= maskpattern(j,i+1);
            if isfinite(l3) & isfinite(l2) & isfinite(l1) & isfinite(c) & isfinite(r1)
                hmask(j,i)= 1;
            end
        end
    end
end

% next: check remainder for good points with two NaN neighbours left
for i=1:nx-4
    for j=1:ny
        if isnan (hmask(j,i))
            c= maskpattern(j,i);
            r1= maskpattern(j,i+1);
            r2= maskpattern(j,i+2);
            r3= maskpattern(j,i+3);
            r4= maskpattern(j,i+4);
            if isfinite(c) & isfinite(r1) & isfinite(r2) & isfinite(r3) & isfinite(r4)
                hmask(j,i)= -2;
            end
        end
    end
end

% next: check remainder for good points with two NaN neighbours right
for i=5:nx
    for j=1:ny
        if isnan(hmask(j,i))        
            l4= maskpattern(j,i-4);
            l3= maskpattern(j,i-3);
            l2= maskpattern(j,i-2);
            l1= maskpattern(j,i-1);
            c= maskpattern(j,i);
            if isfinite(l4) & isfinite(l3) & isfinite(l2) & isfinite(l1) & isfinite(c)
                hmask(j,i)= 2;
            end
        end
    end
end

% vertical processing, same philosophy ...

for j=3:ny-2
    for i=1:nx
        u2= maskpattern(j-2,i);
        u1= maskpattern(j-1,i);
        c= maskpattern(j,i);
        d1= maskpattern(j+1,i);
        d2= maskpattern(j+2,i);
        if isfinite(u2) & isfinite(u1) & isfinite(c) & isfinite(d1) & isfinite(d2)
            vmask(j,i)= 0;
        end
    end
end

for j=2:ny-3
    for i=1:nx
        if isnan(vmask(j,i))
            u1= maskpattern(j-1,i);
            c= maskpattern(j,i);
            d1= maskpattern(j+1,i);
            d2= maskpattern(j+2,i);
            d3= maskpattern(j+3,i);
            if isfinite(u1) & isfinite(c) & isfinite(d1) & isfinite(d2) & isfinite(d3)
                vmask(j,i)= -1;
            end
        end
    end
end

for j=4:ny-1
    for i=1:nx
        if isnan(vmask(j,i))
            u3= maskpattern(j-3,i);
            u2= maskpattern(j-2,i);
            u1= maskpattern(j-1,i);
            c= maskpattern(j,i);
            d1= maskpattern(j+1,i);
            if isfinite(u3) & isfinite(u2) & isfinite(u1) & isfinite(c) & isfinite(d1)
                vmask(j,i)= 1;
            end
        end
    end
end

for j=1:ny-4
    for i=1:nx
        if isnan(vmask(j,i))
            c= maskpattern(j,i);
            d1= maskpattern(j+1,i);
            d2= maskpattern(j+2,i);
            d3= maskpattern(j+3,i);
            d4= maskpattern(j+4,i);
            if isfinite(c) & isfinite(d1) & isfinite(d2) & isfinite(d3) & isfinite(d4)
                vmask(j,i)= -2;
            end
        end
    end
end

for j=5:ny
    for i=1:nx
        if isnan(vmask(j,i))
            u4= maskpattern(j-4,i);
            u3= maskpattern(j-3,i);
            u2= maskpattern(j-2,i);
            u1= maskpattern(j-1,i);
            c= maskpattern(j,i);
            if isfinite(u4) & isfinite(u3) & isfinite(u2) & isfinite(u1) & isfinite(c)
                vmask(j,i)= 2;
            end
        end
    end
end

return
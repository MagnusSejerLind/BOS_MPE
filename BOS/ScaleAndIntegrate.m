function [nxy,alpha] = ScaleAndIntegrate(u,v,BC,opt)

options= optimset ('Display','iter');
alpha0= -1.0;

% BCs
if opt.BC == 1
    % Top/bottom BC
    left= NaN(size(u(:,1)));
    right= NaN(size(u(:,end)));
    top = BC.A;
    bottom = BC.B;
elseif opt.BC == 0
    % Left/right BC
    top= NaN(size(u(1,:)));
    bottom= NaN(size(u(end,:)));
    left = BC.A;
    right = BC.B;
else
    disp('Error in BC selection')
end


% nonlinear least squares optimization
% the geometry scale factor alpha is constrained to be negative
% to account for the inverted sign of the gradients resulting from the
% PIV analysis output
alpha = lsqnonlin(@(alpha)errfun(alpha,u,v,left,right,top,bottom), ...
                    alpha0,[],[0],options)

[nxy,err]= IntegrateDisplacements (u,v,alpha,left,right,top,bottom);

return

function err= errfun (alpha,u,v,left,right,top,bottom)

[nxy,err]= IntegrateDisplacements (u,v,alpha,left,right,top,bottom);

return
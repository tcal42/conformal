function images = bcfun_eval(tt, bndry, W, Wprime, S, ctr)
%%BCFUN_EVAL
%   This code is an interface to kerzstein.m.
%
%   Input:
%     tt = vector of points in [0,2*pi]
%     bndry = complex chebfun on [0,2*pi] defining boundary of domain
%     ctr = complex number defining point in domain to be mapped to 0
%
%   Output:
%     images = vector of images 
    npts = length(tt);
    if npts < 3
        images = kerzstein(bndry, W, Wprime, S, npts, ctr);  
    else        
        images = kerzstein(bndry, W, Wprime, S, npts-1, ctr);  
        images = [images; images(1)];
    end
end

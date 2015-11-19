function win=gausswin(sigma, truncate)
% One-dimensional Gaussian filter.
% 
% Parameters
% ----------
% %(input)s
% sigma : scalar
%     standard deviation for Gaussian kernel
% %(axis)s
% truncate : float, optional
%     Truncate the filter at this many standard deviations.
%     Default is 4.0.
% 
% Returns
% -------
% gaussian_filter1d : ndarray
%From scipy.ndimage.filters library
if nargin<1
    disp('At least two parameters needed')
elseif nargin<2
    truncate = 4;
end

sd              = sigma;
lw              = floor(truncate * sd + 0.5);
weights         = zeros(2*lw + 1, 1);
weights(lw +1)  = 1.0;
sum_g           = 1.0;
sd              = sd * sd;

%calculate the kernel:
for ii = 1 : lw 
    tmp = exp(-0.5 * (ii * ii) / sd);
    weights(lw+1 + ii) = tmp;
    weights(lw+1 - ii) = tmp;
    sum_g = sum_g + 2.0 * tmp;
end

win =  weights./(length(weights));



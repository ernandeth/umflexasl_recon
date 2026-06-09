function [out, comps] = GLR_cleanup(ims, Ncomps)
%function out = GLR_cleanup(ims, Ncomps)
% Global low rank clean up of a time series with varying undersampling
%
% the clean version of each pixel's time course
% is a weighted sum of the main temporal components extracted from the
% whole image.
% Estimate the coefficients of the components (using least squares)
% and add them together as a weighted sum
%
dims = size(ims);
NpixTotal = dims(1)*dims(2)*dims(3);
Nframes = dims(4);

Ncomps = 5;
out = zeros(size(ims));
casorati = reshape(ims, NpixTotal, Nframes);
out = zeros(size(casorati));

% decompose into its components
[u, s, v] = svd(casorati,'econ');
% keep only the biggest temporal components
comps = v(:,1:Ncomps);
plot(real(comps));
pr = 0;

for p = 1 :size(casorati,1)

    tseries = casorati(p,:)';
    out(p,:)= sum (comps*pinv(comps)*real(tseries),2);
    if mod(pr,1000)==0
        fprintf('\rProgress ... %0.2f', 100*pr/size(casorati,1) )
    end
    pr=pr+1;

end

out = reshape(out, dims);


return
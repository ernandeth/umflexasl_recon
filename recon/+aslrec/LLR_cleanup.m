function out = LLR_cleanup(ims, nbrs, Ncomps)
%function out = LLR_cleanup(ims, nbrs, Ncomps)
% Locally low rank clean up of a time series with varying undersampling
%
% the clean version of each pixel's time course
% is a weighted sum of the main temporal components in the neighborhood patch
% estimate the coefficients of the components (using least squares)
% and add them together as a weighted sum
% 
dims = size(ims);
out = zeros(size(ims));

Npix = (nbrs*2+1)^3;
NpixTotal = dims(1)*dims(2)*dims(3);
Nframes = dims(4);
prog=0;
for m = nbrs+1 : dims(1)-nbrs
    for n = nbrs+1 : dims(2)-nbrs
        for p = nbrs+1 : dims(3)-nbrs

            % time series at center f the patch
            tseries = squeeze(ims(m,n,p,:));
            patch = ims(m-nbrs:m+nbrs, n-nbrs:n+nbrs, p-nbrs:p+nbrs, :);

            cpatch = reshape(patch, Npix, Nframes);
            % decompose the patch into its components
            [u, s, v] = svd(cpatch,'econ');
            % keep only the biggest temporal components
            comps = v(:,1:Ncomps);

            out(m,n,p,:) = sum (comps*pinv(comps)*tseries,2);
            if mod(prog/NpixTotal, 1000)==0
                fprintf('\\ Progress: %0.2f    ', 100*prog/NpixTotal)
            end

            prog=prog+1;
        end
    end
end

return
function out = LLR_cleanup(ims, nbrs, Ncomps, X)
%function out = LLR_cleanup(ims, nbrs, Ncomps [,X] )
% Locally low rank clean up of a time series with varying undersampling
%
% the clean version of each pixel's time course
% is a weighted sum of the main temporal components in the neighborhood patch
% estimate the coefficients of the components (using least squares)
% and add them together as a weighted sum
% 

if nargin<4
    X=[];
end

dims = size(ims);
out = zeros(size(ims));

Npix = (nbrs*2+1)^3;
NpixTotal = dims(1)*dims(2)*dims(3);
Nframes = dims(4);
progr=0;
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

            if ~isempty(X)
                % de-correlate the cokmponents from the design matrix in FMRI
                for j=1:Ncomps
                    comps(:,j)= comps(:,j) - X*pinv(X)*comps(:,j);
                end
            end

            % zero mean all the junk regressors
            comps = comps - repmat(mean(comps,1),size(tseries,1),1);

            % estimate and remove components from the signal
            bhat = pinv(comps)*tseries;
            out(m,n,p,:) = tseries - comps*bhat;  

        end
        
    end
    fprintf('\n Progress: %0.2f    ', 100*progr/NpixTotal)
    progr = progr+dims(2)*dims(1);
end

return
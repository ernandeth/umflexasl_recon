function [clean junkcomps] = compcor23(dirty, Ncomps, X)
% function [clean junkcomps] = compcor12(dirty, [,Ncomps] [,DesMat])
%
% written by Luis Hernandez-Garcia at UM (c) 2023
%
% the dirty input data should be a matrix : Tframes by Npixels
% The hdr is just to let us know the dimensions of the images
%
% sample usage:
%   [dirty hdr] = readnii('dirtyFiles.nii');
%   [clean junkcoms] = compcor23(dirty);
%   writenii('cleanFile.nii', clean, 'hdr', hdr);
%
% Preferred method: you can put the junkcomps into the
% design matrix of a GLM.  In that case, it may be good
% to orthogonalize the junk regressors (J) from the
% original matrix (X), like this:
%
%   for n=1:Ncomps
%       J2(:,n) = J(:,n) - X*pinv(X)*J(:,n)
%   end
%

if nargin<2
    Ncomps = 10;
end
if nargin<3
    X=[];
end

showFigs = 0;

TopPercentSigma = 30;
%TopPercentSigma = 10;

fprintf('\nFinding the top %d components in the noisiest pixels \n(pixels with top %d percent of the variance)\n...',Ncomps, TopPercentSigma)
dims = size(dirty);
Npix = dims(1)*dims(2)*dims(3);
Nframes = dims(4);
dirty = reshape(dirty,Npix, Nframes)';

sigma = std(dirty,1);
msk = ccvarmask(sigma, TopPercentSigma);
%msk = ones(size(sigma));

if(showFigs)
    figure;
    subplot(211)
    lightbox(reshape(sigma,dims(1), dims(2), dims(3)));
    title('std. dev. map')
    subplot(212)
    lightbox(reshape(msk,dims(1), dims(2), dims(3)));
    title('pixels for compcor');
    set(gcf,'Name', 'Compcor Uses the Noisiest voxels')
end

mdirty = dirty .* repmat(msk, size(dirty,1),1);
mdirty(isinf(mdirty))=0;
mdirty(isnan(mdirty))=0;

[u, s,v]=svd(mdirty','econ');

junkcomps = v(:,1:Ncomps);

% mean center the components:
junkcomps = junkcomps - repmat(mean(junkcomps,1),size(dirty,1),1);

if ~isempty(X)
    fprintf('\nDecorrelating components from design matrix\n...')
    badinds = [];  % preserve the highest component - usually the mean
    for n=1:size(junkcomps,2)
        % de-correlate the junk components from the design matrix (desired effects)
        junkcomps(:,n)= junkcomps(:,n) - X*pinv(X)*junkcomps(:,n);

%{
        % test the correlation between junkcomps and the design matrix.
        % if they are correlated do not use them for clean up.
        design =  X*pinv(X)*junkcomps(:,n);
        rho = corrcoef( junkcomps(:,n),design);
        % fprintf('\nCorrelation of %d-th component with X = %f', n, rho(1,2));
        if abs(rho(1,2)) > 0.4;
            fprintf('... will not remove noise Component %d');
        else
            badinds = [badinds ; n];
        end
%}
    end
    junkcomps(:,badinds) = [];

end

if showFigs
    figure
    subplot(221), plot(diag(s)); title('Eigenvalues')

    subplot(222)
    plot((junkcomps)), title (sprintf('First %d components',Ncomps));
    set(gcf,'Name', 'SVD identified Noise components')
end

bhat = pinv(junkcomps)*dirty;
clean = dirty - junkcomps*bhat;
cc_sigma = std(clean,1);

% calculate the BIC for this noise model
RSS = sum(clean.^2,1);
n = size(clean,1);
k = size(junkcomps,2);
BIC = n*log(RSS/n) + k*log(n);
BIC = BIC .*msk;
BIC(BIC==0) = nan;

% reshape the output to its original shape
clean = clean';
clean=reshape(clean, dims(1), dims(2), dims(3), dims(4));


if (showFigs)
    subplot(224)
    title('BIC at each voxel');
    subplot(223)
    hist((BIC(:)),50); title(sprintf('BIC histogram. Mean %f', mean(BIC(~isnan(BIC)))));

    figure;
    
    subplot(211)
    lightbox(reshape(sigma,dims(1), dims(2), dims(3)));
    title('std. dev. map BEFORE')

    subplot(212)
    lightbox(reshape(cc_sigma,dims(1), dims(2), dims(3)));
    title('std. dev. map AFTER');
    set(gcf,'Name', 'Compcor Reduction in Noise')
end

return

function msk = ccvarmask(varimage , th)
% makes a mask that preserves the data with the top TH percentage of the
% values

ordered = sort(varimage(:));
Nintgrl = cumsum(ordered)/sum(ordered(:)) * 100;
thval = find(Nintgrl>100-th);
%subplot(211), plot(ordered); title('Ordered Std. Deviations')
%subplot(212), plot(Nintgrl); title('Intrageted , Normalized Std. Deviations')

thval = ordered(thval(1))
msk = varimage;
msk(msk < thval) = 0;
msk(msk>=thval) = 1;
return

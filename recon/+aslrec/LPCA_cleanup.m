function [out, Ncc] = LPCA_cleanup(ims, nbrs, Ncomps, X)
%function [out, Ncc] = LPCA_cleanup(ims, nbrs, Ncomps [,X] )
% 
% Local PCA clean up of a time series with heavy noise
%
% the clean version of each pixel's time course
% is a weighted sum of the main temporal components in the neighborhood patch
% find components in each patch thru SVD
% test the components for correlation with design matrix (X)
% estimate the coefficients of the components's contribution to the timecourse
% (using least squares)
% and then remove them from the time course
% 

if nargin<4
    X=[];
end

dims = size(ims);

Npix = (nbrs*2+1)^3;
NpixTotal = dims(1)*dims(2)*dims(3);
Nframes = dims(4);
progr=0;
msg = [];
fprintf('\nBegin LPCA ...\n')
Nc = Ncomps;

% Make casorati matrix from 4d series
% output is nt x npixels
im_txp = reshape(permute(abs(ims),[4 1:3]),[],NpixTotal);
out = zeros(size(im_txp));

% create a mask for the patches
vmask = zeros((2*nbrs+1)^3, 3);
r=1;
for m = -nbrs:nbrs
    for n = -nbrs:nbrs
        for p = -nbrs:nbrs
            vmask(r,:) = [m n p];
            r=r+1;
        end
    end
end


% go through all voxels in casorati matrix and extract the patch in the
% neightbood
progr = 0;
Ncc=zeros(NpixTotal,1);
parfor r = 1 : NpixTotal
        
    % determine locations of voxels in a patch around each position
    [x,y,z] = ind2sub(dims(1:3), r) ;
    pos = [x y z]+ vmask;

    % make sure we're within the boundaries
    if sum(find(pos(:) >  dims(1))) == 0 && ...
            sum(find(pos(:) < 1 )) == 0

        % translate into indices for casorati matrix
        locs = sub2ind(dims(1:3), pos(:,1), pos(:,2), pos(:,3));

        % extract the patch from the casorati matrix
        cpatch = (im_txp(:,locs))';

        % decompose the patch into its components
        [u, s, v] = svd(cpatch,'econ');

        % consider only the biggest temporal components
        comps = v(:,1:Ncomps);

        % Test the components
        if ~isempty(X)

            % find if components are correlated with design matrix.
            % if that it's the case, we will not remove them.
            % ie - remove them from the junk component pile

            coefs = corr(comps,X);
            coefs = max(abs(coefs),[],2);
            badcomps = find(abs(coefs) > 0.1);
            comps(:,badcomps) = [];

            % diagnostic: track the number of components left
            Ncc(r) = size(comps,2);

            % fprintf('\n%d components correlated with design', length(badcomps))
            % de-correlate the components from the design matrix in FMRI
            % This seems like a bad idea.  It biases the data introducing
            % correlations.
            % for j=1:Nc-1
            %     comps(:,j)= comps(:,j) - X*pinv(X)*comps(:,j);
            % end
        end

        % IMPORTANT: zero mean all the junk regressors
        comps = comps - repmat(mean(comps,1), size(comps,1),1);

        % estimate and remove components from the signal
        tseries = im_txp(:,r);
        bhat = pinv(comps)*tseries;
        out(:,r) = tseries - comps*bhat;

        if mod(r,floor(NpixTotal/1000))==0
            % msg =sprintf('\n Progress: %0.2f %%  \n ', 100*r/NpixTotal);
            % fprintf(repmat('\b', 1, length(msg))),
            % fprintf(msg);
            fprintf('*')
        end
    end
end
out = reshape(out',dims);
Ncc = reshape(Ncc,dims(1:3));
%lbview(Ncc); title('N. Components'); colorbar
%pause(1)

% the edges may be wrong, so remove them
out(1:nbrs,:,:,:) =0;
out(end-nbrs:end,:,:) =0;
out(:, 1:nbrs,:,:) =0;
out(:, end-nbrs:end,:) =0;
out(:,:, 1:nbrs,:) =0;
out(:,:, end-nbrs:end,:) =0;

fprintf('\n...LPCA complete\n')
return
function [dat hdr] = loadScanArchive(fname)
% function [dat hdr] = loadScanArchive(fname)
%
% Output:
%  dat  [ndat nechoes nframes ncoils]
 
if isempty(which('GErecon'))
    fprintf('\nAdding Orchestra to the matlab path ...\n')
    addpath ~/matlab/orchestra-sdk-2.1-1.matlab
end
if(isstruct(fname))
    fname = fname.name;
end
archive = GERecon('Archive.Load', fname);   %'ScanArchive_7347633TMRFIX_20201010_170753212.h5');

etl = archive.DownloadData.rdb_hdr_image.echo_trn_len;
hdr =archive.DownloadData.rdb_hdr_image;

if nargin < 2
    frameStart = 1;
    frameStop = archive.FrameCount;
end
 
% initialize data output matrix
currentControl = GERecon('Archive.Next', archive);

tmp = currentControl.Data;   % size = [ndat ncoils]
[ndat ncoils] = size(tmp);
fprintf('Initializing output data matrix, size [%d %d %d]\n', ndat, ncoils, archive.FrameCount);
dat = zeros(ndat, ncoils, archive.FrameCount);
 
% load data
dat(:,:,1) = currentControl.Data; 
for ii = 2:archive.FrameCount
    if ~mod(ii,1000)
        for ibs = 1:30; fprintf('\b'); end;
        fprintf('Frame %d of %d', ii, archive.FrameCount);
    end
    currentControl = GERecon('Archive.Next', archive);
    dat(:,:,ii) = currentControl.Data;   % size = [ndat ncoils nframes]
end

% split it into echos per frame and reorder as 
% ndat x nechoes x nframes x ncoils
sz = size(dat);
dat = reshape(dat,[sz(1), sz(2), etl, sz(3)/etl]);
dat = permute(dat,[1,3,4,2]);

% Not clear why, but need to clip one from each echo
dat = dat(1:end-1,:,:,:);
return

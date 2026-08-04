function [kdata,klocs,N,fov] = read_data(filename, mrf_mode)
% function [kdata,klocs,N,fov] = read_data(filename, mrf_mode)
% 
% Function to read in the filename and suppporting .txt files 
% ktraj*.txt and kviews*.txt  
% Function formats raw data and kspace locations for recon
%
% filename:  name of the filename containing the raw data
%       6.25.2026: adding support for scanarchives
% mrf_mode: binary flag.  1 means that we the rotation matrices will change from
%       frame to frame.  In that case, there will be different klocs for every
%       frame in the time series.
%

    filetype = 0;  % default filetype =0 mean P file, 
                   % use 1 for Scan Archives
    
    if nargin < 1 || isempty(filename)
        filename = './P*.7'; % default: use first filename on current path
        tmp = dir(filename);
  
        if isempty(tmp)
            fprintf('no filenames found from search string: %s, will try scanArchive\n', filename);
            filetype=1;
            filename = './*.h5'; % default: use first filename on current path
            tmp = dir(filename);
            if isempty(tmp)
                error('no scan Archives found from search string: % either!\n', filename);
            end
        end

    end
    

    if ~exist('mrf_mode')
        mrf_mode=0;
    end

    if filetype==0
        % reading the P file
        filename = tmp(1).name;
        pdir = tmp(1).folder;
        [raw,hdr] = aslrec.ge.read_pfile([pdir,'/',filename]);
        raw(:,all(raw == 0, [1,3:5]),:,:,:) = []; % remove empty views
        raw(:,:,all(raw == 0, [1:2,4:5]),:,:) = []; % remove empty frames
        nviews = size(raw,2);
        nframes = size(raw,3);
        ncoils = size(raw,5);
        kdata = reshape(raw,[],nviews,nframes,ncoils);

        % save matrix size (N) and fov
        N = hdr.image.dim_X * ones(1,3);
        fov = hdr.image.dfov/10 * ones(1,3);

    else
        % Reading a Scan Archive
        filename = tmp(1).name;
        pdir = tmp(1).folder;
        [kdata hdr] = aslrec.loadScanArchive(tmp);
        nviews = size(kdata,2);
        nframes = size(kdata,3);
        ncoils = size(kdata,4);

        fov = hdr.dfov/10 * ones(1,3);
        N = hdr.dim_X * ones(1,3);
    end

    % find and read the ktraj file
    tmp = dir([pdir,'/ktraj*.txt']);
    ktrajfile = tmp(1).name;
    klocs0 = load(ktrajfile);
    
    % find and read the kviews file
    tmp = dir([pdir,'/kviews*.txt']);
    kviewsfile = tmp(1).name;
    kviews = load(kviewsfile);
 
     % make modifications to the sizes if we're in MRF mode.
    if (mrf_mode==1)
        nviews = size(kviews,1);
    end

    % transform kspace locations using rotation matrices
    klocs = zeros(size(klocs0,1),3,nviews); % klocs = [N x 3 x nviews]
    for viewn = 1:nviews
        R = reshape(kviews(viewn,end-8:end)',3,3)';
        klocs(:,:,viewn) = klocs0*R';
    end
    klocs = permute(klocs,[1,3,2]); % klocs = [N x nviews x 3]
    
    
end

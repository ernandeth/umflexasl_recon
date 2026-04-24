flags.ResFac = 256/64;
flags.mip_threshold = 1e-1;
flags.ccf = 8; 
flags.narms = 1;
flags.skipNechoes = 0;  % 5
flags.skipNechoesEnd = 0;
flags.showEchoes = 0;
flags.sv = [];
flags.rm_navs = 1;
flags.doCompCorN = 0;
flags.doRecon = 3 ;
flags.doRegression = 0;
flags.doMoCo = 1;
flags.niter = 5;
flags.hp_filter = 0;
flags.outDir = 'results_niter05_compcor/';
flags.mrf_mode = 1;
flags.doK0correct = 0;
flags.doSmooth = 1;
flags.doCompCorN = 25;


mms = 0;
tseries = [];
mseries = [];

close all
curDir = pwd;
% Names of the directories of the data we want to recon (not the sense
% maps)
d = dir(['asldata_vsas*'])

% directory containing the sense map data
smap_dir = 'asldata_smap'

if flags.doRecon>0
    if isempty(gcp('nocreate'))
        parpool(6)
    end
end

if flags.doRecon>0

    cd(smap_dir)
    if exist('smap.mat')
        % If the maps are already calculated, then we just load them.
        fprintf('Found smap.mat file.  Loading sense maps ...\n')
        load smap.mat
    else
        % calculate the sense maps
        fprintf('Calculating sense maps ... \n')
        % first recon the coil images

        imc = recon3dflex('k0correct',flags.doK0correct,...
            'selectviews', [],...
            'niter', 25, 'Reg',100 , 'coilwise',1 );

        writenii('coils.nii', abs(imc));

        % calculate the maps from the coil images
        smap_mirt = mri_sensemap_denoise(imc, 'niter', 30, ...
            'bodycoil_default','ssos-angle-pca1');

        %  The BART alternative:
        % kc = fftc(imc, 1:3); % bart requires coil-wise kspace - so FT imc along dimensions 1-3
        % smap = bart('ecalib -b0 -m1', kc);

        save smap.mat smap_mirt imc
        writenii('smaps.nii', abs(smap_mirt));

        save smap.mat smap_mirt imc

    end

    % it's helpful to smooth the sense maps
    Sres = size(smap_mirt,1);
    ncoils = size(smap_mirt,4);

    fprintf('\nsmoothing SENSE maps  ... \n')
    for n=1:size(smap_mirt,4)
        smap_mirt(:,:,:,n) = smooth3(smap_mirt(:,:,:,n),'gaussian');
    end
    % display the sense maps
    figure
    for n=1:ncoils
        orthoview(abs(smap_mirt(:,:,:,n))), pause(0.1)
        title(['Coil ' num2str(n) ])
        drawnow
    end


    cd(curDir)
end

%%
if flags.ResFac>1 & flags.doRecon==3
    % interpolate B1 maps
    fprintf("\nInterpolating the sense maps to match the resolution of the images ...")
    ncoils = size(smap_mirt,4);
    Sres = size(smap_mirt,1);
    mtx = ceil(Sres*flags.ResFac);
    x = linspace(0,1, Sres);
    [X,Y,Z] = meshgrid(x,x,x);
    xq = linspace(0,1,mtx);
    [Xq,Yq,Zq] = meshgrid(xq,xq,xq);

    tmp = zeros(mtx, mtx, mtx,ncoils);

    figure
    parfor n=1:size(smap_mirt,4)

        fprintf('interpolating B1 map for coil %d\n', n);
        tmp(:,:,:,n) = interp3(X,Y,Z,smap_mirt(:,:,:,n), Xq, Yq, Zq);

        subplot(211)
        orthoview(smap_mirt(:,:,:,n))
        title(['Coil ' num2str(n) ])

        subplot(212)
        orthoview(tmp(:,:,:,n)), pause(0.1)

        drawnow
    end
    smap_mirt = tmp;
end
%}
%%
% Now do the reconstructions for each of the folders:
for i = 1:length(d)
    save proc_flags.mat flags

    cd(d(i).folder);
    cd(d(i).name);

    % figure out which views to select for the recon (flags.sv)
    [kdata,klocs,N,fov] = aslrec.read_data([],0);
    ndat = size(kdata,1);
    nechoes = size(kdata,2);
    ncoils = size(kdata,3);
    etl =  nechoes/flags.narms;
    
    switch (flags.doRecon)
        case 0  % assume it's already reconned.  Just load the images

            str = ['!mv ' flags.outDir '/* .']
            eval(str)

            if ~isempty(dir('im3d.mat'))
                load im3d.mat
                ims = im3d;
            else
                ims =readnii('im_mag.nii').* exp(1i*readnii('im_ang.nii'));
            end

        case 3  % use the sense map from a separate acquisition

            ims = recon3dflex( ...
                'k0correct',flags.doK0correct, ...
                'niter',flags.niter , ...
                'Reg',1e2, ...
                'ccfac', flags.ccf, ...
                'selectviews', flags.sv,...
                'smap', smap_mirt, ...
                'mrf_mode',flags.mrf_mode, ...
                'hp_filter', flags.hp_filter, ...
                'nonavs', flags.rm_navs);

            writenii('im_mag', abs(ims),'precision','float32','doscl',0);
            writenii('im_ang', angle(ims),'precision','float32','doscl',0);

    end
    save proc_options.mat flags

    fprintf('\nRecon step is Done. Begin  Post Processing ...')

    if flags.doCompCorN
        fprintf('\nCleaning Up time series with Compcor (Real and Imag, separately)...')

        % These the time  course components we want to keep!
        X = ones(size(ims,4), 2);
        X(1:2:end) = -1;

        [tmp_r junkcomps]= aslrec.compcor23(real(ims), flags.doCompCorN, X);
        [tmp_i junkcomps]= aslrec.compcor23(imag(ims), flags.doCompCorN, X);
        ims = complex(tmp_r, tmp_i);
    end

   
    % smooth the data
    if flags.doSmooth
        fprintf('\nGaussian Smoothing time series ...')
        for n=1:size(ims,4)
            ims(:,:,:,n) = smooth3(ims(:,:,:,n),'gaussian',5);
        end
    end
 
    if flags.doMoCo
        fprintf("\nRealign complex  images to the mean volume ...");
        ims = aslrec.complex_realign4D(ims);
    end

    fprintf('\nFlipping Images along z-axis... ');
    ims = ims(:,:,end:-1:1,:);

    fprintf('\Calculating Means and Subtracting contol from label... ');

    m = abs(mean(ims,4));
    s1 = mean(ims(:,:,:,1:2:end), 4);
    s2 = mean(ims(:,:,:,2:2:end), 4);

    s = abs(s1 - s2);
    ms = (mean(s,4));

    fprintf('\nWriting out the NIFTI files Means and Subtracting contol from label... ');

    writenii('mean_im', (m),'precision','float32','doscl',0);
    writenii('mean_sub', (ms),'precision','float32','doscl',0);


    figure
    subplot(411)
    orthoview(ims,'frame',1)
    title(d(i).name, 'Interpreter','none')

    subplot(412)
    orthoview(m);
    title('Mean over time')

    subplot(413)
    orthoview(ms )
    caxis([-1 1]*1)
    title('mean subtraction')

    if flags.showEchoes
        % Show the echo train in the first frame
        % also show the k0 corrected version for comparison... some times
        % it helps
        if ~isempty(flags.sv)
            kdata = kdata(:,flags.sv,:,:);
        end

        etl1 = kdata(:,:,1,1);
        etl1 = etl1(:);

        ndat = size(kdata,1);
        nechoes = size(kdata,2);
        ncoils = size(kdata,3);

        kdatac = aslrec.k0correct(kdata, klocs,1);
        etlc = kdatac(:,:,1,1);
        etlc = etlc(:);

        subplot(414)
        plot(abs(etlc(1:ndat*nechoes)));
        hold on
        plot(abs(etl1(1:ndat*nechoes))); hold off
        legend('corrected', 'uncorrected')
        title(pwd, 'Interpreter','none', 'FontSize',8)

    end
    
    if flags.mip_threshold >0
        subplot(414)

        ms=(ms);
        a =make_mip((ms), flags.mip_threshold, 1,0);
        b =make_mip((ms), flags.mip_threshold, 2,0);
        c =make_mip((ms),flags.mip_threshold, 3,0);
        imagesc([a' b' c'])
        axis image
        axis xy

        title('mean subtraction (MIP)')
    end

    fprintf('\nDone with calcualtations, Moving images and logs to reasonable locations')
    drawnow
    print -dpng summary

    flags
    save proc_flags.mat flags
    
    str = ['!rm proc* results/']; eval(str)
    str = ['!mkdir ' flags.outDir]; eval(str)
    str = ['!mv *.png mean* proc_flags.mat im_* ' flags.outDir]; eval(str)
    disp(str)
end

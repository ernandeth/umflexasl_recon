function [out ] = complex_apply_xform(in, rotmats)
% function [out rotmats] = complex_apply_xform(in, rotmats)
%
% use matlab built in functions to realign a series of 3D images using
% rigid body  transformations, when the 3D images are complex.
% The transfromations have already been calculated!
% - useful for motion correction in segmented data

re = real(in);
im = imag(in);

ref = abs(mean(in,4));



out = zeros(size(in));

parfor n=1:size(in, 4)

    fprintf('\nRealigning ... %d -th frame ', n);
    rtmp = squeeze(re(:,:,:,n));
    itmp = squeeze(im(:,:,:,n));

    tmp = squeeze(im(:,:,:,n));
    

    % calculate the affine transformation matrix
    tform = rigidtform3d(rotmats(:,:,n));

    % Create an output view with the same size/limits as the input
    % The output image will be clipped if the transformation causes content to go outside these limits
    OutView = affineOutputView(size(ref), tform, "BoundsStyle", "SameAsInput");

    rout = imwarp(rtmp, tform, "OutputView", OutView);
    iout = imwarp(itmp, tform, "OutputView", OutView);

    out(:,:,:,n) = rout + 1i*iout;
end


return
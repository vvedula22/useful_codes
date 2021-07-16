clear all; close all; clc; %#ok<CLALL>
t0 = cputime;

%% add paths
srcDir = '/media/vvedula/Data0/Projects/06-Zebrafish/study_02/New/00-ImRegister';
addpath(genpath(srcDir));

projDir = sprintf('%s/../03-rescue/01-4dpf/fish01',srcDir);
cd(projDir);

%% set MIRT options
% Main settings
main.similarity = 'ssd';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI 
main.subdivide  = 4;       % use 3 hierarchical levels
main.okno       = 5;       % mesh window size
main.lambda     = 0.1;     % transformation regularization weight, 0 for none
main.single     = 1;       % show mesh transformation at every iteration

% Optimization settings
optim.maxsteps = 100;      % maximum number of iterations at each hierarchical level
optim.fundif   = 1e-8;     % tolerance (stopping criterion)
optim.gamma    = 1;        % initial optimization step size 
optim.anneal   = 0.8;      % annealing rate on the optimization step    

iso = 0;
ibeg = 1;
iend = 1;
idel = 1;

%% load template stl surface
nsrf = 3;
for isrf=1:nsrf
    fname = sprintf('%s/03-case-setup/03.2-register_seq/%02d_mshsrf_%02d.vertices',projDir,iso,isrf);
    fid = fopen(fname,'r');
    srf(isrf).vertices = fscanf(fid,'%f');
    fclose(fid);
    srf(isrf).vertices = reshape(srf(isrf).vertices, [3, numel(srf(isrf).vertices)/3])';

    fname = sprintf('%s/03-case-setup/03.2-register_seq/%02d_mshsrf_%02d.faces',projDir,iso,isrf);
    srf(isrf).faces = load(fname);
    srf(isrf).faces = srf(isrf).faces(:,2:4) + 1;
    srf(isrf).faces = int32(srf(isrf).faces);
    srf(isrf).disps = zeros(size(srf(isrf).vertices));
    
    srf0(isrf).vertices = srf(isrf).vertices;
    srf0(isrf).faces = srf(isrf).faces;
    srf0(isrf).disps = zeros(size(srf0(isrf).vertices));
end

%% load template (start/source) image
fname = sprintf('%s/02-cropped/cropped_%d.vtk',projDir,iso);
info = vtk_read_header(fname);
nx = info.Dimensions(1);
ny = info.Dimensions(2);
nz = info.Dimensions(3);

orig.x = info.Origin(1);
orig.y = info.Origin(2);
orig.z = info.Origin(3);

spacing.x = info.PixelDimensions(1);
spacing.y = info.PixelDimensions(2);
spacing.z = info.PixelDimensions(3);

im = double(vtk_read_volume(fname));
im = im(:,:,:);

imax = max(im(:));
imin = min(im(:));
im = (im-imin)./(imax-imin);
im = permute(im, [2 1 3]);

%% load ref images and register source with each of them using MIRT
for i=ibeg:idel:iend
    cd(projDir);
    if (i==iso)
        for isrf=1:nsrf
            fname = sprintf('%s/03-case-setup/03.2-register_seq/01-morphed-mshsrf_%02d/%02d_registered',projDir,isrf,i);
%            fprintf('Writing stl file, mesh surface %d..\n', isrf);
%            stlwrite([fname '.stl'],srf(isrf),'mode','binary');
            fprintf('Writing vtk file, mesh surface %d..\n', isrf);
            writeVTK(srf(isrf), [fname '.vtk'], 'a');
        end
    else
        fprintf('**************************************\n');
        fprintf('\nRegistering image %d to image %d..\n\n',i-idel,i);
        fname = sprintf('%s/02-cropped/cropped_%d.vtk',projDir,i);
        refim = double(vtk_read_volume(fname));
        refim = refim(:,:,:);
        refim = permute(refim,[2 1 3]);
        refim = (refim-imin)./(imax-imin);
        refim(refim<0) = 0;
        refim(refim>1) = 1;

        t1 = cputime;
        [res, newim] = mirt3D_register(refim, im, main, optim);
        t2 = cputime;
        fprintf('\nTime for registration: %.1f s\n',t2-t1);
        
        fprintf('\nWriting displacements to file...');
        fname = sprintf('%s/03-case-setup/03.2-register_seq/00-spline-disps/%02d_splineDisp',projDir,i);
        disp = res.X;
        save([fname '.mat'], 'disp');
        
        for isrf=1:nsrf
            fprintf('\nSurface morphing..\n');
            newsrf = morphSurf3D(srf(isrf), res, orig, spacing, newim);
            srf(isrf).vertices = newsrf.vertices;
            srf(isrf).disps = srf0(isrf).vertices - srf(isrf).vertices;

            cd(projDir);
            fname = sprintf('%s/03-case-setup/03.2-register_seq/01-morphed-mshsrf_%02d/%02d_registered',projDir,isrf,i);
%            fprintf('Writing stl file, mesh surface %d..\n', isrf);
%            stlwrite([fname '.stl'],srf(isrf),'mode','binary');
            fprintf('Writing vtk file, mesh surface %d..\n', isrf);
            writeVTK(srf(isrf), [fname '.vtk'], 'ascii');
        end
        im = refim;
    end
end
cd (srcDir);
tn = cputime;
fprintf('\nTotal elapsed time: %.1f s\n', tn-t0);
fprintf('\n**************************************\n');


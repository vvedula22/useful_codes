clear all; close all; clc; %#ok<CLALL>
t0 = cputime;

fishType = 1;   % 1: control; 2: smooth_AG;
fishNo = 1;

%% add paths
srcDir = '/media/vvedula/Data/Projects/06-Zebrafish/study_02/New/00-ImRegister';
addpath(genpath(srcDir));

switch fishType
    case 1
        projDir = sprintf('%s/../01-control/fish%02d',...
            srcDir,fishNo);
    case 2
        projDir = sprintf('%s/../02-AGtype/fish%02d',...
            srcDir,fishNo);
    otherwise
        error('Unknown fish type');
end
cd(projDir);

%% images to be registered
iso = 0;
itar = 17;
slice = 19;
ithresh = 100000;

% moving image (source)
fname = sprintf('%s/02-cropped/cropped_%d.vtk',projDir,iso);
info = vtk_read_header(fname);

orig.x = 1.0;
orig.y = 1.0;

spacing.x = 1.0;
spacing.y = 1.0;

im3d = double(vtk_read_volume(fname));
im3d = permute(im3d, [2 1 3]);
imax = max(im3d(:));
imin = min(im3d(:));
im3d(im3d>ithresh) = imax;
im3d = (im3d-imin)./(imax-imin);
im = im3d(:,:,slice);

% fixed image (reference/target)
fname = sprintf('%s/02-cropped/cropped_%d.vtk',projDir,itar);
refim3d = double(vtk_read_volume(fname));
refim3d = permute(refim3d, [2 1 3]);
refim3d(refim3d>ithresh) = imax;
refim3d = (refim3d-imin)./(imax-imin);
refim3d(refim3d<0) = 0;
refim3d(refim3d>1) = 1;
refim = refim3d(:,:,slice);

%% contour segmentation
figure('units','normalized','outerposition',[0 0 1 1]);
button= 'No';

while length(button)==2
    temp = im;
    imshow(temp); [x,y] = ginput;
    srf.vertices = createCurve(x,y,100);
    hold on; plot(srf.vertices(:,1),srf.vertices(:,2),'ro'); hold off;
    button = questdlg('Is the curve OK?');
end
if (length(button)>3)
    close all;
    return;
end

%% Registration
uiwait(msgbox('Start registration'));

close all;

% Main settings
main.similarity = 'ssd';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI 
main.subdivide  = 4;       % use 3 hierarchical levels
main.okno       = 5;       % mesh window size
main.lambda     = 0.02;   % transformation regularization weight, 0 for none
main.single     = 1;       % show mesh transformation at every iteration
    
% Optimization settings
optim.maxsteps = 100;      % maximum number of iterations at each hierarchical level
optim.fundif   = 1e-8;     % tolerance (stopping criterion)
optim.gamma    = 1;        % initial optimization step size 
optim.anneal   = 0.8;      % annealing rate on the optimization step    

fprintf('**************************************\n');
fprintf('\nRegistering slice %d of image %d to image %d..\n\n',...
    slice,iso,itar);
    
[res, newim] = mirt2D_register(refim, im, main, optim);

fprintf('\nSurface morphing..\n');
newsrf = morphSurf2D(srf, res, orig, spacing);

close all;
%figure('units','normalized','outerposition',[0 0 1 1]);
figure; imshow(im); title('Source Image'); hold on;
plot(srf.vertices(:,1),srf.vertices(:,2),'yo',...
    'MarkerFaceColor','y'); hold off;
figure; imshow(refim); title('Reference Image');
figure; imshow(refim); title('Reference Image'); hold on;
plot(newsrf.vertices(:,1),newsrf.vertices(:,2),'yo',...
    'MarkerFaceColor','y'); hold off;
figure; mirt2D_meshplot(res.X(:,:,1),res.X(:,:,2)); title('Deformation');
figure; imshow(newim); title('Registered Image'); hold on;
plot(newsrf.vertices(:,1),newsrf.vertices(:,2),'yo',...
    'MarkerFaceColor','y'); hold off;

cd (srcDir);
tn = cputime;

fprintf('\nTotal elapsed time: %.1f s\n', tn-t0);
fprintf('\n**************************************\n');


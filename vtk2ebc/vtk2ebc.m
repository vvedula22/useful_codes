clear all; close all; clc;

fname = "lumen_flag_int.vtk";
fid = fopen(fname,'r');

fprintf('   Opened file %s...\n',fname);
dlms = '< =">';
nsd  = 2; % num spatial dimensions

%% POINTS
[tline, iok] = findKwrd(fid, dlms, 'POINTS');
if ( iok > 0 )
    [tokList, ntoks] = parseLine(tline, dlms);
    nNo = sscanf(tokList{2},'%d');
    fprintf('      Number of points: %d\n',nNo);
else
    return;
end

%% POLYGONS
[tline, iok] = findKwrd(fid, dlms, 'POLYGONS');
if ( iok > 0 )
    [tokList, ntoks] = parseLine(tline, dlms);
    nEl = sscanf(tokList{2},'%d');
    eNoN = sscanf(tokList{3},'%d');
    eNoN = (eNoN/nEl) - 1;
    fprintf('      Number of elements: %d\n',nEl);
    fprintf('      Number of nodes per element: %d\n', eNoN);
else
    return;
end

ien  = int32(zeros(nEl, eNoN));
tmpI = zeros(1,eNoN+1);
fprintf('      Reading element connectivity...');
for e=1:nEl
    if (feof(fid))
        fprintf('\tError: EOF reached\n');
        fprintf('\tSome data missing\n');
        return;
    end
    tline = fgetl(fid);
    tmpI(:) = sscanf(tline,'%d %d %d');
    ien(e,:) = tmpI(2:eNoN+1)+1;
end

%% GlobalElementID
[tline, iok] = findKwrd(fid, dlms, 'CELL_DATA');
fprintf('      Reading element IDs...\n');
fgetl(fid);  % SCALARS
fgetl(fid);  % LOOKUP_TABLE
nlines = floor(nEl/9);
nrem = rem(nEl,9);
gE = zeros(nEl,1);
indx = 0;
a = zeros(9,1);
for i=1:nlines
    tline = fgetl(fid);
    a(1:9) = sscanf(tline,'%d');
    for j=1:9
        indx = indx + 1;
        gE(indx) = a(j);
    end
end
a = zeros(nrem,1);
tline = fgetl(fid);
a(1:nrem) = sscanf(tline,'%d');
for j=1:nrem
    indx = indx + 1;
    gE(indx) = a(j);
end

%% GlobalNodeID
[tline, iok] = findKwrd(fid, dlms, 'POINT_DATA');
fprintf('      Reading nodal IDs...\n');
fgetl(fid);  % SCALARS
fgetl(fid);  % LOOKUP_TABLE
nlines = floor(nNo/9);
nrem = rem(nNo,9);
gN = zeros(nNo,1);
indx = 0;
a = zeros(9,1);
for i=1:nlines
    tline = fgetl(fid);
    a(1:9) = sscanf(tline,'%d');
    for j=1:9
        indx = indx + 1;
        gN(indx) = a(j);
    end
end
a = zeros(nrem,1);
tline = fgetl(fid);
a(1:nrem) = sscanf(tline,'%d');
for j=1:nrem
    indx = indx + 1;
    gN(indx) = a(j);
end

fclose(fid);

%% Write ebc file
[~,fname,~] = fileparts(fname);
fname = strcat(fname, '.ebc');
fprintf('   Writing ebc file %s\n', fname);
fid = fopen(fname,'w');
for e=1:nEl
    fprintf(fid,'%d   %d',gE(e), e);
    for a=1:eNoN
        Ac = ien(e,a);
        fprintf(fid,'   %d',gN(Ac));
    end
    if (e < nEl)
        fprintf(fid,'\n');
    end
end
fclose(fid);



function [coords, ien, ptData] = readVTK(fname)

fprintf('----------------------------------------------------------\n');
fid = fopen(fname,'r');
fprintf('   Opened file %s...\n',fname);
dlms = '< =">';
nsd  = 3; % num spatial dimensions
eNoN = 3; % triangles

[tline, iok] = findKwrd(fid, dlms, 'ASCII');
if ~(iok>0)
    fprintf('\tError: only ASCII vtk file can be read\n');
    return
end

%%POINTS
[tline, iok] = findKwrd(fid, dlms, 'POINTS');
if ( iok > 0 )
    [tokList, ntoks] = parseLine(tline, dlms);
    nNo = sscanf(tokList{2},'%d');
    fprintf('      Number of points: %d\n',nNo);
else
    return;
end

coords = zeros(nNo, nsd);
fprintf('      Reading point coordinates...');
for a=1:nNo
    if (feof(fid))
        fprintf('\tError: EOF reached\n');
        fprintf('\tSome data missing\n');
        return;
    end
    tline = fgetl(fid);
    coords(a,:) = sscanf(tline,'%f %f %f');
end
fprintf(' Done!\n');

%% CELLS
[tline, iok] = findKwrd(fid, dlms, 'CELLS');
if ( iok > 0 )
    [tokList, ntoks] = parseLine(tline, dlms);
    nEl = sscanf(tokList{2},'%d');
    fprintf('      Number of elements: %d\n',nEl);
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
    tmpI(:) = sscanf(tline,'%d %d %d %d');
    ien(e,:) = tmpI(2:eNoN+1)+1;
end
fprintf(' Done!\n');

%% POINT_DATA
[tline, iok] = findKwrd(fid, dlms, 'POINT_DATA');
if ~(iok > 0)
    fprintf('\tError: missing POINT_DATA field\n');
    return
end

[tline, iok] = findKwrd(fid, dlms, 'VECTORS');
if ( iok > 0 )
    [tokList, ntoks] = parseLine(tline, dlms);
    stmp = sscanf(tokList{2},'%s');
    if strcmp(stmp,'Displacement')
        ptData = zeros(nNo, nsd);
        fprintf('      Reading point coordinates...');
        for a=1:nNo
            if (feof(fid))
                fprintf('\tError: EOF reached\n');
                fprintf('\tSome data missing\n');
                return;
            end
            tline = fgetl(fid);
            ptData(a,:) = sscanf(tline,'%f %f %f');
        end
        fprintf(' Done!\n');
    else
        ptData = zeros(nNo, nsd);
        return
    end
else
    fprintf('\tError: missing VECTORS field\n');
    return;
end

fclose(fid);

return;
end
clear; clc;

%% parameters to edit

% column spacing (x- spacing)
col_spacing = 0.65;

% row spacing (y- spacing)
row_spacing = 0.65;

% slice thickness (z- spacing)
slice_thickness = 2.0;

%% converting tif to vti
[machine,maxSize,endian]=computer;
if endian=='L'
    endian='ieee-le';
    endianText='LittleEndian';
else
    endian='ieee-be';
    endianText='BigEndian';
end

for iTime = nStart:nFreq:nEnd
    fprintf('%d\n',iTime);
    for k=1:1:slices
        fname = sprintf('%s/%s%d.tif',dir,fhdr,iTime);
        I = imread(fname,'tiff',k);
        [rows,cols] = size(I);
        if (k==1)
            A = uint16(zeros(rows,cols,slices));
        end
        A(:,:,k) = I;
    end
    fname = sprintf('zebrafish3D_%03d.vti',iTime);
    fid = fopen(fname,'w',endian);
    fprintf(fid,'%s%s%s\n',...
        '<VTKFile type="ImageData" version="1.0" byte_order="',...
        endianText,'" header_type="UInt64">');
    fprintf(fid,'%s%d%s%d%s%d%s',...
        '  <ImageData WholeExtent="0 ',(cols-1),...
        ' 0 ',(rows-1),' 0 ',(slices-1),'"');
    fprintf(fid,'%s%f %f %f%s\n',...
        ' Origin="0 0 0" Spacing="',col_spacing,row_spacing,slice_thickness,'">');
 
    fprintf(fid,'%s%d%s%d%s%d%s\n',...
        '    <Piece Extent="0 ',(cols-1),...
        ' 0 ',(rows-1),' 0 ',(slices-1),'">');
    fprintf(fid,'%s\n',...
        '      <PointData Scalars="scalars">');
    fprintf(fid,'%s%s%d%s%d%s\n',...
        '       <DataArray type="UInt16" Name="scalars"',...
        ' format="appended" RangeMin="',min(A(:)),'" RangeMax="',...
        max(A(:)),'" offset="0" />');
    fprintf(fid,'%s\n',...
        '      </PointData>');
    fprintf(fid,'%s\n',...
        '    </Piece>');
    fprintf(fid,'%s\n',...
        '  </ImageData>');

    fprintf(fid,'%s\n%s',...
        '  <AppendedData encoding="raw">','   _');

    fwrite(fid,2*numel(A),'uint64');
    fwrite(fid,reshape(A,1,numel(A)),'uint16','l');
 
    fprintf(fid,'\n%s\n%s',...
        '  </AppendedData>','</VTKFile>');
    fclose(fid);
end




function writeVTK(srf, fname, c)

    if (nargin < 3) 
        c = 'a';
    end
    
	nno = size(srf.vertices,1);
    nsd = size(srf.vertices,2);
	nel = size(srf.faces,1);
    nen = size(srf.faces,2);
    
    if nsd==3
        switch nen
        case {8} % brick
            vtktype = 12;
        case {6} % wedges
            vtktype = 13;
        case {4} % tets
            vtktype = 10;
        case {3}
            vtktype = 5;
        end
    elseif nsd==2
        switch nen
        case {4} % quads
            vtktype = 9;
        case {3} % triangles
            vtktype = 5;
        case {2} % lines
            vtktype = 3;
        end
    else
        error('VTK write error: unknown dimension');
    end
    
    switch c
    case {'b','B','binary'}
        if nsd==2
            fid = fopen(fname,'w');
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, '2D Unstructured Surface\n');
            fprintf(fid, 'BINARY\n');
            fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
            
            fprintf(fid, 'POINTS %d double\n', nno);
            points = double(reshape([srf.vertices zeros(nno,1)]',1,nno+numel(srf.vertices)));
            fwrite(fid, points, 'double', 'b');
            
            fprintf(fid, '\nCELLS %d %d\n',nel, nel*(nen+1));
            cells = int32(reshape([nen+zeros(nel,1), srf.faces-1]',1,nel+numel(srf.faces)));
            fwrite(fid, cells, 'int32', 'b');
            
            cell_types = int32(reshape(vtktype+zeros(nel,1),1,nel));
            fprintf(fid, '\nCELL_TYPES %d\n',nel);
            fwrite(fid, cell_types, 'int32', 'b');
            
            if isfield(srf, 'disps')
                fprintf(fid, 'POINT_DATA %d\n', nno);
                fprintf(fid, 'VECTORS Displacement double\n');
                disps = double(reshape([srf.disps zeros(nno,1)]',1,nno+numel(srf.disps)));
                fwrite(fid, disps, 'double', 'b');
            end
            fclose(fid);
            
        elseif nsd==3
            fid = fopen(fname,'w');
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, '3D Unstructured Surface\n');
            fprintf(fid, 'BINARY\n');
            fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
            
            fprintf(fid, 'POINTS %d double\n', nno);
            points = double(reshape(srf.vertices',1,numel(srf.vertices)));
            fwrite(fid, points, 'double', 'b');

            fprintf(fid, '\nCELLS %d %d\n',nel, nel*(nen+1));
            cells = int32(reshape([nen+zeros(nel,1), srf.faces-1]',1,nel+numel(srf.faces)));
            fwrite(fid, cells, 'int32', 'b');
            
            cell_types = int32(reshape(vtktype+zeros(nel,1),1,nel));
            fprintf(fid, '\nCELL_TYPES %d\n',nel);
            fwrite(fid, cell_types, 'int32', 'b');
            
            if isfield(srf, 'disps')
                fprintf(fid, 'POINT_DATA %d\n', nno);
                fprintf(fid, 'VECTORS Displacement double\n');
                disps = double(reshape(srf.disps',1,numel(srf.disps)));
                fwrite(fid, disps, 'double', 'b');
            end
            
            fclose(fid);
            clear points cells cell_types;
        else
            error('VTK error: unknown problem dimension');
        end
    case {'a','A','ascii'}
        if nsd==2
            fid = fopen(fname,'w');
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, '2D Unstructured Surface\n');
            fprintf(fid, 'ASCII\n');
            fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
            
            fprintf(fid, 'POINTS %d float\n', nno);
            for Ac=1:nno
                fprintf(fid, '%f %f 0.0\n',...
                    srf.vertices(Ac,1), srf.vertices(Ac,2));
            end
            
            fprintf(fid, 'CELLS %d %d\n',nel, nel*(nen+1));
            for e=1:nel
                fprintf(fid,'%d ', nen);
                for a=1:nen
                    fprintf(fid, ' %d',srf.faces(e,a)-1);
                end
                fprintf(fid,'\n');
            end
            
            fprintf(fid, 'CELL_TYPES %d\n',nel);
            for e=1:nel
                fprintf(fid, '%d\n', vtktype);
            end
            
            if isfield(srf, 'disps')
                fprintf(fid, 'POINT_DATA %d\n', nno);
                fprintf(fid, 'VECTORS Displacement float\n');
                for Ac=1:nno
                    fprintf(fid, '%f %f 0.0\n',...
                        srf.disps(Ac,1), srf.disps(Ac,2));
                end
            end
            
            fclose(fid);
            
        elseif nsd==3
            fid = fopen(fname,'w');
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, '3D Unstructured Surface\n');
            fprintf(fid, 'ASCII\n');
            fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
            
            fprintf(fid, 'POINTS %d float\n', nno);
            for Ac=1:nno
                fprintf(fid, '%f %f %f\n',...
                    srf.vertices(Ac,1), srf.vertices(Ac,2), srf.vertices(Ac,3));
            end
            
            fprintf(fid, 'CELLS %d %d\n',nel, nel*(nen+1));
            for e=1:nel
                fprintf(fid,'%d ', nen);
                for a=1:nen
                    fprintf(fid, ' %d',srf.faces(e,a)-1);
                end
                fprintf(fid,'\n');
            end
            
            fprintf(fid, 'CELL_TYPES %d\n',nel);
            for e=1:nel
                fprintf(fid, '%d\n', vtktype);
            end
            
            if isfield(srf, 'disps')
                fprintf(fid, 'POINT_DATA %d\n', nno);
                fprintf(fid, 'VECTORS Displacement float\n');
                for Ac=1:nno
                    fprintf(fid, '%f %f %f\n',...
                        srf.disps(Ac,1), srf.disps(Ac,2), srf.disps(Ac,3));
                end
            end
            
            fclose(fid);
            
        else
            error('VTK error: unknown problem dimension');
        end
    end
end

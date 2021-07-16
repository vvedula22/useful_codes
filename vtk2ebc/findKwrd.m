function [tline,ret] = findKwrd(fid,dlms,literal)
    ret = 0;
    tline = [];
    if ( isempty(fid) || fid<0 )
        fprintf('\tFile status unknown\n');
        return;
    end
    
    fprintf('   Searching for %s...',literal);
    tline = fgetl(fid);
    while ischar(tline)
        [tok,rem] = strtok(tline,dlms);
        if ( strcmp(tok,literal) )
            ret = 1;
            break;
        end
        tline = fgetl(fid);
        if ( feof(fid) )
            fprintf('\tError: end of file reached.\n');
            fprintf('\tKeyword %s not found',literal);
            ret=-1; return;
        end
    end
    fprintf(' Found!\n');
end

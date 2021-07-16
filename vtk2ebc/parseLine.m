function [tlist,cnt] = parseLine(tline,dlms)
    rem = tline;
    cnt = 0;
    while true
        [tok,rem] = strtok(rem,dlms);
        if ( isempty(tok) ), break; end;
        cnt = cnt+1;
        tlist{cnt} = tok;
    end
end

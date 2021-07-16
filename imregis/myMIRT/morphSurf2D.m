function nsrf = morphSurf2D(osrf, res, orig, spacing)
	F = mirt2D_F(res.okno);
    [Xx,Xy] = mirt2D_nodes2grid(res.X, F, res.okno);

    nx = size(Xx,2);
    ny = size(Xy,1);
    [x, y]  = meshgrid(1:nx,1:ny);

    dx = (Xx - x).*spacing.x;
    dy = (Xy - y).*spacing.y;

    x = (x-1).*spacing.x + orig.x;
    y = (y-1).*spacing.y + orig.y;

    disp(:,1) = interp2(x,y,dx,osrf.vertices(:,1),osrf.vertices(:,2),'linear',0);
    disp(:,2) = interp2(x,y,dy,osrf.vertices(:,1),osrf.vertices(:,2),'linear',0);

    if (nnz(isnan(disp))>0)
        fprintf('Warning: some part of surface lies outside of domain\n');
        disp(isnan(disp)) = 0;
    end

    nsrf.vertices = osrf.vertices - disp;
end

function output = createCurve(px, py, Nt)

% px = [22,26,34,43,51,52,46,35,26,19];
% py = [75,79,81,78,69,55,41,41,43,45];

px=px';  py=py';

smoothspline=rscvn([px; py;]);   %perform Hermite interpolation smoothing
points=fnplt(smoothspline);      %get resulting points from spline
px= points(1,:);  py=points(2,:);

% t is the cumulative arclength along the edges of the polygon.
t = cumsum(sqrt([0,diff(px(:)').^2] + [0,diff(py(:)').^2]));

% The total distance around the polygon is t(end)
tmax = t(end);

% create a piecewise linear spline for each of px and py,
% as a function of the cumulative chordwise arclength.
splx = mkpp(t,[diff(px(:))./diff(t'),px(1:(end-1))']);
sply = mkpp(t,[diff(py(:))./diff(t'),py(1:(end-1))']);

% now interpolate the polygon splines, splx and sply.
% Nt is the number of points to generate around the
% polygon. The first and last points should be replicates
% at least to within floating point trash.)
%Nt = 100;
tint = linspace(0,tmax,Nt);

qx = ppval(splx,tint);
qy = ppval(sply,tint);
output= [qx' qy'];

% % plot the polygon itself, as well as the generated points.
% plot(px,py,'k-v',qx,qy,'ro')
% grid on 
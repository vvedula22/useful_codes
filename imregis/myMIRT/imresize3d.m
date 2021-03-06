function Iq = imresize3d(I,mq,nq,pq)
    
    m = size(I,1);
    n = size(I,2);
    p = size(I,3);
    
    x = double(1:n);
    y = double((1:m)');
    z = double(1:p);
    
    xq = double(1:nq);
    yq = double((1:mq)');
    zq = double(linspace(1,p,pq));
    
    Iq = interp3(x,y,z,I,xq,yq,zq,'linear');
    
    return;
end 
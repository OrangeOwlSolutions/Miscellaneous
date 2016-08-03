function mapPolyFunc = mapLagrangeInterpFunction (x,y)
   
    N=length(x);
    mapPolyFunc=zeros(size(x));
    for i=1:N
        mapPolyFunc=mapPolyFunc+poly(x(2:N))./prod((x(1)-x(2:N)))*y(i);
        x=circshift(x.',-1).';
    end

end

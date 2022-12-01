function ret = nf( u, v )
    
    for i = 1:length(u)
        if (u(i) >= v(i))
            if ((u(i)+v(i))/2 > 0)
                ustar(i)=u(i);
            else
                ustar(i)=v(i);
            end
        elseif (u(i)>0)
            ustar(i)=u(i);
        elseif (v(i)<0)
            ustar(i)=v(i);
        else
            ustar(i)=0;
        end
    end
    
    ret =f(ustar);
end

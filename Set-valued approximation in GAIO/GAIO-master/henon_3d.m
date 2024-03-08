function xx = henon_3d(xx,a,b,epsilon,n)
    % 3-D henon extended
    
    if any(size(xx)==1)
        x = xx(1); y = xx(2); theta=xx(3);
    else
        x = xx(1,:); y = xx(2,:); theta=xx(3,:);
    end
    n1 = cos(theta);
    n2 = sin(theta);
    
    for i =1:n
        [x,y,n1,n2] = henon_extended(x,y,n1,n2,a,b,epsilon);
    end
    theta = atan2(n2,n1);
    xx = [x;y;theta];

end
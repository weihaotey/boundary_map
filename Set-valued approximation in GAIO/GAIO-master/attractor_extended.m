%% Attractor approx for extended map
% The H?non map $f(x,y)=(1-ax^+y,bx)$ is a well known example of a map
% which exhibits _complicated dynamics_ . Here is a plot of its attractor
% for \$a=1.4\$ and \$b=0.3\$.
clear all;
a = 0.61; b1 = 0.3; epsilon = 0.0625;                                     % parameters of the map
f = @(x) henon_3d(x,a,b1,epsilon,1);                 % Boundary H?non map

% Preparations
% In order to compute a covering of the attractor, we first choose sample points 
% in the square \$[-1,1]^\$.  Here, we chose these points on the edges of
% the square.
n = 5;x = linspace(-1,1,n);
[XX,YY,ZZ] = meshgrid(x,x,x);
X = [ XX(:) YY(:) ZZ(:) ];                            % sample points in \$[-1,1]^\$

% and initialize the tree using the square [-3,3]^2
%c = [vec3(1), vec3(2), atan2(vec3(4),vec3(3))]; r = [.8 .8 1.6];                         
c = [.5 .1 0]; r = [2,.7,pi];
t = Tree(c, r);


%% Subdivison algorithm
% We can now run the subdivision algorithm for the (relative) global
% attractor. 
viewon = 1;
writevideo = 0;
count = 1;
if writevideo
    vidfile1 = VideoWriter('2henon','MPEG-4');
    vidfile1.FrameRate = 2;
    open(vidfile1);
end
for a = linspace(0.6077,0.6077,1)
    b = 0.3; epsilon = 0.0625;                                     % parameters of the map
    f = @(x) henon_3d(x,a,b,epsilon,2);                 % Boundary H?non map
    tic;
    clear t
    c = [.5 .1 0]; r = [2,.7,pi];
    %c = [-.2 .4 0]; r = [.7,.2,pi];
    t = Tree(c, r);
    dim = t.dim; depth = 40;
    hit = 1; sd = 8;                                 % define flags
    for i = 1:depth                                  % subdivide up to the given depth
        if i < 27
            n = 8;
        else
            n = 5;
        end
        x = linspace(-1,1,n);
        [XX,YY,ZZ] = meshgrid(x,x,x);
        X = [ XX(:) YY(:) ZZ(:) ];                            % sample points in \$[-1,1]^\$
        
        disp(i)
        t.set_flags('all', sd);                      % flag all boxes for subdivision
        t.subdivide(sd);                             % subdivide all flaged boxes
        if i < 1                      
            continue
        end
        b = t.first_box(i);                          % loop over leaves of the tree
        while (~isempty(b))
            c = b(1:dim); r = b(dim+1:2*dim);        % center and radius of current box
            p = X*diag(r) + ones(size(X))*diag(c);   % sample points in current box
            fp = f(p');
            t.set_flags(fp, hit);                 % map points and flag hit boxes
            b = t.next_box(i);
        end
        t.remove(hit);                               % remove all boxes which have *not* been hit
        if i ~= 1
            close
        end
        %boxplot3(t); xlabel('x'); ylabel('y'); zlabel('z'); view(-40,20); grid on;
        %pause(.1)
        if viewon
            close
            boxplot3(t,'edgecolor','none','color','k'); xlabel('x'); ylabel('y'); zlabel('z'); view(0,90); grid on;
            view(0,90);
            axis('equal')
            drawnow;
        end
    end
    close
    %boxplot3(t,'edgecolor','none','color','k'); xlabel('x'); ylabel('y'); zlabel('z'); view(count*3,90); grid on;
    %axis('equal')
    %title(strcat('a =  ',num2str(a)))
    if writevideo, F(count) = getframe(gcf);writeVideo(vidfile1,F(count)); end
    toc;fprintf('\n')
    count = count + 1;
%    pause(.1)
end
if writevideo, close(vidfile1);end

%%
depth = 40;
for i = 36:depth                                  % subdivide up to the given depth
    disp(i)
    t.set_flags('all', sd);                      % flag all boxes for subdivision
    t.subdivide(sd);                             % subdivide all flaged boxes
    if i < 1                      
        continue
    end
    b = t.first_box(i);                          % loop over leaves of the tree
    while (~isempty(b))
        c = b(1:dim); r = b(dim+1:2*dim);        % center and radius of current box
        p = X*diag(r) + ones(size(X))*diag(c);   % sample points in current box
        fp = f(p');
        t.set_flags(fp, hit);                 % map points and flag hit boxes
        b = t.next_box(i);
        disp(length(b(1,:)))
    end
    t.remove(hit);                               % remove all boxes which have *not* been hit
    %boxplot3(t); xlabel('x'); ylabel('y'); zlabel('z'); view(-40,20); grid on;
    %pause(.1)
    if viewon
        if i ~= 1
            close
        end
        close
        boxplot3(t); xlabel('x'); ylabel('y'); zlabel('z'); view(0,90); grid on;
        view(0,90);
        drawnow;
    end
end
boxplot3(t); xlabel('x'); ylabel('y'); zlabel('z'); view(count*3,20); grid on;
title(strcat('a =  ',num2str(a)))

%% Plot of the box collection
boxplot3(t); xlabel('x'); ylabel('y'); zlabel('z'); view(40,20); grid on; 

%% Cleanup
delete(t);

%% faster way?
%n = 20; X1 = linspace(-1,1,n)'; E = ones(size(X1));     % 1d points
%X = [ X1 -E E; X1 E E; X1 -E -E; X1 E -E; -E X1 E; E X1 E; -E X1 -E; E X1 -E;...
%    X1 E -E; X1 E E; X1 -E -E; X1 -E E];                        % sample points in $[-1,1]^2$
viewon = 0;
dim = t.dim; hit = 1; sd = 8; depth = 30; tic
for i=1:depth
    disp(i)
    t.set_flags('all', sd);                             % flag all leaves for subdivision
    t.subdivide;                                        % subdivide flagged leaves
    b = t.boxes(-1); n = size(b,2);                     % get boxes from the leaves
    c = b(1:dim,:); r = b(dim+1:2*dim,1);               % centers and radii of the boxes
    E = ones(n,1);                                      
    P = kron(E,X)*diag(r) + kron(c',ones(size(X,1),1)); % sample points in all boxes
    t.set_flags(f(P'), hit);                            % map points and flag hitted boxes
    t.remove(hit);                                      % remove boxes which have not been hit
    if viewon
        clf;boxplot3(t); xlabel('x'); ylabel('y'); zlabel('z'); view(0,20); grid on;
        pause(.1)
    end
end
toc
boxplot3(t,'edgecolor','none','color','k'); xlabel('x'); ylabel('y'); zlabel('z'); view(0,20); grid on;
title(strcat('a =  ',num2str(a)))



%% Clean up
delete(t);




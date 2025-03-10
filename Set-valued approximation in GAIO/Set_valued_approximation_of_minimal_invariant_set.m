%% Search for minimal invariant set by first determining the recurrent set
%primary assumption: there exist minimal invariant set conatined wholly
%inside the first box
%Furthermore, the domain of attraction is points which goes to minimal
%invariant set without going out of the first big box

clear all ; tic;
stored_work = 0; fine = 0; % whether to use stored transition matrix
savework = 0; % whether to save the tree structure for plotting
saveimage = 0; % whether to save the final approximation figure
extended = 0; % use boundary mapping to find the boundary
finddual = 0; % whether to search for domain of attraction
minimal_invariant = 0; % whether to use saddle point to find dual repelelr of the smaller invariant set
only_boundary = 1; % whether to use only boundary for the epsilon ball
only_center = 1; % only use the center of the box as one sample point
multiplicative = 0; % whether to use multiplicative noise
plotdeter = 0; % whether to plot deterministic point
mapType = 'Hen'; % type of map in map.m
viewon = 1; % whether to plot when computing
t_all = {};

% sample points (uniform grid)
n = 5; x = linspace(-1,1,n)'; [XX,YY] = meshgrid(x); X = [ XX(:) YY(:)];
if only_center == 1, X = [0,0]; end
e_x = 0.6; e_y =e_x; epsilon = [e_x ;e_y]; %epsilon in all direction

% Set box centres and radii and depth 
d = 2; 
% cc = [0 0]; rr = [2 1]; 
% cc = [0.5 0.1]; rr = [1.2 0.5];
cc = [1 0.5]; rr = [3 2]; 
depth1 = 15 ; %depth of deterministic map
depth = 21 ; %depth need to be sufficiently deep to have radius of boxes less than epsilon
%depth = ceil(2*(log2(max(rr)/min(epsilon)))+1); % determined by epsilon and radius
%depth = ceil(2*(log2(max(rr)*2/max(epsilon)))+1);

% Parameter values of the map and initialise the data structures
a_start = 0.18; a_end = 0.18;a_gap = 0.01; % parameter starting and ending parameter to run in the loop
pl = zeros(length(a_start:a_gap:a_end),1);
haus = zeros(length(a_start:a_gap:a_end),1);
min_dist = zeros(length(a_start:a_gap:a_end),1);
b_para = 0.3; 
a_para = 0.60;

% For writing figures into movies or animations
writevideo = 0;
if writevideo
    vidfile1 = VideoWriter('henon a050-057','MPEG-4');
    vidfile1.FrameRate = 2;
    open(vidfile1);
end

% Start the computation looping over desired parameter range
a_all = a_start:a_gap:a_end;   
for l = 1:length(a_all) %for looping around parameters
    clear M_set
    a_para = a_all(l);
    clear M_set M tree tree2 tree3;
    t = Tree(cc,rr);    % the box collection
    f = @(x)map(x,'mapType',mapType,'param',[a_para,b_para]); %single value transformation
    sd = 8; %flag to cover all the boxes
    for i=1:depth
        t.set_flags('all', sd); 
        t.subdivide;
    end
    
    %sample grid for a ball (Xe)
    dist1 = min(rr)/2^(depth/2-1); num = ceil(2*epsilon(1)*pi/dist1);
    if only_boundary
        Xe = n_ball(d,round(num/2)); %only use boundary
    else
        b=t.boxes(depth); dist1=min(b(d+1:2*d,1));
        num = ceil(2*epsilon(1)*pi/dist1);
        if num<10, num=10; end
        num1 = ceil(1/(2*pi/num)+1); tot_num = 0;
        for j = 1:num1
            tot_num = tot_num + ceil(j*num/num1);
        end
        Xe = zeros(tot_num+1,d); Xe(1,:) = [0 0]; i = 1;
        for j = 1:num1
            for k = 1:ceil(j*num/num1)
                x1 = j/num1*cos(2*k*pi/ceil(j*num/num1)); %points on epsilon ball
                x2 = j/num1*sin(2*k*pi/ceil(j*num/num1));
                Xe(i+1,:) = [x1 x2]; i = i+1;
            end
        end
    end
    

    %approximate transpose transition matrix (transpose!!)
    disp('calculating deterministic transition matrix and initial point')
    tic;
    M = tpmatrix(t, f, X, depth1, 1); %single valued

    %approximate invariant set using eigenvector w.r.t. eigenvalue 1 (stationary measure)
    [Vl, ~] = eigs(M,1,'lr');

    %choose one of the most invariant boxes to be the start point
    diver = 1;
    t1 = t.boxes(depth1);
    ccount = 1;
    while diver == 1
        M1 = find(abs(Vl)  == max(abs(Vl)));
        % M1 = find(abs(Vl)  > 1e-17);
        % x = t1(1:2,M1(200));
        x = t1(1:2,M1);
        x0 = x; x = x';
        for trajj = 1:100, x = f(x);end
        x = x'; diver = 0;
        % if any(x>cc+rr) | any(x<cc-rr), diver=1;Vl(M1)=0;end
        ccount = ccount + 1;
        if ccount > 1e2, disp('infinite loop at finding invariant box');break; end
    end
    toc;
    %{
    %for testing the single valued dynamics recurrent points/sets
    v = zeros(length(M),1); v(find(abs(Vl)>eps)) = 1; b = t.boxes(depth); 
    hh = b(1:d,find(abs(Vl)>eps));
    tree = Tree(cc,rr); tree.insert(hh,depth);
    all = [1:size(M,1)]'; I = all; I_comp = I;
    while ~isempty(I_comp)
        M_temp = M(I,I);
        [~,I1] = find(M_temp);
        I1 = unique(I1); I1 = I(I1);
        I_comp = setdiff(I,I1);
        I = I1;
    end
    %{
        all = [1:size(M,1)]'; I = all; 
        [~,J1] = find(M);
    J1 = unique(J1); I = I(J1);
    Jcomp = setdiff(all,I);
    while ~isempty(Jcomp)
        M_temp = M(Jcomp,I);
        [~,Jcomp] = find(M_temp);
        Jcomp = unique(Jcomp);
        Jcomp = I(Jcomp);
        I = setdiff(I,Jcomp);
    end
    %}
    b = t.boxes(depth); dual = b(1:d,I); tree1 = Tree(cc,rr);  
    tree1.insert(dual,depth); %attraction zone
    clf; boxplot2(tree1);
    hold onmap
    boxplot2(tree,'color','y');
    hold off
    title(['a = ' num2str(a_para)])
    axis tight
    %drawnow;
    end
    %}

    %Find the invariant set starting from the starting point
    X1 = [0 0]; %to find which boxes hit by one box
    disp('calculating bounded noise transition')
    if stored_work == 1 % use stored matrix from file
        if fine == 1
            M_set = load(['fineM',num2str(a_para),'.mat']).M_set;
        else
            M_set = load(['M',num2str(a_para),'.mat']).M_set;
        end
    else
        M_set = tpmatrix_set_valued(t, f, epsilon, X, Xe, depth, 1, multiplicative); %set valued 
    end
    %depth = 15; 
    tree = Tree(cc,rr);
    toc;

    %subdivide the boxes in one axis for depth/2 times
    b = zeros(depth,1); b(ceil(depth/2)+1:depth) = 1;
    tree.sd = b; %making boxes numbered such that easier to track
    tree.insert(x, depth); %dimension of problem

    %start to find invariant set
    none = 0; ins = 2; expd = 4; 
    n0 = 0; n1 = tree.count(depth); count = 0;
    eP = [];
    while n1 > n0                                   % while boxes are being added
          ePold = eP;
          tree.change_flags('all', ins, expd);          % flag inserted boxes for expansion
          b = tree.boxes(depth);                           
          flags = b(2*d+1,:); 
          I = find(bitand(flags,expd));                 % find boxes to be expanded  
          c = b(1:d,I); r = b(d+1:2*d,1);               % get center and radii
          t.set_flags(c,ins,depth); 
          b = t.boxes(depth); I = find(b(2*d+1,:));
          [i,~] = find(M_set(:,I)); i = unique(i); %find which box it hit (can do it without the M_set too)
          eP = b(1:d,i);
          t.unset_flags('all',ins);
          tree.insert(eP, depth, ins, none);         % insert boxes which are hit by f(P)
          tree.unset_flags('all', expd);                % unflag recently expanded boxes
          n0 = n1; n1 = tree.count(depth); 
          if viewon
              %disp(sprintf('%d boxes', n1));
              clf; boxplot2(tree,'alpha',0.9); axis([cc(1)-rr(1) cc(1)+rr(1) cc(2)-rr(2) cc(2)+rr(2)]);
              drawnow; 
          end
          count = count + 1;
    end

    disp(['The number of steps for invariant set: ' num2str(count)])
    pl(l) = count; % for tracking count
    b = tree.boxes(depth); inv_set = b(1:d,:); % invariant set points

    % check if the mapping image is not contained in found invariant set
    misinf = 0;
    if any(tree.search(f(inv_set')')==-1), disp('some points in MIS mapped outside'); misinf = 1;end
    if misinf == 0
        % Find dual domain of attraction
        if finddual == 1
            tree3 = Tree(cc,rr);
            tree3.insert(inv_set, depth);
            none = 0; ins = 2; expd = 4; 
            n0 = 0; n1 = tree3.count(depth);
            while n1 > n0                                   % while boxes are being added
              tree3.change_flags('all', ins, expd);          % flag inserted boxes for expansion
              b = tree3.boxes(depth);                           
              flags = b(2*d+1,:); 
              I = find(bitand(flags,expd));                 % find boxes to be expanded  
              c = b(1:d,I); r = b(d+1:2*d,1);               % get center and radii
              t.set_flags(c,ins,depth); b = t.boxes(depth); I = find(b(2*d+1,:));
              [~,i] = find(M_set(I,:)); i = unique(i); % find boxes that will hit
              eP = b(1:d,i);
              t.unset_flags('all',ins);
              tree3.insert(eP, depth, ins, none);         % insert boxes which are hit by f(P)
              tree3.unset_flags('all', expd);                % unflag recently expanded boxes
              n0 = n1; n1 = tree3.count(depth); 
              if viewon
              %disp(sprintf('%d boxes', n1));
              clf; boxplot2(tree3,'alpha',0.9); axis([cc(1)-rr(1) cc(1)+rr(1) cc(2)-rr(2) cc(2)+rr(2)]);
              %light('position',[-1,2,2]); view(20,30); 
              drawnow;
              end
              disp(['dual number of box added',num2str(n0)])
            end
        
    
            %subdivide the boxes in one axis for depth/2 times
            if minimal_invariant
                b = zeros(depth,1); b(ceil(depth/2)+1:depth) = 1;
                %depth = 15; % change the depth of the approximation of dual repeller
                tree1 = Tree(cc,rr); tree1.sd = b; 
        
                %saddle point for henon map
                x1=(b_para-1+sqrt(4*a_para+(b_para-1)^2))/(2*a_para); 
                y1=(b_para^2-b_para+b_para*sqrt(4*a_para+(b_para-1)^2))/(2*a_para);
                x1 = [x1;y1]; %unstable point 
                if (mapType == 'Hen' | mapType == 'Los' ) & tree.search(x1,depth) == -1
                    % find dual repeller/domain of attraction from minimal invariant set from
                    % the saddle point of henon map
                    tree1.insert(x1, depth);
                    %tic
                    none = 0; hit = 1; ins = 2; expd = 4; 
                    n0 = 0; n1 = tree1.count(depth); 
                    while n1 > n0                                   % while boxes are being added
                      tree1.change_flags('all', ins, expd);          % flag inserted boxes for expansion
                      b = tree1.boxes(depth);                           
                      flags = b(2*d+1,:); 
                      I = find(bitand(flags,expd));                 % find boxes to be expanded  
                      c = b(1:d,I); r = b(d+1:2*d,1);               % get center and radii
                      t.set_flags(c,ins,depth); b = t.boxes(depth); I = find(b(2*d+1,:));
                      [~,i] = find(M_set(I,:)); i = unique(i);
                      eP = b(1:d,i);
                      t.unset_flags('all',ins);
                      tree1.insert(eP, depth, ins, none);         % insert boxes which are hit by f(P)
                      tree1.unset_flags('all', expd);                % unflag recently expanded boxes
                      n0 = n1; n1 = tree1.count(depth); 
                      if viewon
                          %disp(sprintf('%d boxes', n1));
                          clf; boxplot2(tree1,'alpha',0.9); axis([cc(1)-rr(1) cc(1)+rr(1) cc(2)-rr(2) cc(2)+rr(2)]);
                      drawnow;
                      end
                      disp(['minimal number of box added ',num2str(n0)])
                    end
            
                    %domain of attraction
                    tree2 = Tree(cc,rr);
                    sd = 8; %flag to cover all the boxes
                    for i=1:depth
                        tree2.set_flags('all', sd); 
                        tree2.subdivide;
                    end
                    b = tree3.boxes(depth); attract = b(1:d,:);
                    tree2.set_flags(attract,2,depth);
                    b = tree2.boxes(depth); repel = b(1:d,find(b(2*d+1,:)==0));
                    tree1.insert(repel,depth);
            
                    b = tree1.boxes(depth); wanted = b(1:d,:);
                    tree2 = Tree(cc,rr);
                    sd = 8; %flag to cover all the boxes
                    for i=1:depth
                        tree2.set_flags('all', sd); 
                        tree2.subdivide;
                    end
                    tree2.set_flags(wanted,2,depth);
                    b = tree2.boxes(depth); wanted = b(1:d,find(b(2*d+1,:)==0));
                    tree2 = Tree(cc,rr); tree2.insert(wanted,depth);
                    %The above calculate either union of domain of attraction with dual
                    %repeller or the larger invariant set, because either saddle point is in
                    %invaraint set or not(assuming it is in the domain)
                    if tree2.count(depth) == 0
                        tree2 = tree3;
                    end
                else
                    tree2 = tree3;
                end
            else
                tree2 = tree3;
            end
        else
            tree2 = t;
        end
        
        if viewon
            clf; boxplot2(tree,'edgecolor','none'); hold on
            boxplot2(tree2,'color','k','edgecolor','none');
            boxplot2(tree,'edgecolor','none'); hold off; axis tight
            if plotdeter
                hold on; scatter(x(1),x(2));hold off;
            end
            title(['Minimal invariant set for Henon map with a = ' num2str(a_para) ', b = ' num2str(b_para)])
            legend('domain of attraction','minimal invariant set','Location','southwest')
            xlabel('x');ylabel('y');
            axis equal
            axis([cc(1)-rr(1) cc(1)+rr(1) cc(2)-rr(2) cc(2)+rr(2)]);
            % axis([-0.7+1.65*l/10 1.7-0.25...
            %     *l/10 -0.4+0.3*l/10 0.6-0.5*l/10])
            axis([-2+1.3-1.3*l/10 2-0.3+0.3*l/10 -1+0.6-0.6*l/10 1-0.4+0.4*l/10])
            drawnow;
        end
    end
    
    
    a = a_para; b = b_para; epsilon1 =epsilon(1);
    %deterministic
    det_record = henon_det_n([0,0],a,b,1000);vec_det = henon_det_n(det_record,a,b,1);
    iter = 0;
    while all(vecnorm((vec_det-det_record)')>1e-3)&&iter<1000
        det_record = [det_record;vec_det];
        vec_det = henon_det_n(vec_det,a,b,1);
        iter = iter +1;
    end
    hold on;
    if plotdeter
        scatter(det_record(:,1),det_record(:,2),50,'filled','b');drawnow;
    end
    jj = 1;
    traj_all_all = {};

    % boundary of minimal invariant set using boundary map
    if extended == 1
        dist = 2; periodicity = 2;
        n = 10; x11 = linspace(-dist,dist,n)';
        [XX,YY,ZZ] = meshgrid(x11,x11,x11/dist*2*pi);
        clear x;
        vec_all = []; 
        % find the periodic points of extended mapping
        for i =1:n^3
            traj_all = [];
            options = optimset('Display','off');
            vec = fsolve(@(x) henon_eq_n(x,a,b,epsilon1,periodicity) ,[XX(i),YY(i),ZZ(i)],options);
            if isempty(vec_all) % compared with other found points
                vecc = 1;
            else
                vecc = all(vecnorm(vec'-vec_all)>1e-2);
            end
            if vec(1)<-1
                continue
            end
            if norm(henon_eq_n(vec,a,b,epsilon1,2))<1e-8 && vecc 
                vec_all = [vec_all,vec'];
                syms x y z
                jac = jacobian(henon_3d([x,y,z],a,b,epsilon1,2),[x,y,z]);
                [evec,eval] = eig(double(subs(jac,{x,y,z},{vec(1),vec(2),vec(3)})));
                [row,~] = find(abs(eval)>1);
                
                if length(row)==1 %one unstable direction
                    unstable = 1;
                    target = evec(:,row);
                    traj = vec'+kron(1e-4*linspace(-1,1,2000),target);
                    for j = 1:400
                        traj_all = [traj_all,traj];
                        traj = henon_3d(traj,a,b,epsilon1,2);
                        if any(traj>1e1)
                            break
                        end
                    end
                    if ~isempty(traj_all)
                        kk = find(abs(vecnorm(traj_all))<10);
                        traj_all = traj_all(:,unique(kk));
                        traj_all_all{floor(jj)} = traj_all;
                        if viewon
                            plot(traj_all(1,:),traj_all(2,:),'.','color','#1f77b4') % plot the unstable manifold
                        end
                        jj=jj+1;
                    end
                    scatter(vec(1),vec(2),40,'filled','c');drawnow;
                elseif length(row)==2 %two unstable directions
                    unstable = 0;
                    [row,~] = find(abs(eval)<1&abs(eval)>0);
                    target = evec(:,row);
                    traj = vec'+kron(1e-4*linspace(-1,1,30000),target);
                    for j = 1:40
                        traj_all = [traj_all,traj];
                        traj = henon_inverse_3d(traj,a,b,epsilon1,2);
                        if any(traj>1e1)
                            break
                        end
                    end
                    if ~isempty(traj_all)
                        kk = find(abs(vecnorm(traj_all))<10);
                        traj_all = traj_all(:,unique(kk));
                        traj_all_all{floor(jj)} = traj_all;
                        if viewon
                            plot(traj_all(1,:),traj_all(2,:),'.','color','#FFA500') % plot the stable manifold
                        end
                        jj=jj+1;
                    end
                    scatter(vec(1),vec(2),40,'filled','MarkerFaceColor','#D7A1F9'); hold on;drawnow;
                else %all stable 
                    scatter(vec(1),vec(2),40,'filled','g');drawnow;
                end
            end
        end
    end

    % xlim([-1,2])
    % ylim([-1,1])
    if viewon
        if mapType == 'Hen'
            title(['Minimal invariant set for Henon map with a = ' num2str(a_para) ', b = ' num2str(b_para)])
        end
        drawnow;
    end
    if viewon
        if finddual == 1
            if plotdeter
                legend('minimal invariant set','domain of attraction','deterministic equilibrium points','Location','SouthWest')
            else
                legend('minimal invariant set','domain of attraction','Location','SouthWest')
            end
        else
            if plotdeter
                legend('minimal invariant set','deterministic equilibrium points','Location','SouthWest')
            else
                legend('minimal invariant set','Location','SouthWest')
            end
        end
    end

    %make video for both tight axis or constant axis
    if writevideo, F(count) = getframe(gcf);writeVideo(vidfile1,F(count)); end
    toc;fprintf('\n')

    % save image or the tree structure for the MIS and dual repeller
    if saveimage
        exportgraphics(gcf,strcat('zoom_henon',num2str(round(a,4)),'.png'),'Resolution',300)
    end
    if savework
        save(['M',num2str(a_para),'.mat'],'M','M_set','-v7.3')
        % plotty = tree.boxes(depth); plotty1 = tree2.boxes(depth);
        % save(['plot',num2str(a_para),'.mat'],'plotty','plotty1','-v7.3')
    end
    if exist('tree2','var')
        t_all{l,1} = {tree,tree2};
    else
        t_all{l,1} = {tree};
    end
    disp(a_para)
    figure
    boxplot2(tree,'alpha',0.9,'color','k'); 
    axis([cc(1)-rr(1) cc(1)+rr(1) cc(2)-rr(2) cc(2)+rr(2)]);
    axis('equal')
    xlabel('x')
    ylabel('y')
    title(['$a = $',num2str(a_para),', $b = $',num2str(b_para),', $\varepsilon = $',num2str(epsilon(1))],'interpreter','latex')
end
if writevideo, close(vidfile1);end


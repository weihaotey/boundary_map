%% global unstable manifold of an equilibrium for extended system
a = 0.1; b = 0.3; epsilon = 0.0625; period = 2;                       % parameters
f = @(x) henon_inverse_3d(x,a,b,epsilon,period);        % extended mapping

% Approximate periodic points
ff = @(x) henon_n([x(1),x(2),x(3),x(4)]',a,b,epsilon,period)-[x(1),x(2),x(3),x(4)]';
n = 10; x = linspace(-1,1,n);
[XX,YY,ZZ] = meshgrid(x*2,x*2,x*pi);
X = [ XX(:) YY(:) cos(ZZ(:)) sin(ZZ(:))];   
target = [];
for i = 1:size(X,1)
    targetx = fsolve(ff,X(i,:)');
    if sqrt(sum(ff(targetx).^2))>1e-2 || targetx(1) < -1 %remove nonroot and dual repeller points
        continue
    end
    if ~isempty(target)
        K = targetx - target;
        if any(sqrt(sum(K.^2,1))<1e-2)
            continue
        end
    end
    target = [target,targetx];
end
sample1 = [];
%%
b = 0.3; % redefine parameter b
% import data
b_vec = target';
%b_vec = readNPY(strcat('tomatlab_vec_',num2str(a),'.npy'));
%b_ = readNPY(strcat('tomatlab_',num2str(a),'.npy'));
syms x y n1 n2
hh = henon_n([x,y,n1,n2],a,b,epsilon,2);
jac1 = jacobian(hh, [x; y; n1; n2]);
disp(b_vec)
vec1 = b_vec(1,:);
disp(vec1)
vec2 = b_vec(2,:);
disp(vec2)
vec3 = b_vec(3,:);
disp(vec3)
vec4 = b_vec(4,:);
disp(vec4)
targetvec = vec4;
disp(targetvec)
[eig1,eig2]=eig(double(subs(jac1,[x,y,n1,n2],targetvec)));
eigenvec1 = eig1(:,1);
eigenvec1_ = [eigenvec1(1),eigenvec1(2),atan2(eigenvec1(4),eigenvec1(3))];
eigenvec2 = eig1(:,2);
eigenvec2_ = [eigenvec2(1),eigenvec2(2),atan2(eigenvec2(4),eigenvec2(3))];
eigenvec3 = eig1(:,3);
eigenvec3_ = [eigenvec3(1),eigenvec3(2),atan2(eigenvec3(4),eigenvec3(3))];
eigenvec4 = eig1(:,4);
eigenvec4_ = [eigenvec4(1),eigenvec4(2),atan2(eigenvec4(4),eigenvec4(3))];
disp(eig2)

%% unstable manifold
%sample1 = [];
sample_dist = 1e-4; sample_no = 5000;
for eigentarget = [eigenvec1,-eigenvec1]
    vec0 = targetvec'+ sample_dist*eigentarget;
    %vec = wedge_func(vec0,a,b,epsilon,2);
    vec = henon_n(vec0,a,b,epsilon,2); % one double henon step
    dist0 = vec-vec0;
    sample = vec0+linspace(1/sample_no,1,sample_no).*dist0; % sample points from two consecutive points

    for i =1:200
        %plot3(sample(1,:),sample(2,:),atan2(sample(4,:),sample(3,:)),'.')
        %hold on
        sample1 = [sample1,sample];
        %sample = wedge_func(sample,a,b,epsilon,period);
        
        sample = henon_n(sample,a,b,epsilon,2);
        %plot3(sample(1,:),sample(2,:),atan2(sample(4,:),sample(3,:)),'.')
        %hold on
        %plot3(b_vec(:,1),b_vec(:,2),atan2(b_vec(:,4),b_vec(:,3)),'.','markersize',20)
        %hold off
        %grid on
        %xlim([targetvec(1)-.1,targetvec(1)+.1]);
        %ylim([targetvec(2)-.1,targetvec(2)+.1])
        %zlim([atan2(targetvec(4),targetvec(3))-0.2,atan2(targetvec(4),targetvec(3))+.2]);
        %pause
        %if any(max(abs(sample))>10)
        %    disp(i)
        %    break
        %end
    end
end
figure(round(a*1000))
plot3(sample1(1,:),sample1(2,:),atan2(sample1(4,:),sample1(3,:)),'.','markersize',.1)
hold on; grid on;
plot3(b_vec(:,1),b_vec(:,2),atan2(b_vec(:,4),b_vec(:,3)),'.','markersize',20)
hold off
%%
sample2 = sample1;
sample2(1:2,1:end) = sample1(1:2,1:end) - epsilon*sample1(3:4,1:end);
hold on
plot3(sample2(1,:),sample2(2,:),atan2(sample2(4,:),sample2(3,:)),'.')
hold off
%% box painting
d = 3; depth = 25; 
n = 100; x = linspace(-1,1,n)';
h = eigenvec2_.*cos(x*pi)+eigenvec3_.*sin(x*pi);
h = kron(linspace(0.01,1,100)',h); 

n = 7; x = linspace(-1,1,n)';
[XX,YY,ZZ] = meshgrid(x,x,x);
X = [ XX(:) YY(:) ZZ(:) ];                      % sample points (uniform grid)
%X = [0 0 0];
%c = [0 0 0]; r = [2 2 4];
c1 = [targetvec(1),targetvec(2),atan2(targetvec(4),targetvec(3))]; r1 = [.1,.1,1];    % the box collection
c = c1; r = r1;
tree = Tree(c,r);

% continuation algorithm
x0 = [vec1(1),vec1(2),atan2(vec1(4),vec1(3))];     % equilibrium  
x1 = [vec2(1),vec2(2),atan2(vec2(4),vec2(3))];     % equilibrium 
x2 = [vec3(1),vec3(2),atan2(vec3(4),vec3(3))];     % equilibrium  
x3 = [vec4(1),vec4(2),atan2(vec4(4),vec4(3))];     % equilibrium 
targetx = [targetvec(1),targetvec(2),atan2(targetvec(4),targetvec(3))];
tree.insert(targetx', depth);
tree.insert((targetx+h*1e-3)',depth);
%tree.insert(x1', depth);
%tree.insert((x1+h*1e-4)',depth);
%tree.insert(x2', depth);
%tree.insert(x3', depth);


none = 0; hit = 1; ins = 2; expd = 4; 
n0 = 0; n1 = tree.count(depth); tic
while n1 > n0                                   % while boxes are being added
  tree.change_flags('all', ins, expd);             % flag inserted boxes for expansion
  b = tree.boxes(depth);                           
  flags = b(2*d+1,:); 
  I = find(bitand(flags,expd));                 % find boxes to be expanded  
  c = b(1:d,I); r = b(d+1:2*d,1);               % get center and radii
  n = length(I); E = ones(n,1);
  P = kron(E,X)*diag(r) + kron(c',ones(size(X,1),1)); % sample points
  tree.insert(f(P'), depth, ins, none);            % insert boxes which are hit by f(P)
  tree.unset_flags('all', expd);                   % unflag recently expanded boxes
  n0 = n1; n1 = tree.count(depth); 
  fprintf('%d boxes\n', n1);
  %boxplot3(tree); view(20,30); axis square; axis tight; drawnow; pause(0)
  %grid on
end
toc

% plot
boxplot3(tree,'edgecolor','none','color','black','alpha',0.1); view(20,30); axis square; axis tight;
xlabel('x'); ylabel('y'); zlabel('z'); 
hold on; grid on;
plot3(b_vec(:,1),b_vec(:,2),atan2(b_vec(:,4),b_vec(:,3)),'.','markersize',20)
hold on
plot3(sample1(1,:),sample1(2,:),atan2(sample1(4,:),sample1(3,:)),'.')
xlim([c1(1)-r1(1),c1(1)+r1(1)])
ylim([c1(2)-r1(2),c1(2)+r1(2)])
zlim([c1(3)-r1(3),c1(3)+r1(3)])
%% cleanup
delete(tree);

%% look at local stable manifold

% weakest direction
sample_dist = .1;
sample_no = 100;
vec0 = targetvec'+ sample_dist*eigenvec1;
dist0 = targetvec'-sample_dist*eigenvec1-vec0;
k = vec0+linspace(1/sample_no,1,sample_no).*dist0; % sample points from two consecutive points

% invariant strong stable submanifold
n = 100; x = linspace(-1,1,n)';
h = eigenvec2'.*cos(x*pi)+eigenvec3'.*sin(x*pi);
h = kron(linspace(0.002,1,500)',h); 
kk = targetvec + h*sample_dist;
hold on
plot3(k(1,:),k(2,:),atan2(k(4,:),k(3,:)))
plot3(kk(:,1),kk(:,2),atan2(kk(:,4),kk(:,3)))

xlim([targetvec(1)-sample_dist,targetvec(1)+sample_dist]);
ylim([targetvec(2)-sample_dist,targetvec(2)+sample_dist])
zlim([atan2(targetvec(4),targetvec(3))-2*sample_dist,atan2(targetvec(4),targetvec(3))+2*sample_dist]);
title(strcat('a = ',num2str(a)))
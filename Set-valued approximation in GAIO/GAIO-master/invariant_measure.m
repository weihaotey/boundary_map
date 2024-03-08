%% Natural invariant measure of some mapping

%% find the boundary of MIS for 1-D map
figure
epsilon = .1;
f = @(x) 3*((1-exp(-x))./(1+exp(-x)));
xmin = -4; xmax = 4; % limit of test sample points
nnn = 400; % number of sample points
test = linspace(xmin,xmax,nnn+1); % initial guess for root finding
a = -0.5;a_end = 1;a_start = a; % parameters
n_n = 10000; % numbers of parameters to run
x_alll = {};
a_all = linspace(a_start,a_end,n_n);
tic;
for ii = 1:n_n
    a = a_all(ii);
    initmin = zeros(nnn,1); initmax = zeros(nnn,1); init = zeros(nnn,1); % initialise for recording data
    for i = 1:nnn+1
        initmin(i)=fzero(@(x)(f(x)-a-epsilon-x),test(i));
        initmax(i)=fzero(@(x)(f(x)-a+epsilon-x),test(i));
        init(i)=fzero(@(x)(f(x)-a-x),test(i));
    end
    x_alll{ii} = {unique(round(initmin,5)),unique(round(initmax,5)),unique(round(init,5))};
    if mod(ii,100)==0
        disp(ii)
    end
end
for ii =1:n_n
    plot(a_all(ii),x_alll{ii}{1},'r.'); hold on
    plot(a_all(ii),x_alll{ii}{3},'k.');
    plot(a_all(ii),x_alll{ii}{2},'b.'); 
end
hold off
toc;

%% invariant measure calculation

d = 1; onepoint = 1; filled = 1; plot_on = 1;
sd = 8; depth = 10; nn = n_n;
data_all = {};
epsilon0 = 0.1; epsilon_add = 0.001; epsilon = epsilon0;
if d == 1
    b = 3;      % the map
else
    mapType = 'Hen'; epsilon = epsilon0; a_add  = 0.;
end
a_para = 0.2; b_para = 3; a_add = 0.05;
cr_all = ones(2,nn); % record the box center and width
for j = 1:nn
    clear t
    nstart = 1; nincrease= 1;
    a_para = a_all(j);
    % a_para = a_start + (nstart+(j-2)*nincrease)*(a_end-a_start)/n_n;
    % epsilon = epsilon0+(j-1)*epsilon_add;
    % a_para = a_para + a_add;
    if d == 2
        f = @(x)map(x,'mapType',mapType,'param',[a_para,b_para]); %single value transformation
    else
        f = @(x) b*((1-exp(-x))./(1+exp(-x)))-a_para; % the function
    end

    % Define the tree
    if d == 1
        n = 10; X = linspace(-1,1,n)';X1=X; Xe = X;    % uniform grid of sample points 
%         c = [2.2]; r = [0.3]; t = Tree(c, r);  % the tree
        % if length(x_alll{nstart+(j-1)*nincrease}{1})==3
        %     rr = x_alll{nstart+(j-1)*nincrease}{2}(end)-x_alll{nstart+(j-1)*nincrease}{1}(end);
        %     cc = x_alll{nstart+(j-1)*nincrease}{2}(end)+x_alll{nstart+(j-1)*nincrease}{1}(end);
        % else
        %     rr = x_alll{nstart+(j-1)*nincrease}{2}(1)-x_alll{nstart+(j-1)*nincrease}{1}(end);
        %     cc = x_alll{nstart+(j-1)*nincrease}{2}(1)+x_alll{nstart+(j-1)*nincrease}{1}(end);
        % end
        if length(x_alll{j}{1})==3
            rr = x_alll{j}{2}(end)-x_alll{j}{1}(end);
            cc = x_alll{j}{2}(end)+x_alll{j}{1}(end);
        else
            rr = x_alll{j}{2}(1)-x_alll{j}{1}(end);
            cc = x_alll{j}{2}(1)+x_alll{j}{1}(end);
        end
        c = cc/2; r = abs(rr/2); t = Tree(c, r);  % the tree

    elseif d == 2
%         c = [0.5 0.15]; r = [1.4 0.4]; t = Tree(c, r);
        c = [0.2 0.15]; r = [1.4 0.6]; t = Tree(c, r);
%         c = [0.2 0.]; r = [1.2 0.5]; t = Tree(c, r);
    else
        error('dimension more than 2')
    end
    
    % Construct full subdivison
    for i=1:depth
        t.set_flags('all', sd);
        t.subdivide;
    end
    
    % Compute invariant vector
    if d == 2
        if onepoint == 1
            X1 = [0 0]; % to find which boxes hit by one box
        else
            n = 4; x = linspace(-1,1,n)'; [XX,YY] = meshgrid(x); X1 = [ XX(:) YY(:)];
        end
    end

    % Sample grid for a ball (Xe)
    if d == 2
        if filled == 1 % use filled epsilon ball
            bb=t.boxes(depth); dist1=min(bb(d+1:2*d,1))*2;
            num = ceil(2*epsilon(1)*pi/dist1); %number of points on boundary of ball
            if num<10, num=10; end % minimimum of points is 10
            num1 = ceil(1/(2*pi/num)+1); tot_num = 0;
            for jj = 1:num1
                tot_num = tot_num + ceil(jj*num/num1);
            end
            Xe = zeros(tot_num+1,d); Xe(1,:) = [0 0]; i = 1;
            for jjj = 1:num1
              for k = 1:ceil(jjj*num/num1)
                x1 = jjj/num1*cos(2*k*pi/ceil(jjj*num/num1)); %points on epsilon ball
                x2 = jjj/num1*sin(2*k*pi/ceil(jjj*num/num1));
                Xe(i+1,:) = [x1 x2]; i = i+1;
              end
            end
            clear bb
        else % only the boundary
            dist1 = min(r)/2^(depth/2-1); num = ceil(2*epsilon(1)*pi/dist1);
            Xe = n_ball(d,round(num/2)); %only use boundary
        end
    elseif d==1
        bb=t.boxes(depth); dist1=min(abs(bb(d+1:2*d,1)))*2; num = ceil(2*epsilon/dist1);
        if filled == 1 % use filled epsilon ball
            Xe = linspace(-1,1,num)';
        else
            Xe = linspace(-1,1,10)';
        end
    end

    disp('calculating matrix')

    % Compute transition matrix
    P = tpmatrix_set_valued(t, f, epsilon, X1, Xe, depth, 1);
    [v,lambda] = eigs(P,1);
    
    % Plot invariant density
    if plot_on == 1
        if d == 1
            n = t.count(depth); x = linspace(c(1)-r(1),c(1)+r(1),n);
            bbb = t.boxes(depth); rrr = bbb(2,1);
            h = abs(v)/(2*rrr.*norm(v,1)); % so that area under the histogram is one
            bar(x,h,1); 
            axis([c(1)-r(1) c(1)+r(1) 0 5]);
%             axis([ 0 5])  
            xlabel('x'); ylabel('density');
            title(strcat('a = ',num2str(a_para),', b = ',num2str(b_para),', epsilon = ',num2str(epsilon)))
            pause(.01)
            data_all{j} = [x;h';ones(1,n)*a_para;ones(1,n)*epsilon];
            cr_all(:,j) = [c,r]; % record the box center and width
        else
            figure(1)
            boxplot3(t,'depth',depth,'density', v ,'alpha',0.7); axis equal
            title(strcat('a = ',num2str(a_para),', b = ',num2str(b_para),', epsilon = ',num2str(epsilon)))
            figure(2)
            boxplot3(t,'depth',depth,'density', log(abs(v)) ,'alpha',0.7); axis equal
            title(strcat('a = ',num2str(a_para),', b = ',num2str(b_para),', epsilon = ',num2str(epsilon)))
            data_all{j} = [v';ones(1,length(v))*a_para,ones(1,length(v))*epsilon];
            cr_all(:,j) = [c,r]; % record the box center and width
        end
    end
    hold off
    fprintf('Step %d \n',j) 
    pause(.01)
end

%%
plotx = ones(1,nn);
ploty = ones(1,nn);
writevideo = 0;
if writevideo
    vidfile1 = VideoWriter('invariant_measure_a07-03_epsilon','MPEG-4');
    vidfile1.FrameRate = 2;
    open(vidfile1);
end
figure('units','pixels','position',[0 0 1440 1080])
for i = 1:nn
    data = data_all{i};
    if d == 1
        n = length(data);
        bar(data_all{i}(1,:),data_all{i}(2,:),1)
        c = cr_all(1,i); r = cr_all(2,i);
        axis([c(1)-r(1) c(1)+r(1) 0 10])
        title(strcat('epsilon=',num2str(epsilon0+epsilon_add*i)))
        % axis([-4 4 0 10])
        pause(.01)
        for j = 1:n
            if sum(data(2,1:j+1)) > n*0.05        
    %         if sum(data(2,n/2-(j-1):n/2+j)) > n*0.05
    %             xx = data(1,n/2-(j-1):n/2+j);
    %             hh = data(2,n/2-(j-1):n/2+j);
    %             N = sum(data(2,n/2-(j-1):n/2+j));
                xx = data(1,1:j+1);
                hh = data(2,1:j+1);
                N = sum(data(2,1:j+1));
                plotx(i) = epsilon0+i*epsilon_add;
                ploty(i) = sum(xx.^2.*hh)/N-(sum(xx.*hh)/N)^2;
                disp(j)
                break
            end
        end
    else
        a_para = 1.05 + i*a_add;
        v = data(1,:);
%         boxplot3(t,'depth',depth,'density', abs(v) ,'alpha',0.7); axis equal
        boxplot3(t,'depth',depth,'density', log(abs(v)) ,'alpha',0.7); axis equal
        title(strcat('a = ',num2str(a_para),', b = ',num2str(b_para),', epsilon = ',num2str(epsilon)))
        
        if writevideo, F(i) = getframe(gcf);writeVideo(vidfile1,F(i)); end
        pause
    end
    
end
if writevideo, close(vidfile1);end



%% find the optimal tail distribution from the stationary distribution
% using least square model
% close all
tic;
writevideo = 0;
if writevideo
    vidfile1 = VideoWriter('brute_tilfit','MPEG-4');
    vidfile1.FrameRate = 5;
    open(vidfile1);
end
brute = 1; linearmethod = 0; extrapolate = 0; logonce = 0;% which method to use
findtol = 1e-6; % the loglog end point line estimate
extol = 1e-10; % exclude data where h value is small
fval_tol = 1e-4; % tolerance of function evaluation
disp_on = 0; % whether to display result
testboun = 0; % whether to test some value away from boundary
% boundata = 0.01; % test some data away from boundary
data_portion = 0.0001; % portion of data point to fit
lambda_all = ones(1,n_n); % record lambda estimate
real_lambda_all = ones(1,n_n); % record real lambda
param = ones(1,n_n); % record parameter
if linearmethod || extrapolate || logonce
    Aloglog = ones(2,n_n); % record resulting param estimate
else
    Aloglog = ones(3,n_n);
end
for j = 1:1
    disp(j)
    fval_tol1 = fval_tol;
    a_para = data_all{j}(3,1); % retrieve a parameter value
    param(j) = a_para; % record parameter
    epsilon = data_all{j}(4,1); % retrieve epsilon parameter value  
    c = cr_all(1,j); r = cr_all(2,j); % retrieve center and radius of the box
    xx = data_all{j}(1,:); % retrieve the x value
    h = data_all{j}(2,:); % retrieve bins
    width = (xx(2)-xx(1)); % width of histogram Bin
    h = h/norm(h,1)/width; % so that area is 1
    bounmin = c-r; % the lower boundary
    bounmax = c+r; % upper boundary
    syms yy; fprime = matlabFunction(diff(f(yy))); % derivative of the map
    % extol = max(h)*data_portion/000; % exclude data where h value is small

    % only take data near left boundary
    if testboun
        w = find(x<(c-r+boundata)); % only test up to epsilon away from boundary
    else
        w = find(cumsum(h)*width<data_portion & abs(h)>extol); % only test some data portion
    end

    % % only take data near right boundary
    % if testboun
    %     w = find(x>(c+r-boundata)); % only test up to epsilon away from boundary
    % else
    %     w = find(cumsum(h)*width>1-data_portion); % only test some data portion
    % end
    w = w(2:end); % get rid of end point 
    % w = w(1:end-1); % get rid of end point 
    x1 = xx(w); % very tail of the distribtuion to test
    x = xx; boun = bounmin;
    
    real_lambda = fprime(boun); % real lambda
    real_lambda_all(1,j) = real_lambda;
    
    opts=optimset('Display','off');
    figure(1)
    % brute or log-log linear regression
    % hyper case
    if brute 
        % hyper case
        opts=optimset('Display','off');
        figure(2)
        plot(xx,h,"Color","blue"); hold on; % plot stationary distribution
        repeat = 0; % to know whether need to recalculate
        fval1 = 1; % initialise large fval
        count = 1; % track iterations in while loop
        while fval1 > fval_tol1 % repeat if tolerance not met
            A1 = optimvar('A1',2);
            fun = log(abs(x1-bounmin));
            fun1 = A1(1)*exp(A1(2)*fun.*fun); % hyper tail
            obj1 = sum((fun1-h(w)).^2); % objective function
            lsqproblem1 = optimproblem("Objective",obj1); % optimization problem
            if j==1
                x0.A1 = [1,-1/10]; % initial point of a
            else
                x0.A1 = sol1.A1';
            end
            if repeat  == 1 % repeat the calculation with random initial
                x0.A1 = [unifrnd(0,1),unifrnd(-10,0)];
            end
            % show(lsqproblem)
            [sol1,fval1] = solve(lsqproblem1,x0,'Options',opts);
            count = count + 1; repeat = 1;
            if mod(count,10) == 0, disp('doubled fval tol'); fval_tol1 = fval_tol1*2; end
        end
        fitted_hyper = evaluate(fun1,sol1); % evaluate function at solution
        plot(x1,fitted_hyper,"Color",'red') % plot hyper best fit
        title(strcat('a = ',num2str(a_para)))
        xlim([c-r,max(x1)])
        lambda = round(exp(1/2/sol1.A1(2)),3); % lambda estimate
        hold off
        disp(fval1)
    elseif linearmethod % approximate linear equation
        A = optimvar('A',2);
        fun = A(1)*log(-log(abs(xx(w)-boun))) + A(2); % linear equation
        obj1 = sum((fun-log(-log(h(w)))).^2); % objective function
        lsqproblem0 = optimproblem("Objective",obj1); % optimization problem
        x0.A = [2,-0.8]; % initial point of a
        % show(lsqproblem)
        [sol0,fval0] = solve(lsqproblem0,x0,'Options',opts);
        Aloglog(:,j) = sol0.A; % record fit parameter
        % plot the result
        loglogx = log(-log(findtol)); % the value of intersection chosen for log-log|x-x0|
        p = log(-log(abs(x(w)-boun)));
        plot(p,log(-log(h(w))),'.-'); hold on;
        plot(p,evaluate(fun,sol0),'.-');
        plot(loglogx,sol0.A(1)*loglogx+sol0.A(2),'.');
        plot(p,2*p+log(-1/(2*log(real_lambda))),'.-');
        plot(loglogx,2*loglogx+log(-1/(2*log(real_lambda))),'.-');
        hold off
        legend('data','best fit','actual line','Location','best')
        aa = (sol0.A(1)-2)*loglogx + sol0.A(2); % calculate the y-intercept for y=2x+c
    elseif extrapolate % technically this method only use last two data to extrapolate
        % findind = find(abs(h) > findtol);
        loglogx = log(-log(findtol)); % the value of intersection chosen for log-log|x-x0|
        p = log(-log(abs(x(w)-boun))); q = log(-log(h(w))); % data points for log-log scale
        plotx = [loglogx p];
        pq = interp1(p,q,plotx,'linear','extrap');
        plot(p,log(-log(h(w)))); hold on;
        plot(plotx,pq,'.-');
        plot(plotx,2*plotx+log(-1/(2*log(real_lambda))),'.-');
        hold off
        legend('data','best fit','actual line','Location','best')
        aa = pq(1)-2*loglogx;
    elseif logonce
        repeat = 0; % to know whether need to recalculate
        fval0 = 1; % initialise large fval
        count = 1; % track iterations in while loop
        while fval0 > fval_tol1 % repeat if tolerance not met
            A = optimvar('A',2);
            pp = log(abs(x1-boun));
            fun = A(2)*pp.*pp + A(1);
            obj1 = sum((fun-log(h(w))).^2); % objective function
            lsqproblem0 = optimproblem("Objective",obj1); % optimization problem
            if j == 1
                x0.A = [2,-0.8]; % initial point of a
            else
                x0.A = sol0.A';
            end
            if repeat == 1
                x0.A = [unifrnd(0,10),unifrnd(-10,0)];
            end
            % show(lsqproblem)
            [sol0,fval0] = solve(lsqproblem0,x0,'Options',opts);
            count = count + 1; repeat = 1;
            if mod(count,10) == 0, disp('doubled fval tol'); fval_tol1 = fval_tol1*2; end
        end
        Aloglog(:,j) = sol0.A; % record fit parameter
        % plot the result
        % findind = find(abs(h) > findtol);
        % loglogx = log(findtol); % the value of intersection chosen for log-log|x-x0|
        p = log(abs(x(w)-boun));
        plot(p,log(h(w)),'.-'); hold on;
        plot(p,evaluate(fun,sol0),'.-');
        plot(p,1/(2*log(real_lambda)).*p.^2,'.-');
        % scatter(loglogx,2*loglogx+sol0.A(2)+sol0.A(1)*exp(-sol0.A(3)*log(loglogx)));
        % scatter(loglogx,2*loglogx+log(-1/(2*log(real_lambda))));
        grid on
        hold off
        xlabel('log(|x-x_0|)')
        ylabel('log(h)')
        legend('data','best fit','actual line','Location','best')
        % aa = sol0.A(1)*exp(-sol0.A(3)*loglogx)+ sol0.A(2); % calculate the y-intercept for y=2x+c
        % aa = sol0.A(2);
        lambda = round(exp(1/2/sol0.A(2)),3); % lambda estimate
        disp(fval0)
    else % exponential term in the log log scale
        repeat = 0; % to know whether need to recalculate
        fval0 = 1; % initialise large fval
        count = 1; % track iterations in while loop
        while fval0 > fval_tol1 % repeat if tolerance not met
            A = optimvar('A',3);
            pp = log(-log(abs(x1-boun)));
            fun = 2*pp + A(2) + A(1)*exp(-A(3)*log(pp));
            obj1 = sum((fun-log(-log(h(w)))).^2); % objective function
            lsqproblem0 = optimproblem("Objective",obj1); % optimization problem
            if j == 1
                x0.A = [2,-0.8,1]; % initial point of a
            else
                x0.A = sol0.A';
            end
            if repeat == 1
                x0.A = [unifrnd(0,10),unifrnd(-10,0),unifrnd(0,100)];
            end
            % show(lsqproblem)
            [sol0,fval0] = solve(lsqproblem0,x0,'Options',opts);
            count = count + 1; repeat = 1;
            if mod(count,10) == 0, disp('doubled fval tol'); fval_tol1 = fval_tol1*2; end
        end
        Aloglog(:,j) = sol0.A; % record fit parameter
        % plot the result
        % findind = find(abs(h) > findtol);
        loglogx = log(-log(findtol)); % the value of intersection chosen for log-log|x-x0|
        p = log(-log(abs(x(w)-boun)));
        plot(p,log(-log(h(w))),'.-'); hold on;
        plot(p,evaluate(fun,sol0),'.-');
        plot(p,2*p+log(-1/(2*log(real_lambda))),'.-');
        scatter(loglogx,2*loglogx+sol0.A(2)+sol0.A(1)*exp(-sol0.A(3)*log(loglogx)));
        scatter(loglogx,2*loglogx+log(-1/(2*log(real_lambda))));
        grid on
        hold off
        xlabel('log-log(|x-x_0|)')
        ylabel('log-log(h)')
        legend('data','best fit','actual line','Location','best')
        % aa = sol0.A(1)*exp(-sol0.A(3)*loglogx)+ sol0.A(2); % calculate the y-intercept for y=2x+c
        aa = sol0.A(2);
        disp(fval0)
    end
    title(strcat('a = ',num2str(a_para)))
    if writevideo, F(i) = getframe(gcf);writeVideo(vidfile1,F(i)); end
    if not(brute) && not(logonce)
        lambda = round(exp(-1/2*exp(-aa)),3);
        figure(2)
        plot(xx,h,"Color","blue"); hold on; % plot stationary distribution
        hold on
        plot(x1,h(w),'.')
        hold off
    end
    lambda_all(1,j) = lambda; % record lambda estimate
    % 
    % tic;
    
    % 
    % %nonhyper case
    % A2 = optimvar('A2',2);
    % fun2 = A2(1)*exp(A2(2)*log(abs(x1-bounmin))./abs(x1-bounmin)); % nonhyper tail
    % obj2 = sum((fun2-h(w)).^2); % objective function
    % lsqproblem2 = optimproblem("Objective",obj2); % optimization problem
    % x0.A2 = [1,0.0001]; % initial point of a
    % [sol2,fval2] = solve(lsqproblem2,x0);
    % 
    % fitted_nonhyper = evaluate(fun2,sol2); % evaluate function at solution
    % plot(x1,fitted_nonhyper,"Color",'black') % plot nonhyper best fit
    % text(x1(end),fitted_nonhyper(end),num2str(lambda));
    % if disp_on
    %     disp(a_para)
    %     disp(sol1.A1)
    %     disp(sol2.A2)
    % end
    % disp(num2str((fval1-fval2)/fval2))
    % legend('stationary dist','hyper fit','nonhyper fit','location','northwest')
    % title(strcat('a = ',num2str(a_para)))
    % xlim([c-r,c+r])
    % ylim([0 10])
    hold off
    pause
    toc;
end
if writevideo, close(vidfile1);end
% plot lambda estimates and the real lambda
figure
plot(param, lambda_all,'.-'); hold on
plot(param, real_lambda_all,'.-')
legend('lambda estimates','real lambda')
grid on
xlabel('a'); ylabel('lambda')
disp(sum(abs(lambda_all-real_lambda_all)))


%% brute force tail distribution hyper and non hyper, one para (old)
close all
testboun = 1; % whether to test some value away from boundary
boundata = 0.01; % test some data away from boundary
data_portion = 0.00000002; % portion of data point to fit
f = @(x,a)(abs(x).^(a.*log(abs(x)))); % hyperbolic tail distribution formula
g = @(x,a)(abs(x).^(a./abs(x))); % non-hyperbolic tail formula
nnn = 1e5; %number of a to test
for j = 200
    a_para = data_all{j}(3,1); % retrieve a parameter value
    epsilon = data_all{j}(4,1); % retrieve epsilon parameter value
    c = cr_all(1,j); r = cr_all(2,j); % retrieve center and radius of the box
    x = data_all{j}(1,:); % retrieve the x value
    width = (x(end)-x(1))/2^depth; % width of each bin
    h = data_all{j}(2,:)/2; % retrieve the h value (stationary distibution height)
    bounmin = c-r; % the lower boundary

    if testboun
        w = find(x<(c-r+boundata)); % only test up to epsilon away from boundary
    else
        w = find(cumsum(h)*width<data_portion); % only test some data portion
    end
    w = w(2:end); % get rid of end point 
    x1 = x(w); % very tail of the distribtuion to test
    dist1 = ones(1,nnn); % record square distance for hyper
    dist2 = ones(1,nnn); % record square distance for nonhyper
    a1 = linspace(-1/nnn,-10,nnn); % test a for hyper 
    a2 = linspace(1/nnn,10,nnn); % test a for nonhyper
    plot(x1,h(w),"Color","blue"); hold on; % plot stationary distribution

    tic;
    parfor i = 1:nnn
        a = a1(i);
        fun = f(x1-bounmin,a);
        obj = sum((fun-h(w)).^2); % objective function
        dist1(i) = obj;
    end
    [min1,argmin1] = min(dist1); % find least square solution
    plot(x1,f(x1-bounmin,a1(argmin1)),"Color",'red') % plot hyper best fit
    
    parfor i = 1:nnn
        a = a2(i);
        fun = g(x1-bounmin,a);
        obj = sum((fun-h(w)).^2); % objective function
        dist2(i) = obj;
    end
    [min2,argmin2] = min(dist2); % find least square solution
    plot(x1,g(x1-bounmin,a2(argmin2)),"Color",'black') % plot nonhyper best fit
    legend('stationary dist','hyper fit','nonhyper fit')
    title(strcat('a = ',num2str(a_para),', b = ',num2str(b_para),', epsilon = '...
        ,num2str(epsilon),' nonhyper deg: ',num2str((min1-min2)/min2)))
    disp([a_para,a1(argmin1),a2(argmin2),(min1-min2)/min2])
    % xlim([-4,4])
    % ylim([0 10])
    shg
    pause
    hold off
    toc;
end

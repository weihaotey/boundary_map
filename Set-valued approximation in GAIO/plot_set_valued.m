% For writing figures into movies or animations
b_para = 0.3;
writevideo = 0;
saveimage = 1;
if writevideo
    vidfile1 = VideoWriter('henon_test','MPEG-4');
    vidfile1.FrameRate = 2;
    open(vidfile1);
end
% Plotting of existing data
d = 2; 
% cc = [0 0]; rr = [2 1];
cc = [0.5 0.1]; rr = [1.5 0.5];
t = Tree(cc,rr);
a_start = 0.41; a_end = 0.49;a_gap = 0.01; % parameter starting and ending parameter to run in the loop
a_all = a_start:a_gap:a_end;   
i = 1;
for a_para = a_all
    % t1 = load(['tree',num2str(a_para),'.mat']);
    % boxplot2(t1.tree,'edgecolor','none'); hold on; 
    % boxplot2(t1.tree2,'color','k','edgecolor','none'); axis tight
    % boxplot2(t1.tree,'edgecolor','none'); hold off; 
    % t1 = load(['plot',num2str(a_para),'.mat']);
    % plotty = t1.plotty; plotty1 = t1.plotty1; 
    % scatter(plotty(1,:),plotty(2,:),10,'r','filled','square'); hold on;
    % scatter(plotty1(1,:),plotty1(2,:),10,'k','filled','square');
    % scatter(plotty(1,:),plotty(2,:),10,'r','filled','square'); hold off;
    tree = t_all{i}{1}; tree2 = t_all{i}{2};
    clf; boxplot2(tree,'edgecolor','none'); hold on
    boxplot2(tree2,'color','k','edgecolor','none');
    boxplot2(tree,'edgecolor','none'); hold off; axis tight
    
    

    
    title(['Minimal invariant set for Henon map with a = ' num2str(a_para) ', b = ' num2str(b_para)])
    legend('minimal invariant set','domain of attraction','Location','SouthWest')
    axis('equal')
    axis([-2+1.3*i/9 2-0.3*i/9 -1+0.6*i/9 1-0.4*i/9])
    drawnow

    %make video for both tight axis or constant axis
    if writevideo, F(count) = getframe(gcf);writeVideo(vidfile1,F(count)); end

    % save image or the tree structure for the MIS and dual repeller
    if saveimage
        exportgraphics(gcf,strcat('zoom_henon',num2str(round(a_para,4)),'.png'),'Resolution',300)
    end
    pause(.1)
    disp(a_para)
    i = i + 1;
end
if writevideo, close(vidfile1);end

%%
for i = 1:10
    a_para = 0.01*i;
    tree = t_all{i}{1}; tree2 = t_all{i}{2};
    plotty = tree.boxes(depth); plotty1 = tree2.boxes(depth);
    save(['plot',num2str(a_para),'.mat'],'plotty','plotty1','-v7.3')
end

%%
% For writing figures into movies or animations for boundary map method
b_para = 0.3;
writevideo = 0;
saveimage = 0;
if writevideo
    vidfile1 = VideoWriter('henon_reverse','MPEG-4');
    vidfile1.FrameRate = 50;
    open(vidfile1);
end
% Plotting of existing data
d = 2; 
% cc = [0 0]; rr = [2 1];
cc = [0.5 0.1]; rr = [1.5 0.5];
a_start = 0.607; a_end = 0.003;a_gap = -0.001; % parameter starting and ending parameter to run in the loop
a_all = a_start:a_gap:a_end;   
i = 1;
for a_para = a_all
    filename = ['a',num2str(a_para),'.npy'];
    if exist(filename,'file')~= 2
        continue
    end
    data = readNPY(filename);
    if isempty(data)
        continue
    end
    plot(data(1,1),data(1,2),'Markersize',1,'Color','#1f77b4'); hold on
    plot(data(:,1),data(:,2),'.','Markersize',1,'Color','#1f77b4'); hold on
    
    filename = ['unstable_a',num2str(a_para),'.npy'];
    data = readNPY(filename);
    scatter(data(:,1),data(:,2),10,'red','filled')
    filename = ['stable_a',num2str(a_para),'.npy'];
    data = readNPY(filename);
    scatter(data(:,1),data(:,2),10,'g','filled')
    filename = ['dual_a',num2str(a_para),'.npy'];
    if exist(filename,'file')== 2
        data = readNPY(filename);
        if ~isempty(data)
            plot(data(1,1),data(1,2),'Markersize',1,'Color','#FFB52E');
            plot(data(:,1),data(:,2),'.','Markersize',1,'Color','#FFB52E')
        end
        filename = ['dual_point_a',num2str(a_para),'.npy'];
        if exist(filename,'file')== 2
            data = readNPY(filename);
            scatter(data(:,1),data(:,2),10,'magenta','filled')
        end
    end
    hold off
    filename = ['dual_point_a',num2str(a_para),'.npy'];
    if exist(filename,'file')== 2 && a_para>= 0.595
        legend('boundary approximation','','saddle with 1-d unstable','stable point','dual repeller boundary','','saddle with 1-d stable','Location','SouthWest')
    else
        legend('boundary approximation','','saddle with 1-d unstable','stable point','Location','SouthWest')
    end
    title(['Boundary of minimal invariant set for Henon map with a = ' num2str(a_para) ', b = ' num2str(b_para)])
    % legend('minimal invariant set','domain of attraction','Location','SouthWest')
    axis('equal')
    axis([-1 2 -1 1])
    xlabel('x');ylabel('y');
    
    drawnow

    %make video for both tight axis or constant axis
    if writevideo, F(i) = getframe(gcf);writeVideo(vidfile1,F(i)); end

    % save image or the tree structure for the MIS and dual repeller
    if saveimage
        exportgraphics(gcf,strcat('boundary_henon',num2str(round(a_para,4)),'.png'),'Resolution',300)
    end
    % pause(.1)
    disp(a_para)
    i = i + 1;
end
if writevideo, close(vidfile1);end

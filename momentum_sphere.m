function momentum_sphere(px1,px2,px3,py1,py2,py3,pz1,pz2,pz3,binsize_p_sphere,p_plot_range)

edge_max_sphere=ceil(max([max([px1;px2;px3;py1;py2;py3;pz1;pz2;pz3]); -min([px1;px2;px3;py1;py2;py3;pz1;pz2;pz3])]));
edges_p_sphere=[-edge_max_sphere:binsize_p_sphere:edge_max_sphere];
  
figure
count_px1_py1 = histcounts2(px1,py1,edges_p_sphere,edges_p_sphere);
count_px1_py1_mod = max(count_px1_py1,1); 

sp_hand1 = subplot(3,3,1)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(edges_p_sphere,edges_p_sphere, ((count_px1_py1_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-p_plot_range p_plot_range -p_plot_range p_plot_range]);
xticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
yticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
set(gca,'FontSize',15);
pbaspect([1 1 1]);
xlabel('px1', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('py1', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10);
set(gca,'colorscale','log');
% pos1 = get(sp_hand1, 'Position') % gives the position of current sub-plot
% new_pos1 = pos1 +[0 0 0 0.05]
% set(sp_hand1, 'Position',new_pos1 ) % set new position of current sub - plot

count_py1_pz1 = histcounts2(py1,pz1,edges_p_sphere,edges_p_sphere);
count_py1_pz1_mod = max(count_py1_pz1,1); 
subplot(3,3,2)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(edges_p_sphere,edges_p_sphere, ((count_py1_pz1_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-p_plot_range p_plot_range -p_plot_range p_plot_range]);
xticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
yticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('py1', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('pz1', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');

count_pz1_px1 = histcounts2(pz1,px1,edges_p_sphere,edges_p_sphere);
count_pz1_px1_mod = max(count_pz1_px1,1); 
subplot(3,3,3)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( edges_p_sphere,edges_p_sphere, ((count_pz1_px1_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-p_plot_range p_plot_range -p_plot_range p_plot_range]);
xticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
yticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('pz1', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('px1', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');

count_px2_py2 = histcounts2(px2,py2,edges_p_sphere,edges_p_sphere);
count_px2_py2_mod = max(count_px2_py2,1); 
hold on;
subplot(3,3,4)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(edges_p_sphere,edges_p_sphere, ((count_px2_py2_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-p_plot_range p_plot_range -p_plot_range p_plot_range]);
xticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
yticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('px2', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('py2', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');


count_py2_pz2 = histcounts2(py2,pz2,edges_p_sphere,edges_p_sphere);
count_py2_pz2_mod = max(count_py2_pz2,1); 
subplot(3,3,5)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( edges_p_sphere,edges_p_sphere, ((count_py2_pz2_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-p_plot_range p_plot_range -p_plot_range p_plot_range]);
xticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
yticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('py2', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('pz2', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');


count_pz2_px2 = histcounts2(pz2,px2,edges_p_sphere,edges_p_sphere);
count_pz2_px2_mod = max(count_pz2_px2,1); 
subplot(3,3,6)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( edges_p_sphere,edges_p_sphere, ((count_pz2_px2_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-p_plot_range p_plot_range -p_plot_range p_plot_range]);
xticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
yticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('pz2', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('px2', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');

count_px3_py3 = histcounts2(px3,py3,edges_p_sphere,edges_p_sphere);
count_px3_py3_mod = max(count_px3_py3,1); 
hold on;
subplot(3,3,7)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(edges_p_sphere,edges_p_sphere, ((count_px3_py3_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-p_plot_range p_plot_range -p_plot_range p_plot_range]);
xticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
yticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('px3', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('py3', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');


count_py3_pz3 = histcounts2(py3,pz3,edges_p_sphere,edges_p_sphere);
count_py3_pz3_mod = max(count_py3_pz3,1); 
subplot(3,3,8)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( edges_p_sphere,edges_p_sphere, ((count_py3_pz3_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-p_plot_range p_plot_range -p_plot_range p_plot_range]);
xticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
yticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('py3', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('pz3', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');


count_pz3_px3 = histcounts2(pz3,px3,edges_p_sphere,edges_p_sphere);
count_pz3_px3_mod = max(count_pz3_px3,1); 
subplot(3,3,9)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( edges_p_sphere,edges_p_sphere, ((count_pz3_px3_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-p_plot_range p_plot_range -p_plot_range p_plot_range]);
xticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
yticks([-p_plot_range -p_plot_range/2 0 p_plot_range/2 p_plot_range]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('pz3', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('px3', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');

end
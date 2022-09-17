function [t1_gate_momentum,t2_gate_momentum,t3_gate_momentum,t23_gate_momentum,...
x1_gate_momentum,x2_gate_momentum,x3_gate_momentum,y1_gate_momentum,y2_gate_momentum,y3_gate_momentum]...
= survived_events(t1_gate_rot,t2_gate_rot,t3_gate_rot,t23_gate_rot,x1_gate_rot,x2_gate_rot,x3_gate_rot,...
                        y1_gate_rot, y2_gate_rot, y3_gate_rot,j_momentum_gate,frag_m_z,frag_m_z_str,t0,k0)

close all;
t1_gate_momentum=t1_gate_rot(j_momentum_gate);
t2_gate_momentum=t2_gate_rot(j_momentum_gate);
t3_gate_momentum=t3_gate_rot(j_momentum_gate);
t23_gate_momentum=t23_gate_rot(j_momentum_gate);

x1_gate_momentum=x1_gate_rot(j_momentum_gate);
x2_gate_momentum=x2_gate_rot(j_momentum_gate);
x3_gate_momentum=x3_gate_rot(j_momentum_gate);

y1_gate_momentum=y1_gate_rot(j_momentum_gate);
y2_gate_momentum=y2_gate_rot(j_momentum_gate);
y3_gate_momentum=y3_gate_rot(j_momentum_gate);

Xedges = min(t1_gate_momentum):1:max(t1_gate_momentum);
Yedges = min(t23_gate_momentum):1:max(t23_gate_momentum);

figure
subplot(2,2,1)
count = histcounts2(t1_gate_momentum,t23_gate_momentum,Xedges,Yedges);
%sum(sum(count))
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( Xedges,Yedges, (log(count')));
colorbar('FontSize', 10);
colormap(myColorMap);
colorbar 
axis xy;
% axis([0 12000 0 12000])
axis equal
hold on;
set(gca,'FontSize',10)
%pbaspect([1 1 1]);

% x=x(j_gate);
tof_range = 1000;
x_range = 2;
y_range = 3;
tedges_gate_momentum = min([t1_gate_momentum;t2_gate_momentum;t3_gate_momentum])-tof_range:4:max([t1_gate_momentum;t2_gate_momentum;t3_gate_momentum])+tof_range;
xedges_gate_momentum = min([x1_gate_momentum;x2_gate_momentum;x3_gate_momentum])-x_range:1:max([x1_gate_momentum;x2_gate_momentum;x3_gate_momentum])+x_range;
yedges_gate_momentum = min([y1_gate_momentum;y2_gate_momentum;y3_gate_momentum])-y_range:1:max([y1_gate_momentum;y2_gate_momentum;y3_gate_momentum])+y_range;


i_t_gate_momentum=1:length(tedges_gate_momentum)-1;
bincent_t_gate_momentum=[];
bincent_t_gate_momentum(i_t_gate_momentum)=(tedges_gate_momentum(i_t_gate_momentum)+tedges_gate_momentum(i_t_gate_momentum+1))/2;
 
[t1_counts_gate_momentum,t_edges_gate_momentum]=histcounts(t1_gate_momentum,tedges_gate_momentum);
[t2_counts_gate_momentum,t_edges_gate_momentum]=histcounts(t2_gate_momentum,tedges_gate_momentum);
[t3_counts_gate_momentum,t_edges_gate_momentum]=histcounts(t3_gate_momentum,tedges_gate_momentum);

%t_hist=[bincent_t' t_counts'];
%dlmwrite('t_hist.csv',t_hist);
subplot(2,2,2)
plot(bincent_t_gate_momentum,t1_counts_gate_momentum,'b',bincent_t_gate_momentum,t2_counts_gate_momentum,'r', bincent_t_gate_momentum,t3_counts_gate_momentum,'g','LineWidth',1)   ;
% grid on
xlabel('TOF/\mus','FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
xlim([min(tedges_gate_momentum) max(tedges_gate_momentum)]);
ylim([0 max(t1_counts_gate_momentum)*1.1]);
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
set(gca,'FontSize',10)
%pbaspect([1 1 1]);
hold on
for i=1:length(frag_m_z);
    xline(  k0*(frag_m_z(i))^0.5+t0, '--k', frag_m_z_str(i),'LineWidth',2,'FontSize', 16);
end


%plot xtof
subplot(2,2,3)
count_xtof_gate_momentum = histcounts2([t1_gate_momentum;t2_gate_momentum;t3_gate_momentum],[x1_gate_momentum;x2_gate_momentum;x3_gate_momentum],tedges_gate_momentum,xedges_gate_momentum);
count_mod_xtof_gate_momentum = max(count_xtof_gate_momentum,1); 
%count_mod = count; 
myColorMap=flipud(hot);
myColorMap(1,:) = 1;
imagesc( tedges_gate_momentum,xedges_gate_momentum, ((count_mod_xtof_gate_momentum')));
colorbar('FontSize', 10);
colormap(myColorMap);
axis xy;
xlabel('TOF /\mus', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('ion position(X) /mm', 'FontWeight', 'normal','FontName', 'Arial');
hold on;
set(gca,'FontSize',10)
set(gca,'colorscale','log');
%pbaspect([1 1 1]);
hold on
for i=1:length(frag_m_z);
    xline(  k0*(frag_m_z(i))^0.5+t0, '--k', frag_m_z_str(i),'LineWidth',2,'FontSize', 16);
end


%plot ytof
subplot(2,2,4)
count_ytof_gate_momentum = histcounts2([t1_gate_momentum;t2_gate_momentum;t3_gate_momentum],[y1_gate_momentum;y2_gate_momentum;y3_gate_momentum],tedges_gate_momentum,yedges_gate_momentum);
count_mod_ytof_gate_momentum = max(count_ytof_gate_momentum,1); 
%count_mod = count; 
myColorMap=flipud(hot);
myColorMap(1,:) = 1;
imagesc( tedges_gate_momentum,yedges_gate_momentum, ((count_mod_ytof_gate_momentum')));
colorbar('FontSize', 10);
colormap(myColorMap);
axis xy;
xlabel('TOF /\mus', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('ion position(Y) /mm', 'FontWeight', 'normal','FontName', 'Arial');
hold on;
set(gca,'FontSize',10)
set(gca,'colorscale','log');
%pbaspect([1 1 1]);
hold on
for i=1:length(frag_m_z);
    xline(  k0*(frag_m_z(i))^0.5+t0, '--k', frag_m_z_str(i),'LineWidth',2,'FontSize', 16);

end
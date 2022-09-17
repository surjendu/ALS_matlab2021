function detector_image(t1_gate_momentum,t2_gate_momentum,t3_gate_momentum,t23_gate_momentum,...
x1_gate_momentum,x2_gate_momentum,x3_gate_momentum,y1_gate_momentum,y2_gate_momentum,y3_gate_momentum,...
x01, y01, x02, y02, x03, y03, rcirc1, rcirc2, rcirc3)

tof_range = 20;
x_range = 3;
y_range = 3;

t1_edges_gate_momentum = min(t1_gate_momentum)-tof_range:1:max(t1_gate_momentum)+tof_range;
t2_edges_gate_momentum = min(t2_gate_momentum)-tof_range:1:max(t2_gate_momentum)+tof_range;
t3_edges_gate_momentum = min(t3_gate_momentum)-tof_range:1:max(t3_gate_momentum)+tof_range;
t23_edges_gate_momentum = min(t23_gate_momentum)-tof_range:1:max(t23_gate_momentum)+tof_range;

x1_edges_gate_momentum = min(x1_gate_momentum-x01)-x_range:1:max(x1_gate_momentum-x01)+x_range;
x2_edges_gate_momentum = min(x2_gate_momentum-x02)-x_range:1:max(x2_gate_momentum-x02)+x_range;
x3_edges_gate_momentum = min(x3_gate_momentum-x03)-x_range:1:max(x3_gate_momentum-x03)+x_range;

y1_edges_gate_momentum = min(y1_gate_momentum-y01)-y_range:1:max(y1_gate_momentum-y01)+y_range;
y2_edges_gate_momentum = min(y2_gate_momentum-y02)-y_range:1:max(y2_gate_momentum-y02)+y_range;
y3_edges_gate_momentum = min(y3_gate_momentum-y03)-y_range:1:max(y3_gate_momentum-y03)+y_range;

  
figure
count_x1_y1 = histcounts2(x1_gate_momentum-x01,y1_gate_momentum-y01,x1_edges_gate_momentum,y1_edges_gate_momentum);
count_x1_y1_mod = max(count_x1_y1,1); 

subplot(3,3,1)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(x1_edges_gate_momentum,y1_edges_gate_momentum, ((count_x1_y1_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
% axis([-300 300 -300 300]);
% xticks([-300 -150 0 150 300]);
% yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('x_1/mm', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('y_1/mm', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');
xline(  0, '--k', 'Linewidth',2 )
yline(  0, '--k', 'Linewidth',2 )

xcirc = -rcirc1:.001:rcirc1;
ycirc = sqrt(rcirc1.^2 - xcirc.^2);
hold on
plot(xcirc,ycirc,'--k', xcirc,-ycirc,'--k', 'Linewidth',2 )

% pos1 = get(sp_hand1, 'Position') % gives the position of current sub-plot
% new_pos1 = pos1 +[0 0 0 0.05]
% set(sp_hand1, 'Position',new_pos1 ) % set new position of current sub - plot


count_t1_x1 = histcounts2(t1_gate_momentum,x1_gate_momentum-x01,t1_edges_gate_momentum,x1_edges_gate_momentum);
count_t1_x1_mod = max(count_t1_x1,1); 

subplot(3,3,2)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(t1_edges_gate_momentum,x1_edges_gate_momentum, ((count_t1_x1_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
% axis([-300 300 -300 300]);
% xticks([-300 -150 0 150 300]);
% yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('t_1/ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('x_1/mm', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');


count_t1_y1 = histcounts2(t1_gate_momentum,y1_gate_momentum-y01,t1_edges_gate_momentum,y1_edges_gate_momentum);
count_t1_y1_mod = max(count_t1_y1,1); 

subplot(3,3,3)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(t1_edges_gate_momentum,y1_edges_gate_momentum, ((count_t1_y1_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
% axis([-300 300 -300 300]);
% xticks([-300 -150 0 150 300]);
% yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('t_1/ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('y_1/mm', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');



count_x2_y2 = histcounts2(x2_gate_momentum-x02,y2_gate_momentum-y02,x2_edges_gate_momentum,y2_edges_gate_momentum);
count_x2_y2_mod = max(count_x2_y2,1); 

subplot(3,3,4)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(x2_edges_gate_momentum,y2_edges_gate_momentum, ((count_x2_y2_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
% axis([-300 300 -300 300]);
% xticks([-300 -150 0 150 300]);
% yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('x_2/mm', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('y_2/mm', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');

xline(  0, '--k', 'Linewidth',2 )
yline(  0, '--k', 'Linewidth',2 )

xcirc = -rcirc2:.001:rcirc2;
ycirc = sqrt(rcirc2.^2 - xcirc.^2);
hold on
plot(xcirc,ycirc,'--k', xcirc,-ycirc,'--k', 'Linewidth',2 )

% pos1 = get(sp_hand1, 'Position') % gives the position of current sub-plot
% new_pos1 = pos1 +[0 0 0 0.05]
% set(sp_hand1, 'Position',new_pos1 ) % set new position of current sub - plot


count_t2_x2 = histcounts2(t2_gate_momentum,x2_gate_momentum-x02,t2_edges_gate_momentum,x2_edges_gate_momentum);
count_t2_x2_mod = max(count_t2_x2,1); 

subplot(3,3,5)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(t2_edges_gate_momentum,x2_edges_gate_momentum, ((count_t2_x2_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
% axis([-300 300 -300 300]);
% xticks([-300 -150 0 150 300]);
% yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('t_2/ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('x_2/mm', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');


count_t2_y2 = histcounts2(t2_gate_momentum,y2_gate_momentum-y02,t2_edges_gate_momentum,y2_edges_gate_momentum);
count_t2_y2_mod = max(count_t2_y2,1); 

subplot(3,3,6)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(t2_edges_gate_momentum,y2_edges_gate_momentum, ((count_t2_y2_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
% axis([-300 300 -300 300]);
% xticks([-300 -150 0 150 300]);
% yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('t_2/ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('y_2/mm', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');


count_x3_y3 = histcounts2(x3_gate_momentum -x03,y3_gate_momentum-y03,x3_edges_gate_momentum,y3_edges_gate_momentum);
count_x3_y3_mod = max(count_x3_y3,1); 

subplot(3,3,7)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(x3_edges_gate_momentum,y3_edges_gate_momentum, ((count_x3_y3_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
% axis([-300 300 -300 300]);
% xticks([-300 -150 0 150 300]);
% yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('x_3/mm', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('y_3/mm', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');

xline(  0, '--k', 'Linewidth',2 )
yline(  0, '--k', 'Linewidth',2 )
xcirc = -rcirc3:.001:rcirc3;
ycirc = sqrt(rcirc3.^2 - xcirc.^2);
hold on
plot(xcirc,ycirc,'--k', xcirc,-ycirc,'--k', 'Linewidth',2 )
% pos1 = get(sp_hand1, 'Position') % gives the position of current sub-plot
% new_pos1 = pos1 +[0 0 0 0.05]
% set(sp_hand1, 'Position',new_pos1 ) % set new position of current sub - plot


count_t3_x3 = histcounts2(t3_gate_momentum,x3_gate_momentum-x03,t3_edges_gate_momentum,x3_edges_gate_momentum);
count_t3_x3_mod = max(count_t3_x3,1); 

subplot(3,3,8)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(t3_edges_gate_momentum,x3_edges_gate_momentum, ((count_t3_x3_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
% axis([-300 300 -300 300]);
% xticks([-300 -150 0 150 300]);
% yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('t_3/ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('x_3/mm', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');


count_t3_y3 = histcounts2(t3_gate_momentum,y3_gate_momentum-y03,t3_edges_gate_momentum,y3_edges_gate_momentum);
count_t3_y3_mod = max(count_t3_y3,1); 

subplot(3,3,9)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(t3_edges_gate_momentum,y3_edges_gate_momentum, ((count_t3_y3_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
% axis([-300 300 -300 300]);
% xticks([-300 -150 0 150 300]);
% yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('t_3/ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('y_3/mm', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');

end
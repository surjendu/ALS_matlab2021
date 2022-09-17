function [t1_gate_rot,t2_gate_rot,t3_gate_rot,t23_gate_rot,...
            x1_gate_rot,x2_gate_rot,x3_gate_rot,...
             y1_gate_rot,y2_gate_rot,y3_gate_rot]...
            =tripico_rotation_retrieved(t1,t2,t3,x1,x2,x3,y1,y2,y3,t23,j_gate,j_rot_gate,sct_peak_tof_1,sct_peak_tof_23)

% [t1_gate_rot,t2_gate_rot,t3_gate_rot,t23_gate_rot,...
%             x1_gate_rot,x2_gate_rot,x3_gate_rot,...
%              y1_gate_rot,y2_gate_rot,y3_gate_rot,delay_step]...
%             =tripico_rotation_retrieved(t1,t2,t3,t23,j_gate,j_rot_gate,sct_peak_tof_1,sct_peak_tof_23)
%including delay

t1_gate=t1(j_gate);
t1_gate_rot=t1_gate(j_rot_gate);
t2_gate=t2(j_gate);
t2_gate_rot=t2_gate(j_rot_gate);
t3_gate=t3(j_gate);
t3_gate_rot=t3_gate(j_rot_gate);
t23_gate=t23(j_gate);
t23_gate_rot=t23_gate(j_rot_gate);




x1_gate=x1(j_gate);
x1_gate_rot=x1_gate(j_rot_gate);
x2_gate=x2(j_gate);
x2_gate_rot=x2_gate(j_rot_gate);
x3_gate=x3(j_gate);
x3_gate_rot=x3_gate(j_rot_gate);


y1_gate=y1(j_gate);
y1_gate_rot=y1_gate(j_rot_gate);
y2_gate=y2(j_gate);
y2_gate_rot=y2_gate(j_rot_gate);
y3_gate=y3(j_gate);
y3_gate_rot=y3_gate(j_rot_gate);

%   delay_step=measurement.data.raw.delay_position((j_gate));
%   delay_step=delay_step(j_rot_gate);


Xedges = min(t1_gate_rot(:,1)):5:max(t1_gate_rot(:,1));
Yedges = min(t23_gate_rot(:,1)):5:max(t23_gate_rot(:,1));

count_gate_rot = histcounts2(t1_gate_rot(:,1),t23_gate_rot(:,1),Xedges,Yedges);

figure
%sum(sum(count))
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( Xedges,Yedges, (log(count_gate_rot')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
% axis([0 12000 0 12000])
% axis equal
hold on;
set(gca,'FontSize',25)
scatter([sct_peak_tof_1], [sct_peak_tof_23],300,'x','g','LineWidth',3)

end
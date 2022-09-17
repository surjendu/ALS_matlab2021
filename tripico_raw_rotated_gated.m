function [j_rot_gate]= tripico_raw_rotated_gated(t_sum,t_diff,t_diff_bin, t_sum_bin, sct_peak_tof_1,sct_peak_tof_23)

j_rot_gate = t_diff > t_diff_bin(1)  &    t_diff  < t_diff_bin(2) &  t_sum > t_sum_bin(1)  & t_sum < t_sum_bin(2);

t_diff_gate=t_diff(j_rot_gate);
t_sum_gate=t_sum(j_rot_gate);

Xedges = min(t_diff_gate(:,1)):3:max(t_diff_gate(:,1));
Yedges = min(t_sum_gate(:,1)):3:max(t_sum_gate(:,1));

count_rotate_gate = histcounts2(t_diff_gate(:,1),t_sum_gate(:,1),Xedges,Yedges);

figure
%sum(sum(count))
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( Xedges,Yedges, (log(count_rotate_gate')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
% axis([0 12000 0 12000])
axis equal
hold on;
set(gca,'FontSize',25)
scatter([sct_peak_tof_23-sct_peak_tof_1], [sct_peak_tof_23+sct_peak_tof_1],300,'x','g','LineWidth',3)

end
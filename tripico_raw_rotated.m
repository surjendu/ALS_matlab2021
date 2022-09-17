function tripico_raw_rotated(t_sum, t_diff, sct_peak_tof_1,sct_peak_tof_23)


Xedges = min(t_diff(:,1)):5:max(t_diff(:,1));
Yedges = min(t_sum(:,1)):5:max(t_sum(:,1));

count_rotate = histcounts2(t_diff(:,1),t_sum(:,1),Xedges,Yedges);

figure
%sum(sum(count))
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( Xedges,Yedges, (log(count_rotate')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
% axis([0 12000 0 12000])
axis equal
hold on;
set(gca,'FontSize',25)
scatter([sct_peak_tof_23-sct_peak_tof_1*1], [sct_peak_tof_23+sct_peak_tof_1*1],300,'x','g','LineWidth',3)

end
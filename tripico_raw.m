function [t1,t2,t3,t23,t1_gate,t2_gate, t3_gate,t23_gate,j_gate]...
        =tripico_raw(measurement,t1_bin,t23_bin,sct_peak_tof_1,sct_peak_tof_23)

t1 = measurement.data.raw.XYT.frst(:,3);
t2 = measurement.data.raw.XYT.scnd(:,3);
t3 = measurement.data.raw.XYT.thrd(:,3);
t23=t2+t3;

j_gate = t1 > t1_bin(1)  &  t1 < t1_bin(2) &  t23 > t23_bin(1)  & t23 < t23_bin(2);

red=length(j_gate)-sum(j_gate);

t1_gate=t1(j_gate);
t2_gate=t2(j_gate);
t3_gate=t3(j_gate);
t23_gate=t23(j_gate);

Xedges = min(t1_gate(:,1)):1:max(t1_gate(:,1));
Yedges = min(t23_gate(:,1)):1:max(t23_gate(:,1));

count = histcounts2(t1_gate(:,1),t23_gate(:,1),Xedges,Yedges);

%sum(sum(count))
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( Xedges,Yedges, ((count')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis xy;
% axis([0 12000 0 12000])
axis equal
% axis equal;
pbaspect([1 1 1])

% axis([0 8 0 17]);
% xtickformat('%.0f');
% ytickformat('%.0f');
% xticks(0 :2 :10);
% yticks(0 :4 :20);

% axis([5.940 6.49 12.900 13.600]);
% xtickformat('%.1f');
% ytickformat('%.1f');
% xticks(5.940 : 0.2: 6.50);
% yticks(12.900 : 0.2 : 13.600);


xlabel('TOF_1/\mus', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('TOF_2 + TOF_3/\mus', 'FontWeight', 'normal','FontName', 'Arial');
%clims=[1 7];
set(gca,'colorscale','log');
hold on;
set(gca,'FontSize',25)
scatter([sct_peak_tof_1], [sct_peak_tof_23],300,'x','g','LineWidth',3)
%scatter(6250, 13300, 300,'x','g','LineWidth',3)

end
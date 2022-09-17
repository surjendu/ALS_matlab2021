function [t1_gate,t2_gate,t3_gate,t23_gate,Xedges,Yedges,count]=import_single_file(t1,t1,t3,t1_bin,t23_bin)

tic
t23=t2+t3;

j_gate = t1 > t1_bin(1)  &  t1 < t1_bin(2) &  t23 > t23_bin(1)  & t23 < t23_bin(2);


t1_gate=t1(j_gate);
t2_gate=t2(j_gate);
t3_gate=t3(j_gate);
t23_gate=t23(j_gate);

 Xedges = min(t1_gate(:,1)):5:max(t1_gate(:,1));
 Yedges = min(t23_gate(:,1)):5:max(t23_gate(:,1));


count = histcounts2(t1_gate(:,1),t23_gate(:,1),Xedges,Yedges);


end
clear 
clc
close all

% fileID = fopen('Bromoform_140eV_4_DAn.bin');
%fileID = fopen('Bromoform_140eV_7_DAn.bin');
%fileID = fopen('Bromocyclopentane_140eV_DAn.00008.bin');
fileID = fopen('air_140eV_Nov21_8am_Before_Bromoform_DAn.bin');
data_ini = fread(fileID,'int32');
%data_ini=data_ini(1:100,:);
%%
i=1;
nloop=0;

nel=sum(abs(data_ini) <= 9); % approximate logic to test the no of loops
if rem(nel, 2) == 0
    nel=nel;
else nel=nel+1;
end
% 28504504 %actual value
nevents=nel+1;
data=NaN(nevents, 2^6); %26 is maximum
for j=1:nevents-1
   nloop=nloop+1;
   n=3*data_ini(i)+1;
   data(nloop,1:n)=data_ini(i:i+n-1)';
   i=i+n;  
   if i-1==length(data_ini)
   break
   end
       
end
data(end,:) = [];

%%
index=[1:2:nevents-1];
data_i=data(index,:);

%%
x_i=data_i(:,2:3:64)./1000;
y_i=data_i(:,3:3:64)./1000;
t_i=data_i(:,4:3:64)./1000;
tic
nhits=8;
% combos = combntns(1:nhits,5); for 5 hits
combos = combntns(1:nhits,2);
length(combos(:,1))
t1=[];x1=[];y1=[];
t2=[];x2=[];y2=[];

parfor i=1:length(combos(:,1))
    combos_1 = combos(:,1);
    combos_2 = combos(:,2);
%     combos_3 = combos(:,3);
%     combos_4 = combos(:,4);
%     combos_5 = combos(:,5);
%     t1=[t1,eval(['t_',num2str(combos_1(i))])];
%     t2=[t2,eval(['t_',num2str(combos_2(i))])];
    
     t1=[t1;t_i(:,combos_1(i))]; x1=[x1;x_i(:,combos_1(i))]; y1=[y1;y_i(:,combos_1(i))];
     t2=[t2;t_i(:,combos_2(i))]; x2=[x2;x_i(:,combos_2(i))]; y2=[y2;y_i(:,combos_2(i))];

end
toc
%%
measurement = [t1 x1 y1 t2 x2 y2];
%  dlmwrite('measurement.csv',measurement);
save('measurement', 'measurement', '-v7.3')

%%
clear 
clc
close all
load('measurement');
t1=measurement(:,1); x1=measurement(:,2); y1=measurement(:,3);
t2=measurement(:,4); x2=measurement(:,5); y2=measurement(:,6);
%%
close all
Xedges = min(t1):1:max(t1);
Yedges = min(t2):1:max(t2);

%t1=data_i(:,4)./1000;
%t2=data_i(:,7)./1000;
count = histcounts2(t1,t2,Xedges,Yedges);

%sum(sum(count))
% subplot(1,2,1)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( Xedges,Yedges, (count'));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
% axis([0 12000 0 12000])
axis equal
pbaspect([1 1 1]);
hold on;
set(gca,'FontSize',25)
set(gca,'colorscale','log');
% scatter([sct_peak_tof_1], [sct_peak_tof_2],300,'x','g','LineWidth',3)

xlabel('TOF_1 /ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('TOF_2 /ns', 'FontWeight', 'normal','FontName', 'Arial');

%
mmns_to_au_conv=0.4571; %mm/ns to a.u. conversion
har_to_ev=  27.211396;
mpvsme=1822.888486424682;
H = 1.00782503223;
D = 2.01410177812;
He =  4.00260325413;
Li = 7.0160034366;
C = 12;
N = 14.00307400443;
O = 15.994914619;
F = 18.99840316273;
Ne =  19.9924401762;
S = 31.9720711744;
Cl = 34.968852682;
Ar = 39.9623831237;
Br79 = 78.9183376;
Br81 = 80.9162906;
I = 126.9044719;
C13 = 13.00335483521;

%frag={ {'H', 1,1} {'He', 1,1}  {'C', 1,1} {'N',1, 1}  {'C',1,'N',1, 1} {'C', 2,1}  {'C', 3,1}  {'C', 6,1} ...
 %   {'I', 1,1} {'C',10,'H',7, 1} {'C',11,'H',7,'N', 1, 1} };
% 
% frag={ {'Br81',1,1}    {'C',1,'H',1, 'Br81',2,1} };
% frag={ {'Br79',1,1}    {'C',1,'H',1, 'Br79',2,1} };
 frag={ {'N',1,1}    {'N',1,1} };
% frag={ {'O',1,1}    {'O',1,1} };
% 
% frag={ {'Br79',1,1}   {'H',1,'Br79',1,1} {'C',1,'H',1, 'Br79',1,1} { 'Br79',2,1} {'C',1, 'Br79',2,1} {'C',1,'H',1, 'Br79',2,1} ...
%     {'Br81',1,1}   {'H',1,'Br81',1,1} {'C',1,'H',1, 'Br81',1,1} { 'Br81',2,1} {'C',1, 'Br81',2,1} {'C',1,'H',1, 'Br81',2,1}};

k0 = 11.2766;
t0= - 288.0202;


frag_m=[]; 
frag_m_z=[];
charge_z=[];
 frag_m_z_str=string.empty;
 for i=1:length(frag)
     frag_m_z_int_sum=0;
    frag_m_z_int_str=[];
     for j=1:(length(frag{i})-1)/2
          frag_m_z_int=frag{i}{2*j}*eval(frag{i}{2*j-1});
          frag_m_z_int_sum = frag_m_z_int_sum + frag_m_z_int;
          frag_m_z_int=[(frag{i}{2*j-1}),'_',num2str(frag{i}{2*j})];
          frag_m_z_int_str=[frag_m_z_int_str,frag_m_z_int];

     end
     charge_z=[charge_z,frag{i}{2*j+1}];
     frag_m=[frag_m, frag_m_z_int_sum];
     frag_m_z=[frag_m_z, frag_m_z_int_sum/frag{i}{2*j+1}];
     frag_m_z_str(i)=[frag_m_z_int_str '^' num2str(frag{i}{2*j+1}) '^' '+' ];
     
 end
[charge_z;frag_m;frag_m_z];
 frag_m=frag_m*mpvsme-1;%q=1
 frag_m_z=frag_m_z*mpvsme-1;%q=1

 
hold on
for i=1:length(frag_m_z);
    xline(  k0*(frag_m_z(i))^0.5+t0, '--r', frag_m_z_str(i))
    yline(  k0*(frag_m_z(i))^0.5+t0, '--r', frag_m_z_str(i))
end

%  sct_peak_tof_1=3984; % in ns %79
%  sct_peak_tof_2=6015;  % in ns 
% %  
%   sct_peak_tof_1=4047; % in ns %81
%  sct_peak_tof_2=6082;  % in ns 

sct_peak_tof_1=1512; % in ns %N+
sct_peak_tof_2=1512;  % in ns 
scatter([sct_peak_tof_1], [sct_peak_tof_2],300,'x','k','LineWidth',3)
%%
t1_bin = [1400 1650] ;
t2_bin = [1400 1650] ;
j_gate = t1 > t1_bin(1)  &  t1 < t1_bin(2) &  t2 > t2_bin(1)  & t2 < t2_bin(2);

red=length(j_gate)-sum(j_gate);

t1_gate=t1(j_gate);
t2_gate=t2(j_gate);

Xedges = min(t1_gate(:,1)):5:max(t1_gate(:,1));
Yedges = min(t2_gate(:,1)):5:max(t2_gate(:,1));

count = histcounts2(t1_gate(:,1),t2_gate(:,1),Xedges,Yedges);

%sum(sum(count))
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( Xedges,Yedges, (log(count')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
%axis([0 10000 0 10000])
axis equal
hold on;
set(gca,'FontSize',25)
scatter([sct_peak_tof_1], [sct_peak_tof_2],300,'x','g','LineWidth',3)

%%
%rotate
t_sum=t1_gate*1+t2_gate;
t_diff=t2_gate-t1_gate*1;

Xedges = min(t_diff(:,1)):2:max(t_diff(:,1));
Yedges = min(t_sum(:,1)):2:max(t_sum(:,1));

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
scatter([sct_peak_tof_2-sct_peak_tof_1*1], [sct_peak_tof_1*1+sct_peak_tof_2],300,'x','g','LineWidth',3)

%%
t_diff_bin = [-20 80]; t_sum_bin = [3005 3050]; %for 79
%  t_diff_bin = [1900 2150]; t_sum_bin = [10110 10140]; %for 81
% t_diff_bin = [1900 2067]; t_sum_bin = [10110 10140]; %for 81

j_rot_gate = t_diff > t_diff_bin(1)  &    t_diff  < t_diff_bin(2) &  t_sum > t_sum_bin(1)  & t_sum < t_sum_bin(2);

t_diff_gate=t_diff(j_rot_gate);
t_sum_gate=t_sum(j_rot_gate);

Xedges = min(t_diff_gate(:,1)):1:max(t_diff_gate(:,1));
Yedges = min(t_sum_gate(:,1)):1:max(t_sum_gate(:,1));

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
scatter([sct_peak_tof_2-sct_peak_tof_1], [sct_peak_tof_1+sct_peak_tof_2],300,'x','g','LineWidth',3)
%%
t1_gate_rot=t1(j_gate);
t1_gate_rot=t1_gate_rot(j_rot_gate);
t2_gate_rot=t2(j_gate);
t2_gate_rot=t2_gate_rot(j_rot_gate);


x1_gate_rot=x1(j_gate);
x1_gate_rot=x1_gate_rot(j_rot_gate);
x2_gate_rot=x2(j_gate);
x2_gate_rot=x2_gate_rot(j_rot_gate);

y1_gate_rot=y1(j_gate);
y1_gate_rot=y1_gate_rot(j_rot_gate);
y2_gate_rot=y2(j_gate);
y2_gate_rot=y2_gate_rot(j_rot_gate);

Xedges = min(t1_gate_rot(:,1)):1:max(t1_gate_rot(:,1));
Yedges = min(t2_gate_rot(:,1)):1:max(t2_gate_rot(:,1));

count_gate_rot = histcounts2(t1_gate_rot(:,1),t2_gate_rot(:,1),Xedges,Yedges);

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
axis equal
hold on;
set(gca,'FontSize',25)

%% Analysis Constants
 clc
close all
mmns_to_au_conv=0.4571*10^(-3); %mm/us to a.u. conversion
run_name = '';
hatoev=  27.211396; % taken from ksu webpage/ twice of ionization potential of H




% % % x01 =+1+.4; mx1 =5.775525175273884537e-01; cx1=-2.650508597571916633e-02;
% % % y01 = 4.6-3-1-.2;    my1 =mx1; cy1=cx1;
% % % t01 = - 288.0202+2*0;   mt1 =-3.391475337227063402e-01; ct1=6.110716762948280802e+02;
% % % 
% % % 
% % % x02 = 1+0.3; mx2 =5.775525175273884537e-01; cx2=-2.650508597571916633e-02;
% % % y02 =5.15-2.8-1-0.1; my2 =mx2; cy2=cx2;
% % % t02 = t01; mt2 =-3.391475337227063402e-01; ct2=6.110716762948280802e+02;




x01 =1.5; mx1 =5.130099261334170047e-01; cx1=-2.201362844723278242e-02;
y01 = 0.9;    my1 =mx1; cy1=cx1;
t01 = -288.35;   mt1 =-3.812251057839449175e-01; ct1=6.869644326205083189e+02;


x02 = 1.3; mx2 =5.130099261334170047e-01; cx2=-2.201362844723278242e-02;
y02 =0.9; my2 =mx2; cy2=cx2;
t02 = t01; mt2 =-3.812251057839449175e-01; ct2=6.869644326205083189e+02;





binsize_px=2; 
binsize_py=2; 
binsize_pz=2;
binsize_p=2;

%%
momentum_gui

%%
%calculating momenta and parameter tuning
% % % % px1 = cal_fac*frag_m_z(1)*((x1_gate_rot-x01).*mx1 + cx1 )*mmns_to_au_conv;
% % % % py1 = cal_fac*frag_m_z(1)*((y1_gate_rot-y01).*my1 + cy1 )*mmns_to_au_conv;
% % % % pz1 = frag_m_z(1)*((t1_gate_rot-t01).*mt1 + ct1)*mmns_to_au_conv;  
% % % % p1=sqrt(px1.^2+py1.^2+pz1.^2);   
% % % % 
% % % % px2 = cal_fac*frag_m_z(2)*((x2_gate_rot-x02).*mx2 + cx2 )*mmns_to_au_conv;
% % % % py2 = cal_fac*frag_m_z(2)*((y2_gate_rot-y02).*my2 + cy2 )*mmns_to_au_conv;
% % % % pz2 = frag_m_z(2)*((t2_gate_rot-t02).*mt2 + ct2)*mmns_to_au_conv;  
% % % % p2=sqrt(px2.^2+py2.^2+pz2.^2);


% % % % % % % % % % % % random 1st and 2nd
% % % % % % % % % % % p_raw=[px1 py1 pz1 px2 py2 pz2 ];
% % % % % % % % % % % ind_rand = zeros(length(p_raw(:,1)),6);
% % % % % % % % % % % 
% % % % % % % % % % % for i=1:length(p_raw(:,1))
% % % % % % % % % % % for j=1:6
% % % % % % % % % % % 
% % % % % % % % % % % ind_rand_12 = (randperm(2))*3-2;
% % % % % % % % % % % 
% % % % % % % % % % % ind_rand_12 = [ind_rand_12(1) ind_rand_12(1)+1  ind_rand_12(1)+2 ...
% % % % % % % % % % %     ind_rand_12(2) ind_rand_12(2)+1  ind_rand_12(2)+2 ];
% % % % % % % % % % % 
% % % % % % % % % % % ind_rand(i,:) = ind_rand_12;
% % % % % % % % % % % 
% % % % % % % % % % % end
% % % % % % % % % % % end
% % % % % % % % % % % 
% % % % % % % % % % % p_raw_new = zeros(length(p_raw(:,1)),6);
% % % % % % % % % % % 
% % % % % % % % % % % for i=1:length(p_raw(:,1))
% % % % % % % % % % % for j=1:6
% % % % % % % % % % %     
% % % % % % % % % % %     p_raw_new(i,j)=p_raw(i,ind_rand(i,j)) ;
% % % % % % % % % % % end
% % % % % % % % % % % end
% % % % % % % % % % % 
% % % % % % % % % % % px1=p_raw_new(:,1); py1=p_raw_new(:,2); pz1=p_raw_new(:,3);
% % % % % % % % % % % px2=p_raw_new(:,4); py2=p_raw_new(:,5); pz2=p_raw_new(:,6);
% % % % % % % % % % % %
%

% 
% % % % % edge_max=ceil(max([max([px1;px2;py1;py2;pz1;pz2]); -min([px1;px2;py1;py2;pz1;pz2])]));
% % % % % 
% % % % % edges_px=[-edge_max:binsize_px:edge_max];
% % % % % edges_py=[-edge_max:binsize_py:edge_max]; 
% % % % % edges_pz=[-edge_max:binsize_pz:edge_max];
% % % % % edges_p=[min([p1;p2]):binsize_p:max([p1;p2])];
% % % % % 
% % % % % i_px=1:length(edges_px)-1;
% % % % % i_py=1:length(edges_py)-1;
% % % % % i_pz=1:length(edges_pz)-1;
% % % % % i_p=1:length(edges_p)-1;
% % % % % 
% % % % % bincent_px=[];
% % % % % bincent_py=[];
% % % % % bincent_pz=[];
% % % % % bincent_p=[];
% % % % % 
% % % % % 
% % % % % bincent_px(i_px)=(edges_px(i_px)+edges_px(i_px+1))/2;
% % % % % bincent_py(i_py)=(edges_py(i_py)+edges_py(i_py+1))/2;
% % % % % bincent_pz(i_pz)=(edges_pz(i_pz)+edges_pz(i_pz+1))/2;
% % % % % bincent_p(i_p)=(edges_p(i_p)+edges_p(i_p+1))/2;
% % % % % 
% % % % % 
% % % % % [px1_counts,px1_edges]=histcounts(px1,edges_px);
% % % % % [py1_counts,py1_edges]=histcounts(py1,edges_py);
% % % % % [pz1_counts,pz1_edges]=histcounts(pz1,edges_pz);
% % % % % [p1_counts,p1_edges]=histcounts(p1,edges_p);
% % % % % 
% % % % % 
% % % % % [px2_counts,px2_edges]=histcounts(px2,edges_px);
% % % % % [py2_counts,py2_edges]=histcounts(py2,edges_py);
% % % % % [pz2_counts,pz2_edges]=histcounts(pz2,edges_pz);
% % % % % [p2_counts,p2_edges]=histcounts(p2,edges_p);
% % % % % 
% % % % % [px_counts,px_edges]=histcounts(px1+px2,edges_px);
% % % % % [py_counts,py_edges]=histcounts(py1+py2,edges_py);
% % % % % [pz_counts,pz_edges]=histcounts(pz1+pz2,edges_pz);
% % % % % 
% % % % % % p_allhits_hist={[bincent_px', px_counts'], [bincent_py', py_counts'], [bincent_pz' pz_counts'], [bincent_p' p_counts']};
% % % % % %  dlmwrite('p_allhits_hist.csv',p_allhits_hist);
% % % % % 
% % % % % close all;
% % % % % %figure to test calibration
% % % % % subplot(2,2,1)
% % % % % plot(bincent_px,px1_counts,'r',bincent_px,fliplr(px1_counts),'b',bincent_px,(px1_counts+fliplr(px1_counts))./2,'k','LineWidth',2);
% % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % grid on;
% % % % % xlim([-300 300]);
% % % % % ylim([0 max(px1_counts)*1.1]);
% % % % % xlabel('px1/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % legend({'px1','-px1','px1\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
% % % % % 
% % % % % 
% % % % % subplot(2,2,2)
% % % % % %figure
% % % % % plot(bincent_px,px2_counts,'r',bincent_px,fliplr(px2_counts),'b',bincent_px,(px2_counts+fliplr(px2_counts))./2,'k','LineWidth',2);
% % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % grid on;
% % % % % xlim([-300 300]);
% % % % % ylim([0 max(px2_counts)*1.1]);
% % % % % xlabel('px2/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % % ytickformat('%0.0e')
% % % % % legend({'px2','-px2','px2\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
% % % % % 
% % % % % subplot(2,2,3) 
% % % % % %figure
% % % % % plot(bincent_px,px_counts,'r','LineWidth',2);
% % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % grid on;
% % % % % xlim([-100 100]);
% % % % % %ylim([0 max(px_counts)*1.1]);
% % % % % xlabel('\Sigmapx/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % legend({'\Sigmapx'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
% % % % % 
% % % % % 
% % % % % 
% % % % % figure
% % % % % subplot(2,2,1)
% % % % % plot(bincent_py,py1_counts,'r',bincent_py,fliplr(py1_counts),'b',bincent_py,(py1_counts+fliplr(py1_counts))./2,'k','LineWidth',2);
% % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % grid on;
% % % % % xlim([-300 300]);
% % % % % ylim([0 max(py1_counts)*1.1]);
% % % % % xlabel('px1/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % legend({'py1','-py1','py1\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
% % % % % 
% % % % % 
% % % % % subplot(2,2,2)
% % % % % %figure
% % % % % plot(bincent_py,py2_counts,'r',bincent_py,fliplr(py2_counts),'b',bincent_py,(py2_counts+fliplr(py2_counts))./2,'k','LineWidth',2);
% % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % grid on;
% % % % % xlim([-300 300]);
% % % % % ylim([0 max(py2_counts)*1.1]);
% % % % % xlabel('py2/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % % ytickformat('%0.0e')
% % % % % legend({'py2','-py2','py2\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
% % % % % 
% % % % % subplot(2,2,3) 
% % % % % %figure
% % % % % plot(bincent_py,py_counts,'r','LineWidth',2);
% % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % grid on;
% % % % % xlim([-100 100]);
% % % % % %ylim([0 max(py_counts)*1.1]);
% % % % % xlabel('\Sigmapy/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % legend({'\Sigmapy'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
% % % % % 
% % % % % 
% % % % % 
% % % % % figure
% % % % % subplot(2,2,1)
% % % % % plot(bincent_pz,pz1_counts,'r',bincent_pz,fliplr(pz1_counts),'b',bincent_pz,(pz1_counts+fliplr(pz1_counts))./2,'k','LineWidth',2);
% % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % grid on;
% % % % % xlim([-300 300]);
% % % % % ylim([0 max(pz1_counts)*1.1]);
% % % % % xlabel('pz1/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % legend({'pz1','-pz1','pz1\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
% % % % % 
% % % % % 
% % % % % subplot(2,2,2)
% % % % % %figure
% % % % % plot(bincent_pz,pz2_counts,'r',bincent_pz,fliplr(pz2_counts),'b',bincent_pz,(pz2_counts+fliplr(pz2_counts))./2,'k','LineWidth',2);
% % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % grid on;
% % % % % xlim([-300 300]);
% % % % % ylim([0 max(pz2_counts)*1.1]);
% % % % % xlabel('pz2/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % % ytickformat('%0.0e')
% % % % % legend({'pz2','-pz2','pz2\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
% % % % % 
% % % % % subplot(2,2,3) 
% % % % % %figure
% % % % % plot(bincent_pz,pz_counts,'r','LineWidth',2);
% % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % grid on;
% % % % % xlim([-100 100]);
% % % % % %ylim([0 max(pz_counts)*1.1]);
% % % % % xlabel('\Sigmapz/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % legend({'\Sigmapz'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
%%
close all
px1_rng = [-600.0 600.0]; px2_rng = [-600 600];     px_sum_rng = [-11 11];
py1_rng = [-600.0 600.0]; py2_rng = [-600 600];     py_sum_rng = [-30 30];
pz1_rng = [-600.0 600.0]; pz2_rng = [-600 600];     pz_sum_rng = [-14 14];

j_momentum_gate = px1 > px1_rng(1)  &    px1  < px1_rng(2) &  py1 > py1_rng(1)  &    py1  < py1_rng(2) &  pz1 > pz1_rng(1)  &    pz1  < pz1_rng(2) & px2 > px2_rng(1)  &    px2  < px2_rng(2) &  py2 > py2_rng(1)...
    &    py2  < py2_rng(2) &  pz2 > pz2_rng(1)  &    pz2  < pz2_rng(2) & px1+px2 > px_sum_rng(1)  &    px1+px2  < px_sum_rng(2) &  py1+py2 > py_sum_rng(1)  &    py1+py2  < py_sum_rng(2) &  pz1+pz2 >  pz_sum_rng(1)  &   pz1+pz2< pz_sum_rng(2);
 
px1 = px1(j_momentum_gate);
py1 = py1(j_momentum_gate);
pz1 = pz1(j_momentum_gate);
p1=p1(j_momentum_gate);   


px2 = px2(j_momentum_gate);
py2 = py2(j_momentum_gate);
pz2 = pz2(j_momentum_gate);
p2= p2(j_momentum_gate);

    
%plot tof in ns
binsize_px=1; 
binsize_py=1; 
binsize_pz=1;
binsize_p=1;

%%
momenta_equalization
%%
% % % % % % % edge_max=ceil(max([max([px1;px2;py1;py2;pz1;pz2]); -min([px1;px2;py1;py2;pz1;pz2])]));
% % % % % % % 
% % % % % % % edges_px=[-edge_max:binsize_px:edge_max];
% % % % % % % edges_py=[-edge_max:binsize_py:edge_max]; 
% % % % % % % edges_pz=[-edge_max:binsize_pz:edge_max];
% % % % % % % edges_p=[min([p1;p2]):binsize_p:max([p1;p2])];
% % % % % % % 
% % % % % % % i_px=1:length(edges_px)-1;
% % % % % % % i_py=1:length(edges_py)-1;
% % % % % % % i_pz=1:length(edges_pz)-1;
% % % % % % % i_p=1:length(edges_p)-1;
% % % % % % % 
% % % % % % % bincent_px=[];
% % % % % % % bincent_py=[];
% % % % % % % bincent_pz=[];
% % % % % % % bincent_p=[];
% % % % % % % 
% % % % % % % 
% % % % % % % bincent_px(i_px)=(edges_px(i_px)+edges_px(i_px+1))/2;
% % % % % % % bincent_py(i_py)=(edges_py(i_py)+edges_py(i_py+1))/2;
% % % % % % % bincent_pz(i_pz)=(edges_pz(i_pz)+edges_pz(i_pz+1))/2;
% % % % % % % bincent_p(i_p)=(edges_p(i_p)+edges_p(i_p+1))/2;
% % % % % % % 
% % % % % % % 
% % % % % % % [px1_counts,px1_edges]=histcounts(px1,edges_px);
% % % % % % % [py1_counts,py1_edges]=histcounts(py1,edges_py);
% % % % % % % [pz1_counts,pz1_edges]=histcounts(pz1,edges_pz);
% % % % % % % [p1_counts,p1_edges]=histcounts(p1,edges_p);
% % % % % % % 
% % % % % % % 
% % % % % % % [px2_counts,px2_edges]=histcounts(px2,edges_px);
% % % % % % % [py2_counts,py2_edges]=histcounts(py2,edges_py);
% % % % % % % [pz2_counts,pz2_edges]=histcounts(pz2,edges_pz);
% % % % % % % [p2_counts,p2_edges]=histcounts(p2,edges_p);
% % % % % % % 
% % % % % % % [px_counts,px_edges]=histcounts(px1+px2,edges_px);
% % % % % % % [py_counts,py_edges]=histcounts(py1+py2,edges_py);
% % % % % % % [pz_counts,pz_edges]=histcounts(pz1+pz2,edges_pz);
% % % % % % % 
% % % % % % % % p_allhits_hist={[bincent_px', px_counts'], [bincent_py', py_counts'], [bincent_pz' pz_counts'], [bincent_p' p_counts']};
% % % % % % % %  dlmwrite('p_allhits_hist.csv',p_allhits_hist);
% % % % % % % % test calibration
% % % % % % % close all;
% % % % % % % %figure 
% % % % % % % subplot(2,2,1)
% % % % % % % plot(bincent_px,px1_counts,'r',bincent_py,py1_counts,'b',bincent_pz,pz1_counts,'k','LineWidth',2);
% % % % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % % % grid on;
% % % % % % % xlim([-300 300]);
% % % % % % % %ylim([0 max(px1_counts)*1.1]);
% % % % % % % xlabel('px1/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % % % legend({'px1','py1','pz1'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
% % % % % % % 
% % % % % % % subplot(2,2,3)
% % % % % % % plot(bincent_px,px2_counts,'r',bincent_py,py2_counts,'b',bincent_pz,pz2_counts,'k','LineWidth',2);
% % % % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % % % grid on;
% % % % % % % xlim([-300 300]);
% % % % % % % %ylim([0 max(px1_counts)*1.1]);
% % % % % % % xlabel('px1/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % % % legend({'px2','py2','pz2'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
% % % % % % % 
% % % % % % % 
% % % % % % % subplot(2,2,2)
% % % % % % % plot(bincent_px,px1_counts,'r',bincent_py,py2_counts,'b','LineWidth',2);
% % % % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % % % grid on;
% % % % % % % xlim([-300 300]);
% % % % % % % %ylim([0 max(px1_counts)*1.1]);
% % % % % % % xlabel('px1/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % % % legend({'px1','py2'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
% % % % % % % 
% % % % % % % 
% % % % % % % subplot(2,2,4)
% % % % % % % plot(bincent_px,px2_counts,'r',bincent_py,py1_counts,'b','LineWidth',2);
% % % % % % % set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
% % % % % % % grid on;
% % % % % % % xlim([-300 300]);
% % % % % % % %ylim([0 max(px1_counts)*1.1]);
% % % % % % % xlabel('px1/au',  'FontWeight', 'normal','FontName', 'Arial');
% % % % % % % ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% % % % % % % legend({'px2','py1'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
%%
close all;
px1=px1_cal; px2=px2_cal;py1=py1_cal; py2=py2_cal; pz1=pz1_cal; pz2=pz2_cal;

py1_pz1 = [py1 pz1];
py2_pz2 = [py2 pz2];

save('py1_pz1','py1_pz1','-v7.3');
save('py2_pz2','py2_pz2','-v7.3');

edge_max=ceil(max([max([px1;px2;py1;py2;pz1;pz2]); -min([px1;px2;py1;py2;pz1;pz2])]));

edges_px=[-edge_max:binsize_px:edge_max];
edges_py=[-edge_max:binsize_py:edge_max]; 
edges_pz=[-edge_max:binsize_pz:edge_max];
edges_p=[min([p1;p2]):binsize_p:max([p1;p2])];

i_px=1:length(edges_px)-1;
i_py=1:length(edges_py)-1;
i_pz=1:length(edges_pz)-1;
i_p=1:length(edges_p)-1;

bincent_px=[];
bincent_py=[];
bincent_pz=[];
bincent_p=[];


bincent_px(i_px)=(edges_px(i_px)+edges_px(i_px+1))/2;
bincent_py(i_py)=(edges_py(i_py)+edges_py(i_py+1))/2;
bincent_pz(i_pz)=(edges_pz(i_pz)+edges_pz(i_pz+1))/2;
bincent_p(i_p)=(edges_p(i_p)+edges_p(i_p+1))/2;


[px1_counts,px1_edges]=histcounts(px1,edges_px);
[py1_counts,py1_edges]=histcounts(py1,edges_py);
[pz1_counts,pz1_edges]=histcounts(pz1,edges_pz);
[p1_counts,p1_edges]=histcounts(p1,edges_p);


[px2_counts,px2_edges]=histcounts(px2,edges_px);
[py2_counts,py2_edges]=histcounts(py2,edges_py);
[pz2_counts,pz2_edges]=histcounts(pz2,edges_pz);
[p2_counts,p2_edges]=histcounts(p2,edges_p);

[px_counts,px_edges]=histcounts(px1+px2,edges_px);
[py_counts,py_edges]=histcounts(py1+py2,edges_py);
[pz_counts,pz_edges]=histcounts(pz1+pz2,edges_pz);


%figure 
subplot(2,2,1)
plot(bincent_px,px1_counts,'r',bincent_px,fliplr(px1_counts),'b',bincent_px,(px1_counts+fliplr(px1_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-400 400]);
ylim([0 max(px1_counts)*1.1]);
xlabel('px1/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'px1','-px1','px1\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')


subplot(2,2,2)
%figure
plot(bincent_px,px2_counts,'r',bincent_px,fliplr(px2_counts),'b',bincent_px,(px2_counts+fliplr(px2_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-400 400]);
ylim([0 max(px2_counts)*1.1]);
xlabel('px2/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% ytickformat('%0.0e')
legend({'px2','-px2','px2\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

subplot(2,2,3) 
%figure
plot(bincent_px,px_counts,'r','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-400 400]);
ylim([0 max(px_counts)*1.1]);
xlabel('\Sigmapx/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'\Sigmapx'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')



figure
subplot(2,2,1)
plot(bincent_py,py1_counts,'r',bincent_py,fliplr(py1_counts),'b',bincent_py,(py1_counts+fliplr(py1_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-400 400]);
ylim([0 max(py1_counts)*1.1]);
xlabel('px1/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'py1','-py1','py1\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')


subplot(2,2,2)
%figure
plot(bincent_py,py2_counts,'r',bincent_py,fliplr(py2_counts),'b',bincent_py,(py2_counts+fliplr(py2_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-400 400]);
ylim([0 max(py2_counts)*1.1]);
xlabel('py2/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% ytickformat('%0.0e')
legend({'py2','-py2','py2\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

subplot(2,2,3) 
%figure
plot(bincent_py,py_counts,'r','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-400 400]);
ylim([0 max(py_counts)*1.1]);
xlabel('\Sigmapy/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'\Sigmapy'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')



figure
subplot(2,2,1)
plot(bincent_pz,pz1_counts,'r',bincent_pz,fliplr(pz1_counts),'b',bincent_pz,(pz1_counts+fliplr(pz1_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-400 400]);
ylim([0 max(pz1_counts)*1.1]);
xlabel('pz1/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'pz1','-pz1','pz1\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')


subplot(2,2,2)
%figure
plot(bincent_pz,pz2_counts,'r',bincent_pz,fliplr(pz2_counts),'b',bincent_pz,(pz2_counts+fliplr(pz2_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-400 400]);
ylim([0 max(pz2_counts)*1.1]);
xlabel('pz2/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% ytickformat('%0.0e')
legend({'pz2','-pz2','pz2\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

subplot(2,2,3) 
%figure
plot(bincent_pz,pz_counts,'r','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([-50 50]);
ylim([0 max(pz_counts)*1.1]);
xlabel('\Sigmapz/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'\Sigmapz'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
%%
close all;

binsize_p_sphere=4; 
edge_max_sphere=ceil(max([max([px1;px2;py1;py2;pz1;pz2]); -min([px1;px2;py1;py2;pz1;pz2])]));
edges_p_sphere=[-edge_max_sphere:binsize_p_sphere:edge_max_sphere];
  
figure
count_px1_py1 = histcounts2(px1,py1,edges_p_sphere,edges_p_sphere);
count_px1_py1_mod = max(count_px1_py1,1); 

subplot(2,3,1)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(edges_p_sphere,edges_p_sphere, ((count_px1_py1_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-300 300 -300 300]);
xticks([-300 -150 0 150 300]);
yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('p1x', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('p1y', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');
% pos1 = get(sp_hand1, 'Position') % gives the position of current sub-plot
% new_pos1 = pos1 +[0 0 0 0.05]
% set(sp_hand1, 'Position',new_pos1 ) % set new position of current sub - plot

count_py1_pz1 = histcounts2(py1,pz1,edges_p_sphere,edges_p_sphere);
save('count_py1_pz1.mat', 'count_py1_pz1','-v7.3');
count_py1_pz1_mod = max(count_py1_pz1,1); 
subplot(2,3,2)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(edges_p_sphere,edges_p_sphere, ((count_py1_pz1_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-300 300 -300 300]);
xticks([-300 -150 0 150 300]);
yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('py1', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('pz1', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');


count_px2_py2 = histcounts2(px2,py2,edges_p_sphere,edges_p_sphere);
count_px2_py2_mod = max(count_px2_py2,1); 
hold on;
subplot(2,3,3)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(edges_p_sphere,edges_p_sphere, ((count_px2_py2_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-300 300 -300 300]);
xticks([-300 -150 0 150 300]);
yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('px2', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('py2', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');


count_pz1_px1 = histcounts2(pz1,px1,edges_p_sphere,edges_p_sphere);
count_pz1_px1_mod = max(count_pz1_px1,1); 
subplot(2,3,4)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( edges_p_sphere,edges_p_sphere, ((count_pz1_px1_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-300 300 -300 300]);
xticks([-300 -150 0 150 300]);
yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('pz1', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('px1', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');


count_py2_pz2 = histcounts2(py2,pz2,edges_p_sphere,edges_p_sphere);
save('count_py2_pz2', 'count_py2_pz2','-v7.3');
count_py2_pz2_mod = max(count_py2_pz2,1); 
subplot(2,3,5)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( edges_p_sphere,edges_p_sphere, ((count_py2_pz2_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-300 300 -300 300]);
xticks([-300 -150 0 150 300]);
yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('py2', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('pz2', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');


count_pz2_px2 = histcounts2(pz2,px2,edges_p_sphere,edges_p_sphere);
count_pz2_px2_mod = max(count_pz2_px2,1); 
subplot(2,3,6)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc( edges_p_sphere,edges_p_sphere, ((count_pz2_px2_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar off
axis xy;
axis([-300 300 -300 300]);
xticks([-300 -150 0 150 300]);
yticks([-300 -150 0 150 300]);
set(gca,'FontSize',15)
pbaspect([1 1 1]);
xlabel('pz2', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('px2', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'FontSize',10)
set(gca,'colorscale','log');

%%
% things survived
close all;
t1_gate_momentum=t1_gate_rot(j_momentum_gate);
t2_gate_momentum=t2_gate_rot(j_momentum_gate);

x1_gate_momentum=x1_gate_rot(j_momentum_gate);
x2_gate_momentum=x2_gate_rot(j_momentum_gate);

y1_gate_momentum=y1_gate_rot(j_momentum_gate);
y2_gate_momentum=y2_gate_rot(j_momentum_gate);


Xedges = min(t1_gate_momentum):1:max(t1_gate_momentum);
Yedges = min(t2_gate_momentum):1:max(t2_gate_momentum);

figure
subplot(2,2,1)
count = histcounts2(t1_gate_momentum,t2_gate_momentum,Xedges,Yedges);
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
tedges_gate_momentum = min([t1_gate_momentum;t2_gate_momentum])-tof_range:1:max([t1_gate_momentum;t2_gate_momentum])+tof_range;
xedges_gate_momentum = min([x1_gate_momentum;x2_gate_momentum])-x_range:1:max([x1_gate_momentum;x2_gate_momentum])+x_range;
yedges_gate_momentum = min([y1_gate_momentum;y2_gate_momentum])-y_range:1:max([y1_gate_momentum;y2_gate_momentum])+y_range;


i_t_gate_momentum=1:length(tedges_gate_momentum)-1;
bincent_t_gate_momentum=[];
bincent_t_gate_momentum(i_t_gate_momentum)=(tedges_gate_momentum(i_t_gate_momentum)+tedges_gate_momentum(i_t_gate_momentum+1))/2;
 
[t_counts_gate_momentum,t_edges_gate_momentum]=histcounts([t1_gate_momentum;t2_gate_momentum],tedges_gate_momentum);
[t1_counts_gate_momentum,t_edges_gate_momentum]=histcounts([t1_gate_momentum],tedges_gate_momentum);
[t2_counts_gate_momentum,t_edges_gate_momentum]=histcounts([t2_gate_momentum],tedges_gate_momentum);

%t_hist=[bincent_t' t_counts'];
%dlmwrite('t_hist.csv',t_hist);

subplot(2,2,2)

plot(bincent_t_gate_momentum,t1_counts_gate_momentum,'LineWidth',1,'Color','r')  
hold on
plot(bincent_t_gate_momentum,t2_counts_gate_momentum,'LineWidth',1,'Color','b') 
% grid on
xlabel('TOF/\mus','FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
xlim([min(tedges_gate_momentum) max(tedges_gate_momentum)]);
ylim([0 max(t_counts_gate_momentum)*1.1]);
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
set(gca,'FontSize',10)
%pbaspect([1 1 1]);

%plot xtof
subplot(2,2,3)
count_xtof_gate_momentum = histcounts2([t1_gate_momentum;t2_gate_momentum],[x1_gate_momentum;x2_gate_momentum],tedges_gate_momentum,xedges_gate_momentum);
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


%plot ytof
subplot(2,2,4)
count_ytof_gate_momentum = histcounts2([t1_gate_momentum;t2_gate_momentum],[y1_gate_momentum;y2_gate_momentum],tedges_gate_momentum,yedges_gate_momentum);
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

%%

close all;

tof_range = 20;
x_range = 3;
y_range = 3;

t1_edges_gate_momentum = min(t1_gate_momentum)-tof_range:1:max(t1_gate_momentum)+tof_range;
t2_edges_gate_momentum = min(t2_gate_momentum)-tof_range:1:max(t2_gate_momentum)+tof_range;

x1_edges_gate_momentum = min(x1_gate_momentum)-x_range:1:max(x1_gate_momentum)+x_range;
x2_edges_gate_momentum = min(x2_gate_momentum)-x_range:1:max(x2_gate_momentum)+x_range;

y1_edges_gate_momentum = min(y1_gate_momentum)-y_range:1:max(y1_gate_momentum)+y_range;
y2_edges_gate_momentum = min(y2_gate_momentum)-y_range:1:max(y2_gate_momentum)+y_range;

  
figure
count_x1_y1 = histcounts2(x1_gate_momentum,y1_gate_momentum,x1_edges_gate_momentum,y1_edges_gate_momentum);
count_x1_y1_mod = max(count_x1_y1,1); 

subplot(2,3,1)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(x1_gate_momentum,y1_gate_momentum, ((count_x1_y1_mod')));
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
% pos1 = get(sp_hand1, 'Position') % gives the position of current sub-plot
% new_pos1 = pos1 +[0 0 0 0.05]
% set(sp_hand1, 'Position',new_pos1 ) % set new position of current sub - plot


count_t1_x1 = histcounts2(t1_gate_momentum,x1_gate_momentum,t1_edges_gate_momentum,x1_edges_gate_momentum);
count_t1_x1_mod = max(count_t1_x1,1); 

subplot(2,3,2)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(t1_gate_momentum,x1_gate_momentum, ((count_t1_x1_mod')));
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


count_t1_y1 = histcounts2(t1_gate_momentum,y1_gate_momentum,t1_edges_gate_momentum,y1_edges_gate_momentum);
count_t1_y1_mod = max(count_t1_y1,1); 

subplot(2,3,3)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(t1_gate_momentum,y1_gate_momentum, ((count_t1_y1_mod')));
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



count_x2_y2 = histcounts2(x2_gate_momentum,y2_gate_momentum,x2_edges_gate_momentum,y2_edges_gate_momentum);
count_x2_y2_mod = max(count_x2_y2,1); 

subplot(2,3,4)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(x2_gate_momentum,y2_gate_momentum, ((count_x2_y2_mod')));
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
% pos1 = get(sp_hand1, 'Position') % gives the position of current sub-plot
% new_pos1 = pos1 +[0 0 0 0.05]
% set(sp_hand1, 'Position',new_pos1 ) % set new position of current sub - plot


count_t2_x2 = histcounts2(t2_gate_momentum,x2_gate_momentum,t2_edges_gate_momentum,x2_edges_gate_momentum);
count_t2_x2_mod = max(count_t2_x2,1); 

subplot(2,3,5)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(t2_gate_momentum,x2_gate_momentum, ((count_t2_x2_mod')));
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


count_t2_y2 = histcounts2(t2_gate_momentum,y2_gate_momentum,t2_edges_gate_momentum,y2_edges_gate_momentum);
count_t2_y2_mod = max(count_t2_y2,1); 

subplot(2,3,6)
myColorMap = jet;
myColorMap(1,:) = 1;
imagesc(t2_gate_momentum,y2_gate_momentum, ((count_t2_y2_mod')));
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


%%

% ke and angular distribution
    
    calfac=0.9349;
%     calfac=1;
    KE1 = calfac*(px1.*px1 + py1.*py1 + pz1.*pz1)/(2*frag_m(1))*hatoev;
    KE2 = calfac*(px2.*px2 + py2.*py2 + pz2.*pz2)/(2*frag_m(2))*hatoev;
    KER=KE1+KE2;

   binsize_ke=0.1; %eV   
   edges_ke=[0:binsize_ke:35];
   i_ke=1:length(edges_ke)-1;
   bincent_ke=[];
   bincent_ke(i_ke)=(edges_ke(i_ke)+edges_ke(i_ke+1))/2;
  
   % ke_raw=[k1 k2 k3 k];
   % dlmwrite('ke_raw.csv',ke_raw);
    
   [KE1_counts,KE1_edges]=histcounts(KE1,edges_ke);
   [KE2_counts,KE2_edges]=histcounts(KE2,edges_ke);
   [KER_counts,KER_edges]=histcounts(KER,edges_ke);
  
   
   ke_hist=[bincent_ke' KE1_counts' KE2_counts' KER_counts'];
   dlmwrite('ke_hist.csv',ke_hist);
  
   figure 
   close all;
%    plot(bincent_ke,KE1_counts,'b', bincent_ke,KE2_counts,'r',bincent_ke,KER_counts,'k','LineWidth',2)   
 plot(bincent_ke,KER_counts,'k','LineWidth',2)  
   set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30 );
   grid on
   xlim([0  edges_ke(end)]);
%    pbaspect([1 1 1]);
   xlabel('KE /eV', 'FontWeight', 'normal','FontName', 'Arial');
   ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
   %legend({'First Hit','Second Hit','KER'}, 'FontSize', 20, 'FontWeight', 'bold')
   %legend({frag_m_z_str(1),frag_m_z_str(2),'KER'}, 'FontSize', 30, 'FontWeight', 'normal')
legend({'KER'}, 'FontSize', 30, 'FontWeight', 'normal')
    %sb
%     ke=[KE1_bin'  KE1_counts KE2_bin'  KE2_counts KE3_bin'  KE3_counts  KER_bin'  KER_counts];
%     dlmwrite('ke.csv',ke);
    %sb

clear 
clc
close all

% fileID = fopen('Bromoform_140eV_4_DAn.bin');
fileID = fopen('Bromoform_140eV_DAn.00000.bin');
%fileID = fopen('Bromocyclopentane_140eV_DAn.00008.bin');
% fileID = fopen('air_140eV_Nov21_8am_Before_Bromoform_DAn.bin');
data_ini = fread(fileID,'int32');
% data_ini=data_ini(1:1000,:);
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
combos = nchoosek(1:nhits,3);
length(combos(:,1))
t1=[];x1=[];y1=[];
t2=[];x2=[];y2=[];
t3=[];x3=[];y3=[];

parfor i=1:length(combos(:,1))
    combos_1 = combos(:,1);
    combos_2 = combos(:,2);
    combos_3 = combos(:,3);
%     combos_4 = combos(:,4);
%     combos_5 = combos(:,5);
%     t1=[t1,eval(['t_',num2str(combos_1(i))])];
%     t2=[t2,eval(['t_',num2str(combos_2(i))])];
    
     t1=[t1;t_i(:,combos_1(i))]; x1=[x1;x_i(:,combos_1(i))]; y1=[y1;y_i(:,combos_1(i))];
     t2=[t2;t_i(:,combos_2(i))]; x2=[x2;x_i(:,combos_2(i))]; y2=[y2;y_i(:,combos_2(i))];
     t3=[t3;t_i(:,combos_3(i))]; x3=[x3;x_i(:,combos_3(i))]; y3=[y3;y_i(:,combos_3(i))];

end
toc
%%
measurement = [t1 x1 y1 t2 x2 y2 t3 x3 y3];
%  dlmwrite('measurement.csv',measurement);
save('measurement', 'measurement', '-v7.3')

%%
clear 
clc
close all
load('measurement');
t1=measurement(:,1); x1=measurement(:,2); y1=measurement(:,3);
t2=measurement(:,4); x2=measurement(:,5); y2=measurement(:,6);
t3=measurement(:,7); x3=measurement(:,8); y3=measurement(:,9);

%%
close all
Xedges = 0:5:max(t1);
Yedges = 0:5:max(t2+t3);

%t1=data_i(:,4)./1000;
%t2=data_i(:,7)./1000;
count = histcounts2(t1,t2+t3,Xedges,Yedges);

%sum(sum(count))
% subplot(1,2,1)
myColorMap = jet;
% myColorMap = flipud(hot);
myColorMap(1,:) = 1;
imagesc( Xedges,Yedges, (count'));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
% axis([0 12000 0 12000])
%axis equal
pbaspect([1 1 1]);
hold on;
set(gca,'FontSize',25)
set(gca,'colorscale','log');
% scatter([sct_peak_tof_1], [sct_peak_tof_2],300,'x','g','LineWidth',3)

xlabel('TOF_1 /ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('TOF_2 + TOF_3 /ns', 'FontWeight', 'normal','FontName', 'Arial');

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

% frag={ {'Br79',1,1} {'Br79',1,1}   {'C',1,'H',1, 'Br79',1,1}  ...
%        {'Br81',1,1} {'Br81',1,1}   {'C',1,'H',1, 'Br81',1,1} };
%    
   frag={  {'Br79',1,1} {'Br79',1,1}   {'C',1,'H',1, 'Br79',1,1} };

   

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
for i=1:length(frag_m_z)
    xline(  k0*(frag_m_z(i))^0.5+t0, '--b', frag_m_z_str(i),'fontweight','bold','fontsize',30 )
end


ii=2;
jj=3;
yline(  k0*( (frag_m_z(ii) ))^0.5+t0 + k0*( (frag_m_z(jj) ))^0.5+t0, '--b', strcat(frag_m_z_str(ii), '+', frag_m_z_str(jj)),'fontweight','bold','fontsize',30 );

 
%%
addpath('D:\amo\als\2021\CHBr3\three_body\three_body_fn');

t1_bin = [3850 4150]; t23_bin = [8250 8550];

sct_peak_tof_1=3987; % in ns %79
sct_peak_tof_23=8321;  % in ns 

close all
clc

[j_gate,t1_gate,t2_gate,t3_gate,t12_gate,t23_gate,Xedges,Yedges,count,t23]=j_gate_tof(t1,t2,t3,t1_bin,t23_bin);


myplot(Xedges,Yedges,t1_gate,t23_gate,t0,k0,frag,sct_peak_tof_1,sct_peak_tof_23)


%% plotting rotated triipico
close all
t_sum=t1_gate*1+t23_gate;
t_diff=t23_gate-t1_gate*1;

tripico_raw_rotated(t_sum, t_diff, sct_peak_tof_1,sct_peak_tof_23)

%% gate on roated tripico
t_diff_bin = [4250 4465]; t_sum_bin = [12290 12325];
[j_rot_gate]=tripico_raw_rotated_gated(t_sum,t_diff,t_diff_bin, t_sum_bin, sct_peak_tof_1,sct_peak_tof_23);

%% tripico rotation retrieved

 [t1_gate_rot,t2_gate_rot,t3_gate_rot,t23_gate_rot,...
            x1_gate_rot,x2_gate_rot,x3_gate_rot,...
             y1_gate_rot,y2_gate_rot,y3_gate_rot]...
            =tripico_rotation_retrieved(t1,t2,t3,x1,x2,x3,y1,y2,y3,t23,j_gate,j_rot_gate,sct_peak_tof_1,sct_peak_tof_23);

%% Analysis Constants
 clc
close all
mmns_to_au_conv=0.4571*10^(-3); %mm/us to a.u. conversion
run_name = '';
hatoev=  27.211396; % taken from ksu webpage/ twice of ionization potential of H



x01 =1.5;        mx1 =2.160968127030588604e-01;   cx1=-9.272382920574947021e-03;
y01 = 0.9;       my1 =mx1;                        cy1=cx1;
t01 = -288.35;   mt1 =-6.764717828664892907e-02;  ct1=2.893873970012056134e+02;


x02 =1.5;        mx2 =2.160968127030588604e-01;   cx2=-9.272382920574947021e-03;
y02 = 0.9;       my2 =mx2;                        cy2=cx2;
t02 = t01;       mt2 =-6.764717828664892907e-02;  ct2=2.893873970012056134e+02;

x03 = 1.3;       mx3 =2.002244236516522724e-01;   cx3=-8.591126955530558468e-03;
y03 =0.9;        my3 =mx3;                        cy3=cx3;
t03 = t01;       mt3 =-5.807525138290989958e-02;  ct3=2.681341306391510102e+02;

binsize_px=2;  binsize_py=2; binsize_pz=2; binsize_p=2;

%%
momentum_gui


%% event 1 and 2 swapping better way of swapping

[px1,px2,px3,py1,py2,py3,pz1,pz2,pz3]=eventswap_momenta_12(px1,px2,px3,py1,py2,py3,pz1,pz2,pz3);


%% plotting ungated momenta
binsize_px=1    ;   binsize_py=4    ;    binsize_pz=0.25   ;    binsize_p=0.5    ; %p=psum vin

px_plot_rng=[-600 600]; px_sum_plot_rng = [-50 50];
py_plot_rng=[-600 600]; py_sum_plot_rng = [-50 50];
pz_plot_rng=[-600 600]; pz_sum_plot_rng = [-50 50];

ungated_momenta(px1,px2,px3,py1,py2,py3,pz1,pz2,pz3,binsize_px,binsize_py,binsize_pz,binsize_p,p1,p2,p3,...
                            px_plot_rng,py_plot_rng,pz_plot_rng,px_sum_plot_rng,py_sum_plot_rng,pz_sum_plot_rng);


%%
close all
px1_rng = [-600.0 600.0]; px2_rng = [-600.0 600.0];  px3_rng = [-600.0 600.0];   px_sum_rng = [-60 60];
py1_rng = [-600.0 600.0]; py2_rng = [-600.0 600.0];  py3_rng = [-600.0 600.0];   py_sum_rng = [-60 60];
pz1_rng = [-600.0 600.0]; pz2_rng = [-600.0 600.0];  pz3_rng = [-600.0 600.0];   pz_sum_rng = [-15.0 15.0];

binsize_px=4       ;  binsize_py=4     ; binsize_pz=4     ; binsize_p=4  ;

[j_momentum_gate,px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,p1,p2,p3]...
            =gate_1d_momenta(px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,p1,p2,p3,...
                    px1_rng, px2_rng, px3_rng, px_sum_rng,py1_rng, py2_rng, py3_rng, py_sum_rng,...
                        pz1_rng, pz2_rng, pz3_rng, pz_sum_rng,binsize_px,binsize_py,binsize_pz,binsize_p);
                    
                    
%% momentum correlation plot

close all;
binsize_px=2   ; binsize_py=2  ; binsize_pz=2  ;    binsize_p=2 ;

momentum_correlation(px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,binsize_px,binsize_py,binsize_pz,binsize_p);                  

%% gate on momentunm correlation plot and replotting the momentum correaltion       use this gate in worst case not in normal situation
% % % % clc;
% % % % close all;
% % % % px1_px2_px3_dif_rng = [-1000.0 1000.0]; px1_px2_px3_sum_rng = [-1000.0 1000.0];  
% % % % py1_py2_py3_dif_rng = [-1000.0 1000.0]; py1_py2_py3_sum_rng = [-1000.0 1000.0]; 
% % % % pz1_pz2_pz3_dif_rng = [-1000.0 1000.0]; pz1_pz2_pz3_sum_rng = [-1000.0 1000.0]; 
% % % % 
% % % % binsize_px=2   ; binsize_py=2  ;  binsize_pz=2   ;   binsize_p=2  ;
% % % % 
% % % % [j_momentum_square_gate,px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,p1,p2,p3]...
% % % %              =gate_momentunm_correlation(px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,p1,p2,p3,...
% % % %              px1_px2_px3_dif_rng, py1_py2_py3_dif_rng, pz1_pz2_pz3_dif_rng,...
% % % %              px1_px2_px3_sum_rng, py1_py2_py3_sum_rng,pz1_pz2_pz3_sum_rng,...
% % % %                binsize_px,binsize_py,binsize_pz,binsize_p);
% % % % 
%% Replotting the 1d momenta again
% % % % 
% % % % binsize_px=1   ;   binsize_py=4   ;  binsize_pz=0.25   ; binsize_p=0.5  ;
% % % % 
% % % % px_plot_rng=[-1000 1000]; px_sum_plot_rng = [-1000 1000];
% % % % py_plot_rng=[-1000 1000]; py_sum_plot_rng = [-1000 1000];
% % % % pz_plot_rng=[-1000 1000]; pz_sum_plot_rng = [-1000 1000];
% % % % 
% % % % 
% % % % replot_1d_momenta(px1,px2,px3,py1,py2,py3,pz1,pz2,pz3,binsize_px,binsize_py,binsize_pz,binsize_p,p1,p2,p3,...
% % % %                                 px_plot_rng,py_plot_rng,pz_plot_rng,px_sum_plot_rng,py_sum_plot_rng,pz_sum_plot_rng)
         
%% plot momentum sphere
close all;

binsize_p_sphere=10;  p_plot_range=400;

momentum_sphere(px1,px2,px3,py1,py2,py3,pz1,pz2,pz3,binsize_p_sphere,p_plot_range);

%% plot survived events
close all;
[t1_gate_momentum,t2_gate_momentum,t3_gate_momentum,t23_gate_momentum,...
x1_gate_momentum,x2_gate_momentum,x3_gate_momentum,y1_gate_momentum,y2_gate_momentum,y3_gate_momentum]...
= survived_events(t1_gate_rot,t2_gate_rot,t3_gate_rot,t23_gate_rot,x1_gate_rot,x2_gate_rot,x3_gate_rot,...
                        y1_gate_rot, y2_gate_rot, y3_gate_rot,j_momentum_gate,frag_m_z,frag_m_z_str,t0,k0);

%% additional gate on survived events can skip this unless necessary

clc
t1_rng = [3900  4100]; t2_rng = [3900  4100];  t3_rng = [4215 4450];   

 [j_momentum_gate_tof,t1_gate_momentum,t2_gate_momentum,t3_gate_momentum,t23_gate_momentum,...
x1_gate_momentum,x2_gate_momentum,x3_gate_momentum,y1_gate_momentum,y2_gate_momentum,y3_gate_momentum]...
= gate_survived_events(t1_gate_momentum,t2_gate_momentum,t3_gate_momentum,t23_gate_momentum,...
x1_gate_momentum,x2_gate_momentum,x3_gate_momentum,y1_gate_momentum,y2_gate_momentum,y3_gate_momentum,...
frag_m_z,frag_m_z_str,t0,k0,t1_rng,t2_rng,t3_rng);


%% plotting detector image and xtof and ytof
close all;
rcirc1=18; rcirc2=18; rcirc3=15;

detector_image(t1_gate_momentum,t2_gate_momentum,t3_gate_momentum,t23_gate_momentum,...
x1_gate_momentum,x2_gate_momentum,x3_gate_momentum,y1_gate_momentum,y2_gate_momentum,y3_gate_momentum,...
x01, y01, x02, y02, x03, y03,rcirc1,rcirc2,rcirc3);

%%
% ke and angular distribution
clc
px1 = px1(j_momentum_gate_tof); py1 = py1(j_momentum_gate_tof); pz1 = pz1(j_momentum_gate_tof); 
px2 = px2(j_momentum_gate_tof); py2 = py2(j_momentum_gate_tof); pz2 = pz2(j_momentum_gate_tof); 
px3 = px3(j_momentum_gate_tof); py3 = py3(j_momentum_gate_tof); pz3 = pz3(j_momentum_gate_tof); 


%%
% ke and angular distribution
clc
    KE1 = (px1.*px1 + py1.*py1 + pz1.*pz1)/(2*frag_m(1))*hatoev;
    KE2 = (px2.*px2 + py2.*py2 + pz2.*pz2)/(2*frag_m(2))*hatoev;
    KE3 = (px3.*px3 + py3.*py3 + pz3.*pz3)/(2*frag_m(3))*hatoev;
    KER=KE1+KE2+KE3;

   %binsize_ke=0.05; %eV  
   binsize_ke=0.2; %eV   
   edges_ke=[0:binsize_ke:20];
   i_ke=1:length(edges_ke)-1;
   bincent_ke=[];
   bincent_ke(i_ke)=(edges_ke(i_ke)+edges_ke(i_ke+1))/2;
  
   % ke_raw=[k1 k2 k3 k];
   % dlmwrite('ke_raw.csv',ke_raw);
    
   [KE1_counts,KE1_edges]=histcounts(KE1,edges_ke);
   [KE2_counts,KE2_edges]=histcounts(KE2,edges_ke);
   [KE3_counts,KE3_edges]=histcounts(KE3,edges_ke);
   [KER_counts,KER_edges]=histcounts(KER,edges_ke);
  
   
%    ke_hist=[bincent_ke' KE1_counts' KE2_counts' KE3_counts' KER_counts'];
%    dlmwrite('ke_hist.csv',ke_hist);
  
   figure 
   close all;
   plot(bincent_ke,KE1_counts,'b', bincent_ke,KE2_counts,'r', bincent_ke,KE3_counts,'g', bincent_ke,KER_counts,'k','LineWidth',2)   
   set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30 );
   grid on
   xlim([0 edges_ke(end)]);
%    pbaspect([1 1 1]);
   xlabel('KER /eV', 'FontWeight', 'normal','FontName', 'Arial');
   ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
   legend({frag_m_z_str(1),frag_m_z_str(2),frag_m_z_str(3),'KER'}, 'FontSize', 20, 'FontWeight', 'normal')


    %sb
%     ke=[KE1_bin'  KE1_counts KE2_bin'  KE2_counts KE3_bin'  KE3_counts  KER_bin'  KER_counts];
%     dlmwrite('ke.csv',ke);
    %sb
    

%%
%gate on KE
% this gate is not so important
close all;
 KE1_rng = [0 30]; KE2_rng = [0 30]; KE3_rng = [0 30]; KER_rng = [0 30];
%     KE1_rng = [0 5.1]; KE2_rng = [0 5.1]; KE3_rng = [0 25]; KER_rng = [0 25];%use this for gating
%    KE1_rng = [5.1 25]; KE2_rng = [0 25]; KE3_rng = [0 25]; KER_rng = [0 25];%use this for gating
j_KE1 = KE1 > KE1_rng(1) &  KE1 < KE1_rng(2); 
KE1_gated=KE1(j_KE1);
    
j_KE2 = KE2 > KE2_rng(1) &  KE2 < KE2_rng(2); 
KE2_gated=KE2(j_KE2);

j_KE3 = KE3 > KE3_rng(1) &  KE3 < KE3_rng(2); 
KE3_gated=KE3(j_KE3);

j_KER = KER > KER_rng(1) &  KER < KER_rng(2); 
KER_gated=KER(j_KER);

% this plot is not so important
edges_ke_gated= edges_ke;
bincent_ke_gated=bincent_ke;
    
[KE1_counts_gated,KE1_edges_gated]=histcounts(KE1_gated,edges_ke_gated);
[KE2_counts_gated,KE2_edges_gated]=histcounts(KE2_gated,edges_ke_gated);
[KE3_counts_gated,KE3_edges_gated]=histcounts(KE3_gated,edges_ke_gated);
[KER_counts_gated,KER_edges_gated]=histcounts(KER_gated,edges_ke_gated);
  
   
% ke_hist=[bincent_gated_ke' KE1_gated_counts' KE2_gated_counts' KE3_gated_counts' KER_gated_counts'];
% dlmwrite('ke_gated_hist.csv',ke_gated_hist);
  
figure 
close all;
plot(bincent_ke_gated,KE1_counts_gated,'b',bincent_ke_gated,KE2_counts_gated,'r',bincent_ke_gated,KE3_counts_gated,'g',bincent_ke_gated,KER_counts_gated,'k','LineWidth',2);   
    
grid on
xlim([0 20]);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30 );
xlabel('KE\_gated / eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
 legend({frag_m_z_str(1),frag_m_z_str(2),frag_m_z_str(3),'KER'}, 'FontSize', 20, 'FontWeight', 'normal')
    
 
%%
% KE all gated
% this gate is important
j_KE_all = j_KE1 & j_KE2 & j_KE3 & j_KER;

% KE1_gated_all=KE1(~j_KE_all);
% KE2_gated_all=KE2(~j_KE_all);
% KE3_gated_all=KE3(~j_KE_all);
% KER_gated_all=KER(~j_KE_all);

KE1_gated_all=KE1(j_KE_all);
KE2_gated_all=KE2(j_KE_all);
KE3_gated_all=KE3(j_KE_all);
KER_gated_all=KER(j_KE_all);

KE1=KE1_gated_all;
KE2=KE2_gated_all;
KE3=KE3_gated_all;
KER=KER_gated_all;


% this plot is important
edges_ke_gated_all= edges_ke;
bincent_ke_gated_all=bincent_ke;
    
[KE1_counts_gated_all,KE1_edges_gated_all]=histcounts(KE1_gated_all,edges_ke_gated_all);
[KE2_counts_gated_all,KE2_edges_gated_all]=histcounts(KE2_gated_all,edges_ke_gated_all);
[KE3_counts_gated_all,KE3_edges_gated_all]=histcounts(KE3_gated_all,edges_ke_gated_all);
[KER_counts_gated_all,KER_edges_gated_all]=histcounts(KER_gated_all,edges_ke_gated_all);
  
   
% ke_gated_hist=[bincent_ke_gated_all' KE1_counts_gated_all' KE2_counts_gated_all' KE3_counts_gated_all' KER_counts_gated_all'];
% dlmwrite('ke_gated_hist.csv',ke_gated_hist);
  
figure 
close all;
plot(bincent_ke_gated_all,KE1_counts_gated_all,'b',bincent_ke_gated_all,KE2_counts_gated_all,'r',bincent_ke_gated_all,KE3_counts_gated_all,'g',bincent_ke_gated_all,KER_counts_gated_all,'K','LineWidth',2);   
    
grid on
xlim([0 edges_ke(end)]);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30 );
xlabel('KE\_gated\_all /eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
 legend({frag_m_z_str(1),frag_m_z_str(2),frag_m_z_str(3),'KER'}, 'FontSize', 20, 'FontWeight', 'normal')
%%
% if do not want to gate, then just not run these three lines for gating
% px1=px1(~j_KE_all); py1=py1(~j_KE_all);pz1=pz1(~j_KE_all);
% px2=px2(~j_KE_all); py2=py2(~j_KE_all);pz2=pz2(~j_KE_all);
% px3=px3(~j_KE_all); py3=py3(~j_KE_all);pz3=pz3(~j_KE_all);

px1=px1(j_KE_all); py1=py1(j_KE_all);pz1=pz1(j_KE_all);
px2=px2(j_KE_all); py2=py2(j_KE_all);pz2=pz2(j_KE_all);
px3=px3(j_KE_all); py3=py3(j_KE_all);pz3=pz3(j_KE_all);

% p_gated=[px1 py1 pz1 px2 py2 pz2 px3 py3 pz3 ];
% dlmwrite('p_gated.csv',p_gated);

%%
%Dalitz plot conventional
%Dalitz plot in a different way
close all;

KE1_n=KE1./KER;
KE2_n=KE2./KER;
KE3_n=KE3./KER;

% %p_sq_sum = 2*frag_m(1)*KE1 + 2*frag_m(2)*KE2+2*frag_m(3)*KE3;
% % KE1_n=2*frag_m(1)*KE1./p_sq_sum;
% % KE2_n=2*frag_m(2)*KE2./p_sq_sum;
% % KE3_n=2*frag_m(3)*KE3./p_sq_sum;

daliz_x = (KE2_n-KE3_n)./sqrt(3);
daliz_y =  KE1_n - 1/3 ;

%making 2d histogram
% Xedges = min(daliz_x(:,1)):(max(daliz_x(:,1)) - min(daliz_x(:,1)) )/500:max(daliz_x(:,1));
% Yedges = min(daliz_y(:,1)):(max(daliz_y(:,1)) - min(daliz_y(:,1)) )/500:max(daliz_y(:,1));

Xedges_d =-0.4:0.8*2/500:0.4;
Yedges_d =Xedges_d;
% count = histcounts2(daliz_x(:,1),daliz_y(:,1),Xedges,Yedges);
count_d = histcounts2(daliz_x(:,1),daliz_y(:,1),Xedges_d,Yedges_d);

count_d_mod = max(count_d,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges_d,Yedges_d, ((count_d')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
%axis([-0.3 0.3 -0.3 0.3]);
axis([-0.4 0.4 -0.4 0.4]);
hold on;

% xtickformat('%.2f');
% ytickformat('%.2f');
% xticks([-0.3 -0.15 0.00 0.15 0.3]);
% yticks([-0.3 -0.15 0.00 0.15 0.3]);

set(gca,'FontSize',30)
%axis equal;
pbaspect([1 1 1]);

xlabel('(\epsilon_{Br^+}_{(2)}-\epsilon_{Br^+}_{(3)})/\surd3', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\epsilon_{CHBr^{2+}}_{(1)}-1/3', 'FontWeight', 'normal','FontName', 'Arial');
hold on;
set(gca,'FontSize',30)
set(gca,'colorscale','log');

% scatter(0.103750179,-0.082863107,300,'x','g','LineWidth',3) %concerted CHBr3
% scatter(0.099948691, -0.085178471,300,'x','w','LineWidth',3) %concerted CHBr3 TS Br migrated
% scatter(0.112864312,-0.106455783,300,'x','b','LineWidth',3) %concerted CHBr3 Br migrated


%%
%Dalitz plot in a different way 1
close all;

KE1_n=KE1./KER;
KE2_n=KE2./KER;
KE3_n=KE3./KER;

% % KE1_n=2*frag_m(1)*KE1./p_sq_sum;
% % KE2_n=2*frag_m(2)*KE2./p_sq_sum;
% % KE3_n=2*frag_m(3)*KE3./p_sq_sum;

daliz_x = (KE1_n-KE2_n)./sqrt(3);
daliz_y =  KE3_n - 1/3 ;

count_d = histcounts2(daliz_x(:,1),daliz_y(:,1),Xedges_d,Yedges_d);
mycolor=[0 1 1];
mycoloredge=[0 0 1];

count_d_mod = max(count_d,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges_d,Yedges_d, ((count_d_mod')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
%axis([-0.3 0.3 -0.3 0.3]);
axis([-0.4 0.4 -0.4 0.4]);
hold on;

% xtickformat('%.2f');
% ytickformat('%.2f');
% xticks([-0.3 -0.15 0.00 0.15 0.3]);
% yticks([-0.3 -0.15 0.00 0.15 0.3]);
xticks([-0.4 :0.2:0.4])
yticks(xticks)
%  set(gca,'YTick',[])
% set(gca,'Yticklabel',[]) 
 set(gca,'FontSize',30)
%axis equal;
pbaspect([1 1 1]);

xlabel('(\epsilon_{Br^+}_{(1)}-\epsilon_{Br^+}_{(2)})/\surd3', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\epsilon_{CHBr^+}_{(3)}-1/3', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');

%scatter(-1.29587E-05,-0.027944246,600,'p','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3) %concerted
%caxis([1 27])
scatter(-1.29587E-05,-0.027944246,500,'x','g','LineWidth',5) %concerted for CHBr3
% scatter(0.041616041,-0.028129597,300,'x','b','LineWidth',3) %concerted for Br migrated
%%
%Dalitz plot in a different way 2
close all;
KE1_n=KE1./KER;
KE2_n=KE2./KER;
KE3_n=KE3./KER;

daliz_x = (KE1_n-KE3_n)./sqrt(3);
daliz_y =  KE2_n - 1/3 ;

count_d = histcounts2(daliz_x(:,1),daliz_y(:,1),Xedges_d,Yedges_d);

count_d_mod = max(count_d,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges_d,Yedges_d, ((count_d_mod')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
%axis([-0.3 0.3 -0.3 0.3]);
axis([-0.4 0.4 -0.4 0.4]);
xticks([-0.4 :0.2:0.4])
yticks(xticks)
%  set(gca,'YTick',[])
hold on;

% xtickformat('%.2f');
% ytickformat('%.2f');
% xticks([-0.3 -0.15 0.00 0.15 0.3]);
% yticks([-0.3 -0.15 0.00 0.15 0.3]);

set(gca,'FontSize',30)
%axis equal;
pbaspect([1 1 1]);

xlabel('(\epsilon_{Br^+}_{(1)}-\epsilon_{Br^+}_{(2)})/\surd3', 'FontWeight', 'normal','FontName', 'Arial');
% ylabel('\epsilon_{CHBr^+}_{(3)}-1/3', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',30)
set(gca,'colorscale','log');
%scatter(0.024193948,0.013983346,600,'p','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3) %concerted

scatter(0.024193948,0.013983346,600,'x','g','LineWidth',3) %concerted

caxis([1 30])

%%


p1_m=sqrt(px1.*px1 + py1.*py1 + pz1.*pz1);
p2_m=sqrt(px2.*px2 + py2.*py2 + pz2.*pz2);
p3_m=sqrt(px3.*px3 + py3.*py3 + pz3.*pz3);


p1_dot_p2 = (px1.*px2 + py1.*py2 + pz1.*pz2);
p1_dot_p3 = (px1.*px3 + py1.*py3 + pz1.*pz3); 
p2_dot_p3 = (px2.*px3 + py2.*py3 + pz2.*pz3); 

cos_theta_12=p1_dot_p2./(p1_m.*p2_m);
cos_theta_13=p1_dot_p3./(p1_m.*p3_m);
cos_theta_23=p2_dot_p3./(p2_m.*p3_m);

close all;
acute_theta_12=(pi - acos(cos_theta_12))*180/pi;
acute_theta_13=(pi - acos(cos_theta_13))*180/pi;
acute_theta_23=(pi - acos(cos_theta_23))*180/pi;

    
%%
%newton_diagram with respect to 1st hit
close all
p2_m_np1=p2_m./p1_m;
p3_m_np1=p3_m./p1_m;

p2_m_np1_x = -(p2_m_np1).*cos(acute_theta_12.*pi/180); %momentum coordinate x
p2_m_np1_y =  (p2_m_np1).*sin(acute_theta_12.*pi/180); %momentum coordinate y

p3_m_np1_x = -(p3_m_np1).*cos(acute_theta_13.*pi/180);
p3_m_np1_y = -(p3_m_np1).*sin(acute_theta_13.*pi/180);

p23_x=[p2_m_np1_x;p3_m_np1_x];
p23_y=[p2_m_np1_y;p3_m_np1_y];

%making 2d histogram
% Xedges = min(p23_x(:,1)):(max(p23_x(:,1)) - min(p23_x(:,1)) )/5000:max(p23_x(:,1));
% Yedges = min(p23_y(:,1)):(max(p23_y(:,1)) - min(p23_y(:,1)) )/5000:max(p23_y(:,1));

Xedges = -20: 2*40/5000: 20;
Yedges = Yedges ;

count = histcounts2(p23_x(:,1),p23_y(:,1),Xedges,Yedges);
count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([-2 2 -2 2])
xticks([-2:1:2])
% yticks([-2:1:2])
%  set(gca,'YTick',[])
%  set(gca,'Yticklabel',[]) 
hold on;
quiver(0,0,1.,0,'-r','LineWidth',2,'MaxHeadSize',3);
set(gca,'FontSize',30)
axis equal;
pbaspect([1 1 1]);

xlabel('rel. P_x', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('rel. P_y', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',50)
set(gca,'colorscale','log');
% 
% annotation('textbox',[0.51 0.61 0 0.04],'EdgeColor','w','String', strcat(frag_m_z_str(1),'(',num2str(1),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.322 0.84 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(2),'(',num2str(2),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.322 0.24 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(3),'(',num2str(3),')'),'FontSize',20,'FontWeight','normal','Color','k');

mycolor=[0 1 1];
mycoloredge=[0 0 1];
% caxis([ 1 120]);

% markers are acccording to computed geometry numbering
%scatter([-0.4440 -0.5560],[0.7841 -0.7841] ,600,'o','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3); %isomer_iso_1 ref - Br(1)   Br(2) on top  CHBr(3) at bottom
%scatter([-0.5468 -0.4531],[0.9657 -0.9656] ,600,'s','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3); %isomer_iso_1  ref - Br(2)   Br(1) on top  CHBr(3) at bottom

% scatter([-0.5558 -0.4442 ],[0.7560 -0.7560 ] ,600,'v','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3); %isomer_iso_2  ref - Br(1)   Br(3) on top  CHBr(2) at bottom
% scatter([-0.6313 -0.3687],[0.8587 -0.8587 ] ,600,'d','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3); %isomer_iso_2  ref - Br(3)   Br(1) on top  CHBr(2) at bottom


% markers are acccording to four and five body numbering
%  scatter([-0.6313 -0.3687],[0.8587 -0.8587 ] ,600,'d','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3); %isomer_iso_2 ref - Br(1)   Br(3) on top  CHBr(2) at bottom
% scatter([-0.5558 -0.4442 ],[0.7560 -0.7560 ] ,600,'v','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3); %isomer_iso_2  ref - Br(3)   Br(1) on top  CHBr(2) at bottom

%scatter([-0.4440 -0.5560],[0.7841 -0.7841] ,600,'o','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3); %isomer_iso_1 ref - Br(3)   Br(2) on top  CHBr(1) at bottom
% scatter([-0.5468 -0.4531],[0.9657 -0.9656] ,600,'s','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3); %isomer_iso_1  ref - Br(2)   Br(3) on top  CHBr(1) at bottom

% scatter([-0.489675374 -0.510307908 ],[0.871919914 -0.871939817] ,600,'p','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3) %concerted for CHBr3

% scatter(-0.489675374,0.871919914,300,'x','g','LineWidth',3) %concerted for CHBr3
% scatter(-0.510307908,-0.871939817,300,'x','g','LineWidth',3) %concerted for CHBr3

% scatter(-0.457415952, 0.684028189,300,'x','b','LineWidth',3) %sequential anti-clockwise
% scatter(-0.542578369,-0.684028187,300,'x','b','LineWidth',3) %sequential anti-clockwise
%%
scatter([-0.489675374 -0.510307908 ],[0.871919914 -0.871939817],300,'x','g','LineWidth',3) %concerted
p_wrt1_seq=dlmread('p_wrt1_seq.csv');
scatter(p_wrt1_seq(:,1),p_wrt1_seq(:,2),'.','b','LineWidth',3)

%%
% x_axis = [0.5 0.5 0.0327380944398187 0.0353319050476945] ;
% upper =  [0.25 0.75 0.0327380944398187 0.0353319050476945]; 
% lower = [0.25 0.4 0.0327380944398187 0.0353319050476945];
% annotation('textbox',x_axis,'String','^{81}Br^+','FontSize',30,'FontWeight','normal','Color','k','FitBoxToText','on');
% annotation('textbox',upper,'String','^{81}Br^+','FontSize',30,'FontWeight','norma','Color','k','FitBoxToText','on');
% annotation('textbox',lower,'String','CH^{81}Br^+','FontSize',30,'FontWeight','norma','Color','k','FitBoxToText','on');
%%
% 
% %newton_diagram with respect to 2nd hit
 
p3_m_np2=p3_m./p2_m;
p1_m_np2=p1_m./p2_m;

p3_m_np2_x = -(p3_m_np2).*cos(acute_theta_23.*pi/180); %momentum coordinate x
p3_m_np2_y =  (p3_m_np2).*sin(acute_theta_23.*pi/180); %momentum coordinate y

p1_m_np2_x = -(p1_m_np2).*cos(acute_theta_12.*pi/180);
p1_m_np2_y = -(p1_m_np2).*sin(acute_theta_12.*pi/180);


p31_x=[p3_m_np2_x;p1_m_np2_x];
p31_y=[p3_m_np2_y;p1_m_np2_y];


%making 2d histogram
% Xedges = min(p31_x(:,1)):(max(p31_x(:,1)) - min(p31_x(:,1)) )/5000:max(p31_x(:,1));
% Yedges = min(p31_y(:,1)):(max(p31_y(:,1)) - min(p31_y(:,1)) )/5000:max(p31_y(:,1));

Xedges = -20: 40/5000: 20;
Yedges = -20: 40/5000: 20;

count = histcounts2(p31_x(:,1),p31_y(:,1),Xedges,Yedges);

count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([-2 2 -2 2])
% xticks([-2:1:2])
% yticks([-2:1:2])
hold on;
quiver(0,0,1.,0,'-r','LineWidth',2,'MaxHeadSize',3);
set(gca,'FontSize',50)
axis equal;
pbaspect([1 1 1]);

xlabel('rel. P_x', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('rel. P_y', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',30)
set(gca,'colorscale','log');


% scatter(-0.489662492,-0.871896976,300,'x','g','LineWidth',3) %concerted
% scatter(-0.510363048,0.871892145,300,'x','g','LineWidth',3) %concerted

scatter(-0.489662492,-0.871896976,300,'x','g','LineWidth',3) %concerted
scatter(-0.510363048,0.871892145,300,'x','g','LineWidth',3) %concerted


% annotation('textbox',[0.322 0.24 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(1),'(',num2str(1),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.51 0.61 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(2),'(',num2str(2),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.322 0.84 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(3),'(',num2str(3),')'),'FontSize',20,'FontWeight','normal','Color','k');

% scatter(-0.324476149,1.010190729,300,'x','b','LineWidth',3) %sequential anti-clockwise
% scatter(-0.675527685,-1.010196468,300,'x','b','LineWidth',3) %sequential anti-clockwise


%%
% %newton_diagram with respect to 3rd hit
close all
p1_m_np3=p1_m./p3_m;
p2_m_np3=p2_m./p3_m;


p1_m_np3_x = -(p1_m_np3).*cos(acute_theta_13.*pi/180); %momentum coordinate x
p1_m_np3_y =  (p1_m_np3).*sin(acute_theta_13.*pi/180); %momentum coordinate y

p2_m_np3_x = -(p2_m_np3).*cos(acute_theta_23.*pi/180);
p2_m_np3_y = -(p2_m_np3).*sin(acute_theta_23.*pi/180);


p12_x=[-100;p1_m_np3_x;p2_m_np3_x;100];
p12_y=[-100;p1_m_np3_y;p2_m_np3_y;100];




%making 2d histogram
% Xedges = min(p12_x(:,1)):(max(p12_x(:,1)) - min(p12_x(:,1)) )/25000:max(p12_x(:,1));
% Yedges = min(p12_y(:,1)):(max(p12_y(:,1)) - min(p12_y(:,1)) )/25000:max(p12_y(:,1));

Xedges = -10: 200/25000: 10;
Yedges = -10: 200/25000: 10;

count = histcounts2(p12_x(:,1),p12_y(:,1),Xedges,Yedges);
count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([-2 2 -2 2])
xticks([-2:1:2])
yticks([-2:1:2])
hold on;
quiver(0,0,1.,0,'-r','LineWidth',2,'MaxHeadSize',3);
set(gca,'FontSize',30)
axis equal;
pbaspect([1 1 1]);

xlabel('rel. P_x', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('rel. P_y', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',50)
set(gca,'colorscale','log');

scatter([-0.49996209 -0.500029267 ],[0.85426239 -0.854238158] ,600,'p','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3) %concerted for CHBr3

%  scatter(-0.49996209,0.85426239,300,'x','g','LineWidth',3) %concerted
%  scatter(-0.500029267,-0.854238158,300,'x','g','LineWidth',3) %concerted

% annotation('textbox',[0.322 0.86 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(1),'(',num2str(1),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.322 0.24 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(2),'(',num2str(2),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.51 0.61 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(3),'(',num2str(3),')'),'FontSize',20,'FontWeight','normal','Color','k');
% scatter(-0.711778096,0.897338171,300,'x','b','LineWidth',3) %sequential anti-clockwise
% scatter(-0.288225948,-0.897333076,300,'x','b','LineWidth',3) %sequential anti-clockwise
%%
%Unnormalized Newton diagram
%newton_diagram with respect to 1st hit
close all;
p2_m_np1=p2_m;
p3_m_np1=p3_m;

p2_m_np1_x = -(p2_m_np1).*cos(acute_theta_12.*pi/180); %momentum coordinate x
p2_m_np1_y =  (p2_m_np1).*sin(acute_theta_12.*pi/180); %momentum coordinate y

p3_m_np1_x = -(p3_m_np1).*cos(acute_theta_13.*pi/180);
p3_m_np1_y = -(p3_m_np1).*sin(acute_theta_13.*pi/180);

p23_x=[p1_m; p2_m_np1_x;p3_m_np1_x];
p23_y=[p1_m.*0; p2_m_np1_y;p3_m_np1_y]; % because no y component

%making 2d histogram
% Xedges = min(p23_x(:,1)):(max(p23_x(:,1)) - min(p23_x(:,1)) )/5000:max(p23_x(:,1));
% Yedges = min(p23_y(:,1)):(max(p23_y(:,1)) - min(p23_y(:,1)) )/5000:max(p23_y(:,1));

Xedges = -400: 40*20/500: 400;
Yedges = -400: 40*20/500: 400;

count = histcounts2(p23_x(:,1),p23_y(:,1),Xedges,Yedges);
count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count')));
  
xlabel('P_x', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('P_y', 'FontWeight', 'normal','FontName', 'Arial');


% annotation('textbox',[0.51 0.61 0 0.04],'EdgeColor','w','String', strcat(frag_m_z_str(1),'(',num2str(1),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.322 0.84 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(2),'(',num2str(2),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.322 0.24 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(3),'(',num2str(3),')'),'FontSize',20,'FontWeight','normal','Color','k');
% 
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([-400 400 -400 400])
xticks([-400:200:400])
yticks([-400:200:400])

set(gca,'FontSize',30)
axis equal;
pbaspect([1 1 1]);

xlabel('P_x / au', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('P_y / au', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
scatter([224.3295238 -109.8486435 -114.4771299 ],[0 195.5973791 -195.6018439] ,600,'p','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3) %concerted for CHBr3


% % % %scatter([-0.49996209 -0.500029267 ],[0.85426239 -0.854238158] ,600,'p','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3) %concerted for CHBr3
% % % 
% % % 
% % % 
% % % scatter(-109.8486435,195.5973791,300,'x','g','LineWidth',3) %concerted for CHBr3
% % % scatter(-114.4771299,-195.6018439,300,'x','g','LineWidth',3) %concerted for CHBr3
% % % % 
% % % % 
% % % % scatter(-101.5321169, 179.2951395,300,'x','b','LineWidth',3) %concerted for Br migration
% % % % scatter(-127.1247373, -179.2892877,300,'x','b','LineWidth',3) %concerted for Br migration
% % % 
% % % % scatter(-0.457415952, 0.684028189,300,'x','b','LineWidth',3) %sequential anti-clockwise
% % % % scatter(-0.542578369,-0.684028187,300,'x','b','LineWidth',3) %sequential anti-clockwise
%%
%plotting many scattered dots
%  scatter(p3_m_np1_x(:,1),p3_m_np1_y(:,1),'.','r','LineWidth',1) 
p_wrt1_seq=dlmread('p_wrt1_seq_abs.csv');
scatter(p_wrt1_seq(:,1),p_wrt1_seq(:,2),'.','b','LineWidth',3) 


%%
% 
% %newton_diagram with respect to 2nd hit
 
p3_m_np2=p3_m;
p1_m_np2=p1_m;

p3_m_np2_x = -(p3_m_np2).*cos(acute_theta_23.*pi/180); %momentum coordinate x
p3_m_np2_y =  (p3_m_np2).*sin(acute_theta_23.*pi/180); %momentum coordinate y

p1_m_np2_x = -(p1_m_np2).*cos(acute_theta_12.*pi/180);
p1_m_np2_y = -(p1_m_np2).*sin(acute_theta_12.*pi/180);


p31_x=[p3_m_np2_x;p1_m_np2_x];
p31_y=[p3_m_np2_y;p1_m_np2_y];


%making 2d histogram
% Xedges = min(p31_x(:,1)):(max(p31_x(:,1)) - min(p31_x(:,1)) )/5000:max(p31_x(:,1));
% Yedges = min(p31_y(:,1)):(max(p31_y(:,1)) - min(p31_y(:,1)) )/5000:max(p31_y(:,1));

Xedges = -400: 40*20/500: 400;
Yedges = -400: 40*20/500: 400;

count = histcounts2(p31_x(:,1),p31_y(:,1),Xedges,Yedges);

count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
% axis([-2 2 -2 2])
hold on;
quiver(0,0,mean(p2_m),0,'-r','LineWidth',2,'MaxHeadSize',3);
set(gca,'FontSize',30)
axis equal;
pbaspect([1 1 1]);

xlabel('P_x', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('P_y', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',30)
set(gca,'colorscale','log');

annotation('textbox',[0.322 0.24 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(1),'(',num2str(1),')'),'FontSize',20,'FontWeight','normal','Color','k');
annotation('textbox',[0.51 0.61 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(2),'(',num2str(2),')'),'FontSize',20,'FontWeight','normal','Color','k');
annotation('textbox',[0.322 0.84 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(3),'(',num2str(3),')'),'FontSize',20,'FontWeight','normal','Color','k');


% scatter(-109.8471985,-195.5948062,300,'x','g','LineWidth',3) %concerted for CHBr3
% scatter(-114.4910055,195.5937225,300,'x','g','LineWidth',3) %concerted for CHBr3
% 
% scatter(-112.6743583,-198.9711768,300,'x','b','LineWidth',3) %concerted for Br migration
% scatter(-93.3690847,198.9662322,300,'x','b','LineWidth',3) %concerted for Br migration
% 

% scatter(-0.324476149,1.010190729,300,'x','b','LineWidth',3) %sequential anti-clockwise
% scatter(-0.675527685,-1.010196468,300,'x','b','LineWidth',3) %sequential anti-clockwise


%%
% %newton_diagram with respect to 3rd hit

p1_m_np3=p1_m;
p2_m_np3=p2_m;


p1_m_np3_x = -(p1_m_np3).*cos(acute_theta_13.*pi/180); %momentum coordinate x
p1_m_np3_y =  (p1_m_np3).*sin(acute_theta_13.*pi/180); %momentum coordinate y

p2_m_np3_x = -(p2_m_np3).*cos(acute_theta_23.*pi/180);
p2_m_np3_y = -(p2_m_np3).*sin(acute_theta_23.*pi/180);


p12_x=[-100;p1_m_np3_x;p2_m_np3_x;100];
p12_y=[-100;p1_m_np3_y;p2_m_np3_y;100];




%making 2d histogram
% Xedges = min(p12_x(:,1)):(max(p12_x(:,1)) - min(p12_x(:,1)) )/25000:max(p12_x(:,1));
% Yedges = min(p12_y(:,1)):(max(p12_y(:,1)) - min(p12_y(:,1)) )/25000:max(p12_y(:,1));

Xedges = -400: 40*20/500: 400;
Yedges = -400: 40*20/500: 400;

count = histcounts2(p12_x(:,1),p12_y(:,1),Xedges,Yedges);
count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([-400 400 -400 400])
hold on;
quiver(0,0,mean(p3_m),0,'-r','LineWidth',2,'MaxHeadSize',3);
set(gca,'FontSize',30)
axis equal;
pbaspect([1 1 1]);

xlabel('P_x', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('P_y', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',30)
set(gca,'colorscale','log');

annotation('textbox',[0.322 0.86 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(1),'(',num2str(1),')'),'FontSize',20,'FontWeight','normal','Color','k');
annotation('textbox',[0.322 0.24 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(2),'(',num2str(2),')'),'FontSize',20,'FontWeight','normal','Color','k');
annotation('textbox',[0.51 0.61 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(3),'(',num2str(3),')'),'FontSize',20,'FontWeight','normal','Color','k');



scatter(-113.3107517,193.6089068,300,'x','g','LineWidth',3) %concerted
scatter(-113.3259767,-193.6034148,300,'x','g','LineWidth',3) %concerted

scatter(-132.2577671,186.5286124,300,'x','b','LineWidth',3) %concerted for Br migration
scatter(-87.53315193,-186.5300649,300,'x','b','LineWidth',3) %concerted for Br migration

% scatter(-0.711778096,0.897338171,300,'x','b','LineWidth',3) %sequential anti-clockwise
% scatter(-0.288225948,-0.897333076,300,'x','b','LineWidth',3) %sequential anti-clockwise

%%
%angular ditribution
costheta1z= pz1./p1_m;
costheta2z= pz2./p2_m;
costheta3z= pz3./p3_m;

% costheta1z=costheta1z(j_KE_all);
% costheta2z=costheta2z(j_KE_all);
% costheta3z=costheta3z(j_KE_all);

binsize_ke=0.2; %eV   
edges_ke=[0:binsize_ke:28];

i_ke=1:length(edges_ke)-1;
bincent_ke=[];
bincent_ke(i_ke)=(edges_ke(i_ke)+edges_ke(i_ke+1))/2;



binsize_costheta=0.02; %eV   
edges_costheta=[-1:binsize_costheta:1];

i_costheta=1:length(edges_costheta)-1;
bincent_costheta=[];
bincent_costheta(i_costheta)=(edges_costheta(i_costheta)+edges_costheta(i_costheta+1))/2;


%%
%plot angular distribution first hit
close all;
count_1z = histcounts2(KE1(:,1),costheta1z(:,1),edges_ke,edges_costheta);
count_mod_1z = max(count_1z,1); 
figure
subplot(2,3,1)
% count_mod = count; 
 myColorMap = jet;
% myColorMap=flipud(hot);
myColorMap(1,:) = 1;
imagesc( edges_ke,edges_costheta, ((count_mod_1z')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([0 15 -1 1])
set(gca,'FontSize',30)
set(gca,'colorscale','log');
pbaspect([1 1 1]);
xlabel('KE1_{Br(1)^+}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('cos\theta\_zaxis', 'FontWeight', 'normal','FontName', 'Arial');

% figure
% plot(bincent_ke,sum(count_mod_1z'));

 
subplot(2,3,4)
plot(bincent_costheta,sum(count_mod_1z),'r','LineWidth',2);
xlabel('cos\theta\_zaxis','FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
set(gca,'FontSize',30)
% xlim([-1 1]);
ylim([0 max(sum(count_mod_1z))*1.1]);
pbaspect([1 1 1]);


subplot(2,3,2)
%plot angular distribution second hit
count_2z = histcounts2(KE2(:,1),costheta2z(:,1),edges_ke,edges_costheta);
count_mod_2z = max(count_2z,1); 
% count_mod = count; 
% myColorMap = jet;
% figure;
imagesc( edges_ke,edges_costheta, ((count_mod_2z')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([0 15 -1 1])
set(gca,'FontSize',30)
set(gca,'colorscale','log');
pbaspect([1 1 1]);
xlabel('KE2_{Br(2)^+}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('cos\theta\_zaxis', 'FontWeight', 'normal','FontName', 'Arial');

subplot(2,3,5)
plot(bincent_costheta,sum(count_mod_2z),'r','LineWidth',2);
xlabel('cos\theta\_zaxis','FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
set(gca,'FontSize',30)
% xlim([-1 1]);
ylim([0 max(sum(count_mod_2z))*1.1]);
pbaspect([1 1 1]);

% close all;
% pos1 = [0.15 0.15 0.5 0.5];
% pos2 = [0.15 0.65 0.5 0.25];
% pos3 = [0.65 0.15 0.25 0.5];
% subplot('Position',pos1)
%plot angular distribution third hit
count_3z = histcounts2(KE3(:,1),costheta3z(:,1),edges_ke,edges_costheta);
count_mod_3z = max(count_3z,1); 
% figure
subplot(2,3,3)
imagesc( edges_ke,edges_costheta, ((count_mod_3z')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([0 15 -1 1])
set(gca,'FontSize',30)
set(gca,'colorscale','log');
pbaspect([1 1 1]);
xlabel('KE3_{CHBr^+}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('cos\theta\_zaxis', 'FontWeight', 'normal','FontName', 'Arial');

% figure
subplot(2,3,6)
plot(bincent_costheta,sum(count_mod_3z),'r','LineWidth',2);
xlabel('cos\theta\_zaxis','FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
set(gca,'FontSize',30)
% xlim([-1 1]);
ylim([0 max(sum(count_mod_3z))*1.1]);
pbaspect([1 1 1]);
% subplot('Position',pos2)
% hold on
% plot(bincent_ke,sum(count_mod_3z'));
% subplot('Position',pos3)
% hold on
% plot(bincent_costheta,sum(count_mod_3z));
%%
%plot angular distribution first and second hits average this does not mean
%anything
close all
figure;
subplot(1,2,1)
imagesc( edges_ke,edges_costheta, ((count_mod_1z'+count_mod_2z')./2));  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([0 15 -1 1])
set(gca,'FontSize',30)
set(gca,'colorscale','log');
pbaspect([1 1 1]);
xlabel('(KE1+KE2)\_avg/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('cos\theta\_zaxis', 'FontWeight', 'normal','FontName', 'Arial');

% figure
subplot(1,2,2)
plot(bincent_costheta,sum(count_mod_1z+count_mod_2z)/2,'r','LineWidth',2);
xlabel('cos\theta\_zaxis','FontWeight', 'normal','FontName', 'Arial');
ylabel('counts\_1+2\_avg', 'FontWeight', 'normal','FontName', 'Arial');
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
set(gca,'FontSize',30)
% xlim([-1 1]);
ylim([0 max(sum(count_mod_1z+count_mod_2z)/2)*1.1]);
pbaspect([1 1 1]);

%%
%plot KE1+KE2 vs ThetaP1P2  if does or does not run commment or uncomment  line 2248 % cos_theta_12=cos_theta_12(j_KE_all);
close all
figure;
% subplot(2,2,1)
KE12=KE1+KE2;
% cos_theta_12=cos_theta_12(j_KE_all);
count_12_rel = histcounts2(KE12(:,1),cos_theta_12(:,1),edges_ke,edges_costheta);
count_12_rel_mod = max(count_12_rel,1); 

myColorMap=flipud(hot);
% myColorMap(1,:) = 1;
imagesc( edges_ke,edges_costheta, ((count_12_rel')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([0 20 -1 1])
hold on;
set(gca,'FontSize',30)
set(gca,'colorscale','log');
pbaspect([1 1 1]);
xlabel('KE_{Br(1)^+} + KE_{Br(2)^+}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('cos\theta_{Br(1)^+- Br(2)^+}', 'FontWeight', 'normal','FontName', 'Arial');
scatter(+ 9.2841,-0.48967,300,'x','g','LineWidth',3) %concerted CES resutls 119.3191 4.6419 4.6422 
scatter(+8.7388,-0.49273,300,'x','b','LineWidth',3) %concerted CES results 119.52 4.8227 3.9161

%%
%plot KE1 vs KE2 vs ThetaP1P2  if does or does not run commment or uncomment  line 2248 % cos_theta_12=cos_theta_12(j_KE_all);
close all
binsize_ke=0.05; %eV   
edges_ke=[0:binsize_ke:28];
figure;

% cos_theta_12=cos_theta_12(j_KE_all);
count_12_ke = histcounts2(KE1,KE2,edges_ke,edges_ke);
count_12_ke_mod = max(count_12_ke,1); 

myColorMap=flipud(hot);
% myColorMap(1,:) = 1;
imagesc( edges_ke,edges_ke, ((count_12_ke_mod')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([0 8 0 8])
hold on;
xticks([0:2:8]);
yticks([0:2:8]);
set(gca,'colorscale','log');
pbaspect([1 1 1]);
xlabel('KE_{Br(1)^+}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('KE_{Br(2)^+}/eV', 'FontWeight', 'normal','FontName', 'Arial');
scatter(4.6419,4.6422,1000,'p','MarkerEdgeColor',mycoloredge,'MarkerFaceColor',mycolor,'LineWidth',3) %concerted CES resutls 119.3191 4.6419 4.6422 
set(gca,'FontSize',40)
%%
%native frame intermediate 1st and 2nd
close all;
m(1)=frag_m(1);
m(2)=frag_m(2);
m(3)=frag_m(3);

m_123=(m(1)+m(2)+m(3));
m_12=(m(1)+m(2));
mu_12=1/(1/m(1)+1/m(2));


p12_x=mu_12*(px2./m(2)-px1./m(1));
p12_y=mu_12*(py2./m(2)-py1./m(1));
p12_z=mu_12*(pz2./m(2)-pz1./m(1));
%second step
p12_m=sqrt(p12_x.^2 + p12_y.^2 + p12_z.^2);


p12_3_x=(m_12/m_123)*px3 - (m(3)/m_123)*(px1+px2);
p12_3_y=(m_12/m_123)*py3 - (m(3)/m_123)*(py1+py2);
p12_3_z=(m_12/m_123)*pz3 - (m(3)/m_123)*(pz1+pz2);
%first step
p12_3_m=sqrt(p12_3_x.^2 + p12_3_y.^2 + p12_3_z.^2);

p12_3_dot_p12=(p12_3_x.*p12_x + p12_3_y.*p12_y +p12_3_z.*p12_z );

ke_12 = p12_m.^2/(2*mu_12)*hatoev;
theta_12_3=acos(p12_3_dot_p12./(p12_3_m.*p12_m))*180/pi;

Xedges = min(ke_12(:,1)):(max(ke_12(:,1)) - min(ke_12(:,1)) )/200:max(ke_12(:,1));
Yedges = min(theta_12_3(:,1)):(max(theta_12_3(:,1)) - min(theta_12_3(:,1)) )/250:max(theta_12_3(:,1));

count = histcounts2(ke_12(:,1),theta_12_3(:,1),Xedges,Yedges);

count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([0 15 0 180])
xticks([0:5:15]);
yticks([0:60:180]);
hold on;


% axis equal;
pbaspect([1 1 1]);

xlabel(strcat('KE','[',frag_m_z_str(1),'(',num2str(1),')','-',frag_m_z_str(2),'(',num2str(2),')', ']'  ,'/eV'),'FontWeight', 'normal','FontName', 'Arial');
ylabel(strcat('\theta','[',frag_m_z_str(1),'(',num2str(1),')','-',frag_m_z_str(2),'(',num2str(2),')',',',frag_m_z_str(3),'(',num2str(3),')', ']','/^{\circ}'), 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');


%%
%native frame intermediate 2nd and 3rd
close all; 
m(1)=frag_m(1);
m(2)=frag_m(2);
m(3)=frag_m(3);

m_123=(m(1)+m(2)+m(3));
m_23=(m(2)+m(3));
mu_23=1/(1/m(2)+1/m(3));

p23_x=mu_23*(px3./m(3)-px2./m(2));
p23_y=mu_23*(py3./m(3)-py2./m(2));
p23_z=mu_23*(pz3./m(3)-pz2./m(2));
p23_m=sqrt(p23_x.^2 + p23_y.^2 + p23_z.^2);

p23_1_x=(m_23/m_123)*px1 - (m(1)/m_123)*(px2+px3);
p23_1_y=(m_23/m_123)*py1 - (m(1)/m_123)*(py2+py3);
p23_1_z=(m_23/m_123)*pz1 - (m(1)/m_123)*(pz2+pz3);
p23_1_m=sqrt(p23_1_x.^2 + p23_1_y.^2 + p23_1_z.^2);

p23_1_dot_p23=(p23_1_x.*p23_x + p23_1_y.*p23_y +p23_1_z.*p23_z );

ke_23 = p23_m.^2/(2*mu_23)*hatoev;
% ke_23_1 = p23_1_m.^2/(2*(m(2)+m(3)))*hatoev;
% ke123_1=ke_23+ke_23_1;
theta_23_1=acos(p23_1_dot_p23./(p23_1_m.*p23_m))*180/pi;

Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));

count = histcounts2(ke_23(:,1),theta_23_1(:,1),Xedges,Yedges);

count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis([0 15 0 180])
xticks([0:5:15]);
yticks([0:60:180]);
hold on;

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel(strcat('KE','[',frag_m_z_str(2),'(',num2str(2),')','-',frag_m_z_str(3),'(',num2str(3),')', ']'  ,'/eV'),'FontWeight', 'normal','FontName', 'Arial');
ylabel(strcat('\theta','[',frag_m_z_str(2),'(',num2str(2),')','-',frag_m_z_str(3),'(',num2str(3),')',',',frag_m_z_str(1),'(',num2str(1),')', ']','/^{\circ}'), 'FontWeight', 'normal','FontName', 'Arial');


hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
%% gate on native frame
close all
ke_23 = p23_m.^2/(2*mu_23)*hatoev;
j_nat_ke23 = ke_23 < 4.15;

theta_23_1=acos(p23_1_dot_p23./(p23_1_m.*p23_m))*180/pi;

ke_23 = ke_23(j_nat_ke23 );
theta_23_1 = theta_23_1(j_nat_ke23);
Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));

count_23_cor = histcounts2(ke_23(:,1),theta_23_1(:,1),Xedges,Yedges);

count_mod = max(count_23_cor,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_23_cor')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis([0 15 0 180])
xticks([0:5:15]);
yticks([0:60:180]);
hold on;

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel(strcat('KE','[',frag_m_z_str(2),'(',num2str(2),')','-',frag_m_z_str(3),'(',num2str(3),')', ']'  ,'/eV'),'FontWeight', 'normal','FontName', 'Arial');
ylabel(strcat('\theta','[',frag_m_z_str(2),'(',num2str(2),')','-',frag_m_z_str(3),'(',num2str(3),')',',',frag_m_z_str(1),'(',num2str(1),')', ']','/^{\circ}'), 'FontWeight', 'normal','FontName', 'Arial');


hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');

%%
%native frame intermediate 3rd and 1st
% close all;
figure
m(1)=frag_m(1);
m(2)=frag_m(2);
m(3)=frag_m(3);

m_123=(m(1)+m(2)+m(3));
m_31=(m(3)+m(1));
mu_31=1/(1/m(3)+1/m(1));

p31_x=mu_31*(px1./m(1)-px3./m(3));
p31_y=mu_31*(py1./m(1)-py3./m(3));
p31_z=mu_31*(pz1./m(1)-pz3./m(3));

p31_m=sqrt(p31_x.^2 + p31_y.^2 + p31_z.^2);

p31_2_x=(m_31/m_123)*px2 - (m(2)/m_123)*(px3+px1);
p31_2_y=(m_31/m_123)*py2 - (m(2)/m_123)*(py3+py1);
p31_2_z=(m_31/m_123)*pz2 - (m(2)/m_123)*(pz3+pz1);
p31_2_m=sqrt(p31_2_x.^2 + p31_2_y.^2 + p31_2_z.^2);

p31_2_dot_p31=(p31_2_x.*p31_x + p31_2_y.*p31_y +p31_2_z.*p31_z );

ke_31 = p31_m.^2/(2*mu_31)*hatoev;
theta_31_2=acos(p31_2_dot_p31./(p31_2_m.*p31_m))*180/pi;

Xedges = min(ke_31(:,1)):(max(ke_31(:,1)) - min(ke_31(:,1)) )/200:max(ke_31(:,1));
Yedges = min(theta_31_2(:,1)):(max(theta_31_2(:,1)) - min(theta_31_2(:,1)) )/250:max(theta_31_2(:,1));

count = histcounts2(ke_31(:,1),theta_31_2(:,1),Xedges,Yedges);

count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis([0 15 0 180])
xticks([0:5:15]);
yticks([0:60:180]);
hold on;


set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);


xlabel(strcat('KE','[',frag_m_z_str(3),'(',num2str(3),')','-',frag_m_z_str(1),'(',num2str(1),')', ']'  ,'/eV'),'FontWeight', 'normal','FontName', 'Arial');
ylabel(strcat('\theta','[',frag_m_z_str(3),'(',num2str(3),')','-',frag_m_z_str(1),'(',num2str(1),')',',',frag_m_z_str(2),'(',num2str(2),')', ']','/^{\circ}'), 'FontWeight', 'normal','FontName', 'Arial');


hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');

%%
close all
ke_31 = p31_m.^2/(2*mu_31)*hatoev;
j_nat_ke31 = ke_31 < 4.25;

theta_31_2=acos(p31_2_dot_p31./(p31_2_m.*p31_m))*180/pi;

ke_31 = ke_31(j_nat_ke31 );
theta_31_2 = theta_31_2(j_nat_ke31);
Xedges = min(ke_31(:,1)):(max(ke_31(:,1)) - min(ke_31(:,1)) )/200:max(ke_31(:,1));
Yedges = min(theta_31_2(:,1)):(max(theta_31_2(:,1)) - min(theta_31_2(:,1)) )/250:max(theta_31_2(:,1));

count_31_cor = histcounts2(ke_31(:,1),theta_31_2(:,1),Xedges,Yedges);

count_mod = max(count_31_cor,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis([0 15 0 180])
xticks([0:5:15]);
yticks([0:60:180]);
hold on;

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);


xlabel(strcat('KE','[',frag_m_z_str(3),'(',num2str(3),')','-',frag_m_z_str(1),'(',num2str(1),')', ']'  ,'/eV'),'FontWeight', 'normal','FontName', 'Arial');
ylabel(strcat('\theta','[',frag_m_z_str(3),'(',num2str(3),')','-',frag_m_z_str(1),'(',num2str(1),')',',',frag_m_z_str(2),'(',num2str(2),')', ']','/^{\circ}'), 'FontWeight', 'normal','FontName', 'Arial');


hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');

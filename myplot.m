function myplot(Xedges,Yedges,t1_gate,t23_gate,t0,k0,frag,sct_peak_tof_1,sct_peak_tof_23)

Xedges = min(t1_gate(:,1)):5:max(t1_gate(:,1));
Yedges = min(t23_gate(:,1)):5:max(t23_gate(:,1));

count = histcounts2(t1_gate(:,1),t23_gate(:,1),Xedges,Yedges);
count_mod=max(count,1);
myColorMap = jet;
myColorMap = flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% caxis([-100 1000])
axis xy;
%  axis([4000 7000 12000 15000])
% % axis([0 8000 0 17000])
% %axis([0 10000 0 20000])
%  axis([0 3200 0 12000])
%axis([3200 7200 0 12000])
% axis([0 3200 12000 17000])
% axis([3200 7200 12000 17000])
% 
% axis([3200 7200 0 12000])
% axis([3200 7200 11000 17000])
%axis equal
% pbaspect([1 1 1])

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


xlabel('TOF_1/ns', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('TOF_2 + TOF_3/ns', 'FontWeight', 'normal','FontName', 'Arial');
%clims=[1 7];
set(gca,'colorscale','log');
hold on;
set(gca,'FontSize',25)
% clearvars -except fileno nfiles t1_bin t23_bin;
scatter([sct_peak_tof_1], [sct_peak_tof_23],2000,'x','g','LineWidth',3)

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

% m=(C+H+Br79*3)*mpvsme-q;


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

 

 for i=1:length(frag_m_z)
         xline(  k0*(frag_m_z(i))^0.5+t0, '--b', frag_m_z_str(i),'LabelVerticalAlignment' ,'bottom','fontweight','bold','fontsize',30);
 end

ii=2;
jj=3;
yline(  k0*( (frag_m_z(ii) ))^0.5+t0 + k0*( (frag_m_z(jj) ))^0.5+t0, '--b', strcat(frag_m_z_str(ii), '+', frag_m_z_str(jj)),'fontweight','bold','fontsize',30 );

% ii=2;
% 
% yline(  k0*( (frag_m_z(ii) ))^0.5+t0, '--b', strcat(frag_m_z_str(ii)),'fontweight','bold','fontsize',30);
% 
% ii=3;
% 
% yline(  k0*( (frag_m_z(ii) ))^0.5+t0, '--b', strcat(frag_m_z_str(ii)),'fontweight','bold','fontsize',30);

end
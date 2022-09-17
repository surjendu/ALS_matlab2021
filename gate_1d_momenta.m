function [j_momentum_gate,px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,p1,p2,p3]...
            =gate_1d_momenta(px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,p1,p2,p3,...
                    px1_rng, px2_rng, px3_rng, px_sum_rng,py1_rng, py2_rng, py3_rng, py_sum_rng,...
                        pz1_rng, pz2_rng, pz3_rng, pz_sum_rng,binsize_px,binsize_py,binsize_pz,binsize_p)


j_momentum_gate = px1 > px1_rng(1) &  px1 < px1_rng(2) &  py1 > py1_rng(1) &  py1 <  py1_rng(2) &  pz1 > pz1_rng(1)  &    pz1  < pz1_rng(2)...
                & px2 > px2_rng(1) &  px2 < px2_rng(2) &  py2 > py2_rng(1) &  py2 <  py2_rng(2) &  pz2 > pz2_rng(1)  &    pz2  < pz2_rng(2)...
                & px3 > px3_rng(1) &  px3 < px3_rng(2) &  py3 > py3_rng(1) &  py3 <  py3_rng(2) &  pz3 > pz3_rng(1)  &    pz3  < pz3_rng(2)...
                & px1+px2+px3 > px_sum_rng(1) & px1+px2+px3 < px_sum_rng(2) &  py1+py2+py3 > py_sum_rng(1) &  py1+py2+py3 < py_sum_rng(2) &  pz1+pz2+pz3 > pz_sum_rng(1)  &   pz1+pz2+pz3 < pz_sum_rng(2);

            
px1 = px1(j_momentum_gate);
py1 = py1(j_momentum_gate);
pz1 = pz1(j_momentum_gate);
p1  = p1 (j_momentum_gate);   


px2 = px2(j_momentum_gate);
py2 = py2(j_momentum_gate);
pz2 = pz2(j_momentum_gate);
p2  = p2 (j_momentum_gate);

px3 = px3(j_momentum_gate);
py3 = py3(j_momentum_gate);
pz3 = pz3(j_momentum_gate);
p3  = p3 (j_momentum_gate);

%    delay_step=delay_step(j_momentum_gate);
    
%plot tof in ns

edge_max=ceil(max([max([px1;px2;px3;py1;py2;py3;pz1;pz2;pz3]); -min([px1;px2;px3;py1;py2;py3;pz1;pz2;pz3])]));

edges_px=[-edge_max:binsize_px:edge_max];
edges_py=[-edge_max:binsize_py:edge_max]; 
edges_pz=[-edge_max:binsize_pz:edge_max];
edges_p=[min([p1;p2;p3]):binsize_p:max([p1;p2;p3])];

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

[px3_counts,px3_edges]=histcounts(px3,edges_px);
[py3_counts,py3_edges]=histcounts(py3,edges_py);
[pz3_counts,pz3_edges]=histcounts(pz3,edges_pz);
[p3_counts,p3_edges]=histcounts(p2,edges_p);

[px_counts,px_edges]=histcounts(px1+px2+px3,edges_px);
[py_counts,py_edges]=histcounts(py1+py2+py3,edges_py);
[pz_counts,pz_edges]=histcounts(pz1+pz2+pz3,edges_pz);

% p_allhits_hist={[bincent_px', px_counts'], [bincent_py', py_counts'], [bincent_pz' pz_counts'], [bincent_p' p_counts']};
%  dlmwrite('p_allhits_hist.csv',p_allhits_hist);

close all;
%figure 
subplot(2,2,1)
plot(bincent_px,px1_counts,'r',bincent_px,fliplr(px1_counts),'b',bincent_px,(px1_counts+fliplr(px1_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([px1_rng(1)*1.25 px1_rng(2)*1.25]);
ylim([0 max(px1_counts)*1.1]);
xlabel('px1/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'px1','-px1','px1\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')


subplot(2,2,2)
%figure
plot(bincent_px,px2_counts,'r',bincent_px,fliplr(px2_counts),'b',bincent_px,(px2_counts+fliplr(px2_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([px2_rng(1)*1.25 px2_rng(2)*1.25]);
ylim([0 max(px2_counts)*1.1]);
xlabel('px2/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% ytickformat('%0.0e')
legend({'px2','-px2','px2\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

subplot(2,2,3)
%figure
plot(bincent_px,px3_counts,'r',bincent_px,fliplr(px3_counts),'b',bincent_px,(px3_counts+fliplr(px3_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([px3_rng(1)*1.25 px3_rng(2)*1.25]);
ylim([0 max(px3_counts)*1.1]);
xlabel('px2/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% ytickformat('%0.0e')
legend({'px3','-px3','px3\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

subplot(2,2,4) 
%figure
plot(bincent_px,px_counts,'r','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([px_sum_rng(1)*1.5 px_sum_rng(2)*1.5]);
ylim([0 max(px_counts)*1.1]);
xlabel('\Sigmapx/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'\Sigmapx'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')



figure
subplot(2,2,1)
plot(bincent_py,py1_counts,'r',bincent_py,fliplr(py1_counts),'b',bincent_py,(py1_counts+fliplr(py1_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([py1_rng(1)*1.25 py1_rng(2)*1.25]);
ylim([0 max(py1_counts)*1.1]);
xlabel('px1/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'py1','-py1','py1\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')


subplot(2,2,2)
%figure
plot(bincent_py,py2_counts,'r',bincent_py,fliplr(py2_counts),'b',bincent_py,(py2_counts+fliplr(py2_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([py2_rng(1)*1.25 py2_rng(2)*1.25]);
ylim([0 max(py2_counts)*1.1]);
xlabel('py2/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% ytickformat('%0.0e')
legend({'py2','-py2','py2\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

subplot(2,2,3)
%figure
plot(bincent_py,py3_counts,'r',bincent_py,fliplr(py3_counts),'b',bincent_py,(py3_counts+fliplr(py3_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([py3_rng(1)*1.25 py3_rng(2)*1.25]);
ylim([0 max(py3_counts)*1.1]);
xlabel('py3/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% ytickformat('%0.0e')
legend({'py3','-py3','py3\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

subplot(2,2,4) 
%figure
plot(bincent_py,py_counts,'r','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([py_sum_rng(1)*1.5 py_sum_rng(2)*1.5]);
ylim([0 max(py_counts)*1.1]);
xlabel('\Sigmapy/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'\Sigmapy'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')



figure
subplot(2,2,1)
plot(bincent_pz,pz1_counts,'r',bincent_pz,fliplr(pz1_counts),'b',bincent_pz,(pz1_counts+fliplr(pz1_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([pz1_rng(1)*1.25 pz1_rng(2)*1.25]);
ylim([0 max(pz1_counts)*1.1]);
xlabel('pz1/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'pz1','-pz1','pz1\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')


subplot(2,2,2)
%figure
plot(bincent_pz,pz2_counts,'r',bincent_pz,fliplr(pz2_counts),'b',bincent_pz,(pz2_counts+fliplr(pz2_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([pz2_rng(1)*1.25 pz2_rng(2)*1.25]);
ylim([0 max(pz2_counts)*1.1]);
xlabel('pz2/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% ytickformat('%0.0e')
legend({'pz2','-pz2','pz2\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

subplot(2,2,3)
%figure
plot(bincent_pz,pz3_counts,'r',bincent_pz,fliplr(pz3_counts),'b',bincent_pz,(pz3_counts+fliplr(pz3_counts))./2,'k','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([pz3_rng(1)*1.25 pz3_rng(2)*1.25]);
ylim([0 max(pz3_counts)*1.1]);
xlabel('pz3/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
% ytickformat('%0.0e')
legend({'pz3','-pz3','pz3\_avg'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')

subplot(2,2,4) 
%figure
plot(bincent_pz,pz_counts,'r','LineWidth',2);
set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30,'GridColor', [0.1, 0.7, 0.2],'GridAlpha',0.95,'LineWidth',1 );
grid on;
xlim([pz_sum_rng(1)*1.5 pz_sum_rng(2)*1.5]);
ylim([0 max(pz_counts)*1.1]);
xlabel('\Sigmapz/au',  'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
legend({'\Sigmapz'}, 'FontSize', 20, 'FontWeight', 'normal','FontName', 'Arial')
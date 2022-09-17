function momentum_sum_correlation(px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,binsize_px,binsize_py,binsize_pz,binsize_p);

edge_max=ceil(max([max([px1;px2+px3;py1;py2+py3;pz1;pz2+pz3]); -min([px1;px2+px3;py1;py2+py3;pz1;pz2+pz3])]));

edges_px=[-edge_max:binsize_px:edge_max];
edges_py=[-edge_max:binsize_py:edge_max]; 
edges_pz=[-edge_max:binsize_pz:edge_max];

px_sum=px1+px2+px3;  py_sum=py1+py2+py3;  pz_sum=pz1+pz2+pz3;

edge_max_sum=ceil(max([max([px_sum;py_sum;pz_sum]); -min([px_sum;py_sum;pz_sum])]));

edges_px_sum=[-edge_max_sum:binsize_px:edge_max_sum];
edges_py_sum=[-edge_max_sum:binsize_py:edge_max_sum]; 
edges_pz_sum=[-edge_max_sum:binsize_pz:edge_max_sum];

figure
subplot(1,3,1)
myColorMap = jet;
myColorMap(1,:) = 0;
% subplot(1,2,1)
count = histcounts2(px_sum,py_sum,edges_px_sum,edges_py_sum);
count_mod=max(count,1);
imagesc( edges_px_sum,edges_py_sum, ((count_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
xlabel('p_xsum / au', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('p_ysum  / au', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'colorscale','log');
set(gca,'FontSize',25)
axis xy;
pbaspect([1 1 1])

subplot(1,3,2)
myColorMap = jet;
myColorMap(1,:) = 0;
% subplot(1,2,1)
count = histcounts2(pz_sum,px_sum,edges_pz_sum,edges_px_sum);
count_mod=max(count,1);
imagesc( edges_pz_sum,edges_px_sum, ((count_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
xlabel('p_zsum / au', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('p_xsum  / au', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'colorscale','log');
set(gca,'FontSize',25)
axis xy;
pbaspect([1 1 1])

subplot(1,3,3)
myColorMap = jet;
myColorMap(1,:) = 0;
% subplot(1,2,1)
count = histcounts2(pz_sum,py_sum,edges_pz_sum,edges_py_sum);
count_mod=max(count,1);
imagesc( edges_pz_sum,edges_py_sum, ((count_mod')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
xlabel('p_zsum / au', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('p_ysum  / au', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'colorscale','log');
set(gca,'FontSize',25)
axis xy;
pbaspect([1 1 1])


end
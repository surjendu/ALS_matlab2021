function momentum_sum_correlation(px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,binsize_px,binsize_py,binsize_pz,binsize_p);

edge_max=ceil(max([max([px1;px2+px3;py1;py2+py3;pz1;pz2+pz3]); -min([px1;px2+px3;py1;py2+py3;pz1;pz2+pz3])]));

edges_px=[-edge_max:binsize_px:edge_max];
edges_py=[-edge_max:binsize_py:edge_max]; 
edges_pz=[-edge_max:binsize_pz:edge_max];

px_sum=px1+px2+px3;  py_sum=py1+py2+py3;  pz_sum=pz1+pz2+pz3;
px_diff=px1-px2-px3; py_diff=py1-py2-py3; pz_diff=pz1-pz2-pz3;

edge_max_sum=ceil(max([max([px_sum;py_sum;pz_sum;px_diff;py_diff;pz_diff]); -min([px_sum;py_sum;pz_sum;px_diff;py_diff;pz_diff])]));

edges_px_sum=[-edge_max_sum:binsize_px:edge_max_sum];
edges_py_sum=[-edge_max_sum:binsize_py:edge_max_sum]; 
edges_pz_sum=[-edge_max_sum:binsize_pz:edge_max_sum];

myColorMap = jet;
myColorMap(1,:) = 0;
subplot(1,2,1)
count = histcounts2(px1,px2+px3,edges_px,edges_px);
imagesc( edges_px,edges_px, ((count')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
xlabel('px_1 / au', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('px_2 + px_3 / au', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'colorscale','log');
set(gca,'FontSize',25)
axis xy;



subplot(1,2,2)
count = histcounts2((px1-px2-px3),(px1+px2+px3),edges_px_sum,edges_px_sum);
imagesc( edges_px_sum,edges_px_sum, ((count')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
xlabel('px_1 - px_2 - px_3 / au', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('px_1 + px_2 + px_3 / au', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'colorscale','log');
set(gca,'FontSize',25)
axis xy;

figure
subplot(1,2,1)
count = histcounts2(py1,py2+py3,edges_py,edges_py);
imagesc( edges_py,edges_py, ((count')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
xlabel('py_1 / au', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('py_2 + py_3 / au', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'colorscale','log');
set(gca,'FontSize',25)
axis xy;

subplot(1,2,2)
count = histcounts2((py1-py2-py3),(py1+py2+py3),edges_px_sum,edges_px_sum);
imagesc( edges_py_sum,edges_py_sum, ((count')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
xlabel('py_1 - py_2 - py_3 / au', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('py_1 + py_2 + py_3 / au', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'colorscale','log');
set(gca,'FontSize',25)
axis xy;

figure
subplot(1,2,1)
count = histcounts2(pz1,pz2+pz3,edges_pz,edges_pz);
imagesc( edges_pz,edges_pz, ((count')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
xlabel('pz_1 / au', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('pz_2 + pz_3 / au', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'colorscale','log');
set(gca,'FontSize',25)
axis xy;

subplot(1,2,2)
count = histcounts2((pz1-pz2-pz3),(pz1+pz2+pz3),edges_pz_sum,edges_pz_sum);
imagesc( edges_pz_sum,edges_pz_sum, ((count')));
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
xlabel('pz_1 - pz_2 - pz_3 / au', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('pz_1 + pz_2 + pz_3 / au', 'FontWeight', 'normal','FontName', 'Arial');
set(gca,'colorscale','log');
set(gca,'FontSize',25)
axis xy;
end
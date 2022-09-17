function [px1,px2,px3,py1,py2,py3,pz1,pz2,pz3]=eventswap_momenta_23(px1,px2,px3,py1,py2,py3,pz1,pz2,pz3)

p1_swap=[px1 py1 pz1];
p2_swap=[px2 py2 pz2];
p3_swap=[px3 py3 pz3];

l=length(p2_swap(:,1));
nel = round (l/2);   
index = randperm(l,nel);
index=index';

    
scnd_1=p2_swap;
scnd_1(index,:)=[];

thrd_2=p3_swap;
thrd_2(index,:)=[];

scnd_swap=[p2_swap(index,:) ;thrd_2];
thrd_swap=[p3_swap(index,:) ;scnd_1];

frst_2=p1_swap;
frst_1=frst_2(index,:);
frst_2(index,:)=[];
frst=[frst_1;frst_2];

p1_swap=[];
p2_swap=[];
p3_swap=[];

p1_swap=frst;
p2_swap=scnd_swap;
p3_swap=thrd_swap;

px1=p1_swap(:,1);  py1=p1_swap(:,2);  pz1=p1_swap(:,3); 
px2=p2_swap(:,1);  py2=p2_swap(:,2);  pz2=p2_swap(:,3); 
px3=p3_swap(:,1);  py3=p3_swap(:,2);  pz3=p3_swap(:,3); 

end
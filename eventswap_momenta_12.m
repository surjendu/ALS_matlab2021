function [px1,px2,px3,py1,py2,py3,pz1,pz2,pz3]=eventswap_momenta_12(px1,px2,px3,py1,py2,py3,pz1,pz2,pz3)

p1_swap=[px1 py1 pz1];
p2_swap=[px2 py2 pz2];
p3_swap=[px3 py3 pz3];

l=length(p1_swap(:,1));
nel = round (l/2);   
index = randperm(l,nel);
index=index';

frst_1=p1_swap;
frst_1(index,:)=[];
    
scnd_2=p2_swap;
scnd_2(index,:)=[];
  
frst_swap=[p1_swap(index,:) ; scnd_2  ];
scnd_swap=[p2_swap(index,:) ;frst_1];
  
 thrd_2=p3_swap;
 thrd_1=thrd_2(index,:);
 thrd_2(index,:)=[];
 thrd=[thrd_1;thrd_2];

 p1_swap=[];
 p2_swap=[];
 p3_swap=[];
   
 p1_swap=frst_swap;
 p2_swap=scnd_swap;
 p3_swap=thrd;
 
 px1=p1_swap(:,1);  py1=p1_swap(:,2);  pz1=p1_swap(:,3); 
 px2=p2_swap(:,1);  py2=p2_swap(:,2);  pz2=p2_swap(:,3); 
 px3=p3_swap(:,1);  py3=p3_swap(:,2);  pz3=p3_swap(:,3); 

end
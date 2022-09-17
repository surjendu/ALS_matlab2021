function [measurement]=eventswap_data_23(measurement)

l=length(measurement.data.raw.XYT.scnd(:,1));
nel = round (l/2);   
index = randperm(l,nel);
index=index';
  
scnd_1=measurement.data.raw.XYT.scnd;
scnd_1(index,:)=[];
    
thrd_2=measurement.data.raw.XYT.thrd;
thrd_2(index,:)=[];
  
scnd_swap=[measurement.data.raw.XYT.scnd(index,:) ;  thrd_2 ];
thrd_swap=[measurement.data.raw.XYT.thrd(index,:) ;  scnd_1 ];
  
frst_2=measurement.data.raw.XYT.frst;
frst_1=frst_2(index,:);
frst_2(index,:)=[];
frst=[frst_1;frst_2];
   
measurement.data.raw.XYT.frst=[];
measurement.data.raw.XYT.scnd=[];
measurement.data.raw.XYT.thrd=[];
  
measurement.data.raw.XYT.frst=frst;
measurement.data.raw.XYT.scnd=scnd_swap;
measurement.data.raw.XYT.thrd=thrd_swap;

end
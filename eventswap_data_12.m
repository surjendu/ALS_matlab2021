function [measurement]=eventswap_data_12(measurement)

  l=length(measurement.data.raw.XYT.frst(:,1));
  nel = round (l/2);   
  index = randperm(l,nel);
  index=index';

  frst_1=measurement.data.raw.XYT.frst;
  frst_1(index,:)=[];
    
  scnd_2=measurement.data.raw.XYT.scnd;
  scnd_2(index,:)=[];
  
  frst_swap=[measurement.data.raw.XYT.frst(index,:) ;  scnd_2   ];
  scnd_swap=[measurement.data.raw.XYT.scnd(index,:) ; frst_1];
  
  thrd_2=measurement.data.raw.XYT.thrd;
  thrd_1=thrd_2(index,:);
  thrd_2(index,:)=[];
  thrd=[thrd_1;thrd_2];
   
  measurement.data.raw.XYT.frst=[];
  measurement.data.raw.XYT.scnd=[];
  measurement.data.raw.XYT.thrd=[];
  
  measurement.data.raw.XYT.frst=frst_swap;
  measurement.data.raw.XYT.scnd=scnd_swap;
  measurement.data.raw.XYT.thrd=thrd;

end

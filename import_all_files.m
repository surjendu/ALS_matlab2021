function [measurement]=import_all_files(nfiles,fileno,t1_bin,t23_bin)

tic
measurement_data_delay_position = [];
measurement_data_XYT_frst=[];
measurement_data_XYT_scnd=[];
measurement_data_XYT_thrd=[];
disp('Importing data now')

for i=1:nfiles 
filename=num2str(i);
files= dir(['Run',num2str(fileno),'_XYT3bin',filename,'.out']);
file = files';
fileID_real=fopen(file.name);
x=fread(fileID_real,[11 Inf],'*float');
fclose(fileID_real);
x = x';

id = find(x(:,1));
delay_position = x(id,11);
XYT_frst = x(id, 1:3); 
XYT_scnd = x(id, 4:6);
XYT_thrd = x(id, 7:9); 
% laser_shot = x(id,10);  




clearvars x

t1 = XYT_frst(:,3);
t2 = XYT_scnd(:,3);
t3 = XYT_thrd(:,3);
t23=t2+t3;
j_gate = t1 > t1_bin(1)  &  t1 < t1_bin(2) &  t23 > t23_bin(1)  & t23 < t23_bin(2);

delay_position=delay_position(j_gate,:);
XYT_frst=XYT_frst(j_gate,:);
XYT_scnd=XYT_scnd(j_gate,:);
XYT_thrd=XYT_thrd(j_gate,:);
% laser_shot = laser_shot(j_gate,:);

measurement_data_delay_position = [measurement_data_delay_position;delay_position];
measurement_data_XYT_frst = [measurement_data_XYT_frst; XYT_frst];
measurement_data_XYT_scnd = [measurement_data_XYT_scnd; XYT_scnd];
measurement_data_XYT_thrd = [measurement_data_XYT_thrd; XYT_thrd];


disp(['Done Importing',' ',filename]);
pause(0.1);
end
measurement.files = {files.name};
measurement.data.raw.delay_position = measurement_data_delay_position;
measurement.data.raw.XYT.frst =  measurement_data_XYT_frst; 
measurement.data.raw.XYT.scnd =  measurement_data_XYT_scnd;
measurement.data.raw.XYT.thrd =  measurement_data_XYT_thrd;

clearvars measurement_data_XYT_frst measurement_data_XYT_scnd measurement_data_XYT_thrd XYT_frst XYT_scnd XYT_thrd

save('measurement', 'measurement', '-v7.3'); % the measurement is the file name second is the variable itself
toc

end
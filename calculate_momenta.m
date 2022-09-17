function [px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,p1,p2,p3]...
                = calculate_momenta(frag_m,frag_m_z,x1_gate_rot, x0 ,v0x, y1_gate_rot,y0,v0y,z0, t0, k0,...
                   x2_gate_rot,y2_gate_rot,x3_gate_rot,y3_gate_rot, t1_gate_rot,t2_gate_rot,t3_gate_rot,v0z,...
                   mmns_to_au_conv,t0_off_1,t0_off_2,t0_off_3,z_cal_1,z_cal_2,z_cal_3,image_cal_xy)


px1 = frag_m(1)*(image_cal_xy)*((x1_gate_rot-x0)./(t1_gate_rot-t0)-v0x)*mmns_to_au_conv;
py1 = frag_m(1)*(image_cal_xy)*((y1_gate_rot-y0)./(t1_gate_rot-t0)-v0y)*mmns_to_au_conv;
pz1 = z_cal_1*(t1_gate_rot-t0-t0_off_1)*mmns_to_au_conv;
%pz1 = frag_m(1)*(  (Lion-z0)./(t1_gate_rot-t0)  - 0.5* (t1_gate_rot-t0)*(1 * 2 * (Lion - z0) / (frag_m_z(1)*k0^2) )-v0z)*mmns_to_au_conv;
p1=sqrt(px1.^2+py1.^2+pz1.^2);   

px2 = frag_m(2)*(image_cal_xy)*((x2_gate_rot-x0)./(t2_gate_rot-t0)-v0x)*mmns_to_au_conv;
py2 = frag_m(2)*(image_cal_xy)*((y2_gate_rot-y0)./(t2_gate_rot-t0)-v0y)*mmns_to_au_conv;
pz2 = z_cal_2*(t2_gate_rot-t0-t0_off_2)*mmns_to_au_conv;
%pz2 = frag_m(2)*(  (Lion-z0)./(t2_gate_rot-t0)  - 0.5* (t2_gate_rot-t0)*(1 * 2 * (Lion - z0) / (frag_m_z(2)*k0^2) )-v0z)*mmns_to_au_conv;
p2=sqrt(px2.^2+py2.^2+pz2.^2);  

px3 = frag_m(3)*(image_cal_xy)*((x3_gate_rot-x0)./(t3_gate_rot-t0)-v0x)*mmns_to_au_conv;
py3 = frag_m(3)*(image_cal_xy)*((y3_gate_rot-y0)./(t3_gate_rot-t0)-v0y)*mmns_to_au_conv;
pz3 = z_cal_3*(t3_gate_rot-t0-t0_off_3)*mmns_to_au_conv;
%pz3 = frag_m(3)*(  (Lion-z0)./(t3_gate_rot-t0)  - 0.5* (t3_gate_rot-t0)*(1 * 2 * (Lion - z0) / (frag_m_z(3)*k0^2) )-v0z)*mmns_to_au_conv;
p3=sqrt(px3.^2+py3.^2+pz3.^2); 

end
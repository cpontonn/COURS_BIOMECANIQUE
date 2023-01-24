function R=transfo_cap_to_opensim(q1,q2,q3)
% q1 rotation autour de x en ° 
% q2 rotation autour de y en °
% q3 rotation autour de z en °
R=[1 0 0;0 cosd(q1) -sind(q1);0 sind(q1) cosd(q1)]*...
    [cosd(q2) 0 sind(q2);0 1 0;-sind(q2) 0 cosd(q2)]*...
    [cosd(q3) -sind(q3) 0;sind(q3) cosd(q3) 0; 0  0 1];

end
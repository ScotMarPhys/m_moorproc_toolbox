function [ur,vr]=uvrot_test(u,v,rot)
% function [ur,vr]=uvrot(u,v,rot)
% rotate velocities
% by angle rot in deg
% M. Visbeck
rot = -rot*pi/180;
cr=cos(rot);
sr=sin(rot);
ur=u.*cr-v.*sr;
vr=u.*sr+v.*cr;
return
end
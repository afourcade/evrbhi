function vc_new=vc2vcdeform(vc,vals);

[np ndum]=size(vc.vc);
rads=sqrt(sum(vc.vc.^2,2));
vc_new=vc;
vc_new.vc=vc.vc.*repmat((rads+vals)./rads,1,ndum);

return;
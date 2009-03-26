function [phi,psi] = run_le(nxy,curr_dim)

NX = size(nxy,1);
NY = size(nxy,2);
pxy0 = nxy/sum(nxy(:));

xv0 = rand(NX,curr_dim);
yv0 = rand(NY,curr_dim);

x0 = [xv0(:);yv0(:)];

cg_params.NX = NX;
cg_params.NY = NY;
cg_params.dim = curr_dim;
cg_params.pxy0 = pxy0;

opts = -1*zeros(1,9);
x = conj_grad('logemb_grad',cg_params,x0,opts);
[phi,psi] = read_legrad_x(x,cg_params);

if 0 
clf ; hold on ;
for(i_x=1:NX)
  h=plot(phi(i_x,1),phi(i_x,2),'.');  
  cm = uicontextmenu;
  set(h,'UIContextMenu', cm);
  item1 =uimenu(cm,'Label',sprintf('%d',i_x));
end

 for(i_x=1:NY)
  h=plot(psi(i_x,1),psi(i_x,2),'r.');    
  cm = uicontextmenu;
  set(h,'UIContextMenu', cm);
  item1 =uimenu(cm,'Label',sprintf('%d',i_x));
end
end
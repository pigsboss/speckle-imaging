function g=shift(f,x,y)
%SHIFT Shift data in f by x along x-axis and y along y-axis, so that
%g(x',y') = f(x'-x,y'-y).

siz=size(f);
h=[rot90(f(2:siz(1),2:siz(2)),2),flipud(f(2:siz(1),:)),rot90(f(2:siz(1),1:siz(2)-1),2);...
  fliplr(f(:,2:siz(2))),f,fliplr(f(:,1:siz(2)-1));...
  rot90(f(1:siz(1)-1,2:siz(2)),2),flipud(f(1:siz(1)-1,:)),rot90(f(2:siz(1),2:siz(2)),2)];
g=h(siz(1)-y:siz(1)+siz(1)-1-y,siz(2)-x:siz(2)+siz(2)-1-x);
return
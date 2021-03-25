function [he] = FoldKe(Cood, List)
rkj = [Cood(List(2),:)-Cood(List(1),:)]'; 
rij = [Cood(List(3),:)-Cood(List(1),:)]'; 
rkl = [Cood(List(2),:)-Cood(List(4),:)]'; 
rmj = cross(rij,rkj); rnk = cross(rkj,rkl);
sgn = ((abs(rnk'*rij)>1e-8)*sign(rnk'*rij)+(abs(rnk'*rij)<=1e-8)*1);
he = real(acos(rmj'*rnk/(norm(rmj)*norm(rnk)))); 
he = real(sgn*he);
if he<0 
    he = 2*pi+he; 
end;
function [F] = LoadFun(Node,U,icrm)
%% Compute current configuration 
Nodenw = Node;
Nodenw(:,1) = Node(:,1)+U(1:3:end);
Nodenw(:,2) = Node(:,2)+U(2:3:end);
Nodenw(:,3) = Node(:,3)+U(3:3:end);

%% Define Load here based on icrm and Nodenw
if icrm<=0 
    error('Wrong increment!'); 
else
    DIRC = cross(Nodenw(109,:)-Nodenw(98,:),Nodenw(120,:)-Nodenw(109,:));
    Load = [109 0.3*DIRC/norm(DIRC)];
end

%% Wrap up Load info to force/displacement vector F
m = size(Node,1);
F = zeros(3*m,1);
indp = Load(:,1);
F(3*indp-2) = Load(:,2); 
F(3*indp-1) = Load(:,3); 
F(3*indp) = Load(:,4);


function [flow_x,fitness_flow]=initialization(alpha,dim,ub,lb)

for i=1:alpha
    flow_x(i,:)=lb+rand(1,dim).*(ub-lb);%position of each flow
end
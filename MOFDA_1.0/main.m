%___________________________________________________________________%
%  Multi-Objective Flow Direction Algorithm (FDA): source codes version 1.0%
% To run MOFDA: [Best_fitness,BestX,FDA_cg_curve]=MOFDA(Max_iteration,lb,ub,dim,fobj,alpha,beta);
%__________________________________________


alpha=50; % Number of flows
beta=8; %Number of neighborhood
Max_iteration=100; % Maximum number of iterations

switch MultiObjFnc
    case 'Schaffer'         % Schaffer
        num_of_objectives = 2;
        MultiObj.fun = @(x) [x(:).^2, (x(:)-2).^2];
        MultiObj.nVar = 1;
        MultiObj.var_min = -5;
        MultiObj.var_max = 5;
        load('Schaffer.mat');
        MultiObj.truePF = PF;
    case 'Kursawe'          % Kursawe 
        num_of_objectives = 2;
        MultiObj.fun = @(x) [-10.*(exp(-0.2.*sqrt(x(:,1).^2+x(:,2).^2)) + exp(-0.2.*sqrt(x(:,2).^2+x(:,3).^2))), ...
                             sum(abs(x).^0.8 + 5.*sin(x.^3),2)];
        MultiObj.nVar = 3;
        MultiObj.var_min = -5.*ones(1,MultiObj.nVar);
        MultiObj.var_max = 5.*ones(1,MultiObj.nVar);
        load('Kursawe.mat');
        MultiObj.truePF = PF;
    case 'Poloni'           % Poloni's two-objective
        num_of_objectives = 2;
        A1 = 0.5*sin(1)-2*cos(1)+sin(2)-1.5*cos(2);
        A2 = 1.5*sin(1)-cos(1)+2*sin(2)-0.5*cos(2);
        B1 = @(x,y) 0.5.*sin(x)-2.*cos(x)+sin(y)-1.5.*cos(y);
        B2 = @(x,y) 1.5.*sin(x)-cos(x)+2.*sin(y)-0.5.*cos(y);
        f1 = @(x,y) 1+(A1-B1(x,y)).^2+(A2-B2(x,y)).^2;
        f2 = @(x,y) (x+3).^2+(y+1).^2;
        MultiObj.fun = @(x) [f1(x(:,1),x(:,2)), f2(x(:,1),x(:,2))];
        MultiObj.nVar = 2;
        MultiObj.var_min = -pi.*ones(1,MultiObj.nVar);
        MultiObj.var_max = pi.*ones(1,MultiObj.nVar);
    case 'Viennet2'         % Viennet2
        num_of_objectives = 3;
        f1 = @(x,y) 0.5.*(x-2).^2+(1/13).*(y+1).^2+3;
        f2 = @(x,y) (1/36).*(x+y-3).^2+(1/8).*(-x+y+2).^2-17;
        f3 = @(x,y) (1/175).*(x+2.*y-1).^2+(1/17).*(2.*y-x).^2-13;
        MultiObj.fun = @(x) [f1(x(:,1),x(:,2)), f2(x(:,1),x(:,2)), f3(x(:,1),x(:,2))];
        MultiObj.nVar = 2;
        MultiObj.var_min = [-4, -4];
        MultiObj.var_max = [4, 4];
        load('Viennet2.mat');
        MultiObj.truePF = PF;
        
   
    case 'Viennet3'         % Viennet3
        num_of_objectives = 3;
        f1 = @(x,y) 0.5.*(x.^2+y.^2)+sin(x.^2+y.^2);
        f2 = @(x,y) (1/8).*(3.*x-2.*y+4).^2 + (1/27).*(x-y+1).^2 +15;
        f3 = @(x,y) (1./(x.^2+y.^2+1))-1.1.*exp(-(x.^2+y.^2));
        MultiObj.fun = @(x) [f1(x(:,1),x(:,2)), f2(x(:,1),x(:,2)), f3(x(:,1),x(:,2))];
        MultiObj.nVar = 2;
        MultiObj.var_min = [-3, -10];
        MultiObj.var_max = [10, 3];
        load('Viennet3.mat');
        MultiObj.truePF = PF;
    case 'ZDT1'             % ZDT1 (convex)
        num_of_objectives = 2;
        g = @(x) 1+9.*sum(x(:,2:end),2)./(size(x,2)-1);
        MultiObj.fun = @(x) [x(:,1), g(x).*(1-sqrt(x(:,1)./g(x)))];
        MultiObj.nVar = 30; 
        MultiObj.var_min = zeros(1,MultiObj.nVar);
        MultiObj.var_max = ones(1,MultiObj.nVar);
        load('ZDT1.mat');
        MultiObj.truePF = PF;
    case 'ZDT2'             % ZDT2 (non-convex)
        num_of_objectives = 2;
        f = @(x) x(:,1);
        g = @(x) 1+9.*sum(x(:,2:end),2)./(size(x,2)-1);
        h = @(x) 1-(f(x)./g(x)).^2;
        MultiObj.fun = @(x) [f(x), g(x).*h(x)];
        MultiObj.nVar = 30; 
        MultiObj.var_min = zeros(1,MultiObj.nVar);
        MultiObj.var_max = ones(1,MultiObj.nVar);
        load('ZDT2.mat');
        MultiObj.truePF = PF;
    case 'ZDT3'             % ZDT3 (discrete)
        num_of_objectives = 2;
        f = @(x) x(:,1);
        g  = @(x) 1+(9/size(x,2)-1).*sum(x(:,2:end),2);
        h  = @(x) 1 - sqrt(f(x)./g(x)) - (f(x)./g(x)).*sin(10.*pi.*f(x));
        MultiObj.fun = @(x) [f(x), g(x).*h(x)];
        MultiObj.nVar = 30;
        MultiObj.var_min = 0.*ones(1,MultiObj.nVar);
        MultiObj.var_max = 1.*ones(1,MultiObj.nVar);
        load('ZDT3.mat');
        MultiObj.truePF = PF;
    case 'ZDT6'             % ZDT6 (non-uniform)
        num_of_objectives = 2;
        f = @(x) 1 - exp(-4.*x(:,1)).*sin(6.*pi.*x(:,1));
        g = @(x) 1 + 9.*(sum(x(:,2:end),2)./(size(x,2)-1)).^0.25;
        h = @(x) 1 - (f(x)./g(x)).^2;
        MultiObj.fun = @(x) [f(x), g(x).*h(x)];
        MultiObj.nVar = 10;
        MultiObj.var_min = 0.*ones(1,MultiObj.nVar);
        MultiObj.var_max = 1.*ones(1,MultiObj.nVar);
        load('ZDT6.mat');
        MultiObj.truePF = PF;
end
if num_of_objectives > 2
    figure(1)
    scatter3((MultiObj.truePF(:,1).'),(MultiObj.truePF(:,2).'),(PF(:,3).'));
else
    figure(1)
    scatter((MultiObj.truePF(:,1).'),(MultiObj.truePF(:,2).'));
end
    


weights = ones(1,num_of_objectives);


[Best_fitness,BestX,Convergence_curve]=MOFDA(Max_iteration,MultiObj.fun,MultiObj.var_min,MultiObj.var_max,MultiObj.nVar,num_of_objectives,weights,alpha,beta);

Pareto_Sum_Fitness = ones(1,Max_iteration);
for i=1:Max_iteration
    Pareto_Sum_Fitness(1,i) = sum(weights.*(Convergence_curve(:,i).'));
end


figure(2)
%hold on;
%Draw objective space
plot(Pareto_Sum_Fitness)
title('Convergence Curve')
xlabel('Iterations');
ylabel('Pareto weigted fitness sum');
grid on
box on
legend('MOFDA')

display(['The best solution obtained by DFA is : ', num2str(BestX)]);
display(['The best optimal value of the objective funciton found by DFA is F1:']);
display(Convergence_curve(:,Max_iteration));


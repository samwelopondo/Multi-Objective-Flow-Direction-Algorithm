 function [f1,f2,f3] = Schaffer(x,params)         % Schaffer
        c1 = @(x,y) 0.5.*(x-2).^2+(1/13).*(y+1).^2+3;
        c2 = @(x,y) (1/36).*(x+y-3).^2+(1/8).*(-x+y+2).^2-17;
        c3 = @(x,y) (1/175).*(x+2.*y-1).^2+(1/17).*(2.*y-x).^2-13;
        f1 = c1(x(:,1),x(:,2));
        f2 = c2(x(:,1),x(:,2));
        f3 = c3(x(:,1),x(:,2));
        MultiObj.nVar = 2;
        MultiObj.var_min = [-4, -4];
        MultiObj.var_max = [4, 4];
        load('Viennet2.mat');
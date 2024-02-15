clc;
clear all;
close all;


% Multi-Objective algorithm
%Alg = 'MOFDA';
%Alg = 'MOPSO';
%Alg = 'GA-MUL';

% Multi-objective function
MultiObjFnc = 'Schaffer';
%MultiObjFnc = 'Kursawe';
%MultiObjFnc = 'Poloni';
%MultiObjFnc = 'Viennet2';
%MultiObjFnc = 'Viennet3';
%MultiObjFnc = 'ZDT1';
%MultiObjFnc = 'ZDT2';
%MultiObjFnc = 'ZDT3';
MultiObjFnc = 'ZDT6'; 
%all;

results1 = [];
results2 = [];
results3 = [];
        
for g = 1:1

    switch Alg
    case 'MOFDA'         % Multi-obj Flow Directional Algorithm
        main;
        results1(g) = Convergence_curve(1,Max_iteration);
        results2(g) = Convergence_curve(2,Max_iteration);
        results3(g) = Convergence_curve(3,Max_iteration);
        pause(30)

    case 'MOPSO'         % Multi-obj Particle Swarm Optimizer
        example;
        results1(g) = REP.pos_fit(indexes(1),1);
        results2(g) = REP.pos_fit(indexes(1),2);
        results3(g) = REP.pos_fit(indexes(1),3);
    case 'GA-MUL'        %Multi-obj Genetic Algorithm
        GAM;
        results1(g) = Fval(indexes(1),1);
        results2(g) = Fval(indexes(1),2);
        results3(g) = Fval(indexes(1),3);
 
    end
end
display(results1.');
display(results2.');
display(results3.');
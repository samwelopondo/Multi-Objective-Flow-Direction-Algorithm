% multi_area_grid - multiobjective function for a multiarea grid
% Inputs:
%     * x: vector of generation dispatch for each generating station
%     * area_loads: a vector of length num_areas containing the loads in each area
%     * gen_costs: vector of size num_areas containing the costs of generating power in each area.
%     * transmission_loss_coeff: a matrix of size num_areas x num_areas containing the transmission loss coefficiets between each pair of areas
% Outputs:
%     * f1: the cost of the grid
%     * f2: the losses in the grid
function [f1, f2] = multi_area_grid(x,params)
    %Load grid parameters
    area_loads = params.area_loads;
    gen_costs = params.gen_costs;
    loss_coefficient = params.loss_coefficient;
    f1 = inf;
    f2 = inf;

    %Ensure generation dispatch meets grid load
    if sum(x) >= sum(area_loads)
        % Calculate the total generation cost for each area
        gen_costs = x.*gen_costs;
        total_generation_cost = sum(gen_costs);

        %calculate available remaining supply and demand in each area
        demand = zeros(1,length(area_loads));
        supply = zeros(1,length(x));
        diff = x - area_loads;
        for i=1:length(diff)
            if diff(i) > 0
                supply(i) = diff(i);
            end
            if diff(i) < 0
                demand(i) = abs(diff(i));
            end
        end
            
        % Initialize the power flow matrix
        power_flow = zeros(length(demand), length(supply));

        % Offset the demand from all available power
        for i = 1:length(demand)
            for j = 1:length(supply)
                if supply(j) >= demand(i)
                    power_flow(i, j) = demand(i);
                    supply(j) = supply(j) - demand(i);
                    demand(i) = 0;
                else
                    power_flow(i, j) = supply(j);
                    demand(i) = demand(i) - supply(j);
                    supply(j) = 0;
                end
            end
        end
        % Compute the total loss
        total_loss = 0;
        for i = 1:length(demand)
            for j = 1:length(supply)
            total_loss = total_loss + power_flow(i, j) * loss_coefficient(i, j);
            end
        end
        % Return the two objective values
        f1 = total_generation_cost;
        f2 = total_loss;
end
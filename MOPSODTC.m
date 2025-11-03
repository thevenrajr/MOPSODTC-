function [archive, execution_time] = MOPSODTC(problem, MaxIt, N, archive_size, gamma, c, w, c1, c2, pCrossover, pMutation, eta_c, eta_m)
    tic; % Start timer

    %% Problem Details
    CostFunction = problem.CostFunction;
    nVar = problem.nVar;
    VarMin = problem.VarMin;
    VarMax = problem.VarMax;
    nObj = problem.nObj;

    %% Initialization
    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.pbest.Position = [];
    empty_particle.pbest.Cost = [];
    
    pop = repmat(empty_particle, N, 1);
    
    for i = 1:N
        pop(i).Position = unifrnd(VarMin, VarMax, [1, nVar]);
        pop(i).Velocity = zeros(1, nVar);
        pop(i).Cost = CostFunction(pop(i).Position);
        pop(i).pbest.Position = pop(i).Position;
        pop(i).pbest.Cost = pop(i).Cost;
    end
    
    % Initialize Archive
    archive = UpdateArchive([], pop, archive_size);

    %% Main Loop
    for it = 1:MaxIt
        
        % === 1. Perturbation (SBX and PM) ===
        % Non-dominated sort to find ranks 1 and 2
        all_costs = vertcat(pop.Cost);
        [ranks, ~] = sort(all_costs, 1);
        % This is a simplified ranking, a full non-dominated sort is more robust
        % For simplicity, we just take a pool of good candidates
        [~, sorted_indices] = sort(sum(all_costs, 2));
        parent_pool_indices = sorted_indices(1:ceil(0.5*N)); % Select top 50% as parent pool

        num_offspring = round(pCrossover * numel(parent_pool_indices) / 2) * 2;
        offspring = repmat(empty_particle, num_offspring, 1);
        
        for k = 1:2:num_offspring
            p_indices = randsample(parent_pool_indices, 2);
            p1 = pop(p_indices(1));
            p2 = pop(p_indices(2));
            
            [o1_pos, o2_pos] = SBX(p1.Position, p2.Position, eta_c, VarMin, VarMax);
            
            o1_pos = PM(o1_pos, pMutation, eta_m, VarMin, VarMax);
            o2_pos = PM(o2_pos, pMutation, eta_m, VarMin, VarMax);
            
            offspring(k).Position = o1_pos;
            offspring(k+1).Position = o2_pos;
            
            offspring(k).Cost = CostFunction(offspring(k).Position);
            offspring(k+1).Cost = CostFunction(offspring(k+1).Position);
        end
        
        % Replace worst particles with offspring
        [~, worst_indices] = sort(sum(all_costs, 2), 'descend');
        pop(worst_indices(1:num_offspring)) = offspring;
        
        % === 2. Dynamic Population Sizing (Eqs. 4 & 5) ===
        NE = round(N - gamma - (N - 2*gamma) * (it/MaxIt));
        ND = N - NE; % The rest go to DS
        
        % === 3. Population Division ===
        % --- Diversity Metric for ES (Eq. 6) ---
        alpha = 4; % Crowding threshold
        all_costs = vertcat(pop.Cost);
        
        % Create adaptive grid
        min_costs = min(all_costs, [], 1);
        max_costs = max(all_costs, [], 1);
        grid_div = 10;
        grid_locs = zeros(N, nObj);
        for j = 1:nObj
            grid_locs(:,j) = floor((all_costs(:,j) - min_costs(j)) ./ ((max_costs(j) - min_costs(j) + eps) / grid_div));
        end
        
        diversity_metric = zeros(N, 1);
        for i = 1:N
            pn = sum(all(grid_locs == grid_locs(i,:), 2));
            if pn <= alpha
                diversity_metric(i) = 1 - pn/alpha;
            else
                diversity_metric(i) = 0;
            end
        end
        
        % --- Convergence Metric for DS (ASF, Eqs. 7 & 8) ---
        zmin = min(all_costs, [], 1); % Ideal point
        asf_metric = zeros(N, 1);
        for i = 1:N
            w_vec = (pop(i).Cost - zmin) ./ (sum(pop(i).Cost - zmin) + eps);
            w_vec(w_vec == 0) = 1e-6;
            asf_metric(i) = max((pop(i).Cost - zmin) ./ w_vec);
        end

        % Assign particles
        [~, sorted_diversity_idx] = sort(diversity_metric, 'descend');
        [~, sorted_convergence_idx] = sort(asf_metric, 'ascend');
        
        es_indices = sorted_diversity_idx(1:NE);
        ds_indices = sorted_convergence_idx(1:ND);
        
        % === 4. Update ES and DS Populations ===
        % Select gbest for all particles from archive using roulette wheel (density based)
        archive_costs = vertcat(archive.Cost);
        if size(archive_costs, 1) > 1
             % Inversely proportional to density. We use crowding distance as a proxy.
             dists = pdist2(archive_costs, archive_costs);
             dists(logical(eye(size(dists)))) = inf;
             sorted_dists = sort(dists, 1);
             crowding_dist = sum(sorted_dists(1,:), 1)'; % A simple density estimator
             selection_probs = crowding_dist / sum(crowding_dist);
        else
            selection_probs = 1;
        end
       
        % --- Update ES (Exploration Swarm) ---
        for i = es_indices'
            if ~isempty(archive)
                gbest_idx = randsample(1:numel(archive), 1, true, selection_probs);
                gbest = archive(gbest_idx);
            else
                gbest = pop(randsample(1:N,1)); % Fallback
            end

            r = rand(1, nVar);
            % Velocity update (Eq. 9)
            pop(i).Velocity = c * r .* (gbest.Position - pop(i).Position);
            % Position update (Eq. 3)
            pop(i).Position = pop(i).Position + pop(i).Velocity;
            
            pop(i).Position = max(pop(i).Position, VarMin);
            pop(i).Position = min(pop(i).Position, VarMax);
            
            pop(i).Cost = CostFunction(pop(i).Position);
        end
        
        % --- Update DS (Development Swarm) ---
        for i = ds_indices'
            if ~isempty(archive)
                gbest_idx = randsample(1:numel(archive), 1, true, selection_probs);
                gbest = archive(gbest_idx);
            else
                gbest = pop(randsample(1:N,1)); % Fallback
            end

            r1 = rand(1, nVar);
            r2 = rand(1, nVar);
            % Velocity update (Eq. 2)
            pop(i).Velocity = w*pop(i).Velocity ...
                            + c1*r1.*(pop(i).pbest.Position - pop(i).Position) ...
                            + c2*r2.*(gbest.Position - pop(i).Position);
            % Position update (Eq. 3)
            pop(i).Position = pop(i).Position + pop(i).Velocity;

            pop(i).Position = max(pop(i).Position, VarMin);
            pop(i).Position = min(pop(i).Position, VarMax);

            pop(i).Cost = CostFunction(pop(i).Position);
        end
        
        % === 5. Update pbest ===
        for i = 1:N
            if Dominates(pop(i).Cost, pop(i).pbest.Cost)
                pop(i).pbest.Position = pop(i).Position;
                pop(i).pbest.Cost = pop(i).Cost;
            elseif ~Dominates(pop(i).pbest.Cost, pop(i).Cost) && rand < 0.5
                pop(i).pbest.Position = pop(i).Position;
                pop(i).pbest.Cost = pop(i).Cost;
            end
        end

        % === 6. Update Archive ===
        archive = UpdateArchive(archive, pop, archive_size);
        
        % Display Iteration Information
        disp(['Iteration ' num2str(it) ': Archive Size = ' num2str(numel(archive)) ...
              ', ES Size = ' num2str(NE) ', DS Size = ' num2str(ND)]);
    end
    execution_time = toc; % End timer
end

function archive = UpdateArchive(archive, pop, archive_size)
    % Combine current population and archive
    combined = [archive; pop];
    
    % Find non-dominated solutions
    costs = vertcat(combined.Cost);
    is_dominated = false(size(costs, 1), 1);
    for i = 1:size(costs, 1)
        for j = 1:size(costs, 1)
            if i ~= j && Dominates(costs(j,:), costs(i,:))
                is_dominated(i) = true;
                break;
            end
        end
    end
    
    archive = combined(~is_dominated);
    
    % Prune archive if it exceeds the size limit
    while numel(archive) > archive_size
        archive_costs = vertcat(archive.Cost);
        n_archive = numel(archive);
        
        % === Adaptive Grid (Eq. 10) ===
        beta = rand();
        Kd = 50 + round(sin(beta * pi / 2) * 10);
        
        min_costs = min(archive_costs, [], 1);
        max_costs = max(archive_costs, [], 1);
        
        grid_locs = zeros(n_archive, size(archive_costs, 2));
        for j = 1:size(archive_costs, 2)
            grid_locs(:,j) = floor((archive_costs(:,j) - min_costs(j)) ./ ((max_costs(j) - min_costs(j) + eps) / Kd));
        end
        
        % Find densest grid cell
        [unique_cells, ~, ic] = unique(grid_locs, 'rows');
        counts = accumarray(ic, 1);
        [~, max_count_idx] = max(counts);
        densest_cell_members_idx = find(ic == max_count_idx);
        
        % === Shift-based Mean (SM) Density Estimation (Eqs. 11-13) ===
        if numel(densest_cell_members_idx) > 1
            cell_costs = archive_costs(densest_cell_members_idx,:);
            n_cell = size(cell_costs, 1);
            sm_density = zeros(n_cell, 1);
            
            for i = 1:n_cell
                dist_sum = 0;
                for q = 1:n_cell
                    if i ~= q
                        % F_iq is a complex comparison, for simplicity we use Euclidean distance
                        dist_sum = dist_sum + norm(cell_costs(i,:) - cell_costs(q,:));
                    end
                end
                sm_density(i) = dist_sum / (n_cell - 1);
            end
            
            [~, min_sm_idx] = min(sm_density); % Smallest SM value means densest region
            prune_idx = densest_cell_members_idx(min_sm_idx);
        else % If only one member, just pick it
            prune_idx = densest_cell_members_idx(1);
        end
        
        archive(prune_idx) = [];
    end
end

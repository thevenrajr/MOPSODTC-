function d = Dominates(x, y)
    % Checks if solution x dominates solution y
    d = all(x <= y) && any(x < y);
end

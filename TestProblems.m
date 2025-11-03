function p = TestProblems(name)
    switch name
        case 'ZDT1'
            p.CostFunction = @(x) ZDT1(x);
            p.nVar = 30;
            p.nObj = 2;
            p.VarMin = 0;
            p.VarMax = 1;
            load('ZDT1_PF.mat', 'PF');
            p.TruePF = PF;
        case 'ZDT2'
            p.CostFunction = @(x) ZDT2(x);
            p.nVar = 30;
            p.nObj = 2;
            p.VarMin = 0;
            p.VarMax = 1;
            load('ZDT2_PF.mat', 'PF');
            p.TruePF = PF;
        case 'ZDT3'
            p.CostFunction = @(x) ZDT3(x);
            p.nVar = 30;
            p.nObj = 2;
            p.VarMin = 0;
            p.VarMax = 1;
            load('ZDT3_PF.mat', 'PF');
            p.TruePF = PF;
        case 'ZDT4'
            p.CostFunction = @(x) ZDT4(x);
            p.nVar = 10;
            p.nObj = 2;
            p.VarMin = [0, -5*ones(1,9)];
            p.VarMax = [1, 5*ones(1,9)];
            load('ZDT4_PF.mat', 'PF');
            p.TruePF = PF;
        case 'ZDT6'
            p.CostFunction = @(x) ZDT6(x);
            p.nVar = 10;
            p.nObj = 2;
            p.VarMin = 0;
            p.VarMax = 1;
            load('ZDT6_PF.mat', 'PF');
            p.TruePF = PF;
    end
end

% ZDT Test Functions
function f = ZDT1(x)
    f1 = x(1);
    g = 1 + 9*sum(x(2:end))/(numel(x)-1);
    f2 = g * (1 - sqrt(f1/g));
    f = [f1, f2];
end

function f = ZDT2(x)
    f1 = x(1);
    g = 1 + 9*sum(x(2:end))/(numel(x)-1);
    f2 = g * (1 - (f1/g)^2);
    f = [f1, f2];
end

function f = ZDT3(x)
    f1 = x(1);
    g = 1 + 9*sum(x(2:end))/(numel(x)-1);
    f2 = g * (1 - sqrt(f1/g) - (f1/g)*sin(10*pi*f1));
    f = [f1, f2];
end

function f = ZDT4(x)
    f1 = x(1);
    g = 1 + 10*(numel(x)-1) + sum(x(2:end).^2 - 10*cos(4*pi*x(2:end)));
    f2 = g * (1 - sqrt(f1/g));
    f = [f1, f2];
end

function f = ZDT6(x)
    f1 = 1 - exp(-4*x(1))*sin(6*pi*x(1))^6;
    g = 1 + 9*(sum(x(2:end))/(numel(x)-1))^0.25;
    f2 = g * (1 - (f1/g)^2);
    f = [f1, f2];
end

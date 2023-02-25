function simOut = extract_sim_data(equ, t, state)
    %% extract_sim_data
    %
    % input: [equ, t, state]
    % equ       @(t,state)      simulation equation
    % t         1D double       time vector
    % state     2D double       state
    %
    % output: simOut
    % simOut    sutructure      simulation data out
    %
    % update:2022/02/10
    % Author:Hóng Jyùn Yaò
    
    %% --------------------------------------
    lenT = length(t);
    for i = 1:lenT
        [~,simOutTemp(i)] = equ(t(i),state(i,:)');
    end
    
    simOut = ArrayStruc2StrucArray(simOutTemp);
end

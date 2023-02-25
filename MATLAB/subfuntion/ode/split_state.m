function [X,stateOut] = split_state(stateIn,length)
    %% split_state
    %
    % input: (stateIn,length)
    % stateIn   1D double
    % length    double      length of X
    %
    % output: [X,stateOut]
    % X         1D double   system state
    % stateOut  1D double
    %
    % update:2022/02/11
    % Author:Hóng Jyùn Yaò
    
    %% --------------------------------------
    X = stateIn(1:length);
    stateIn(1:length) = [];
    if size(stateIn,2)
        stateOut = stateIn;
    else
        stateOut = stateIn';
    end
end


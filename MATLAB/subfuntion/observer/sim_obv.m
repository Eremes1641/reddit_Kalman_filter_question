function [tout,Xout,Yout] = sim_obv(observer,u,Y,IC)
    %% sim_obv
    % simulation discrete observer
    %
    % input: (observer,u,Y,IC)
    % Kalman    sys         discrete observer
    % u         array       input
    % Y         array       measurement
    % IC        array       initial condition
    %
    % output: [tout,Xout,Yout]
    % tout      array       time array
    % Xout      array       estimated state X out
    % Yout      array       estimated measurement Y out
    %
    % update:2021/12/29
    % Author:Hóng Jyùn Yaò
    
    %% --------------------------------------
    arguments
        observer
        u
        Y
        IC = zeros(length(observer.A),1)
    end
    
    %% --------------------------------------
    [Aobv,Bobv,Cobv,Dobv] = ssdata(observer);
    Xobv(:,1) = IC;
    Fs = 1/observer.Ts;
    tout = (0:length(Y)-1)'/Fs;

    % run
    for k = 1:length(tout)
        % substitude
        if ~isempty(u)
            uobv = [u(k,:) Y(k,:)]';
        else
            uobv = Y(k,:)';
        end
        % observer
        Xobv(:,k+1) = Aobv*Xobv(:,k) + Bobv*uobv;
        Yobv(:,k) = Cobv*Xobv(:,k) + Dobv*uobv;
    end
    
    %% return
    numY = size(Y,2);
    Yout = Yobv(1:numY,:)';
    Xout = Yobv(1+numY:end,:)';
end

%     Xout = Yobv(1+numMeasure:end,:)';
%     Yout = Yobv';
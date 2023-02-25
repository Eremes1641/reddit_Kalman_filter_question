function y = LPF(x,bandwidth,dt)
    if bandwidth > 1/dt
        bandwidth = 1/dt;
    end
    y = zeros(length(x),1);
    for i = 2:length(x)
        vel_est_temp = ( x(i) - x(i-1) )/dt;
        y(i) = (1-bandwidth*dt)*y(i-1) + bandwidth*dt*vel_est_temp;
    end
end
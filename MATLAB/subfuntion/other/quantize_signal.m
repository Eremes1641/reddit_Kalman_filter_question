function y = quantize_signal(x,resolution)
    y = zeros(length(x),1);
    for i = 1:length(x)
        y(i) = fix(x(i) / resolution) * resolution;
    end
end
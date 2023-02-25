function fig_narrow(fig,Wscale,Hscale)
    %% fig_narrow
    %
    % input: (fig) / (fig,scale)
    % fig       figure      figure in
    % scale     double      scale
    %
    %
    % update:2022/02/04
    % Author:Hóng Jyùn Yaò
    
    %% --------------------------------------
    if nargin < 2
        Wscale = 1/5*3;    %2.7
        Hscale = 1/5*5;
    end
    fig.Position(3) = fig.Position(3)*Wscale;
    fig.Position(4) = fig.Position(4)*Hscale;
end


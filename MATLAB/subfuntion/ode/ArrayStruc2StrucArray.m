function StrucArray = ArrayStruc2StrucArray(ArrayStruc)
    %% ArrayStruc2StrucArray
    % convert array of structure to structure of array
    %
    % input: ArrayStruc
    % ArrayStruc    array       array of structure
    %
    % output: StrucArray
    % StrucArray    structure   structure of array
    %
    % update:2022/02/01
    % Author:Hóng Jyùn Yaò
    
    %% --------------------------------------
    fieldName = fields(ArrayStruc(1));
    for i = 1:length(fieldName)
        StrucArray.(fieldName{i}) = [ArrayStruc.(fieldName{i})]';
    end
end


function [gradient] = numeric_pderiv1(numeric_pderivs,estdata,mdata)
%x is the state, n is the index of the parameter that the partial is
%respect to
%outputs partial derivative of cost with respect to some parameter, used to
%fill gradient matrix

    max_j = size(numeric_pderivs,2);
    gradient = zeros(max_j, 1);
    for j = 1:max_j
        d = 0;
        for i = 1:length(mdata)
            d = d + numeric_pderivs(i,j)*(mdata(i)-estdata(i));
        end
        d = d * -1;
        gradient(j) = d;
    end
end

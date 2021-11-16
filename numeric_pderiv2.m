function [d] = numeric_pderiv2(numeric_pderivs,numtrials,datapoints)

%x is the state, n is the index of the parameter that the partial is
%respect to
%outputs full second gradient of cost

    max_n = size(numeric_pderivs,2);
    d = zeros(max_n, max_n);
    
     for j = 1:(numtrials*datapoints)
        grad_p = zeros(2, max_n);
         for c = 1:max_n
            %p_deriv = numeric_pderiv(estdata,x,c,tdata,mdata,max_t);
            grad_p(1, c) = numeric_pderivs(2*j-1,c);
            grad_p(2, c) = numeric_pderivs(2*j,c);
        end
        
        d = d + transpose(grad_p)*grad_p;
    end
    

%       grad_p = zeros(2, 5);
%         for c = 1:5
%             p_deriv = numeric_pderiv(x,c,tdata,mdata,max_t);
%              grad_p(1, c) = p_deriv(1,:);
%             grad_p(2, c) = p_deriv(2,:);
%          end
%          
%         %d = d + transpose(grad_p);
%     grad_p_y2 = zeros(2, 5);
%     for c = 1:5
%         p_deriv = numeric_pderiv(x,c,tdata,mdata,max_t);
%         grad_p_y2(1, c) = p_deriv(3,:);
%         grad_p_y2(2, c) = p_deriv(4,:);
%     end
%          d = transpose(grad_p)*grad_p + transpose(grad_p_y2)*grad_p_y2;
end
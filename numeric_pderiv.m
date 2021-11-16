function [p_deriv] = numeric_pderiv(estdata,x,n,tdata,mdata,max_t,model)
%calculates just the partial derivative of estdata with respect to x(n)

    dB = 0.001;
    
    %used to find change in cost over small nudge to nth parameter
    x(n) = x(n) + dB;
            
    %NEED TO CHANGE TO DEVAL INSTEAD OF INTERP1?
    [t_est, s_est] = ode45(model, [0, max_t], x);
    xdata_est = interp1(t_est, s_est(:,1), tdata); 
    ydata_est = interp1(t_est, s_est(:,2), tdata);
    
    %sol = ode45(@Equations, [0, max_t], x);
    %all_data = deval(sol, tdata);
    %xdata_est = transpose(all_data(1,:));
    %ydata_est = transpose(all_data(2,:));

    estdata2 = zeros(size(mdata));
    estdata2(1:2:end) = xdata_est;
    estdata2(2:2:end) = ydata_est;

    p_deriv = (estdata2-estdata)/(dB);
end

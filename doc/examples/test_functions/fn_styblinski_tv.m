function f = fn_styblinski_tv( x_e )
    persistent ctr;
    % Styblinski-Tang function
    % https://en.wikipedia.org/wiki/Test_functions_for_optimization
    % f(x*) = -39.16616 * n at x* = ( -2.903534, -2.903534, ... )   

    if istable( x_e )
        x = table2array( x_e );
    else
        x = x_e;
    end

    if isempty(ctr)
        ctr = 0;
    else
        ctr = ctr+1;
    end

    if length(x) == 2
        theta = ctr/25;
        R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
        x = R*x;
    end

    f = sum( x .^ 4 - 16 * x .^ 2 + 5 * x ) / 2;
end
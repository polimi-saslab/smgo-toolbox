function f = fn_parab( x_e )
    % simple parabolic function 

    if istable( x_e )
        x = table2array( x_e );
    else
        x = x_e;
    end

    f = sum( ( x - [ -5 -5 ]' ) .^ 2 ) - 64;
end
function f = fn_plane( x_e )
    % simple parabolic function 

    if istable( x_e )
        x = table2array( x_e );
    else
        x = x_e;
    end

    f = x(1) - x(2);
end
function [ f, g ] = fn_g10b( x )
    x = table2array( x );
    f = x( 1 ) + x( 2 ) + x( 3 );
    g1 = -1 + 0.0025 * ( x( 4 ) + x( 6 ) );
    g2 = -1 + 0.0025 * ( x( 5 ) + x( 7 ) - x( 4 ) );
    g3 = -1 + 0.01 * ( x( 8 ) - x( 5 ) );
    g4 = -x( 1 ) * x( 6 ) + 833.33252 * x( 4 ) + 100 * x( 1 ) - 83333.333;
    g5 = -x( 2 ) * x( 7 ) + 1250 * x( 5 ) + x( 2 ) * x( 4 ) - 1250 * x( 4 );
    g6 = -x( 3 ) * x( 8 ) + 1250000 + x( 3 ) * x( 5 ) - 2500 * x( 5 );
    g = [ g1 g2 g3 g4 g5 g6 ];
end
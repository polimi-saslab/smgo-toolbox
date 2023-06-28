function [ f, g ] = fn_g24b( x )
    x = table2array( x );
    f = -x( 1 ) - x( 2 );
    g1 = -2 * x( 1 ) ^ 4 + 8 * x( 1 ) ^ 3 - 8 * x( 1 ) ^ 2 + x( 2 ) - 2;
    g2 = -4 * x( 1 ) ^ 4 + 32 * x( 1 ) ^ 3 - 88 * x( 1 ) ^ 2 + 96 * x( 1 ) + x( 2 ) - 36;
    g = [ g1 g2 ];
end
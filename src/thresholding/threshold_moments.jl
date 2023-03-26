function threshold_moments( img ) 

    nh, nb, no = normhisto( img );

    # Zeroth order moment of the normalized histogram
    m0 = 1.0;
   
    # Calculate the first, second, and third order moments of the histogram
    m1 = m2 = m3 = 0.0;
    for ih in 1:length(nh)
      m1 += ih * nh[ih];
      m2 += ih * ih * nh[ih];
      m3 += ih * ih * ih * nh[ih];
    end
   
    #=  First 4 moments of the gray-level image should match the first 4 moments
        of the target binary image. This leads to 4 equalities whose solutions 
        are given in the Appendix of Ref. 1 
    =#

    cd = m0 * m2 - m1 * m1;
    c0 = ( -m2 * m2 + m1 * m3 ) / cd;
    c1 = ( m0 * -m3 + m2 * m1 ) / cd;
    z0 = 0.5 * ( -c1 - sqrt( c1 * c1 - 4.0 * c0 ) );
    z1 = 0.5 * ( -c1 + sqrt( c1 * c1 - 4.0 * c0 ) );

    # Fraction of the object pixels in the target binary image
    p0 = ( z1 - m1 ) / ( z1 - z0 );

    #=  The threshold is the gray-level closest  
        to the p0-tile of the normalized histogram 
    =#

    threshold = 0;
    sum = 0.0;
    for ih in 1:length(nh)
        sum += nh[ih];
        if ( sum > p0 )
            threshold = nb[ih];
            break;
        end
    end
end
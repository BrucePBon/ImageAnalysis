function mean_extrema_pad( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fmin=1, fmax=1, ovp=(0,0) ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    return mean_extrema_pad!( img, zeros(Float64, padsize .+ 1), 
                              rad, th_hard, 
                              zeros(Bool,size(img)), zeros(Bool,size(img)),
                              fmin=fmin, fmax=fmax, ovp=ovp );
end

function mean_extrema_pad!( img::Array{T,2}, intA::Array{<:AbstractFloat,2},
                            rad::Tuple{Int64,Int64}, th_hard, 
                            minima::Array{Bool,2}, maxima::Array{Bool,2};
                            fmin=1, fmax=1, ovp=0 ) where {T<:Real}
    # integral array
    @inbounds @simd for idx in 1:length(intA)
        intA[idx] = 0.0
    end
    integralArray_pad!( img, intA, 3 .* rad .+ 1 )   
    
    # convenient quantities
    lows  = (  1,  1  ) .+ ( 3 .* rad .+ 1 );
    highs = size( img ) .+ ( 3 .* rad .+ 1 ); 
    npix  = prod(  (  2 .* rad .+ 1  )  ); 
    vert  = ( 2*rad[1] + 1 - ovp[1], 0  ); 
    horz  = ( 0, 2*rad[2] + 1  - ovp[2] );
    
    # these will accomodate the means and signs of the 9x9 mean neighbourhood
    sums9 = zeros( Float64, 3, 3 ); 
    signs = zeros(   Int64, 3, 3 ); 
    
    @inbounds for x in lows[2]:highs[2], y in lows[1]:highs[1]
        
        tl = ( y,x ) .- rad .- 1; 
        br = ( y,x ) .+ rad; 
        s0 = ImageAnalysis.integralArea( intA, tl, br );   
             
        ( s0/npix <= th_hard ) && ( continue; ) # s0/npix == mean of the central square

        idx = 1;
        for xoff in -1:1, yoff in -1:1
            off = vert .* yoff .+ horz .* xoff; 
            sums9[idx] = ImageAnalysis.integralArea( intA, tl .+ off, br .+ off ); 
            signs[idx] = sign( s0 - sums9[idx] ); 
            idx += 1;   
        end

        #=
        sums9[1,1] = ImageAnalysis.integralArea( intA, tl .- vert .- horz, br .- vert .- horz );
        sums9[2,1] = ImageAnalysis.integralArea( intA, tl         .- horz, br         .- horz );
        sums9[3,1] = ImageAnalysis.integralArea( intA, tl .+ vert .- horz, br .+ vert .- horz ); 
        sums9[1,2] = ImageAnalysis.integralArea( intA, tl .- vert        , br .- vert         );
        sums9[2,2] = ImageAnalysis.integralArea( intA, tl                , br                 );
        sums9[3,2] = ImageAnalysis.integralArea( intA, tl .+ vert        , br .+ vert         ); 
        sums9[1,3] = ImageAnalysis.integralArea( intA, tl .- vert .+ horz, br .- vert .+ horz );
        sums9[2,3] = ImageAnalysis.integralArea( intA, tl         .+ horz, br         .+ horz );
        sums9[3,3] = ImageAnalysis.integralArea( intA, tl .+ vert .+ horz, br .+ vert .+ horz );
        =#

        score_min = 0
        score_max = 0
        for idx in 1:4
            score_min += ( signs[idx] == -1 && signs[end-idx+1] == -1 )   
            score_max += ( signs[idx] ==  1 && signs[end-idx+1] ==  1 )
        end
        
        minima[y-3*rad[1]-1,x-3*rad[2]-1] = score_min > fmin; 
        maxima[y-3*rad[1]-1,x-3*rad[2]-1] = score_max > fmax; 
    end

    return minima, maxima
end





# sum ( x - mean )^2 = sum( x^2 + mean^2 - 2*x*mean ) = sum x^2 + N*mean^2 - 2*mean*sum x
#                                                     = sum x^2 + mean*sum x - 2*mean*sum x
#                                                     = sum x^2 - mean*sum x

function mean_std_pad( img::Array{<:Real,2}, rad::Tuple{Int64,Int64}; typ=Float64 )

    padsize = size( img ) .+ 2 .* rad;
    return mean_std_pad!( img, zeros( typ, padsize .+ 1), zeros( typ, padsize .+ 1),
                          rad, zeros( typ, size( img ) ) );
end

function mean_std_pad!( img::Array{<:Real,2}, intA::Array{T,2}, intA2::Array{T,2},
                        rad::Tuple{Int64,Int64}, stds::Array{T,2} ) where {T<:AbstractFloat}

    # integral array
    @inbounds @simd for idx in 1:length(intA)
         intA[idx] = 0.0
        intA2[idx] = 0.0
    end
    integralArray_pad!( img,  intA, rad, fun=(x)->(T(x)) ); 
    integralArray_pad!( img, intA2, rad, fun=(x)->(T(x)^2) );   
        
    # convenience variables
    lows  = (  1,  1  ) .+ rad;
    highs = size( img ) .+ rad;
    npix  = prod( ( 2 .* rad .+ 1 ) );
    
    @inbounds for x in lows[2]:highs[2], y in lows[1]:highs[1]
        
        tl  = ( y, x ) .- rad .- 1; 
        br  = ( y, x ) .+ rad; 
        s0  = ImageAnalysis.integralArea(  intA, tl, br ); 
        s02 = ImageAnalysis.integralArea( intA2, tl, br );
        num = s02 - s0*s0/npix;
        
        stds[y-rad[1],x-rad[2]] = sqrt( num / npix )
    end

    return stds
end




function mean_thresh_pad( img::Array{T,2}, rad::Tuple{Int64,Int64}; fun=(x,y)->(x>y) ) where {T<:Real}

    padsize = size( img ) .+ 2 .*rad;
    return mean_thresh_pad!( img, zeros(Float64, padsize .+ 1),
                             rad, zeros(Bool,size(img)), fun=fun );
end

function mean_thresh_pad!(  img::Array{T,2}, 
                            intA::Array{<:AbstractFloat,2},
                            rad::Tuple{Int64,Int64},
                            filt::Array{Bool,2}; fun=(x,y)->(x>y) ) where {T<:Real}
    # integral array
    @inbounds @simd for idx in 1:length(intA)
        intA[idx] = 0.0
    end
    integralArray_pad!( img, intA, rad )   
    
    # convenient quantities
    lows  = (  1,  1  ) .+ rad;
    highs = size( img ) .+ rad;
    npix  = prod( ( 2 .* rad .+ 1 ) ); 
    
    @inbounds for x in lows[2]:highs[2], y in lows[1]:highs[1]
        
        tl = ( y,x ) .- rad .- 1; 
        br = ( y,x ) .+ rad; 
        s0 = ImageAnalysis.integralArea( intA, tl, br );    
        mean = s0/npix;     
        
        filt[y-rad[1],x-rad[2]] = fun( img[ y-rad[1],x-rad[2] ], mean ) ; 
    end

    return filt
end




function mean_extrema_pad2( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fmin=1, fmax=1, ovp=(0,0) ) where {T<:Real}

    return mean_extrema_pad2!( img, zeros(Float64, size(img) .+ 1), rad, th_hard, zeros(Bool,size(img)), zeros(Bool,size(img)),
                              fmin=fmin, fmax=fmax, ovp=ovp );
end

function mean_extrema_pad2!( img::Array{T,2}, intA::Array{<:AbstractFloat,2},
                             rad::Tuple{Int64,Int64}, th_hard, 
                             minima::Array{Bool,2}, maxima::Array{Bool,2};
                             fmin=1, fmax=1, ovp=0 ) where {T<:Real}

    # integral array
    @inbounds @simd for idx in 1:length(intA)
    intA[idx] = 0.0
    end
    integralArray!( img, intA )   

    if th_hard == nothing
    th_hard = intA[ end, end ] / length(img); # == mean of the image
    end

    # convenient quantities
    npix  = prod( ( 2 .* rad .+ 1 ) ); 
    vert  = ( 2*rad[1] + 1 - ovp[1], 0 ); 
    horz  = ( 0, 2*rad[2] + 1  - ovp[2]);
    offs  = ( -1 .* (vert .+ horz), -1 .* horz, vert .- horz, -1 .* vert, (0,0), vert, horz .- vert, horz, horz .+ vert )

    # these will accomodate the means and signs of the 9x9 mean neighbourhood
    sums9 = zeros( Float64, 3, 3 ); 
    signs = zeros(   Int64, 3, 3 ); 

    @inbounds for x in 1:size(img,2), y in 1:size(img,1)

    tl = max.( 1, min.( size(img), ( y,x ) .- rad .- 1 ) ); 
    br = max.( 1, min.( size(img), ( y,x ) .+ rad ) ); 

    s0 = ImageAnalysis.integralArea( intA, tl, br );        
    ( s0/npix <= th_hard ) && ( continue; ) # s0/npix == mean of the central square

    cont = 1; 
    for off in offs
        tl_ = max.( 1, min.( size(img), tl .+ off ) ); 
        br_ = max.( 1, min.( size(img), br .+ off ) ); 
        sums9[cont] = ImageAnalysis.integralArea( intA, tl_, br_ ); 
        cont += 1; 
    end

    @inbounds @simd for idx in 1:9
        signs[idx] = sign( s0 .- sums9[idx] ); 
    end

    score_min = 0
    score_max = 0
    for idx in 1:4
        score_min += ( signs[idx] == -1 && signs[end-idx+1] == -1 )   
        score_max += ( signs[idx] ==  1 && signs[end-idx+1] ==  1 )
    end

    minima[y,x] = score_min > fmin; 
    maxima[y,x] = score_max > fmax; 
    end

    return minima, maxima
end

function std_extrema_pad_filt( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fmin=1, fmax=1, m=1, fun=(x)->(x) ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    return std_extrema_pad_filt!( img,  
                                  zeros(Float64, padsize .+ 1), zeros(Float64, padsize .+ 1), 
                                  zeros(Float64, padsize .+ 1), zeros(Float64, padsize .+ 1), zeros(Float64, padsize .+ 1),
                                  rad, th_hard, zeros(Float32,size(img)), 
                                  fmin=fmin, fmax=fmax, m=m, fun=fun );
end

# sum( ( x - mean )^2 - ( xf - mean )^2 ) = sum( x^2 + mean^2 - 2*x*mean ) - sum( xf^2 + mean^2 - 2*xf*mean )

function std_extrema_pad_filt!( img::Array{T,2}, 
                                intA::Array{<:AbstractFloat,2}, intA2::Array{<:AbstractFloat,2}, 
                                intF::Array{<:AbstractFloat,2}, intF2::Array{<:AbstractFloat,2}, Ns::Array{<:AbstractFloat,2},
                                rad::Tuple{Int64,Int64}, th_hard, stds::Array{Float32,2};
                                fmin=1, fmax=1, m=1, fun=(x)->(x) ) where {T<:Real}
    # integral array
    @inbounds @simd for idx in 1:length(intA)
         intA[idx] = 0.0
        intA2[idx] = 0.0
    end

    integralArray_pad!( img,  intA, 3 .* rad .+ 1 )  
    integralArray_pad!( img, intA2, 3 .* rad .+ 1, fun=(x)->(x*x) ) 
    integralArray_pad!( img,  intF, 3 .* rad .+ 1, fun=fun )  
    integralArray_pad!( img, intF2, 3 .* rad .+ 1, fun=(x)->(fun(x)^2) )  
    integralArray_pad!( img,    Ns, 3 .* rad .+ 1, fun=(x)->(fun(x)>0) );  


    if th_hard == nothing
        th_hard = ( intA[ end, end ] - intF[ end, end ] )/ length(img); # == mean of the image
    end
    
    # convenient quantities
    lows  = (  1,  1  ) .+ ( 3 .* rad .+ 1 );
    highs = size( img ) .+ ( 3 .* rad .+ 1 ); 
    npix  = prod( ( 2 .* rad .+ 1 ) ); 
    vert  = ( 2*rad[1] + 1, 0 ); 
    horz  = ( 0, 2*rad[2] + 1 );
    
    @inbounds for x in lows[2]:highs[2], y in lows[1]:highs[1]
        
        tl   = ( y,x ) .- rad .- 1; 
        br   = ( y,x ) .+ rad; 

        s0   = ImageAnalysis.integralArea(  intA, tl, br );  
        s02  = ImageAnalysis.integralArea( intA2, tl, br );
        sF   = ImageAnalysis.integralArea(  intF, tl, br );  
        sF2  = ImageAnalysis.integralArea( intF2, tl, br );
        nF   = ImageAnalysis.integralArea(    Ns, tl, br );

        mean = ( s0 - sF )/( npix - nF ); 
        std  = ( s02 + npix*mean*mean - 2*s0*mean ) - ( sF2 + nF*mean*mean - 2*sF*mean ); 

        stds[y-3*rad[1]-1,x-3*rad[2]-1] = sqrt( abs(std) / npix ) # img[y-3*rad[1]-1,x-3*rad[2]-1] - mean # sqrt( std / npix ); 
    end

    return stds
end

function mean_thresh_pad_filt( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fun=(x)->(x) ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    return mean_thresh_pad_filt!( img, 
                                  zeros(Float64, padsize .+ 1),
                                  zeros(Float64, padsize .+ 1), zeros(Float64, padsize .+ 1),
                                  rad, th_hard, zeros(Bool,size(img)),
                                  fun=fun );
end

function mean_thresh_pad_filt!(  img::Array{T,2}, 
                                 intA::Array{<:AbstractFloat,2},
                                 intF::Array{<:AbstractFloat,2}, numF::Array{<:AbstractFloat,2},
                                 rad::Tuple{Int64,Int64}, th_hard, 
                                 filt::Array{Bool,2};
                                 fun=(x)->(x) ) where {T<:Real}
    # integral array
    @inbounds @simd for idx in 1:length(intA)
        intA[idx] = 0.0
    end
    integralArray_pad!( img, intA, 3 .* rad .+ 1 )   
    integralArray_pad!( img, intF, 3 .* rad .+ 1, fun=fun )  # fun is a condition, x = cond * x
    integralArray_pad!( img, numF, 3 .* rad .+ 1, fun=(x)->(fun(x)>0) ); # counts the num of pixels that fulfill the condition
    
    if th_hard == nothing
        th_hard = intA[ end, end ] / length(img); # == mean of the image
    end
    
    # convenient quantities
    lows  = (  1,  1  ) .+ ( 3 .* rad .+ 1 );
    highs = size( img ) .+ ( 3 .* rad .+ 1 ); 
    npix  = prod( ( 2 .* rad .+ 1 ) ); 
    
    @inbounds for x in lows[2]:highs[2], y in lows[1]:highs[1]
        
        tl = ( y,x ) .- rad .- 1; 
        br = ( y,x ) .+ rad; 
        s0 = ImageAnalysis.integralArea( intA, tl, br );  
        sF = ImageAnalysis.integralArea( intF, tl, br ); 
        nF = ImageAnalysis.integralArea( numF, tl, br );   
        mean = ( s0 - sF )/( npix - nF );     
        
        filt[y-3*rad[1]-1,x-3*rad[2]-1] = img[y-3*rad[1]-1,x-3*rad[2]-1] > mean; 
    end

    return filt
end




function mean_slope_pad( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fun=(x)->(x>1), ovp=(0,0) ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    return mean_slope_pad!( img, zeros(Float64, padsize .+ 1), rad, th_hard, zeros(Bool,size(img)),
                            fun=fun, ovp=ovp );
end

function mean_slope_pad!( img::Array{T,2}, intA::Array{<:AbstractFloat,2},
                          rad::Tuple{Int64,Int64}, th_hard, 
                          slopes::Array{Bool,2};
                          fun=(x)->(x>1), ovp=0 ) where {T<:Real}
    # integral array
    @inbounds @simd for idx in 1:length(intA)
        intA[idx] = 0.0
    end
    integralArray_pad!( img, intA, 3 .* rad .+ 1 )   
    
    if th_hard == nothing
        th_hard = intA[ end, end ] / length(img); # == mean of the image
    end
    
    # convenient quantities
    lows  = (  1,  1  ) .+ ( 3 .* rad .+ 1 );
    highs = size( img ) .+ ( 3 .* rad .+ 1 ); 
    npix  = prod( ( 2 .* rad .+ 1 ) ); 
    vert  = ( 2*rad[1] + 1 - ovp[1], 0 ); 
    horz  = ( 0, 2*rad[2] + 1  - ovp[2]);
    
    # these will accomodate the means and signs of the 9x9 mean neighbourhood
    sums9 = zeros( Float64, 3, 3 ); 
    signs = zeros(   Int64, 3, 3 ); 
    
    @inbounds for x in lows[2]:highs[2], y in lows[1]:highs[1]
        
        tl = ( y,x ) .- rad .- 1; 
        br = ( y,x ) .+ rad; 
        s0 = ImageAnalysis.integralArea( intA, tl, br );        
        ( s0/npix <= th_hard ) && ( continue; ) # s0/npix == mean of the central square

        sums9[1,1] = ImageAnalysis.integralArea( intA, tl .- vert .- horz, br .- vert .- horz );
        sums9[2,1] = ImageAnalysis.integralArea( intA, tl         .- horz, br         .- horz );
        sums9[3,1] = ImageAnalysis.integralArea( intA, tl .+ vert .- horz, br .+ vert .- horz ); 
        sums9[1,2] = ImageAnalysis.integralArea( intA, tl .- vert        , br .- vert         );
        sums9[2,2] = ImageAnalysis.integralArea( intA, tl                , br                 );
        sums9[3,2] = ImageAnalysis.integralArea( intA, tl .+ vert        , br .+ vert         ); 
        sums9[1,3] = ImageAnalysis.integralArea( intA, tl .- vert .+ horz, br .- vert .+ horz );
        sums9[2,3] = ImageAnalysis.integralArea( intA, tl         .+ horz, br         .+ horz );
        sums9[3,3] = ImageAnalysis.integralArea( intA, tl .+ vert .+ horz, br .+ vert .+ horz );

        @inbounds @simd for idx in 1:9
            signs[idx] = sign( s0 .- sums9[idx] ); 
        end

        score_slope = 0
        for idx in 1:4
            score_slope += ( signs[idx]*signs[end-idx+1] == -1 )   
        end
        
        slopes[y-3*rad[1]-1,x-3*rad[2]-1] = fun( score_slope ); 
    end

    return slopes
end





function mean_thresh_count( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fun=(x,y)->(x<y), f=1.0 ) where {T<:Real}

    padsize = size( img ) .+ 2 .* rad; #( 3 .* rad .+ 1 );
    return mean_thresh_count!( img, 
                               zeros(Float64, padsize .+ 1),
                               rad, th_hard, zeros(Int64,size(img)), fun=fun, f=f );
end

function mean_thresh_count!(  img::Array{T,2}, 
                              intA::Array{<:AbstractFloat,2},
                              rad::Tuple{Int64,Int64}, th_hard, 
                              counts::Array{Int64,2}; fun=(x,y)->(x<y), f=1.0 ) where {T<:Real}
    # Building integral array
    # 1-. Setting intA to zero (redudant when calling out-of-place).
    @inbounds @simd for idx in 1:length(intA)
        intA[idx] = 0.0
    end
    # 2-. Populating integral array
    integralArray_pad!( img, intA, rad ); #3 .* rad .+ 1 )   
    
    # Checking the threshold value
    if th_hard == nothing
        th_hard = intA[ end, end ] / length(img); # == mean of the image
    end
    
    # Convenient quantities
    npix   = prod( ( 2 .* rad .+ 1 ) ); 
    y0, x0 = (  1,  1  ) .+ rad; #( 3 .* rad .+ 1 ); 
    y1, x1 = size( img ) .+ rad; #( 3 .* rad .+ 1 ); 
    ry, rx = rad; 
    
    for x in x0:x1, y in y0:y1
        
        # top left (tl) and bottom right (br) of the central patch
        tl = ( y, x ) .- rad .- 1; 
        br = ( y, x ) .+ rad;      
        s0 = ImageAnalysis.integralArea( intA, tl, br );    
        local_mean = s0*f/npix; 

        if local_mean < th_hard
            counts[ y-y0+1, x-x0+1 ] = 1000000; 
        end
        
        # For each pixel in the central patch that falls inside the image (we added padding to the integral array)
        # add +1 if the pixel is above the local mean. 
        
        for col in max(tl[2]+1,x0):min(br[2],x1)
            for row in max(tl[1]+1,y0):min(br[1],y1)

                counts[ row-y0+1,col-x0+1 ] += fun( img[ row-y0+1,col-x0+1 ], local_mean )
        end end
        
        #filt[y-3*rad[1]-1,x-3*rad[2]-1] = img[ y-3*rad[1]-1,x-3*rad[2]-1 ] > local_mean*0.95 ; 
    end

    for idx in 1:length(counts)
        if counts[idx] > 1000000
            counts[idx] = -1
        end
    end

    return counts
end





function mean_laplacian_pad( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fmin=1, fmax=1, ovp=(0,0) ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    return mean_laplacian_pad!( img, zeros(Float64, padsize .+ 1), rad, th_hard, zeros(Float32,size(img)),
                              fmin=fmin, fmax=fmax, ovp=ovp );
end

function mean_laplacian_pad!( img::Array{T,2}, intA::Array{<:AbstractFloat,2},
                              rad::Tuple{Int64,Int64}, th_hard, 
                              laplacian::Array{Float32,2};
                              fmin=1, fmax=1, ovp=(0,0) ) where {T<:Real}
    # integral array
    @inbounds @simd for idx in 1:length(intA)
        intA[idx] = 0.0
    end
    integralArray_pad!( img, intA, 3 .* rad .+ 1 )   
    
    if th_hard == nothing
        th_hard = intA[ end, end ] / length(img); # == mean of the image
    end
    
    # convenient quantities
    lows  = (  1,  1  ) .+ ( 3 .* rad .+ 1 );
    highs = size( img ) .+ ( 3 .* rad .+ 1 ); 
    npix  = prod( ( 2 .* rad .+ 1 ) ); 
    vert  = ( 2*rad[1] + 1 - ovp[1], 0 ); 
    horz  = ( 0, 2*rad[2] + 1  - ovp[2]);

    sums9 = zeros( Float32, 3, 3 ); 
    
    @inbounds for x in lows[2]:highs[2], y in lows[1]:highs[1]
        
        tl = ( y,x ) .- rad .- 1; 
        br = ( y,x ) .+ rad; 
        s0 = ImageAnalysis.integralArea( intA, tl, br );        
        ( s0/npix <= th_hard ) && ( continue; ) # s0/npix == mean of the central square

        sums9[1,1] = ImageAnalysis.integralArea( intA, tl .- vert .- horz, br .- vert .- horz );
        sums9[2,1] = ImageAnalysis.integralArea( intA, tl         .- horz, br         .- horz );
        sums9[3,1] = ImageAnalysis.integralArea( intA, tl .+ vert .- horz, br .+ vert .- horz ); 
        sums9[1,2] = ImageAnalysis.integralArea( intA, tl .- vert        , br .- vert         );
        sums9[2,2] = ImageAnalysis.integralArea( intA, tl                , br                 );
        sums9[3,2] = ImageAnalysis.integralArea( intA, tl .+ vert        , br .+ vert         ); 
        sums9[1,3] = ImageAnalysis.integralArea( intA, tl .- vert .+ horz, br .- vert .+ horz );
        sums9[2,3] = ImageAnalysis.integralArea( intA, tl         .+ horz, br         .+ horz );
        sums9[3,3] = ImageAnalysis.integralArea( intA, tl .+ vert .+ horz, br .+ vert .+ horz );
        
        laplacian[y-3*rad[1]-1,x-3*rad[2]-1] = abs( 4*s0 - sums9[1,2] - sums9[2,1] - sums9[2,3] - sums9[3,2] );
    end

    return laplacian
end




function mean_laplacian_pad_filt( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fun=(x)->(x), ovp=(0,0) ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    return mean_laplacian_pad_filt!( img, 
                                  zeros(Float64, padsize .+ 1),
                                  zeros(Float64, padsize .+ 1), zeros(Float64, padsize .+ 1),
                                  rad, th_hard, zeros(Float32,size(img)),
                                  fun=fun, ovp=ovp );
end

function mean_laplacian_pad_filt!(  img::Array{T,2}, 
                                 intA::Array{<:AbstractFloat,2},
                                 intF::Array{<:AbstractFloat,2}, numF::Array{<:AbstractFloat,2},
                                 rad::Tuple{Int64,Int64}, th_hard, 
                                 laplacian::Array{Float32,2};
                                 fun=(x)->(x), ovp=(0,0) ) where {T<:Real}
    # integral array
    @inbounds @simd for idx in 1:length(intA)
        intA[idx] = 0.0
    end
    integralArray_pad!( img, intA, 3 .* rad .+ 1 )   
    integralArray_pad!( img, intF, 3 .* rad .+ 1, fun=fun )  
    integralArray_pad!( img, numF, 3 .* rad .+ 1, fun=(x)->(fun(x)>0) );  
    
    if th_hard == nothing
        th_hard = intA[ end, end ] / length(img); # == mean of the image
    end
    
    # convenient quantities
    lows  = (  1,  1  ) .+ ( 3 .* rad .+ 1 );
    highs = size( img ) .+ ( 3 .* rad .+ 1 ); 
    npix  = prod( ( 2 .* rad .+ 1 ) ); 
    vert  = ( 2*rad[1] + 1 - ovp[1], 0 ); 
    horz  = ( 0, 2*rad[2] + 1  - ovp[2]);

    sums9 = zeros( Float32, 3, 3 ); 
    
    
    @inbounds for x in lows[2]:highs[2], y in lows[1]:highs[1]
        
        tl = ( y,x ) .- rad .- 1; 
        br = ( y,x ) .+ rad; 
        s0 = ImageAnalysis.integralArea( intA, tl, br );  
        sF = ImageAnalysis.integralArea( intF, tl, br ); 
        nF = ImageAnalysis.integralArea( numF, tl, br );   
        mean = ( s0 - sF )/( 1 + npix - nF );     

        sums9[1] = ( ImageAnalysis.integralArea( intA, tl .- vert .- horz, br .- vert .- horz ) - ImageAnalysis.integralArea( intF, tl .- vert .- horz, br .- vert .- horz ) )/( 1 + npix - ImageAnalysis.integralArea( numF, tl .- vert .- horz, br .- vert .- horz ) );
        sums9[2] = ( ImageAnalysis.integralArea( intA, tl         .- horz, br         .- horz ) - ImageAnalysis.integralArea( intF, tl         .- horz, br         .- horz ) )/( 1 + npix - ImageAnalysis.integralArea( numF, tl         .- horz, br         .- horz ) );
        sums9[3] = ( ImageAnalysis.integralArea( intA, tl .+ vert .- horz, br .+ vert .- horz ) - ImageAnalysis.integralArea( intF, tl .+ vert .- horz, br .+ vert .- horz ) )/( 1 + npix - ImageAnalysis.integralArea( numF, tl .+ vert .- horz, br .+ vert .- horz ) ); 
        sums9[4] = ( ImageAnalysis.integralArea( intA, tl .+ vert        , br .+ vert         ) - ImageAnalysis.integralArea( intF, tl .+ vert        , br .+ vert         ) )/( 1 + npix - ImageAnalysis.integralArea( numF, tl .+ vert        , br .+ vert         ) );
        sums9[5] = ( ImageAnalysis.integralArea( intA, tl .+ vert .+ horz, br .+ vert .+ horz ) - ImageAnalysis.integralArea( intF, tl .+ vert .+ horz, br .+ vert .+ horz ) )/( 1 + npix - ImageAnalysis.integralArea( numF, tl .+ vert .+ horz, br .+ vert .+ horz ) );
        sums9[6] = ( ImageAnalysis.integralArea( intA, tl         .+ horz, br         .+ horz ) - ImageAnalysis.integralArea( intF, tl         .+ horz, br         .+ horz ) )/( 1 + npix - ImageAnalysis.integralArea( numF, tl         .+ horz, br         .+ horz ) );
        sums9[7] = ( ImageAnalysis.integralArea( intA, tl .- vert .+ horz, br .- vert .+ horz ) - ImageAnalysis.integralArea( intF, tl .- vert .+ horz, br .- vert .+ horz ) )/( 1 + npix - ImageAnalysis.integralArea( numF, tl .- vert .+ horz, br .- vert .+ horz ) );
        sums9[8] = ( ImageAnalysis.integralArea( intA, tl .- vert        , br .- vert         ) - ImageAnalysis.integralArea( intF, tl .- vert        , br .- vert         ) )/( 1 + npix - ImageAnalysis.integralArea( numF, tl .- vert        , br .- vert         ) );
        sums9[9] = sums9[1]; #( ImageAnalysis.integralArea( intA, tl                , br                 ) - ImageAnalysis.integralArea( intF, tl                , br                 ) )/( 1 + npix - ImageAnalysis.integralArea( numF, tl                , br                 ) );
        
        difs = 0.0
        for idx in 1:8
            difs = max( abs( sums9[idx] - sums9[idx+1] ), difs ); 
        end
        laplacian[y-3*rad[1]-1,x-3*rad[2]-1] = difs;
    end

    return laplacian
end

function mean_corr_pad( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fmin=1, fmax=1, ovp=(0,0) ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    return mean_corr_pad!( img, zeros(Float64, padsize .+ 1), zeros(Float64, padsize .+ 1),
                           rad, th_hard, zeros(Float32,size(img)), fmin=fmin, fmax=fmax, ovp=ovp );
end

function mean_corr_pad!( img::Array{T,2}, 
                         intA::Array{<:AbstractFloat,2}, intA2::Array{<:AbstractFloat,2},
                         rad::Tuple{Int64,Int64}, th_hard, 
                         corrs::Array{Float32,2};
                         fmin=1, fmax=1, ovp=(0,0) ) where {T<:Real}

    # integral array
    @inbounds @simd for idx in 1:length(intA)
         intA[idx] = 0.0;
        intA2[idx] = 0.0;
    end
    integralArray!( img, intA );
    integralArray!( img, intA2, fun=(x)->(x^2) );  
    
    if th_hard == nothing
        th_hard = intA[ end, end ] / length(img); # == mean of the image
    end
    
    # convenient quantities
    lows  = (  1,  1  ) .+ rad;
    highs = size( img ) .- rad; 
    npix  = prod( ( 2 .* rad .+ 1 ) ); 
    dists = zeros( Float32, 2 .* rad .+ 1 ); 

    # x = image intensity
    # y = distance from minimum value
    # corr( x, y ) = sum( ( x - meanX )( y - meanY ) )/ ( stdx * stdy ); 
    #              = sum( x*y - x*meanY - y*meanX + meanX*meanY )
    #              = sum( x*y ) - meanY*sum( x ) - meanX*sum( y ) + N*meanX*meanY
    #              =   loop       loop * IA          IA * loop         IA * loop

    # stdx = sum( ( x - meanX )^2 ) = sum( x^2 ) + N*meanX^2 - meanX*sum( x )  -> IA
    # stdy = sum( ( y - meanY )^2 ) = sum( y^2 ) + N*meanY^2 - meanY*sum( y )  -> loop
    
    @inbounds for x in lows[2]:highs[2], y in lows[1]:highs[1]
        
        tl    = ( y,x ) .- rad .- 1; 
        br    = ( y,x ) .+ rad; 
        sumI  = ImageAnalysis.integralArea( intA , tl, br );
        sumI2 = ImageAnalysis.integralArea( intA2, tl, br ); 
        meanI = sumI/npix;  
        stdI  = sqrt( ( sumI2 + npix*meanI^2 - meanI*sumI )/npix );    

        ( meanI <= th_hard ) && ( continue; )

        #println( tlimg, "  ", brimg )

        sumD, sumD2, meanD, sumID = maxdist!( img, dists, tl[1] + 1, br[1], tl[2] + 1, br[2] ); 
        stdD = sqrt( ( sumD2 + npix*meanD^2 - meanD*sumD )/npix ); 

        corr = ( sumID - meanD*sumI - meanI*sumD + npix*meanI*meanD )/( stdI * stdD )
                                        
        corrs[y,x] = corr;    
    end

    return corrs
end

function mindist!( a, dists, y0, y1, x0, x1 ) 
	minpos = ( y0, x0 ); 
	minval = a[minpos...]; 
	N = 1; 
	for col in x0:x1, row in y0:y1
			
        if     ( a[row,col] <  minval )
            minval = a[row,col]; 
            minpos = (row,col); 
            N = 1; 
        elseif ( a[row,col] == minval )
            minpos = minpos .+ (row,col); 
            N += 1; 
        end
	end

    minpos = minpos ./ N

    sumdist2 = 0.0
    sumdist  = 0.0
    sumXY    = 0.0
	@inbounds for col in 1:size( dists, 2 ) 
		d2 = ( x0 + col - 1 - minpos[2] )^2;
		@simd for row in 1:size( dists, 1 )
			d1 = ( y0 + row - 1 - minpos[1] )^2; 
			dists[row,col] = sqrt( d1 + d2 )
            sumdist  += dists[row,col]
            sumdist2 += d1 + d2
            sumXY    += a[ y0+row-1, x0+col-1 ] * dists[row,col]; 
		end
	end
    
	return sumdist, sumdist2, sumdist/length(dists ), sumXY # mean
end

function maxdist!( a, dists, y0, y1, x0, x1 ) 
	maxpos = ( y0, x0 ); 
	maxval = a[maxpos...]; 
	N = 1; 
	for col in x0:x1, row in y0:y1
			
        if     ( a[row,col] > maxval )
            maxval = a[row,col]; 
            maxpos = (row,col); 
            N = 1; 
        elseif ( a[row,col] == maxval )
            maxpos = maxpos .+ (row,col); 
            N += 1; 
        end
	end

    maxpos = maxpos ./ N;

    sumdist2 = 0.0
    sumdist  = 0.0
    sumXY    = 0.0
	@inbounds for col in 1:size( dists, 2 ) 
		d2 = ( x0 + col - 1 - maxpos[2] )^2;
		@simd for row in 1:size( dists, 1 )
			d1 = ( y0 + row - 1 - maxpos[1] )^2; 
			dists[row,col] = sqrt( d1 + d2 )
            sumdist  += dists[row,col]
            sumdist2 += d1 + d2
            sumXY    += a[ y0+row-1, x0+col-1 ] * dists[row,col]; 
		end
	end
    
	return sumdist, sumdist2, sumdist/length(dists ), sumXY # mean
end





























"""  3D """

function mean_extrema_pad3( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fmin=1, fmax=1, ovp=(0,0) ) where {T<:Real}

    return mean_extrema_pad3!( img, zeros(Float64, size(img) .+ 1), rad, th_hard, zeros(Bool,size(img)), zeros(Bool,size(img)),
                              fmin=fmin, fmax=fmax, ovp=ovp );
end


function mean_extrema_pad3!( img::Array{T,2}, intA::Array{<:AbstractFloat,2},
                             rad::Tuple{Int64,Int64}, th_hard, 
                             minima::Array{Bool,2}, maxima::Array{Bool,2};
                             fmin=1, fmax=1, ovp=0 ) where {T<:Real}

    # integral array
    @inbounds @simd for idx in 1:length(intA)
    intA[idx] = 0.0
    end
    integralArray!( img, intA )   

    if th_hard == nothing
    th_hard = intA[ end, end ] / length(img); # == mean of the image
    end

    # convenient quantities
    npix  = prod( ( 2 .* rad .+ 1 ) ); 
    vert  = ( 2*rad[1] + 1 - ovp[1], 0 ); 
    horz  = ( 0, 2*rad[2] + 1  - ovp[2]);

    # these will accomodate the means and signs of the 9x9 mean neighbourhood
    sums9 = zeros( Float64, 3, 3 ); 
    signs = zeros(   Int64, 3, 3 ); 

    @inbounds for x in 1:size(img,2), y in 1:size(img,1)

    tl = max.( 1, min.( size(img), ( y,x ) .- rad .- 1 ) ); 
    br = max.( 1, min.( size(img), ( y,x ) .+ rad ) ); 

    s0 = ImageAnalysis.integralArea( intA, tl, br );        
    ( s0/npix <= th_hard ) && ( continue; ) # s0/npix == mean of the central square

    cont = 1; 
    for off in ( -1 .* (vert .+ horz), -1 .* horz, vert .- horz, -1 .* vert, (0,0), vert, horz .- vert, horz, horz .+ vert )
        tl_ = max.( 1, min.( size(img), tl .+ off ) ); 
        br_ = max.( 1, min.( size(img), br .+ off ) ); 
        sums9[cont] = ImageAnalysis.integralArea( intA, tl_, br_ ); 
        cont += 1; 
    end

    @inbounds @simd for idx in 1:9
        signs[idx] = sign( s0 .- sums9[idx] ); 
    end

    score_min = 0
    score_max = 0
    for idx in 1:4
        score_min += ( signs[idx] == -1 && signs[end-idx+1] == -1 )   
        score_max += ( signs[idx] ==  1 && signs[end-idx+1] ==  1 )
    end

    minima[y,x] = score_min > fmin; 
    maxima[y,x] = score_max > fmax; 
    end

    return minima, maxima
end



function int_maxima_pad!( vol::Array{T,3}, 
                          intA::Array{Float32,3}, 
                          maxima::Array{Bool,3}, minima::Array{Bool,3},
                          rad::Tuple{Int64,Int64,Int64}, th_hard; f=12, m=1 ) where {T<:Real}
   
    h, w, d = size( vol ); 
    npix = prod( ( 2 .* rad .+ 1 ) );
    
    # integral array of padvolume
    
    intA = ImageAnalysis.integralArray( vol, padvol, type=Float32 );
    
    # Convenience vectors to move 1 side in each direction
    
    vert = (2*rad[1] + 1,0,0); 
    horz = (0,2*rad[2] + 1,0); 
    depz = (0,0,2*rad[3] + 1);
    offs = ( -1 .* (vert .+ horz .+ depz ), -1 .* ( horz .+ depz ), vert .- horz .- depz, -1 .* ( vert .+ depz ), -1 .* depz, vert .- depz, horz .- vert .- depz, horz .- depz, vert .+ horz .- depz, 
             -1 .* (vert .+ horz ), -1 .* ( horz ), vert .- horz, -1 .* ( vert ), (0,0,0), vert, horz .- vert, horz, vert .+ horz,   
             -1 .* (vert .+ horz .- depz ), -1 .* ( horz .- depz ), vert .- horz .+ depz, -1 .* ( vert .- depz ), depz, vert .+ depz, horz .- vert .+ depz, horz .+ depz, vert .+ horz .+ depz, 
           )
    
    # array to store the 27-neighbour sums
    
    sums27 = zeros( Float64, 3, 3, 3 ); 
    signs  = zeros(   Int64, 3, 3, 3 ); 
   
    for z in 1:d, x in 1:w, y in 1:h

        # top-left and bottom-right corners
        tl = min.( 1, max.( size(vol), ( y,x,z ) .- rad .- 1 ) ); 
        br = min.( 1, max.( size(vol), ( y,x,z ) .+ rad ) ); 
        
        s0   = ImageAnalysis.integralArea( intA, tl, br );
        mean = s0/npix;
        
        if mean <= th_hard
            continue; 
        end

        cont = 1; 
        for off in offs
            tl_ = min.( 1, max.( size(vol), tl .+ off ) ); 
            br_ = min.( 1, max.( size(vol), br .+ off ) ); 

            sums27[cont] = ImageAnalysis.integralArea( intA, tl_, br_ ); 
            cont += 1; 
        end
        
        signs = sign.( s0*m .- sums27 ); 
        score = 0
        for idx in 1:13
            score += signs[idx] == 1 && signs[end-idx+1] == 1
        end
        
        maxima[y-3*rad[1]-1,x-3*rad[2]-1,z-3*rad[3]-1] = score > f; 
    end

    return maxima;
end

#=   
    sums27[1,1,1] = ImageAnalysis.integralArea( intA, tl .- vert .- horz .- depz, br .- vert .- horz .- depz );
    sums27[2,1,1] = ImageAnalysis.integralArea( intA, tl         .- horz .- depz, br         .- horz .- depz );
    sums27[3,1,1] = ImageAnalysis.integralArea( intA, tl .+ vert .- horz .- depz, br .+ vert .- horz .- depz ); 
    sums27[1,2,1] = ImageAnalysis.integralArea( intA, tl .- vert         .- depz, br .- vert         .- depz );
    sums27[2,2,1] = ImageAnalysis.integralArea( intA, tl                 .- depz, br                 .- depz );
    sums27[3,2,1] = ImageAnalysis.integralArea( intA, tl .+ vert         .- depz, br .+ vert         .- depz ); 
    sums27[1,3,1] = ImageAnalysis.integralArea( intA, tl .- vert .+ horz .- depz, br .- vert .+ horz .- depz );
    sums27[2,3,1] = ImageAnalysis.integralArea( intA, tl         .+ horz .- depz, br         .+ horz .- depz );
    sums27[3,3,1] = ImageAnalysis.integralArea( intA, tl .+ vert .+ horz .- depz, br .+ vert .+ horz .- depz );
    
    sums27[1,1,2] = ImageAnalysis.integralArea( intA, tl .- vert .- horz        , br .- vert .- horz         );
    sums27[2,1,2] = ImageAnalysis.integralArea( intA, tl         .- horz        , br         .- horz         );
    sums27[3,1,2] = ImageAnalysis.integralArea( intA, tl .+ vert .- horz        , br .+ vert .- horz         ); 
    sums27[1,2,2] = ImageAnalysis.integralArea( intA, tl .- vert                , br .- vert                 );
    sums27[2,2,2] = s0 # 
    sums27[3,2,2] = ImageAnalysis.integralArea( intA, tl .+ vert                , br .+ vert                 ); 
    sums27[1,3,2] = ImageAnalysis.integralArea( intA, tl .- vert .+ horz        , br .- vert .+ horz         );
    sums27[2,3,2] = ImageAnalysis.integralArea( intA, tl         .+ horz        , br         .+ horz         );
    sums27[3,3,2] = ImageAnalysis.integralArea( intA, tl .+ vert .+ horz        , br .+ vert .+ horz         );
    
    sums27[1,1,3] = ImageAnalysis.integralArea( intA, tl .- vert .- horz .+ depz, br .- vert .- horz .+ depz );
    sums27[2,1,3] = ImageAnalysis.integralArea( intA, tl         .- horz .+ depz, br         .- horz .+ depz );
    sums27[3,1,3] = ImageAnalysis.integralArea( intA, tl .+ vert .- horz .+ depz, br .+ vert .- horz .+ depz ); 
    sums27[1,2,3] = ImageAnalysis.integralArea( intA, tl .- vert         .+ depz, br .- vert         .+ depz );
    sums27[2,2,3] = ImageAnalysis.integralArea( intA, tl                 .+ depz, br                 .+ depz );
    sums27[3,2,3] = ImageAnalysis.integralArea( intA, tl .+ vert         .+ depz, br .+ vert         .+ depz ); 
    sums27[1,3,3] = ImageAnalysis.integralArea( intA, tl .- vert .+ horz .+ depz, br .- vert .+ horz .+ depz );
    sums27[2,3,3] = ImageAnalysis.integralArea( intA, tl         .+ horz .+ depz, br         .+ horz .+ depz );
    sums27[3,3,3] = ImageAnalysis.integralArea( intA, tl .+ vert .+ horz .+ depz, br .+ vert .+ horz .+ depz );
=#


# gr = grid radius; sr = square radius; big_radius = (2 .* sr .+  1) .* gr .+ sr 
# padvol = sizevol + 2.* big_radius
# you can avoid padding by testing wether coordinates are within range, but that is somewhat slower 

function local_grid_pad( vol, gr, sr )
    big_radius = ( 2 .* sr .+  1 ) .* gr .+ sr;
    padsize    = size( vol ) .+ 2 .* big_radius;
    local_grid_pad!( vol, zeros( padsize .+ 1 ), zeros( size(vol) ), gr, sr )
end

function local_grid_pad!( vol::Array{T,3}, 
                          intA::Array{Float32,3}, 
                          out::Array{Bool,3},
                          gr::Tuple{Int64,Int64,Int64},
                          sr::Tuple{Int64,Int64,Int64} ) where {T<:Real}

    # integral array of padvolume
    big_rad = ( 2 .* sr .+  1 ) .* gr .+ sr;
    intA = ImageAnalysis.integralArray_pad!( vol, intA, ones(3).*big_rad, type=Float32 );

    h, w, d = size( vol ); 
    n = prod( ( 2 .* square_rad .+ 1 ) ); # number of pi/voxels in each grid element

    # Convenience vectors to move 1 side in each direction

    vert = (2*rad[1] + 1,0,0); 
    horz = (0,2*rad[2] + 1,0); 
    depz = (0,0,2*rad[3] + 1);

    offs = ( -1 .* (vert .+ horz .+ depz ), -1 .* ( horz .+ depz ), vert .- horz .- depz, -1 .* ( vert .+ depz ), -1 .* depz, vert .- depz, horz .- vert .- depz, horz .- depz, vert .+ horz .- depz, 
    -1 .* (vert .+ horz ), -1 .* ( horz ), vert .- horz, -1 .* ( vert ), (0,0,0), vert, horz .- vert, horz, vert .+ horz,   
    -1 .* (vert .+ horz .- depz ), -1 .* ( horz .- depz ), vert .- horz .+ depz, -1 .* ( vert .- depz ), depz, vert .+ depz, horz .- vert .+ depz, horz .+ depz, vert .+ horz .+ depz, 
    )

    # array to store the 27-neighbour sums

    sums27 = zeros( Float64, 3, 3, 3 ); 
    signs  = zeros(   Int64, 3, 3, 3 ); 

    for z in 1:d, x in 1:w, y in 1:h

        # top-left and bottom-right corners
        tl = min.( 1, max.( size(vol), ( y,x,z ) .- rad .- 1 ) ); 
        br = min.( 1, max.( size(vol), ( y,x,z ) .+ rad ) ); 

        s0   = ImageAnalysis.integralArea( intA, tl, br );
        mean = s0/npix;

        if mean <= th_hard
        continue; 
        end

        cont = 1; 
        for off in offs
            tl_ = min.( 1, max.( size(vol), tl .+ off ) ); 
            br_ = min.( 1, max.( size(vol), br .+ off ) ); 

            sums27[cont] = ImageAnalysis.integralArea( intA, tl_, br_ ); 
            cont += 1; 
        end

        signs = sign.( s0*m .- sums27 ); 
        score = 0
        for idx in 1:13
            score += signs[idx] == 1 && signs[end-idx+1] == 1
        end

        maxima[y-3*rad[1]-1,x-3*rad[2]-1,z-3*rad[3]-1] = score > f; 
    end

    return maxima;
end
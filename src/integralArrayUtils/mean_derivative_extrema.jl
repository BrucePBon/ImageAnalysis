##### BOTH MINIMA AND MAXIMA

function mean_extrema_pad( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fmin=1, fmax=1, ovp=(0,0), f=1 ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    return mean_extrema_pad!( img, zeros(Float64, padsize .+ 1), rad, th_hard, 
                              zeros(Bool,size(img)), zeros(Bool,size(img)),
                              fmin=fmin, fmax=fmax, ovp=ovp, f=f );
end

function mean_extrema_pad!( img::Array{T,2}, intA::Array{<:AbstractFloat,2},
                            rad::Tuple{Int64,Int64}, th_hard, 
                            minima::Array{Bool,2}, maxima::Array{Bool,2};
                            fmin=1, fmax=1, ovp=0, f=1 ) where {T<:Real}
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
            signs[idx] = sign( s0*f - sums9[idx] ); 
            idx += 1;   
        end

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

function mean_extrema_pad( vol::Array{T,3}, rad::Tuple{Int64,Int64,Int64}, th_hard; fmin=1, fmax=1, ovp=(0,0,0) ) where {T<:Real}

    padsize = size( vol ) .+ 2 .* ( 3 .* rad .+ 1 );
    return mean_extrema_pad!( vol, zeros(Float64, padsize .+ 1), rad, th_hard, 
                              zeros(Bool,size(vol)), zeros(Bool,size(vol)),
                              fmin=fmin, fmax=fmax, ovp=ovp );
end

function mean_extrema_pad!( img::Array{T,3}, intA::Array{<:AbstractFloat,3},
                            rad::Tuple{Int64,Int64,Int64}, th_hard, 
                            minima::Array{Bool,3}, maxima::Array{Bool,3};
                            fmin=1, fmax=1, ovp=(0,0,0) ) where {T<:Real}
    # integral array
    @inbounds @simd for idx in 1:length(intA)
        intA[idx] = 0.0
    end
    integralArray_pad!( img, intA, 3 .* rad .+ 1 )   
    
    # convenient quantities
    lows  = (  1,  1, 1  ) .+ ( 3 .* rad .+ 1 );
    highs = size( img ) .+ ( 3 .* rad .+ 1 ); 
    npix  = prod(  (  2 .* rad .+ 1  )  ); 
    vert  = ( 2*rad[1] + 1 - ovp[1], 0,  0 ); 
    horz  = ( 0, 2*rad[2] + 1  - ovp[2], 0 );
    deep  = ( 0, 0, 2*rad[2] + 1  - ovp[2] );

    
    # these will accomodate the means and signs of the 3x3x3 mean neighbourhood
    sums27 = zeros( Float64, 3, 3, 3 ); 
    signs  = zeros(   Int64, 3, 3, 3 ); 
    
    @inbounds for z in lows[3]:highs[3], x in lows[2]:highs[2], y in lows[1]:highs[1]
        
        tl = ( y,x,z ) .- rad .- 1; 
        br = ( y,x,z ) .+ rad; 
        s0 = ImageAnalysis.integralArea( intA, tl, br );   
             
        ( s0/npix <= th_hard ) && ( continue; ) # s0/npix == mean of the central square

        idx = 1;
        for zoff in -1:1, xoff in -1:1, yoff in -1:1
            off = vert .* yoff .+ horz .* xoff .+ deep .* zoff; 
            sums27[idx] = ImageAnalysis.integralArea( intA, tl .+ off, br .+ off ); 
            signs[ idx] = sign( s0 - sums27[idx] ); 
            idx += 1;   
        end

        score_min = 0
        score_max = 0
        for idx in 1:13
            score_min += ( signs[idx] == -1 && signs[end-idx+1] == -1 )   
            score_max += ( signs[idx] ==  1 && signs[end-idx+1] ==  1 )
        end
        
        minima[y-3*rad[1]-1,x-3*rad[2]-1,z-3*rad[3]-1] = score_min > fmin; 
        maxima[y-3*rad[1]-1,x-3*rad[2]-1,z-3*rad[3]-1] = score_max > fmax; 
    end

    return minima, maxima
end


#### ONLY MAXIMA 

function mean_maxima_pad( vol::Array{T,3}, rad::Tuple{Int64,Int64,Int64}, th_hard; fmax=1, ovp=(0,0,0), f=1 ) where {T<:Real}

    padsize = size( vol ) .+ 2 .* ( 3 .* rad .+ 1 );
    intA    = zeros( Float64, padsize .+ 1 )
    maxima  = zeros( Bool, size(vol) ); 

    integralArray_pad!( vol, intA, 3 .* rad .+ 1 );
    mean_maxima_pad!( vol, intA, rad, th_hard, maxima, fmax=fmax, ovp=ovp, f=f );

    return maxima
end


# for nsqecc struc surface
function mean_maxima_pad!( img::Array{T,3}, intA::Array{<:AbstractFloat,3},
                           rad::Tuple{Int64,Int64,Int64}, th_hard, 
                           maxima::Array{Bool,3}; fmax=1, ovp=(0,0,0), f=1 ) where {T<:Real}

    # convenient quantities
    lows  = (  1,  1, 1  ) .+ ( 3 .* rad .+ 1 );
    highs = size( img ) .+ ( 3 .* rad .+ 1 ); 
    npix  = prod(  (  2 .* rad .+ 1  )  ); 
    vert  = ( 2*rad[1] + 1 - ovp[1], 0,  0 ); 
    horz  = ( 0, 2*rad[2] + 1  - ovp[2], 0 );
    deep  = ( 0, 0, 2*rad[3] + 1  - ovp[3] );

    # these will accomodate the means and signs of the 3x3x3 mean neighbourhood
    sums27 = zeros( Float64, 3, 3, 3 ); 
    signs  = zeros(   Int64, 3, 3, 3 ); 

    @inbounds for z in lows[3]:highs[3], x in lows[2]:highs[2], y in lows[1]:highs[1]

        tl = ( y,x,z ) .- rad .- 1; 
        br = ( y,x,z ) .+ rad; 

        s0 = ImageAnalysis.integralArea( intA, tl, br );   

        ( s0/npix <= th_hard ) && ( continue; ) # s0/npix == mean of the central square

        idx = 1;
        for zoff in -1:1, xoff in -1:1, yoff in -1:1
            off = vert .* yoff .+ horz .* xoff .+ deep .* zoff; 
            sums27[idx] = ImageAnalysis.integralArea( intA, tl .+ off, br .+ off ); 
            signs[ idx] = sign( s0*f - sums27[idx] ); 
            idx += 1;   
        end

        score_max = 0
        for idx in 1:13
            @inbounds score_max += ( signs[idx] ==  1 && signs[end-idx+1] ==  1 )
        end

        maxima[y-3*rad[1]-1,x-3*rad[2]-1,z-3*rad[3]-1] = score_max > fmax; 
    end

    return maxima
end

#=
    # TODO: optimize this. We only need to recompute maxima around the deleted maxima, not the whole
    #       volume

    mxsum = zeros( Bool, size( vol ) ); 
    for rep in 1:10
        print( rep, ", " ); 
        mx1   = ImageAnalysis.mean_maxima_pad( vol, mxrad, 0, fmax=fmax, ovp=mxovp );
        mxsum = Base.or_int.( mxsum, mx1 )
        false && ( VTJK.volume2VTK( UInt8.(mx1), path=pwd(), fn="5ME_$(rep).vtk" ); )
        volth = vol .* Base.not_int.(mx1); 
    end
=#
function mean_maxima_pad_masked!( img::Array{T,3}, mask::Array{Bool,3},
                                  intA::Array{<:AbstractFloat,3},
                                  rad::Tuple{Int64,Int64,Int64}, th_hard, 
                                  maxima::Array{Bool,3}; fmax=1, ovp=(0,0,0) ) where {T<:Real}

    # convenient quantities
    lows  = (  1,  1, 1  ) .+ ( 3 .* rad .+ 1 );
    highs = size( img ) .+ ( 3 .* rad .+ 1 ); 
    npix  = prod(  (  2 .* rad .+ 1  )  ); 
    vert  = ( 2*rad[1] + 1 - ovp[1], 0,  0 ); 
    horz  = ( 0, 2*rad[2] + 1  - ovp[2], 0 );
    deep  = ( 0, 0, 2*rad[3] + 1  - ovp[3] );

    # these will accomodate the means and signs of the 3x3x3 mean neighbourhood
    sums27 = zeros( Float64, 3, 3, 3 ); 
    signs  = zeros(   Int64, 3, 3, 3 ); 

    @inbounds for z in lows[3]:highs[3], x in lows[2]:highs[2], y in lows[1]:highs[1]

        if !mask[y-3*rad[1]-1,x-3*rad[2]-1,z-3*rad[3]-1]
            continue
        end

        tl = ( y,x,z ) .- rad .- 1; 
        br = ( y,x,z ) .+ rad; 

        s0 = ImageAnalysis.integralArea( intA, tl, br );   

        ( s0/npix <= th_hard ) && ( continue; ) # s0/npix == mean of the central square

        idx = 1;
        for zoff in -1:1, xoff in -1:1, yoff in -1:1
            off = vert .* yoff .+ horz .* xoff .+ deep .* zoff; 
            sums27[idx] = ImageAnalysis.integralArea( intA, tl .+ off, br .+ off ); 
            signs[ idx] = sign( s0 - sums27[idx] ); 
            idx += 1;   
        end

        score_max = 0
        for idx in 1:13
            @inbounds score_max += ( signs[idx] ==  1 && signs[end-idx+1] ==  1 )
        end

        maxima[y-3*rad[1]-1,x-3*rad[2]-1,z-3*rad[3]-1] = score_max > fmax; 
    end

    return maxima
end

function iterative_maxima_old( vol1, rad, th_hard; fmax=1, ovp=(0,0,0), iters=5, fmax0=10 )

    vol   = copy(vol1); 
    h,w,d = size(vol); 

    # initial maxima computation
    padsize = size( vol ) .+ 2 .* ( 3 .* rad .+ 1 );
    intA    = zeros(Float64, padsize .+ 1)
    maxima  = zeros( Bool, size(vol) ); 

    integralArray_pad!( vol, intA, 3 .* rad .+ 1 ); 
    mean_maxima_pad!( vol, intA, rad ,th_hard, maxima, fmax=fmax0, ovp=ovp ); 

    mask  = zeros( Bool, size(vol) ); 
    print( "iteration number: ")
    for i in 1:iters
        print( i, "," );
        for z in 2:d-1, x in 2:w-1, y in 2:h-1
            # eliminate maxima and label neighbours for processing
            if maxima[y,x,z] 
                for yoff in -1:1, xoff in -1:1, zoff in -1:1
                    isvol = vol[y+yoff,x+xoff,z+zoff] > 0
                    mask[y+yoff,x+xoff,z+zoff] = true && isvol; 
                end
                mask[y,x,z] = false;
                 vol[y,x,z] = 0;
            end
        end 
        intA = integralArray_pad!( vol, intA, 3 .* rad .+ 1 ); 
        mean_maxima_pad_masked!( vol, mask, intA, rad ,th_hard, maxima, fmax=fmax, ovp=ovp ); 
    end
    println()

    return maxima
end

function iterative_maxima( vol1, rad, th_hard; fmax=1, ovp=(0,0,0), iters=5, fmax0=10 )

    vol   = copy(vol1); 
    h,w,d = size(vol); 

    # initial maxima computation
    padsize = size( vol ) .+ 2 .* ( 3 .* rad .+ 1 );
    intA    = zeros(Float64, padsize .+ 1)
    maxima  = zeros( Bool, size(vol) ); 

    integralArray_pad!( vol, intA, 3 .* rad .+ 1 ); 
    mean_maxima_pad!( vol, intA, rad ,th_hard, maxima, fmax=fmax0, ovp=ovp ); 

    mask  = zeros( Bool, size(vol) ); 
    print( "iteration number: ")
    for i in 1:iters
        print( i, "," );
        for z in 2:d-1, x in 2:w-1, y in 2:h-1
            # eliminate maxima and label neighbours for processing
            if maxima[y,x,z] 
                minval = Inf
                for yoff in -1:1, xoff in -1:1, zoff in -1:1
                    minval = min( vol[y+yoff,x+xoff,z+zoff], minval )
                    mask[y+yoff,x+xoff,z+zoff] = true; 
                end
                 vol[y,x,z] = minval;
            end
        end 
        intA = integralArray_pad!( vol, intA, 3 .* rad .+ 1 ); 
        mean_maxima_pad_masked!( vol, mask, intA, rad ,th_hard, maxima, fmax=fmax, ovp=ovp ); 
        @inbounds @simd for idx in 1:length(mask) 
            mask[idx] = false;
        end
    end
    println()

    return maxima
end


#### ONLY MINIMA 

function mean_minima_pad( img::Array{T,3}, rad::Tuple{Int64,Int64,Int64}, th_hard; fmin=1, ovp=(0,0,0), f=1 ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    return mean_minima_pad!( img, zeros(Float64, padsize .+ 1), rad, th_hard, 
                             zeros(Bool,size(img)), fmin=fmin, ovp=ovp, f=f );
end


# for nsqecc struc surface
function mean_minima_pad!( img::Array{T,3}, intA::Array{<:AbstractFloat,3},
                           rad::Tuple{Int64,Int64,Int64}, th_hard, 
                           minima::Array{Bool,3}; fmin=1, ovp=(0,0,0), f=1 ) where {T<:Real}
    # integral array
    @inbounds @simd for idx in 1:length(intA)
    intA[idx] = 0.0
    end
    integralArray_pad!( img, intA, 3 .* rad .+ 1 )   

    # convenient quantities
    lows  = (  1,  1, 1  ) .+ ( 3 .* rad .+ 1 );
    highs = size( img ) .+ ( 3 .* rad .+ 1 ); 
    npix  = prod(  (  2 .* rad .+ 1  )  ); 
    vert  = ( 2*rad[1] + 1 - ovp[1], 0,  0 ); 
    horz  = ( 0, 2*rad[2] + 1  - ovp[2], 0 );
    deep  = ( 0, 0, 2*rad[3] + 1  - ovp[3] );

    # these will accomodate the means and signs of the 3x3x3 mean neighbourhood
    sums27 = zeros( Float64, 3, 3, 3 ); 
    signs  = zeros(   Int64, 3, 3, 3 ); 

    @inbounds for z in lows[3]:highs[3], x in lows[2]:highs[2], y in lows[1]:highs[1]

        tl = ( y,x,z ) .- rad .- 1; 
        br = ( y,x,z ) .+ rad; 
        s0 = ImageAnalysis.integralArea( intA, tl, br );   

        ( s0/npix <= th_hard ) && ( continue; ) # s0/npix == mean of the central square

        idx = 1;
        for zoff in -1:1, xoff in -1:1, yoff in -1:1
            off = vert .* yoff .+ horz .* xoff .+ deep .* zoff; 
            sums27[idx] = ImageAnalysis.integralArea( intA, tl .+ off, br .+ off ); 
            signs[ idx] = sign( s0*f - sums27[idx] ); 
            idx += 1;   
        end

        score_min = 0
        for idx in 1:13
            score_min += ( signs[idx] ==  -1 && signs[end-idx+1] == -1 )
        end

        minima[y-3*rad[1]-1,x-3*rad[2]-1,z-3*rad[3]-1] = score_min > fmin; 
    end

    return minima
end


#### segmentation of spheroids, adaptations




function mean_minima_pad_onmax( img::Array{T,3}, mask, rad1, rad2, th_hard; f=1 ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( rad1 .+ rad2 );

    return mean_minima_pad_onmax!( img, mask, rad1, rad2, th_hard, 
                                   zeros(Float64, padsize .+ 1), zeros(Float64, padsize .+ 1), zeros(Bool,size(img)), f=f );
end


# for nsqecc struc surface
function mean_minima_pad_onmax!( img::Array{T,3}, mask, rad1, rad2, th_hard,
                                 intN, intA, minima; f=1 ) where {T<:Real}
    # integral array
    @inbounds @simd for idx in 1:length(intA)
        intA[idx] = 0.0;
        intN[idx] = 0.0;
    end
    integralArray_pad!( img .* mask, intA, rad1 .+ rad2 )
    integralArray_pad!(     mask   , intN, rad1 .+ rad2 )

    # convenient quantities
    lows  = ( 1, 1, 1 ) .+ ( rad1 .+ rad2 );
    highs = size( img ) .+ ( rad1 .+ rad2 ); 

    # these will accomodate the means and signs of the 3x3x3 mean neighbourhood
    sums27 = zeros( Float64, 3, 3, 3 ); 
    signs  = zeros(   Int64, 3, 3, 3 ); 

    for z in lows[3]:highs[3], x in lows[2]:highs[2], y in lows[1]:highs[1]

        if !mask[y-lows[1]+1,x-lows[2]+1,z-lows[3]+1]
            continue
        end

        tl = ( y,x,z ) .- rad1 .- 1; 
        br = ( y,x,z ) .+ rad1; 

        sin  = ImageAnalysis.integralArea( intA, tl, br );  
        Nin  = ImageAnalysis.integralArea( intN, tl, br ); 

        sout = ImageAnalysis.integralArea( intA, tl .- rad2, br .+ rad2 ) - sin;
        Nout = ImageAnalysis.integralArea( intN, tl .- rad2, br .+ rad2 ) - Nin;

        minima[y-lows[1]+1,x-lows[2]+1,z-lows[3]+1] = sin/Nin*f < sout/Nout; 
    end

    return minima
end

function mean_minima_pad_fstd( img::Array{<:Real,3}, rad, th_hard; fmin=1, ovp=(0,0,0), f=1 )

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    intA    = zeros(Float64, padsize .+ 1)
    intA2   = zeros(Float64, padsize .+ 1)
    minima  = zeros(Bool,size(img))

    return mean_minima_pad_fstd!( img, intA, intA2, rad, th_hard, minima, fmin=fmin, ovp=ovp, f=f );
end


function mean_minima_pad_fstd!( img::Array{<:Real,3}, intA, intA2, rad, th_hard, minima; fmin=1, ovp=(0,0,0), f=1 )
    # integral array
    @inbounds @simd for idx in 1:length(intA)
         intA[idx] = 0.0
        intA2[idx] = 0.0
    end
    integralArray_pad!( img,  intA, 3 .* rad .+ 1 )   
    integralArray_pad!( img, intA2, 3 .* rad .+ 1, fun=(x)->(x*x) ) 
    
    # std_ = sum( I - M )^2  = sum( I2 + M2 - 2IM ) = sum( I2 ) + sum( M2 ) - 2Msum(I) = sum(I2) + N(M2) - 2/N sum(I)sum(I) = sum(I2) + sum(I)^2/N - 2sum(I)^2/N = sum(I2) - sum(I)^2/N


    # convenient quantities
    lows  = (  1,  1, 1  ) .+ ( 3 .* rad .+ 1 );
    highs = size( img ) .+ ( 3 .* rad .+ 1 ); 
    npix  = prod(  (  2 .* rad .+ 1  )  ); 
    vert  = ( 2*rad[1] + 1 - ovp[1], 0,  0 ); 
    horz  = ( 0, 2*rad[2] + 1  - ovp[2], 0 );
    deep  = ( 0, 0, 2*rad[3] + 1  - ovp[3] );
    stdpix = prod( 2 .* ( 3 .* rad .+ 1 ) .+ 1 ); 

    # these will accomodate the means and signs of the 3x3x3 mean neighbourhood
    sums27 = zeros( Float64, 3, 3, 3 ); 
    signs  = zeros(   Int64, 3, 3, 3 ); 

    @inbounds for z in lows[3]:highs[3], x in lows[2]:highs[2], y in lows[1]:highs[1]

        tl = ( y,x,z ) .- rad .- 1; 
        br = ( y,x,z ) .+ rad; 
        s0 = ImageAnalysis.integralArea( intA, tl, br );   

        ( s0/npix <= th_hard ) && ( continue; ) # s0/npix == mean of the central square

        std_I2 = ImageAnalysis.integralArea( intA2, tl .- horz .- vert .- deep, br .+ horz .+ vert .+ deep ); 
        std_II = ImageAnalysis.integralArea( intA, tl .- horz .- vert .- deep, br .+ horz .+ vert .+ deep )^2/stdpix;
        std    = sqrt( std_I2 - std_II )/stdpix; 


        idx = 1;
        for zoff in -1:1, xoff in -1:1, yoff in -1:1
            off = vert .* yoff .+ horz .* xoff .+ deep .* zoff; 
            sums27[idx] = ImageAnalysis.integralArea( intA, tl .+ off, br .+ off ); 
            signs[ idx] = sign( s0 + f*std - sums27[idx] ); 
            idx += 1;   
        end

        score_min = 0
        for idx in 1:13
            score_min += ( signs[idx] ==  -1 && signs[end-idx+1] == -1 )
        end

        minima[y-3*rad[1]-1,x-3*rad[2]-1,z-3*rad[3]-1] = score_min > fmin; 
    end

    return minima
end

function mean_minima_pad_uneven( img::Array{T,3}, rad1, rad2, th_hard; fmin=1, ovp=(0,0,0), f=1 ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 2 .* rad2 .+ 1 .+ rad1 );
    intA    = zeros(Float64, padsize .+ 1)
    minima  = zeros(Bool,size(img))

    return mean_minima_pad_uneven!( img, intA, rad1, rad2, th_hard, minima, fmin=fmin, ovp=ovp, f=f );
end


# for nsqecc struc surface
function mean_minima_pad_uneven!( img::Array{T,3}, intA, rad1, rad2, th_hard, minima; fmin=1, ovp=(0,0,0), f=1 ) where {T<:Real}
    
    pad   = ( 2 .* rad2 .+ 1 .+ rad1 )

    # integral array
    @inbounds @simd for idx in 1:length(intA)
    intA[idx] = 0.0
    end
    integralArray_pad!( img, intA, pad )   

    # convenient quantities
    lows  = (  1,  1, 1  ) .+ pad;
    highs = size( img ) .+ pad; 
    npix1 = prod(  (  2 .* rad1 .+ 1  )  ); 
    npix2 = prod(  (  2 .* rad2 .+ 1  )  ); 
    vert  = ( 2*rad2[1] + 1 + rad1[1] - ovp[1], 0,  0 ); 
    horz  = ( 0, 2*rad2[2] + 1 + rad1[2] - ovp[2], 0 );
    deep  = ( 0, 0, 2*rad2[3] + 1  + rad1[3] - ovp[3] );

    # these will accomodate the means and signs of the 3x3x3 mean neighbourhood
    sums27 = zeros( Float64, 3, 3, 3 ); 
    signs  = zeros(   Int64, 3, 3, 3 ); 

    @inbounds for z in lows[3]:highs[3], x in lows[2]:highs[2], y in lows[1]:highs[1]

        tl = ( y,x,z ) .- rad1 .- 1; 
        br = ( y,x,z ) .+ rad1; 
        s0 = ImageAnalysis.integralArea( intA, tl, br );   

        ( s0/npix1 <= th_hard ) && ( continue; ) # s0/npix == mean of the central square

        idx = 1;
        for zoff in -1:1, xoff in -1:1, yoff in -1:1
            off = vert .* yoff .+ horz .* xoff .+ deep .* zoff; 
            sums27[idx] = ImageAnalysis.integralArea( intA, tl .+ off, tl .+ off .+ ( 2 .* rad2 ) ); 
            signs[ idx] = sign( f*s0/npix1 - sums27[idx]/npix2 ); 
            idx += 1;   
        end

        score_min = 0
        for idx in 1:13
            score_min += ( signs[idx] ==  -1 && signs[end-idx+1] == -1 )
        end

        minima[y-pad[1],x-pad[2],z-pad[3]] = score_min > fmin; 
    end

    return minima
end













############################## VISUAL TESTING 

function VTmean_extrema_pad( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fmin=1, fmax=1, ovp=(0,0), f=1 ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    return VTmean_extrema_pad!( img, zeros(Float64, padsize .+ 1), rad, th_hard, 
                                zeros(Int64,size(img)), zeros(Int64,size(img)),
                                fmin=fmin, fmax=fmax, ovp=ovp, f=f );
end

function VTmean_extrema_pad!( img::Array{T,2}, intA::Array{<:AbstractFloat,2},
                              rad::Tuple{Int64,Int64}, th_hard, 
                              minima::Array{Int64,2}, maxima::Array{Int64,2};
                              fmin=1, fmax=1, ovp=0, f=1 ) where {T<:Real}
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
            signs[idx] = sign( s0*f - sums9[idx] ); 
            idx += 1;   
        end

        score_min = 0
        score_max = 0
        for idx in 1:4
            score_min += ( signs[idx] == -1 && signs[end-idx+1] == -1 )   
            score_max += ( signs[idx] ==  1 && signs[end-idx+1] ==  1 )
        end
        
        minima[y-3*rad[1]-1,x-3*rad[2]-1] = score_min; 
        maxima[y-3*rad[1]-1,x-3*rad[2]-1] = score_max; 
    end

    return minima, maxima
end
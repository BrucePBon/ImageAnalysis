##### BOTH MINIMA AND MAXIMA

function mean_extrema_pad( img::Array{T,2}, rad::Tuple{Int64,Int64}, th_hard; fmin=1, fmax=1, ovp=(0,0) ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    return mean_extrema_pad!( img, zeros(Float64, padsize .+ 1), rad, th_hard, 
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

function mean_maxima_pad( vol::Array{T,3}, rad::Tuple{Int64,Int64,Int64}, th_hard; fmax=1, ovp=(0,0,0) ) where {T<:Real}

    padsize = size( vol ) .+ 2 .* ( 3 .* rad .+ 1 );
    intA    = zeros( Float64, padsize .+ 1 )
    maxima  = zeros( Bool, size(vol) ); 

    integralArray_pad!( vol, intA, 3 .* rad .+ 1 );
    mean_maxima_pad!( vol, intA, rad, th_hard, maxima, fmax=fmax, ovp=ovp );

    return maxima
end


# for nsqecc struc surface
function mean_maxima_pad!( img::Array{T,3}, intA::Array{<:AbstractFloat,3},
                           rad::Tuple{Int64,Int64,Int64}, th_hard, 
                           maxima::Array{Bool,3}; fmax=1, ovp=(0,0,0) ) where {T<:Real}

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
    deep  = ( 0, 0, 2*rad[2] + 1  - ovp[2] );

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


#### ONLY MINIMA 

function mean_minima_pad( img::Array{T,3}, rad::Tuple{Int64,Int64,Int64}, th_hard; fmin=1, ovp=(0,0,0) ) where {T<:Real}

    padsize = size( img ) .+ 2 .* ( 3 .* rad .+ 1 );
    return mean_minima_pad!( img, zeros(Float64, padsize .+ 1), rad, th_hard, 
                             zeros(Bool,size(img)), fmin=fmin, ovp=ovp );
end


# for nsqecc struc surface
function mean_minima_pad!( img::Array{T,3}, intA::Array{<:AbstractFloat,3},
                           rad::Tuple{Int64,Int64,Int64}, th_hard, 
                           minima::Array{Bool,3}; fmin=1, ovp=(0,0,0) ) where {T<:Real}
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
        for idx in 1:13
            score_min += ( signs[idx] ==  -1 && signs[end-idx+1] == -1 )
        end

        minima[y-3*rad[1]-1,x-3*rad[2]-1,z-3*rad[3]-1] = score_min > fmin; 
    end

    return minima
end
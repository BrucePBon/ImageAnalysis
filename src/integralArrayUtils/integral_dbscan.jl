function integral_dbscan_pad( vol::Array{T,3}, rad::Tuple{Int64,Int64,Int64}; th_hard=0, ovp=(0,0,0) ) where {T<:Real}

    padsize = size( vol ) .+ 2 .* ( 3 .* rad .+ 1 );
    intA    = zeros(Float64, padsize .+ 1); 
    counts  = zeros( Int64, size(vol) ); 
    displ   = zeros( Int64, size(vol) ); # UInt32 is probably more than fine
    return integral_dbscan_pad!( vol, intA, rad, displ, counts, th_hard=th_hard, ovp=ovp );
end

function integral_dbscan_pad!( vol::Array{T,3}, intA::Array{<:AbstractFloat,3},
                               rad::Tuple{Int64,Int64,Int64}, displ::Array{Int64,3}, 
                               counts::Array{Int64,3}; th_hard=0, ovp=0 
                             ) where {T<:Real}  

    # convenient quantities
    lows  = (  1,  1, 1  ) .+ ( 3 .* rad .+ 1 );
    highs = size( vol ) .+ ( 3 .* rad .+ 1 ); 
    npix  = prod(  (  2 .* rad .+ 1  )  ); 
    vert  = ( 2*rad[1] + 1 - ovp[1], 0,  0 ); 
    horz  = ( 0, 2*rad[2] + 1  - ovp[2], 0 );
    deep  = ( 0, 0, 2*rad[2] + 1  - ovp[2] );
    offs  = 1, size(vol,1), size(vol,1)*size(vol,2); 
    
    # these will accomodate the means and signs of the 9x9 mean neighbourhood
    disp27 = [ offs[1]*y + offs[2]*x + offs[3]*z for y in -1:1, x in -1:1, z in -1:1 ]; 
    sums27 = zeros( Float64, 3, 3, 3 ); 
    signs  = zeros(   Int64, 3, 3, 3 ); 
    
    lidx = 1; 
    @inbounds for z in lows[3]:highs[3], x in lows[2]:highs[2], y in lows[1]:highs[1]
        
        if !vol[y-3*rad[1]-1,x-3*rad[2]-1,z-3*rad[3]-1]
            continue;
        end

        tl = ( y,x,z ) .- rad .- 1; 
        br = ( y,x,z ) .+ rad; 
        s0 = ImageAnalysis.integralArea( intA, tl, br );  
             
        ( s0/npix <= th_hard ) && ( continue; ) # s0/npix == mean of the central square

        maxoff = 0;
        maxdif = -Inf; 
        for zoff in -1:1, xoff in -1:1, yoff in -1:1
            off = vert .* yoff .+ horz .* xoff .+ deep .* zoff;
            sN  = ImageAnalysis.integralArea( intA, tl .+ off, br .+ off ); 

            if sN - s0 > maxdif
                maxoff = lidx + ( yoff * offs[1] + xoff * offs[2] + zoff * offs[3] )
                maxdif = sN - s0
            end 
        end
        lidx += 1
        displ[y-3*rad[1]-1,x-3*rad[2]-1,z-3*rad[3]-1] = maxoff
    end

    for idx in 1:length(displ)
        
        id0 = idx
        for rep in 1:10
            id0 = displ[id0]
        end
        counts[id0] += 1
    end

    return counts
end
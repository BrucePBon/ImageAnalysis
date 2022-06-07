
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


function mean_thresh_count( img::Array{T,3}, rad::Tuple{Int64,Int64,Int64}, th_hard; fun=(x,y)->(x<y), f=1.0 ) where {T<:Real}

    padsize = size( img ) .+ 2 .* rad; #( 3 .* rad .+ 1 );
    return mean_thresh_count!( img, 
                               zeros(Float64, padsize .+ 1),
                               rad, th_hard, zeros(Int64,size(img)), fun=fun, f=f );
end

function mean_thresh_count!(  img::Array{T,3}, 
                              intA::Array{<:AbstractFloat,3},
                              rad::Tuple{Int64,Int64,Int64}, th_hard, 
                              counts::Array{Int64,3}; fun=(x,y)->(x<y), f=1.0 ) where {T<:Real}
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
    npix       = prod( ( 2 .* rad .+ 1 ) ); 
    y0, x0, z0 = ( 1, 1, 1 ) .+ rad; #( 3 .* rad .+ 1 ); 
    y1, x1, z1 = size( img ) .+ rad; #( 3 .* rad .+ 1 ); 
    ry, rx, rz = rad; 
    
    for z in z0:z1, x in x0:x1, y in y0:y1
        
        # top left (tl) and bottom right (br) of the central patch
        tl = ( y, x, z ) .- rad .- 1; 
        br = ( y, x, z ) .+ rad;      
        s0 = ImageAnalysis.integralArea( intA, tl, br );    
        local_mean = s0*f/npix; 

        if local_mean < th_hard
            counts[ y-y0+1, x-x0+1, z-z0+1 ] = 1000000; 
        end
        
        # For each pixel in the central patch that falls inside the image (we added padding to the integral array)
        # add +1 if the pixel is above the local mean. 
        for zet in max(tl[3]+1,z0):min(br[3],z1)
            for col in max(tl[2]+1,x0):min(br[2],x1)
                for row in max(tl[1]+1,y0):min(br[1],y1)
                    counts[ row-y0+1,col-x0+1,zet-z0+1 ] += fun( img[ row-y0+1,col-x0+1,zet-z0+1 ], local_mean )
        end end end
        
        #filt[y-3*rad[1]-1,x-3*rad[2]-1] = img[ y-3*rad[1]-1,x-3*rad[2]-1 ] > local_mean*0.95 ; 
    end

    for idx in 1:length(counts)
        if counts[idx] > 1000000
            counts[idx] = -1
        end
    end

    return counts
end
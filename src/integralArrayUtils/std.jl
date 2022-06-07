
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

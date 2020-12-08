function linearBinSearch( array::Array{T1,1}, query::T2 ) where {T1<:Number,T2<:Number}

    idx = 0; 	
    @inbounds for e in 1:length(array)
        query <= array[ e ] && break
        idx = e; 
    end

    return idx
end

function binaryBinSearch( array::Array{T1,1}, query::T2 ) where {T1<:Number,T2<:Number}

    hi, lo, probe = length( array ) + 1, 0, 0; 
    
    @inbounds while 1 < hi - lo
        probe = ( hi + lo ) >>> 1  
        if array[probe] < query
            lo = probe 
        else 
            hi = probe
        end
    end
    
    return lo
end

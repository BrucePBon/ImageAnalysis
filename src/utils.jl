function linearBinSearch( array::Array{T1,1}, querry::T2 ) where {T1<:Number,T2<:Number}

    querry  > array[end] && ( return 0 )
    querry  < array[1]   && ( return 0 )
	querry == array[1]   && ( return 1 )

    idx = 1; 	

    @inbounds for e in 1:length(array)
        idx = e; 
        querry <= array[ e ] && break
    end

    return idx-1
end

function binaryBinSearch( array::Array{T1,1}, querry::T2 ) where {T1<:Number,T2<:Number}

    querry  > array[end] && ( return 0 )
    querry  < array[1]   && ( return 0 )
	querry == array[1]   && ( return 1 )

    hi, lo, probe = length( array ) + 1, 0, 0; 
    
    @inbounds while 1 < hi - lo
        probe = ( hi + lo ) >>> 1  
        if array[probe] < querry
            lo = probe 
        else 
            hi = probe
        end
    end
    
    return lo
end

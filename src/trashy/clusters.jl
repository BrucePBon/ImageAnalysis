# My first attempt at comming up with efficient connected components 

function groupTRUES( mask::Array{Bool,2} ) 

	h, w = size( mask ); 

	finalGroups   = [];
	prevColGroups = []; 

	for col in 1:w

		# Finding connected TRUES ( groups ) in the current column, 1D

	    colGroups = groupsInColumn( view( mask, :, col ), col, length(finalGroups) ); 
	    numGroups = length( colGroups ); 
	       
	    # For each group in colGroups, we have to check if it overlaps with any group in the previous column
	    # Since the groups are ordered from "top" to "bottom", we use the variable "start" to skip all groups
	    # in the previous column that are above the currently considered group in colGroups

	    start = 1; 

	    for idx1 in 1:numGroups

	        thisGroup = colGroups[idx1];  
	        overlaps  = false; 

	        for idx2 in start:length(prevColGroups)

	            prevGroup = prevColGroups[idx2]

	            if     thisGroup[1] > prevGroup[2] # thisGroup is above prevGroup

	                start = idx2+1; 
	            
	            elseif thisGroup[2] < prevGroup[1] # thisGroup is below prevGroup
	                break; 
	            
	            else                               # thisGroup and prevGroup overlap

	                # adding coordinates of thisGroup to prevGroup's parent
	                push!( finalGroups[ prevGroup[4] ], thisGroup[1:3] ); 

	                # changing thisGroup's parent to prevGroup's parent
	                colGroups[idx1][4] = prevGroup[4]; 
	                diditoverlap = true;      
	            end
	        end

	        # If thisGroup does not overlap, it will be added to the list of groups
	        # and have the change to nucleate its own group

	        if !overlaps
	            push!( finalGroups, [ thisGroup, ] )
	        end
	    end
	    prevColGroups = colGroups; 
	end

	return finalGroups
end

function groupsInColumn( column::AbstractArray{Bool,1}, colIdx, last )
	
	groups    = []; 
	prevPixel = false; 
    
    for row in 1:length(column) 
        pixel = column[ row ]
        if pixel 
            if prevPixel
                groups[end][2] = row; 
            else
                push!( groups, [ row, row, colIdx, last+1 ] );                 
            end
        end
        prevPixel = pixel; 
    end

    return groups; 
end

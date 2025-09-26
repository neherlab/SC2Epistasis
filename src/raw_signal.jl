""" 
    Functions to quantify the raw epistatic signal
    contained in the fitness discrepancies data
"""

# Compute the normalized fitness discrepancies over a dataframe
function z_dfit(dfit::DataFrame)

    z = zeros(Float64, size(dfit, 1))

    z .= (dfit.fit1 .- dfit.fit2) ./ sqrt.(dfit.std_fit1 .^ 2 .+ dfit.std_fit2 .^ 2)

    return z

end
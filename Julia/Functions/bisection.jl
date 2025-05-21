# Bisection Function

# Recursive implementation of the bisection algorithm

function bisection(low, high, W, u, eta)
    mid = (low + high) / 2
    flow = rCESS(W, u, low, eta)
    fmid = rCESS(W, u, mid, eta)
    fhigh = rCESS(W, u, high, eta)

    if flow * fhigh > 0
        error("Invalid endpoint for bisection.")
    end

    try 
        if low >= high
        error("bisection overlap")
        end

        if abs(fmid) < 1e-10 || (high - low)/2 < 1e-10
            return mid
        end
    
        if flow * fmid < 0
            return bisection(low, mid, W, u, eta)
        end
    
        if fhigh * fmid < 0
            return bisection(mid, high, W, u, eta)
        end

    catch e
        println("Caught error:", e)
    end

    error("bisection flawed")

end
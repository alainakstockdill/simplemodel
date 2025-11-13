# consumer resource model - want to run different equations in parallel
using Plots, DifferentialEquations, Parameters, BifurcationKit, Statistics, DataFrames, CSV

function rosemac!(du, u, params, t)

    # these should all be uncoupled but when i change d forcing it alters results from other equations...
    @unpack r, K, a, b, e, d, A, p = params
    # forcing = A * sin(2 * pi * p * t);

    # Khopf = -b * (d + a * e) / (d - a * e) + .05
    # dhopf = a * e * (K - b) / (K + b) + .05

    R, C = u 
    #R, C, RK, CK, Rd, Cd = u

    du[1] = R * (r * (1 - R / (K)) - C * a / (b + R)) 
    du[2] = C * ((e * a * R / (b + R)) - d)
    # du[3] = RK * (r * (1 - RK / (Khopf + forcing)) - CK * a / (b + RK)) 
    # du[4] = CK * ((e * a * RK / (b + RK)) - d)
    # du[5] = Rd * ((r) * (1 - Rd / (K)) - Cd * a / (b + Rd)) 
    # du[6] = Cd * ((e * (a) * Rd / (b + Rd)) - (dhopf + forcing))

end

function rosemacK!(du, u, params, t)

    # these should all be uncoupled but when i change d forcing it alters results from other equations...
    @unpack r, K, a, b, e, d, A, p = params
    forcing = A * sin(2 * pi * p * t);

    Khopf = -b * (d + a * e) / (d - a * e) + .05

    RK, CK = u
    du[1] = RK * (r * (1 - RK / (Khopf + forcing)) - CK * a / (b + RK)) 
    du[2] = CK * ((e * a * RK / (b + RK)) - d)

end

function rosemacd!(du, u, params, t)

    # these should all be uncoupled but when i change d forcing it alters results from other equations...
    @unpack r, K, a, b, e, d, A, p = params
    forcing = A * sin(2 * pi * p * t) / 5;


    dhopf = a * e * (K - b) / (K + b) + .05

    Rd, Cd = u
    du[1] = Rd * ((r) * (1 - Rd / (K)) - Cd * a / (b + Rd)) 
    du[2] = Cd * ((e * (a) * Rd / (b + Rd)) - (dhopf + forcing))

end
# redo with separate functions
growth = [1.5, 5, 10, 20];
period = logrange(0.01, 1000.0, length = 100)
u0 = [2.0, 2.0]
cvs = zeros(length(growth), length(period))

for (j, r) in enumerate(growth)
    
    for (i, p) in enumerate(period)

        #println(string("period:", p))
        params = (r = r, K = 2, a = 1.3, b = 1.0, e = .7, d = .2, A = .5, p = p); # these are the defaults
        tspan = (0.0, 50000.0)

        prob = ODEProblem(rosemac!, u0, tspan, params)
        probK = ODEProblem(rosemacK!, u0, tspan, params)
        #probd = ODEProblem(rosemacd!, u0, tspan, params)

        output = solve(prob, Rodas5(); maxiters=1e8, saveat = 1)
        outputK = solve(probK, Rodas5(); maxiters=1e8, saveat = 1)
        #outputd = solve(probd, Rodas5(); maxiters=1e8, saveat = 1)
        
        solution = hcat(output.u...)
        solutionK = hcat(outputK.u...)
        #solutiond = hcat(outputd.u...)
        
        sol = [solution; solutionK]

        #file = string("outputs/solution_r_",r,"_p_", ceil(Int, p), ".csv")
        #solutiondf = DataFrame(sol, :auto) 
        #CSV.write(file, solutiondf)

        cvrange = ceil(Int, (maximum((1000, 4 ./ p))))
        
        ssf = solutionK[1, (end-cvrange):end]
        ss = solution[1, (end - cvrange): end]

        cvs[j, i] = std(ssf) / mean(ssf) - std(ss) / mean(ss)


    end


end

cvsdf = DataFrame(cvs, :auto) 
CSV.write("cvs.csv", cvsdf)
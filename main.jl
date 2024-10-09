using Plots, DelimitedFiles,Roots,Measures
include("functions.jl")

data = readdlm("Shock_Nozzle_Data.tsv", '\t','\n',header = false)

backpressures = data[2:end,end]
chokedindex = 7
x = [-1.0,-.2, .3, 1.0, 1.7, 2.4, 3.1, 3.8, 4.2]
xmax = 4.2
gamma = 1.4
R = 287
T0 = 21.11111
P0 = 87149.732186

w = @. 0.4 + 0.04762*x
w[1:2] = [1e20, 0.4]
wmax = 0.4 + 0.04762*xmax
A = @. w*0.0254 #m^2
At = 0.4*0.0254
Ae = wmax*0.0254

AoverAstar = Ae/At
Me = machfromarea(AoverAstar,gamma,false)
Pe = pressurefrommach(Me,P0,gamma)

datasets::Vector{Vector} = []
prechokedatasets::Vector{Vector} = []
chokeddatasets::Vector{Vector} = []
for i = 2:1:length(data[:,1])
    vector::Vector{Float64} = []
    push!(vector, P0)
    for j = 1:1:length(data[1,:])
        push!(vector, data[i,j]*6894.7572932)
    end
    push!(datasets,vector)
    if i in [6, 7]
        push!(prechokedatasets,vector)
    elseif i in [8,9,10,11]
        push!(chokeddatasets,vector)
    end
end
Pprechokedtheory = []
for i = 1:1:2
    pressures = prechokedatasets[i]
    Pb = pressures[end]
    #println(Pb)
    local Me = machfrompressure(P0,Pb,gamma)
    #println(Me)
    Astar=astarfrommach(Me,Ae,gamma)
    #println(Astar)
    #println(At/Astar)
    #println(@. A/Astar)
    Ms = @. machfromarea(A/Astar,gamma,false)
    #println(Ms)
    Ptheory = @. pressurefrommach(Ms,P0,gamma)
    #Ptheory[end] = Pb
    #Ptheory[1] = P0
    push!(Pprechokedtheory,Ptheory)
end


Pchokedtheory = []
recordedshocklocations = [0.65, 1.35, 2.05, 2.75]
theoryshocklocation = zeros(4)
theoryshockM1 = zeros(4)
shockPbs = zeros(4)
for i = 1:1:4
    pressures = chokeddatasets[i]
    Pb = pressures[end]
    shockPbs[i] = deepcopy(Pb)
    error = 1
    margin = 1e-7
    Ashock = At+ 0.005#0.0254*(0.4 + 0.04762*recordedshocklocations[i]) #guess from experiment
    Ashockold = 0.0
    maxitters = 40000
    itters = 1
    while (error > margin) && (itters <= maxitters)
        #println(Ashock)
        M1 = machfromarea(Ashock/At,gamma,true)
        #I just plugged in gamma to eq 4.15:
        P0ratio = ( ((1.2*M1^2)/(1+.2*M1^2))^(1.4/.4) ) * ( 1/((2.8/2.4)*M1^2-.4/2.4) )^(1/.4)
        AetoAstar2 = (Ae/At)*P0ratio
        local Me = machfromarea(AetoAstar2,gamma,false)
        PetoP02 = (1+.2*Me^2)^(-1.4/.4)
        PetoP01 = PetoP02*P0ratio
        Ashockold = deepcopy(Ashock)
        step = 1e-4*0.999^itters
        if PetoP01 < (Pb/P0)
            Ashock = Ashock - step
        else
            Ashock = Ashock + step
        end
        error = abs(Ashock-Ashockold)
        itters = itters +1
    end
    theoryshockM1[i] = machfromarea(Ashock/At,gamma,true)
    theoryshocklocation[i]=((Ashock/0.0254) - 0.4)/0.04762
    theoryPs = []

    for j = 1:1:length(x)
        Mach = 1.0
        if x[j] < -0.2
            Mach = machfromarea(A[j]/At,gamma,false)
            push!(theoryPs,pressurefrommach(Mach,P0,gamma))
        elseif (x[j] > -0.2) && (x[j] < theoryshocklocation[i])
            Mach = machfromarea(A[j]/At,gamma,true)
            push!(theoryPs,pressurefrommach(Mach,P0,gamma))
        elseif (x[j] >= theoryshocklocation[i])#after shock
            M1 = theoryshockM1[i]
            #println(M1)
            P02 = P0*( ((1.2*M1^2)/(1+.2*M1^2))^(1.4/.4) ) * ( 1/((2.8/2.4)*M1^2-.4/2.4) )^(1/.4)
            Astar2 = At*P0/P02
            Mach = machfromarea(A[j]/Astar2,gamma,false)
            push!(theoryPs,pressurefrommach(Mach,P02,gamma))
        else
            push!(theoryPs,pressurefrommach(Mach,P0,gamma))
        end
        
    end
    push!(Pchokedtheory,theoryPs)# =#
end
shockAs = @. 0.0254*(0.4 + 0.04762*recordedshocklocations)
M1s = @. machfromarea(shockAs/At,gamma,true)
P02s = @. P0*( ((1.2*M1s^2)/(1+.2*M1s^2))^(1.4/.4) ) * ( 1/((2.8/2.4)*M1s^2-.4/2.4) )^(1/.4)
Astar2s = @. At*P0/P02s
Mes = @. machfromarea(Ae/Astar2s,gamma,false)
theoreticalPback = @. pressurefrommach(Mes,P02s,gamma)

#Ideally expanded
Mi = machfromarea(Ae/At,gamma,true)
Pi = pressurefrommach(Mi,P0,gamma)

red = RGB(0.8,0,0)
blue = RGB(0,0.1,0.7)
purple = RGB(0.4,0,0.4)
green = RGB(0,0.7,0.)
p3 = plot(x,(1/P0)*Pchokedtheory, labels = ["Theory 1" "Theory 2" "Theory 3" "Theory 4"],c = [red green blue purple], legend = :outertopright,size = (800,350),leftmargin = 5mm,bottommargin = 5mm) #TODO: fix plot
plot!(x,(1/P0)*chokeddatasets, labels = ["Experimental 1" "Experimental 2" "Experimental 3" "Experimental 4"],c = [red green blue purple],markershape = :circle, linestyle = :dash)
xlabel!("Nozzle location (in)")
ylabel!("P/Pambient")

partA = (Pe,mdotmax(P0,T0,At,gamma,R))
partB = plot(x,(1/P0)*Pprechokedtheory, labels = ["Theory 1" "Theory 2"],c = [red blue])
plot!(x,(1/P0)*prechokedatasets[1:2],labels = ["Experimental 1" "Experimental 2"],c = [red blue],markershape = :circle, linestyle = :dash)
xlabel!("Nozzle location (in)")
ylabel!("P/Pambient")
partC=(theoryshocklocation, recordedshocklocations,p3)
partD = (theoreticalPback,shockPbs)
partE = Pi


#p1 = plot(x,datasets)
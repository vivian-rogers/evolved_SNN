push!(LOAD_PATH, "./src/")
push!(LOAD_PATH, "./")

using Models
using SNN
using Random



maxt = 5*10^-7


layerSizes = [4,20,50,1000,40,20,5,3]
iRes = 4
n_res = layerSizes[iRes]

M = initModel(layerSizes)
addLSM(M,iRes,round.(rand(n_res,n_res)))

#printModel(M)

function stimulus(t)
	#return rand(layerSizes[1]) 
	cutoff = maxt/2
	return max(cutoff - t,0)*cutoff^-1*rand(layerSizes[1]) 
end
test(M,stimulus,stimulus,0.5*10^-9,maxt,zeros(n_res),zeros(n_res),true)

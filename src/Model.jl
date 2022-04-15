using LinearAlgebra
using Random

module SNN

struct Model
	Rₚ::Float64
	n::Int 
	TMR₀::Float64
	Pinhib::Float64
	Vb::Float64
	LayerSizes::Vector{Neuron} # Vector of sizes of each layer
	W::Matrix{Matrix} # matrix of matrices of weights
	Jₙ::Matrix{Matrix} # matrix of ones for each layer
	G₀::Matrix{Matrix} # matrix of matrices of conductances
	Boolmask::Matrix{Bool} # matrix of boolean masks
	ConnectionList::Vector
end


function randWeights(m,n,nweights,Pinhibitory)
	W = Matrix{Rational{Int}}(m,n)
	for i in i = 1:m
		if(rand() > Pinhibitory)
			S = -1
		else
			S = 1
		end
		for j in 1:n
			W[i,j] = S*rand(0:nweights)//nweights
		end
	end
	return W
end

function Ggen(M::Model)
	for C in ConnectionList
		M.G₀[C[2],C[1]] = M.Rₚ^-1*M.W[C[2],C[1]] + (M.TMR₀*M.Rₚ + M.Rₚ)*(M.Jₙ[C[2],C[1]] .- M.W[C[2],C[1]])
	end
end

function initModel(LayerSizes,nweights=1,Pinhib=0.05,Rₚ=1,TMR₀=2.00,Vb = 0.5)
	n = size(LayerSizes)[1]

	# weights for the weight matrix of each layer
	W = Matrix{Matrix}(n,n)
	
	# Conductances for each forward layer
	G₀ = Matrix{Matrix}(n,n)
	# Can this layer be modified?
	Boolmask = Bool.(zeros(n,n))
	
	# list of all of the layers with connections 
	ConnectionList = []
	# set up the forward prop
	for i = 2:n
		mᵢ = LayerSizes[i-1]; nᵢ = LayerSizes[i] 
		push!(ConnectionList,[i-1,i])
		W[i,i-1] = randWeights(mᵢ,nᵢ,nweights,Pinhib)
		G₀[i,i-1] = zeros(mᵢ,nᵢ) # initialize it to 0 for now
		Jₙ[i,i-1] = ones(mᵢ,nᵢ) # initialize it to 0 for now
		Boolmask[i,i-1] = 1
	end
	tempModel = Model(Rₚ, n, TMR₀, Pinhib, Vb, LayerSizes, W, Jₙ, G₀, Boolmask, ConnectionList)
	Ggen(tempModel)
	return tempModel
	

	

	

function update_neuron(n::Neuron,dt::Float64) # updates neuron, sends out spike if saturated
	if(n.integrated < n.threshold)
		n.integrated -= dt*n.decay_rate
		return 0
	elseif
		if(firinglock == true)
			if(n.t_since_sat > n.sat_period) #
				n.lock == true
				n.
			else
				return n.output
		else

		end
	end
end


function +(n::Neuron, in::Float64)
	if(n.firinglock == false)
		integrated = max(in,in+integrated)
	end
end


module Models

export initModel, Ggen, addLSM

push!(LOAD_PATH, "./")
using LinearAlgebra
using Random



mutable struct Model
	Rₚ::Float64
	n::Int 
	TMR₀::Float64
	Pinhib::Float64
	Vb::Float64
	LayerSizes::Vector{Int} # Vector of sizes of each layer
	W::Matrix{Matrix} # matrix of matrices of weights
	Jₙ::Matrix{Matrix} # matrix of ones for each layer
	G₀::Matrix{Matrix} # matrix of matrices of conductances
	Boolmask::Matrix{Bool} # matrix of boolean masks
	ConnectionList::Vector
	iRes::Int # index of the reservoir in the W
end


function printModel(M::Model)
	println("===== Model structure ====")
	println("Misc network parameters: Vb = $(M.Vb), Rₚ = $(M.Rₚ), TMR₀ = $(M.TMR₀*100)")
	println("Layer structure: ")
	show(M.LayerSizes)
	if(M.iRes > 0)
		println("\n$(M.n) x $(M.n) reservoir on layer $(M.iRes)")
	end
	print("\nWeight matrices:")
	for C in M.ConnectionList
		print("\n\nLayer $(C[1]) to $(C[2]) connection:")
		display(round.(M.W[C[2],C[1]]; sigdigits = 4))
	end
	println("")
	println("")
end

function randWeights(m,n,nweights,Pinhibitory)
	W = Matrix{Rational{Int}}(undef,m,n)
	for i = 1:m
		if(rand() < Pinhibitory)
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
	for C in M.ConnectionList
		M.G₀[C[2],C[1]] = M.Rₚ^-1*M.W[C[2],C[1]] + (M.TMR₀*M.Rₚ + M.Rₚ)*(M.Jₙ[C[2],C[1]] .- M.W[C[2],C[1]])
	end
end

function initModel(LayerSizes,nweights=1,Pinhib=0.05,Rₚ=1,TMR₀=2.00,Vb = 0.5)
	n = size(LayerSizes)[1]

	# weights for the weight matrix of each layer
	W = Matrix{Matrix}(undef,n,n)
	
	# Conductances for each forward layer
	G₀ = Matrix{Matrix}(undef,n,n)
	Jₙ = Matrix{Matrix}(undef,n,n)
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
	tempModel = Model(Rₚ, n, TMR₀, Pinhib, Vb, LayerSizes, W, Jₙ, G₀, Boolmask, ConnectionList,0)
	Ggen(tempModel)
	return tempModel
end

function addLSM(M::Model,iRes::Int, W_res::Matrix)
	M.iRes = iRes
	M.W[iRes,iRes] = W_res
	push!(M.ConnectionList,[iRes,iRes])
	# Could train the reservoir, but don't
	return true
end	


end

using LinearAlgebra
using Random

module SNN

struct Model
	name
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
	Neurons::
end



function Error(Output,Expected)
	sum = 0
	n_points = size(Output)[1]
	for i in eachindex(Output)
		sum += √(sum((Output[i] - Expected[i]).^2))
	end
	return sum/n_points
end

function mutate(M::Model,P::Float64 = 0.02)
	M_new = deepcopy(M)
	for i in 1:M.n
		for j in 1:M.n
			# Can this sub-matrix be modified?
			if(Boolmask[i,j])
				W = M.W[i,j]
				mᵢ = size(W)[1]; nᵢ = size(W)[2]
				for iW = 1:nᵢ
					# flip the neuron to inhibit?
					if(rand() < M.Pinhib)
						S = -1
					else
						S = 1
						W[:,iW] = abs.(W[:,iW])
					end
					for jW = 1:mᵢ
						W[nᵢ,mᵢ] = S
					

					end
				end
			end
		end
	end
	return M_new
end

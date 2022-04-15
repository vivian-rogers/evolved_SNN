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
	Neurons::
end

struct Network
	M::Model
	Neurons::Vector{Neuron}
end

function initNeurons()
	# basically make a default constructor for a bunch of neurons
	#return vector of neurons
end

function ForwardProp(M::Model,stimulus,dt,nt,Q::Float64=1)
	net = Network(M,initNeurons())
	AllSpikes = Any[]
	OutSpikes = BitArray[]
	Times = [t for t in 0:dt:(nt*dt)]
	for input in stimulus
		# update neurons
		# Spikesᵢ (Vector of bitarrays) = updateNeurons(net.Neurons)
		run(`clear`)
		display(Spikesᵢ)
		push!(AllSpikes,Spikesᵢ)
		OutSpikesᵢ = Spikesᵢ[M.n]
		push!(OutSpikes,OutSpikesᵢ)
		# encode to pulse train somehow
		# M.Neurons + inSpikes
		for C in net.M.ConnectionList
			L₁ = C[2]; L₀ = C[1];
			net.Neurons[Lᵢ] + Q*G₀[L₁,L₀]*Spikesᵢ
		end
	end
	# decode OutSpikes to data 
	# Output = Decode(dt,OutSpikes)
	return OutSpikes
end

			

using CSV
using DataFrames
using PlotlyBase
using Random
using StatsBase
using Graphs
using Dates
using JLD2
using Logging
using SparseArrays
using Distributions
using Parameters


######################################################################################################################################
# bacteria_na.jl includes the code for the bacterial submodel, a standalone ABM that is an attribute of each infected cattle host
######################################################################################################################################

# Components of the model will be imported as modules. This submodel is loaded as a precompiled module

"""
Module - bact_data
Contains a vector of possible neighbouring positions for each bacterial agent in the 33*33 grid.
"""
module bact_data
    using JLD2
    @load "./data/bact_neighbours.jld2" 
end


"""
Module - randoms
Includes packages for PRNG.
Includes a struct of defined vectors of random numbers. These are sampled for Bernoulli evaluations within the bacterial ABM. This a faster alternative to generating a new random number in this submodel
"""
module randoms
using StatsBase
using Random
export randsamplers

struct RandSamplers
    randoms::Vector{Float64}
    clinical_colonised::Vector{Vector{Int}}
    subclinical_colonised::Vector{Vector{Int}}
    get_fit::Vector{Float64}
    clincal_untreated_calf_mortality::Vector{Float64}
    clinical_treated_calf_mortality::Vector{Float64}
    clinical_untreated_mortality::Vector{Float64}
    clinical_treated_mortality::Vector{Float64}
    carrier_colonised::Vector{Vector{Int}}
end



randsamplers = RandSamplers(
 [rand() for _ in 1:1089],
[sample(1:1089, Int(round(rand((0.15:0.01:0.30))*1089)), replace = false) for _ in 1:100],
[sample(1:1089, Int(round(rand((0.1:0.01:0.2))*1089)), replace = false) for _ in 1:100],
[rand(0.96:0.001:0.99) for _ in 1:1089],
[rand(0.10:0.001:0.12) for _ in 1:10000],
[rand(0.6:0.001:0.72) for _ in 1:10000],
[rand(0.04:0.001:0.05) for _ in 1:10000],
[rand(0.024:0.001:0.03) for _ in 1:10000],
[sample(1:1089, Int(round(rand((0.01:0.001:0.1))*1089)), replace = false) for _ in 1:100]



)

end

"""
BacterialAgent
An agent.
A constructor for defining bacterial agents in each model
    Fields:
    * id : Integer ID of the colony
    * status: 8-bit integer of status (0 = not salmonella, 1 = sensitive salmonella, 2 = resistant salmonella)
    * fitness: arbitrary Float64 describing the firness of the colony in relationship to other colonies.
    * processed: placeholder bool describing whether the colony has been processed in each model step.
"""
mutable struct BacterialAgent 
    id::Int16
    status::Int8
    fitness::Float16
    processed::Bool
end

"""
BacterialModel
A constructor for a bacterial ABM.
    Fields:
    * pop_r: proportion of the population that is resistant
    * pop_p: proportion of the population that is sensitive
    * pop_s: proportion of the population that is other, non-pathogenic bacterial colonies
    * colonies: An array of agents of type BacterialAgent
    * total_status: The total status of the animal hosting the BacteiralModel. Set by events and attributes of the AnimalAgent
    * days_treated: The number of days the AnimalAgent has been treated with antimicrobials
    * days_exposed: The number of days that have elapsed since the AnimalAgent was exposed to salmonella
    * days_recovered: The number of days that the AnimalAgent has been recovering from a salmonella infection
    * days_carrier: The number of days that the AnimalAgent has been a carrier of salmonella
    * seed: Seed for RNG
    * clinical: Does the AnimalAgent have a clinical salmonella infection
"""
mutable struct BacterialModel
    pop_r::Float64
    pop_s::Float64
    pop_p::Float64
    pop_d::Float64
    colonies::Vector{BacterialAgent}
    total_status::Int
    days_treated::Int
    days_exposed::Int
    days_recovered::Int
    days_carrier::Int
    stress::Bool
    seed::Int
    rng::Xoshiro
    clinical::Bool
end

"""
BacterialData
Struct for creating BacterialData objects that hold information about the bacterial model.
    Fields:
    * id: Identification of the AnimalAgent that the BacterialModel is hosted in
    * timestep: Time elapsed in the model
    * pop_r/pop_s/pop_p/pop_d: Proportion of resistant, non-pathogenic, sensitive and dead bacterial colonies in the model
"""
 struct BacterialData
    id::Array{Int}
    timestep::Array{Int}
    pop_r::Array{Float64}
    pop_s::Array{Float64}
    pop_p::Array{Float64}
    pop_d::Array{Float64}
end

"""
bernoulli_less
Function for bernoulli trial with an outcome of less than p
"""
bernoulli_less(p, rng) = rand(rng) < p

"""
bernoulli_more
Function for bernoulli trial with an outcome of greater than p
"""
bernoulli_more(p, rng) = rand(rng) > p
#Utility ====================================================================================

"""
initialiseBacteria()
Function to initialise a bacterialModel of type BacterialModel
kwargs
* total_status: The infection status of the host AnimalAgent
* days_treated: Number of days the AnimalAgent has been treated with antibiotics for
* days_exposed: Number of days since the AnimalAgent was exposed to salmonella
* days_recovered: Number of days that the AnimalAgent has been recovering from Salmonella infection
* stress: Is the AnimalAgent in a stress period?
* seed: Sed for RNG
* rng: RNG
"""
function initialiseBacteria(;
    total_status::Int = Int(0),
    days_treated::Int = Int(0),
    days_exposed::Int = Int(0),
    days_recovered::Int = Int(0),
    stress::Bool = false,
    seed::Int = Int(42),
    rng::Xoshiro = Xoshiro(seed))

       # Agent space =========================================================================================
    # Create the colonies object
    colonies = [BacterialAgent(0,0,0,0) for _ in 1:1089]
    #Create the initial model properties ==========================================================================
    pop_r = Float64(0.0)
    pop_s = Float64(0.0)
    pop_p = Float64(0.0)
    pop_d = Float64(0.0)
    days_carrier = 0
    clinical = false

        for i in eachindex(colonies)
            colonies[i].id = i
            colonies[i].status = 0
            colonies[i].fitness = rand(rng, 0.98:0.001:0.99)
            colonies[i].processed = false
        end

    #Create an object of type BacterialModel
    bacterialModel = BacterialModel( pop_r, pop_s, pop_p, pop_d, colonies, total_status, days_treated, days_exposed, days_recovered, days_carrier, stress, seed, rng, clinical)
        
    bacterialModel.pop_r = 0
    bacterialModel.pop_s = 0
    bacterialModel.pop_p = 0
    bacterialModel.pop_d = 0 
    
    #Enumerate the colonies

    total_pop = length(bacterialModel.colonies)
    colony_status = countmap([colony.status for colony in bacterialModel.colonies])

    bacterialModel.pop_p = get(colony_status, 1, 0)/total_pop
    bacterialModel.pop_r = get(colony_status, 2, 0)/total_pop
    bacterialModel.pop_s = get(colony_status, 0, 0)/total_pop
    

return bacterialModel

end

# Example model
#@time bacterialModel =  initialiseBacteria(total_status = Int(0), days_treated = Int(0), days_exposed = Int(0), days_recovered = Int(0), stress = false, seed = Int(42));


"""
bact_treatment!
Allow treatment of the host AnimalAgent 
* Resistant colonies are not impacted by treatment
* The effect of treatment is mediated by the number of days that the animals have been treated for.
* The status of the colony changes to 10 (dead) when they have been treated
"""
function bact_treatment!(bacterialModel::BacterialModel, colony::BacterialAgent)
    colony.processed == true && return
    bacterialModel.days_treated == 0 && return
    colony.status > 1 && return
    bernoulli_more(1 - ℯ^(-bacterialModel.days_treated/5), bacterialModel.rng) && return
    colony.status = 10
    colony.fitness = 0
    colony.processed = true 
end

"""
bact_repopulate!
Repopulate after treatment.
After treatment, the remaining colonies compete for posititions on the grid, unless they are dead.
The vacant squares are colonised by a neighbour.
"""
function bact_repopulate!(bacterialModel::BacterialModel, colony::BacterialAgent)
    colony.processed == true && return
    colony.status != 10 && return # Dead colonies look for their neighbours
    
    competing_neighbour = bacterialModel.colonies[sample(bacterialModel.rng, bact_data.possible_bact_neighbours[colony.id])]

    
        if bacterialModel.days_treated != 0 #New colonisation with susceptible bacteria does not occur in the face of treatment.
            if competing_neighbour.status == 2
                rand(bacterialModel.rng, Bool) == false && return
                colony.status = 2
                colony.processed = true
                colony.fitness = competing_neighbour.fitness
            end
        elseif (bacterialModel.days_treated == 0 && bacterialModel.total_status ≤ 1)  #Otherwise they simply acquire the status of their neighbour.
            rand(bacterialModel.rng, Bool) == false && return
            colony.status = competing_neighbour.status
            colony.processed = true
            colony.fitness = competing_neighbour.fitness
        end

end



"""
bact_export!
Export bacterial data at the colony level. This is used for troubleshooting, or the amount of data generated will be excessive.
"""
function bact_export!(bacterialModel, bacterialData)
    push!(bacterialData.id, bacterialModel.id)
    push!(bacterialData.timestep, bacterialModel.timestep)
    push!(bacterialData.pop_r, bacterialModel.pop_r)
    push!(bacterialData.pop_s, bacterialModel.pop_s)
    push!(bacterialData.pop_p, bacterialModel.pop_p)
    push!(bacterialData.pop_d, bacterialModel.pop_d)
end
     
"""
bact_carrier!
Set the bacterial population of carrier animals.
This is only active if the host AnimalAgent is in a carrier status (with a total status of 5 or 6. ).
The proportion and quantity of the carrier state is determined by the level of expression by the carrier animal.
This is randomly assigned and subsequently modified by the bacterialModel. It is not an emergeny property of the model.
"""

function bact_carrier!(bacterialModel::BacterialModel)
    bacterialModel.days_carrier == 0 && return
    bacterialModel.total_status != 5 && bacterialModel.total_status != 6 && return
    to_status = bacterialModel.total_status == 5 ? 1 : 2
    [x.status = 0 for x in bacterialModel.colonies]
    colonised = sample(bacterialModel.rng, randoms.randsamplers.carrier_colonised)
    for colony in colonised 
        bacterialModel.colonies[colony].status = to_status 
        bacterialModel.colonies[colony].processed = true
        bacterialModel.colonies[colony].fitness = rand(bacterialModel.rng, 0.98:0.001:0.99)

    end
end

"""
bact_fitness!
Competition between bacterial colonies.
When two colonies meet, they compete. The fittest colony survives this challenge and goes on to colonise the square of the outcompeted neighbour.
"""
function bact_fitness!(bacterialModel::BacterialModel, colony::BacterialAgent)
    colony.processed == true && return
    competing_neighbour = bacterialModel.colonies[sample(bacterialModel.rng, bact_data.possible_bact_neighbours[colony.id])]
    colony.fitness > competing_neighbour.fitness && return
    rand(bacterialModel.rng, Bool) == false && return
    colony.status = competing_neighbour.status
    colony.fitness = competing_neighbour.fitness
end

"""
bact_processed
ensure bacteria do not get processed twice by marking the processed flag.
"""
function bact_processed!(colony::BacterialAgent)
    colony.processed = false
end
"""
bact_timestep!
step bacterialModel time
"""
function bact_timestep!(bacterialModel::BacterialModel)
    bacterialModel.timestep += 1
end

"""
bact_exposed
When an AnimalAgent is exposed, the bacterial type of the infecting animal is introduced into its BacterialModel.
"""
function bact_infection!(bacterialModel::BacterialModel)
    bacterialModel.days_exposed != 1 && return
    bacterialModel.total_status != 3 && bacterialModel.total_status != 4 && return
    to_status = bacterialModel.total_status == 3 ? 1 : 2
    if bacterialModel.clinical == true 
        colonised = sample(bacterialModel.rng, randoms.randsamplers.clinical_colonised)
    else
        colonised = sample(bacterialModel.rng, randoms.randsamplers.subclinical_colonised)
    end
    for colony in colonised 
        bacterialModel.colonies[colony].status = to_status 
        bacterialModel.colonies[colony].processed = true
        if bacterialModel.clinical == false
           bacterialModel.colonies[colony].fitness = randoms.randsamplers.get_fit[colony]

        end
    end
end

"""
find_unprocessed_infected
Search neighbouring indices and find uninfected colonies.
"""
function find_unprocessed_infected(colonies::Vector{BacterialAgent})
    indices = Vector{Int}()
    for i in eachindex(colonies)
        if !colonies[i].processed && (colonies[i].status == 1 || colonies[i].status == 2)
            push!(indices, i)
        end
    end
    return indices
end


"""
bact_latent_period
Period of latent infection after effective transmission occurs for a cattle host.
"""
function bact_latent_period!(bacterialModel::BacterialModel, colony::BacterialAgent)
    colony.processed == true && return
    bacterialModel.days_exposed < 2 && return
    bacterialModel.total_status != 3 && bacterialModel.total_status != 4 && return
    colonies = find_unprocessed_infected(bacterialModel.colonies)
    for neighbour in bact_data.possible_bact_neighbours[colony.id]
                    colony.status != 1 && colony.status != 2 && continue
                    competing_neighbour = bacterialModel.colonies[neighbour]              
                    competing_neighbour.status != 0 && competing_neighbour.status != 10  && continue
                    colony.fitness < competing_neighbour.fitness && continue
                    competing_neighbour.status = colony.status
                    competing_neighbour.processed = true
        end
    end


"""
bact_recovery!(bacterialModel)
Immune response to pathogenic bacteria. Salmonella are gradually removed by the host immune system over a number of days.
"""
function bact_recovery!(bacterialModel::BacterialModel, colony::BacterialAgent)
    bacterialModel.days_recovered == 0 && return
    colony.processed == true && return
    colony.status != 1 && colony.status != 2 && return
    bernoulli_more((1 - ℯ^(-bacterialModel.days_recovered/rand(bacterialModel.rng, 1:5))), bacterialModel.rng) && return
    colony.status = 0
    colony.processed = true
end

"""
bact_step!
Update attributes over time.
Shuffle the agents  between steps.
"""
function bact_step!(bacterialModel::BacterialModel)
    shuffle!(bacterialModel.colonies)
    bact_processed!.(bacterialModel.colonies)#Reset the processing flag
    bact_infection!(bacterialModel) 
    bact_carrier!(bacterialModel) 
  for colony in bacterialModel.colonies
        colony.processed == true && continue
        if bacterialModel.days_exposed > 1 bact_latent_period!(bacterialModel, colony) end
        if bacterialModel.days_treated != 0 bact_treatment!(bacterialModel, colony) end
        bact_repopulate!(bacterialModel, colony)
        bact_recovery!(bacterialModel, colony)
        bact_fitness!(bacterialModel, colony)
   end 
      
   bacterialModel.pop_r = 0
   bacterialModel.pop_s = 0
   bacterialModel.pop_p = 0
   bacterialModel.pop_d = 0

   total_pop = length(bacterialModel.colonies)
   colony_status = countmap([colony.status for colony in bacterialModel.colonies])

   bacterialModel.pop_p = get(colony_status, 1, 0)/total_pop
   bacterialModel.pop_r = get(colony_status, 2, 0)/total_pop
   bacterialModel.pop_s = get(colony_status, 0, 0)/total_pop

end

####################################################################################################################
### old_animal.jl defines the animal (within-herd) model 
# This is precompiled and loaded as a module 
####################################################################################################################

"""
data
A module containing data on potential neighbour positions on the grid and a randomly generated culling curve.
"""
module modata
  using JLD2
  using Distributions
  @load "./data/possible_neighbours.jld2"
  cull_curve = sort(([rand(truncated(Rayleigh(0.025), 0.00, 0.025)) for i in 1:(20*365)]./365))
end

# Load the data module 
import .modata


"""
Pens
A constructor defining each calf pen.
Fields:
  * step: Model timestep
  * pen: An array of pen numbers
  * num: Positions occupied in the pen?
"""
@with_kw mutable struct Pens
  step::Array{Int} = []
  pen::Array{Int} = []
  num::Array{Int} = []
end

"""
Infections
A constructor for defining an object containing information on infections tracked over time, to be exported at the end of each model run.
Fields:
* id: animal id
* status: infection status of animal
* step: model timestep
* stage: animal lifestage
* clin: Whether animal is clinical
* death: Whether the animal died in that step.
* days_inf: Number of days infected.
* days_exposed: Number of days exposed.
* vaccinated: vaccination status
* fpt: whether or not animal has failure of passive transfer
* age: age of animal in days
* cull_reason: reason animal culled from herd if culled.
"""
@with_kw mutable struct Infections
  id::Array{Int}  = []
  status::Array{Int}  = []
  step::Array{Int}  = []
  stage::Array{Int}  = []
  clin::Array{Bool}  = []
  death::Array{Bool}  = []
  days_inf::Array{Int}  = []
  days_exposed::Array{Int}  = [] 
  vaccinated::Array{Bool}  = [] 
  fpt::Array{Bool}  = []
  age::Array{Int}  = []
  cull_reason::Array{Symbol} = []
end


"""
Agent type - AnimalAgent
Define an agent of type AnimalAgent.
Fields:  
id::Int - Animal ID
pos::Tuple{Int, Int, Int} - Position in model space
status::Int - Infection status
stage::Int - Lifestage
days_infected::Int - Num days infected
days_exposed::Int - Num days exposed 
days_carrier::Int - Num days carrier
days_recovered::Int - Num days recovered
days_treated::Int - Num days treated
treatment::Bool - Bool indicating treatment status
pop_p::Float64 - Proportion infected with sensitive salmonella
pop_d::Float64 - Proportion dead
pop_r::Float64 - Proportion infected with resistant salmonella
bacteriaSubmodel::Union{Nothing,BacterialModel} - bacterial population of host animal
dic::Int - Days in calf
dim::Int - Days in milk
stress::Bool - Bool indicating whether animal is stressed
sex::Int - Sex, M/F
calving_season::Int - Calving season animal due to calve in based on birth
age::Int - Age in days
lactation::Int - Lactation number
pregstat::Int - Pregnancy status (0 = Empty, 1 = in calf)
trade_status::Int - Whether animal can be traded
neighbours::NTuple{8, Tuple{Int, Int, Int}} - Neighbouring animals at given step
processed::Bool - Flag - has the animal been processed in this step?
carryover::Bool - Flag - has the animal been carried over between calving seasons?
fpt::Bool - Flag - has the animal had failure of passive transfer
vaccinated::Bool - Flag. Has the animal been vaccinated?
susceptibility::Float64 - Float describing susceptibility to infection. Mediated by vacc, fpt, prior infection status. Proxy for immunity. Emergent.
clinical::Bool - Flag. Is the animal's infection clinical?
pen::Int - Pen calf present in 
times_treated::Int - Number of times the animal has been treated.
"""
mutable struct AnimalAgent
    id::Int
    pos::Tuple{Int, Int, Int}
    status::Int
    stage::Int
    days_infected::Int
    days_exposed::Int
    days_carrier::Int
    days_recovered::Int
    days_treated::Int
    treatment::Bool
    pop_p::Float64
    pop_d::Float64
    pop_r::Float64
    bacteriaSubmodel::Union{Nothing,BacterialModel}
    dic::Int
    dim::Int
    stress::Bool
    sex::Int
    calving_season::Int
    age::Int
    lactation::Int
    pregstat::Int
    trade_status::Int
    neighbours::NTuple{8, Tuple{Int, Int, Int}}
    processed::Bool
    carryover::Bool
    fpt::Bool
    vaccinated::Bool
    susceptibility::Float64
    clinical::Bool
    pen::Int
    times_treated::Int
end

"""
Transmissions: A constructor for containing vectors detailing animal transmissions on each simulation day
"""
@with_kw mutable struct Transmissions
  step::Array{Int} = []
  from_id::Array{Int}= []
  to_id::Array{Int}= []
  stage::Array{Int}= []
  from::Array{Int}= []
  to::Array{Int}= []
  type::Array{Symbol}= []
  effective::Array{Bool}= []
  clinical::Array{Bool}= []
end

"""
AnimalData
Struct for animal data. Contains details on infectons and population dynamics on each simulation day.
This implicitly limits the duration of the model to 2555 days (five years). Can be changed, but changes are also required elsewhere.
"""
 @with_kw mutable struct AnimalData
    id::Array{Int} = Array{Int}(zeros(2556))
    timestep::Array{Int} = Array{Int}(zeros(2556))
    pop_r::Array{Int} = Array{Int}(zeros(2556))
    pop_s::Array{Int} = Array{Int}(zeros(2556))
    pop_p::Array{Int} = Array{Int}(zeros(2556))
    pop_d::Array{Int} = Array{Int}(zeros(2556))
    pop_rec_r::Array{Int} = Array{Int}(zeros(2556))
    pop_rec_p::Array{Int} = Array{Int}(zeros(2556))
    pop_car_p::Array{Int} = Array{Int}(zeros(2556))
    pop_car_r::Array{Int} = Array{Int}(zeros(2556))
    num_calves::Array{Int} = Array{Int}(zeros(2556))
    num_weaned::Array{Int} = Array{Int}(zeros(2556))
    num_dh::Array{Int} = Array{Int}(zeros(2556))
    num_heifers::Array{Int} = Array{Int}(zeros(2556))
    num_lactating::Array{Int} = Array{Int}(zeros(2556))
    num_dry::Array{Int} = Array{Int}(zeros(2556))
    pop_er::Array{Int} = Array{Int}(zeros(2556))
    pop_ep::Array{Int} = Array{Int}(zeros(2556))
    
    inf_calves::Array{Int} = Array{Int}(zeros(2556))
    car_calves::Array{Int} = Array{Int}(zeros(2556))
    ri_calves::Array{Int} = Array{Int}(zeros(2556))
    rc_calves::Array{Int} = Array{Int}(zeros(2556))

    inf_heifers::Array{Int} = Array{Int}(zeros(2556))
    car_heifers::Array{Int} = Array{Int}(zeros(2556))
    ri_heifers::Array{Int} = Array{Int}(zeros(2556))
    rc_heifers::Array{Int} = Array{Int}(zeros(2556))

    inf_weaned::Array{Int} = Array{Int}(zeros(2556))
    car_weaned::Array{Int} = Array{Int}(zeros(2556))
    ri_weaned::Array{Int} = Array{Int}(zeros(2556))
    rc_weaned::Array{Int} = Array{Int}(zeros(2556))

    inf_dh::Array{Int} = Array{Int}(zeros(2556))
    car_dh::Array{Int} = Array{Int}(zeros(2556))
    ri_dh::Array{Int} = Array{Int}(zeros(2556))
    rc_dh::Array{Int} = Array{Int}(zeros(2556))

    inf_dry::Array{Int} = Array{Int}(zeros(2556))
    car_dry::Array{Int} = Array{Int}(zeros(2556))
    ri_dry::Array{Int} = Array{Int}(zeros(2556))
    rc_dry::Array{Int} = Array{Int}(zeros(2556))


    inf_lac::Array{Int} = Array{Int}(zeros(2556))
    car_lac::Array{Int} = Array{Int}(zeros(2556))
    ri_lac::Array{Int} = Array{Int}(zeros(2556))
    rc_lac::Array{Int} = Array{Int}(zeros(2556))

    clinical::Array{Int} = Array{Int}(zeros(2556))
    subclinical::Array{Int} = Array{Int}(zeros(2556))
    current_b1::Array{Int} = Array{Int}(zeros(2556))
    current_b2::Array{Int} = Array{Int}(zeros(2556))
    current_b3::Array{Int} = Array{Int}(zeros(2556))
    current_b4::Array{Int}  = Array{Int}(zeros(2556))

end


"""
AllData is an empty constructor that can be used to record information on every animal on every simulation day.
This is intensive in terms of processing time and data usage. It is turned off by default.
"""
@with_kw mutable struct AllData
  id::Array{Int} = []
  step::Array{Int} = []
  stage::Array{Int} = []
  pregstat::Array{Int} = []
  status::Array{Int} = []
  dic::Array{Int} = []
  dim::Array{Int} = []
  age::Array{Int} = []
  days_infected::Array{Int} = []
  days_exposed::Array{Int} = []
  days_recovered::Array{Int} = []
  days_carrier::Array{Int} = []
  days_treated::Array{Int} = []
  treatment::Array{Bool} = []
  pop_r::Array{Float64} = []
  pop_p::Array{Float64} = []
  susceptibility::Array{Float64} = []
end

"""
Environment is an empty constructor for recording the state of the environment at each model step.
"""
@with_kw mutable struct Environment
  contamination::Vector{SparseMatrixCSC{Float64, Int}} = [spzeros(250,250) for i in 1:25]
  contam_time::Vector{SparseMatrixCSC{Float64, Int}} = [spzeros(250,250) for i in 1:25]
  contam_type::Vector{SparseMatrixCSC{Float64, Int}} = [spzeros(250,250) for i in 1:25]
end

"""
AnimalModel
Type container for animal model.
This struct is "the model" and brings together all components of the model than can then be executed using the animal_model! function
"""
  mutable struct AnimalModel
    farmno::Int
    animals::Array{AnimalAgent}
    timestep::Int
    date::Date
    rng::Xoshiro
    system::Int
    msd::Date
    msd_2::Date
    msd_3::Date
    msd_4::Date
    seed::Int
    optimal_stock::Int
    treatment_prob::Float64
    treatment_length::Int
    carrier_prob::Float64
    current_stock::Int
    current_lactating::Int
    optimal_lactating::Int
    current_heifers::Int
    optimal_heifers::Int
    current_calves::Int
    optimal_calves::Int
    current_weaned::Int
    optimal_weaned::Int
    current_dh::Int
    optimal_dh::Int
    current_dry::Int
    optimal_dry::Int
    tradeable_stock::Int
    sending::Vector{Int}
    receiving::Vector{Int}
    density_lactating::Int
    density_calves::Float64
    density_dry::Int
    positions::Array{Tuple{Int, Int, Int}}
    pop_r::Int
    pop_s::Int
    pop_p::Int
    pop_d::Int
    id_counter::Int
    vacc_rate::Float64
    fpt_rate::Float64
    prev_r::Float64
    prev_p::Float64
    prev_cr::Float64
    prev_cp::Float64
    current_autumn::Int
    optimal_autumn::Int
    current_spring::Int
    optimal_spring::Int
    current_b1::Int
    current_b2::Int
    current_b3::Int
    current_b4::Int
    optimal_b1::Int
    optimal_b2::Int
    optimal_b3::Int
    optimal_b4::Int
    sim::AnimalData
    environment::Environment
    pen_counter::Int
    calf_pen::Int 
    pen_decon::Bool
    transmissions::Transmissions
    infections::Infections
    treat_calf::Float64
    treat_dry::Float64
    treat_lac::Float64
    vacc_protection::Float64
    vacc_shedding::Float64
    vacc_duration::Float64
    p_retreatment::Float64
    cull_curve::Vector{Float64}
    longevity::Float64
    allData::AllData
    pens::Pens
    spare_bact::BacterialModel
    
end

"""
count_animals!(animalModel)
Count stock classes and types at each step of the animalModel
"""
function count_animals!(animalModel::AnimalModel)
  #Stages 
  animalModel.current_calves = count(i->(i.stage == 1), animalModel.animals)
  animalModel.current_weaned = count(i->(i.stage == 2), animalModel.animals)
  animalModel.current_dh = count(i->(i.stage == 4), animalModel.animals)
  animalModel.current_heifers = count(i->(i.stage == 3), animalModel.animals)
  animalModel.current_lactating = count(i->(i.stage == 5), animalModel.animals)
  animalModel.current_dry = Int16(count(i->(i.stage == 6), animalModel.animals))
  animalModel.current_stock = length(animalModel.animals)

  #assigned_statuses
  animalModel.pop_p = Int16(count(i->(i.status == 1 && i.sex == 1), animalModel.animals))
  animalModel.pop_r = count(i->(i.status == 2 && i.sex == 1), animalModel.animals)
  animalModel.pop_s = count(i->(i.status == 0 && i.sex == 1), animalModel.animals)
  animalModel.pop_d = count(i->((i.status == 1 || i.status == 2) && i.sex == 1), animalModel.animals)

  animalModel.current_stock = animalModel.current_weaned + animalModel.current_calves + animalModel.current_dh + animalModel.current_heifers + animalModel.current_lactating + animalModel.current_dry 

  #Split systems

  animalModel.current_spring = count(i-> (i.calving_season == 1 && i.stage == 5), animalModel.animals)
  animalModel.current_autumn = count(i-> (i.calving_season == 2 && i.stage == 5), animalModel.animals)
  
  #Batch systems 
  animalModel.current_b1 = count(i-> (i.calving_season == 1 && i.stage == 5), animalModel.animals)
  animalModel.current_b2 = count(i-> (i.calving_season == 2 && i.stage == 5), animalModel.animals)
  animalModel.current_b3 = count(i-> (i.calving_season == 3 && i.stage == 5), animalModel.animals)
  animalModel.current_b4 = count(i-> (i.calving_season == 4 && i.stage == 5), animalModel.animals)
end





# Initialisation functions ============================================

"""
initialiseSpring
Function for generating system 1 farms (Spring calving, CS1).
Sets the state of the system at init time (which by default is July).
These form the bulk of the code. There are certainly more elegant ways of doing this.
"""
  function initialiseSpring(;
    farmno::Int = FarmAgent.id,
    system::Int,
    msd::Date,
    seed::Int,
    optimal_stock::Int,
    optimal_lactating::Int,
    treatment_prob::Float64,
    treatment_length::Int,
    carrier_prob::Float64,
    timestep::Int,
    density_lactating::Int,
    density_dry::Int,
    density_calves::Float64,
    date::Date,
    vacc_rate::Float64,
    fpt_rate::Float64,
    prev_r::Float64,
    prev_p::Float64,
    prev_cr::Float64,
    prev_cp::Float64,
    pen_decon::Float64,
    treat_calf::Float64,
    treat_dry::Float64,
    treat_lac::Float64,
    vacc_protection::Float64,
    vacc_shedding::Float64,
    vacc_duration::Float64,
    p_retreatment::Float64,
    longevity::Float64
    )

    # Initialise the objects that will be used by the model 
    spare_bact = initialiseBacteria()
    allData = AllData()
    pens = Pens()
    animals = Array{AnimalAgent}[]
    sim = AnimalData()
    environment = Environment()
    transmission = Transmissions()
    infections = Infections()

    cull_curve = modata.cull_curve
    #Create the initial model parameters ===============================
    msd_2 = msd_3 = msd_4 = Date(0)
    current_stock = current_lactating = current_dry = current_heifers = current_dh = current_weaned = current_calves = 0
    optimal_dry = optimal_heifers = optimal_dh = optimal_weaned = optimal_calves = 0
    tradeable_stock = 0
    sending = receiving = Vector{Int}()
    rng = Xoshiro(seed)
    pop_p = pop_r = pop_s = pop_d = 0
    id_counter = 0
    positions = Array{Array{Int}}[]
    processed = false
    current_spring = current_autumn = optimal_spring = optimal_autumn = current_b1 = current_b2 = current_b3 = current_b4 = 0
    optimal_b1 = optimal_b2 = optimal_b3 = optimal_b4 = 0
    #Set up the model ====================================================

    # Specify pen counting functions
    pen_counter = 0
    calf_pen = 8

    # Determine if the herd uses pen decontamination
    pen_decon = ifelse(pen_decon < 0.5, false, true)

    # Specify an instance of the model with the parameters above
    animalModel = AnimalModel(farmno, animals, timestep, date, rng, system, msd, msd_2, msd_3, msd_4, seed,  optimal_stock, treatment_prob, treatment_length, carrier_prob, current_stock, current_lactating, optimal_lactating, current_heifers, optimal_heifers, current_calves, optimal_calves, current_weaned, optimal_weaned, current_dh, optimal_dh, current_dry, optimal_dry, tradeable_stock, sending, receiving, density_lactating, density_calves, density_dry, positions, pop_r, pop_s, pop_p, pop_d, id_counter, vacc_rate, fpt_rate, prev_r, prev_p, prev_cr, prev_cp,   current_autumn, optimal_autumn, current_spring, optimal_spring, current_b1, current_b2, current_b3, current_b4, optimal_b1, optimal_b2, optimal_b3, optimal_b4, sim, environment, pen_counter, calf_pen, pen_decon, transmission, infections, treat_calf, treat_dry, treat_lac, vacc_protection, vacc_shedding, vacc_duration, p_retreatment, cull_curve, longevity, allData, pens, spare_bact)
  
    
    # Set the initial stock parameters
    animalModel.optimal_heifers = animalModel.optimal_dry = animalModel.optimal_weaned = animalModel.optimal_calves = animalModel.optimal_dh = animalModel.optimal_heifers = ceil(0.3*animalModel.optimal_lactating)
    
    # Add the dry cows ---------------------------------------------

    animalModel.id_counter = 0
  # Add the dry cows ---------------------------------------------
    #Dry stage is 6, Dry plane is 6. Model opens on day before psc
    animalModel.id_counter = 0
     for cow in 1:(animalModel.optimal_lactating - 0.9*animalModel.optimal_dh)
        animalModel.id_counter += 1
        id = Int(animalModel.id_counter)
        stage = Int(6)
        range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry))))
        pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
         while pos in animalModel.positions == true
            pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
        end 
        push!(animalModel.positions, pos)
        status = 0
        days_infected = 0
        days_exposed = Int(0)
        days_carrier = 0
        days_recovered = Int(0)
        days_treated = Int(0)
        treatment = false
        pop_d = Float64(0.0)
        bacteriaSubmodel = nothing
        dic =  Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(259), 199, 280))))   
        dim = Int(0)
        pop_p = 0
        pop_r = 0
        stress = false
        sex = 1#Female
        calving_season = 1#Spring
        age = Int(ceil(rand(animalModel.rng, truncated(Rayleigh(animalModel.longevity*365*0.75),(3*365), (10*365))))) 
        lactation= round(age/365) - 1
        pregstat = 1#Pregnant
        trade_status = 0#false
        neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
        carryover = false
        fpt = false
        vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
        susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5 : rand(animalModel.rng, 0.45:0.01:0.55)
        clinical = false
        pen = 0
        times_treated = 0
        animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
        push!(animalModel.animals, animal)
    end

# Add the heifers ---------------------------------------------
#Heifers about to calve, heifer stage 4
    for heifer in 1:animalModel.optimal_dh
        animalModel.id_counter += 1
        id = Int(animalModel.id_counter)
        stage = 4
        range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_heifers)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_heifers))))
        pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
         while pos in animalModel.positions == true
            pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
        end 
        push!(animalModel.positions, pos)
        status = 0
        days_infected = 0
        days_exposed = Int(0)
        days_carrier = 0
        days_recovered = Int(0)
        days_treated = Int(0)
        treatment = false
        pop_p = Float64(0.0)
        bacteriaSubmodel = nothing
        pop_p = 0
        pop_r = 0
        dic =  Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(259), 199, 280))))   
        dim = 0
        stress = false
        sex = 1#Female
        calving_season = 1#Spring
        age = Int(ceil(rand(animalModel.rng, truncated(Rayleigh(2*365),(22*30), (25*30))))) 
        lactation= 0
        pregstat = 1#Pregnant
        trade_status = 0#false
        neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
        carryover = false
        fpt = false
        vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
        susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5 : rand(animalModel.rng, 0.45:0.01:0.55)
        clinical = false
        pen = 0
        times_treated = 0
        animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
        push!(animalModel.animals, animal)
    end

    #Add weaned animals-------------------------------------------------------------------------------------------------------------------

    for weaned in 1:ceil(animalModel.optimal_lactating*0.2)
        animalModel.id_counter += 1
        id = Int(animalModel.id_counter)
        stage = 2
        range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_weaned)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_weaned))))
        pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
         while pos in animalModel.positions == true
            pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
        end 
        push!(animalModel.positions, pos)
        status = 0
        days_infected = status == 1 || status == 2 ? 1 : 0
        days_exposed = Int(0)
        days_carrier = 0
        days_recovered = Int(0)
        days_treated = Int(0)
        treatment = false
        pop_p = Float64(0.0)
        bacteriaSubmodel = nothing
        pop_p = 0
        pop_r = 0
        dic =  Int(0)   
        dim = 0
        stress = false
        sex = 1#Female
        calving_season = 1#Spring
        age = Int(ceil(rand(animalModel.rng, truncated(Rayleigh(365),(295), (350))))) 
        lactation= 0
        pregstat = 0#Pregnant
        trade_status = 0#false
        neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
        carryover = false 
        fpt = false
        vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
        susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
        clinical = false
        pen = 0
        times_treated = 0
        animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
        push!(animalModel.animals, animal)
    end
    #Distribute the initial infections
    num_r = Int(ceil(animalModel.prev_r*length(animalModel.animals)))
    num_p = Int(ceil(animalModel.prev_p*length(animalModel.animals)))
    num_cr = Int(ceil(animalModel.prev_cr*length(animalModel.animals)))
    num_cp = Int(ceil(animalModel.prev_cp*length(animalModel.animals)))
    
    #Set resistant 
    uninfected = findall(x-> x.status == 0, animalModel.animals)
    resistant = sample(animalModel.rng, uninfected, num_r)
    for resist in resistant
      animalModel.animals[resist].status = 4
      animalModel.animals[resist].bacteriaSubmodel = initialiseBacteria()
      animalModel.animals[resist].bacteriaSubmodel.rng = Xoshiro(resist)

      animalModel.animals[resist].bacteriaSubmodel.total_status = 4
      animalModel.animals[resist].days_exposed = 1
      animalModel.animals[resist].clinical = true
      bact_step!(animalModel.animals[resist].bacteriaSubmodel)
    end

    #Set pathogenic 
    uninfected = findall(x-> x.status == 0, animalModel.animals)
    pathogenic = sample(animalModel.rng, uninfected, num_p)
    for pathogen in pathogenic
      animalModel.animals[pathogen].bacteriaSubmodel = initialiseBacteria()
      animalModel.animals[pathogen].bacteriaSubmodel.rng = Xoshiro(pathogen)
      animalModel.animals[pathogen].status = 3
      animalModel.animals[pathogen].bacteriaSubmodel.total_status = 3
      animalModel.animals[pathogen].days_exposed = 1
      animalModel.animals[pathogen].clinical = true
      bact_step!(animalModel.animals[pathogen].bacteriaSubmodel)
    end

    #Set carrier resistant
    uninfected = findall(x-> x.status == 0, animalModel.animals)
    carrier_resistant = sample(animalModel.rng, uninfected, num_cr)
    for resist in carrier_resistant
      animalModel.animals[resist].bacteriaSubmodel = initialiseBacteria( total_status = Int(0), days_treated = Int(0), days_exposed = Int(0), days_recovered = Int(0), stress = false, seed = Int(resist))

      animalModel.animals[resist].status = 6
      animalModel.animals[resist].bacteriaSubmodel.total_status = 6
      animalModel.animals[resist].days_carrier = 1
      bact_step!(animalModel.animals[resist].bacteriaSubmodel)
    end

    #Set carrier pathogenic 
    uninfected = findall(x-> x.status == 0, animalModel.animals)
    carrier_pathogenic = sample(animalModel.rng, uninfected, num_cp)
    for pathogen in carrier_pathogenic
      animalModel.animals[pathogen].bacteriaSubmodel = initialiseBacteria()
      animalModel.animals[pathogen].bacteriaSubmodel.rng = Xoshiro(pathogen)
      animalModel.animals[pathogen].status = 5
      animalModel.animals[pathogen].bacteriaSubmodel.total_status = 5
      animalModel.animals[pathogen].days_carrier = 1
      bact_step!(animalModel.animals[pathogen].bacteriaSubmodel)
    end
    
    for animal in 1:length(animalModel.animals)
      animalModel.animals[animal].bacteriaSubmodel === nothing && continue
      animalModel.animals[animal].bacteriaSubmodel.seed = animalModel.animals[animal].id
    end
    count_animals!(animalModel)

    optimal_stock = length(animalModel.animals)
    animalModel.timestep = 1
    return animalModel

end


"""
initialiseBatch(;kwargs)
Initialise a batch calving farm, as above, but for CS3
"""
function initialiseBatch(;
    farmno::Int = FarmAgent.id,
    system::Int,
    msd::Date,
    seed::Int,
    optimal_stock::Int,
    optimal_lactating::Int,
    treatment_prob::Float64,
    treatment_length::Int,
    carrier_prob::Float64,
    timestep::Int,
    density_lactating::Int,
    density_dry::Int,
    density_calves::Float64,
    date::Date,
    vacc_rate::Float64,
    fpt_rate::Float64,
    prev_r::Float64,
    prev_p::Float64,
    prev_cr::Float64,
    prev_cp::Float64,   
    pen_decon::Float64,
    treat_calf::Float64,
    treat_dry::Float64,
    treat_lac::Float64,
    vacc_protection::Float64,
    vacc_shedding::Float64,
    vacc_duration::Float64,
    p_retreatment::Float64,
    longevity::Float64
    )
    # Initialise the objects that will be used by the model 
    spare_bact = initialiseBacteria()
    allData = AllData()
    pens = Pens()
    animals = Array{AnimalAgent}[]
    sim = AnimalData()
    environment = Environment()
    transmission = Transmissions()
    infections = Infections()
    cull_curve = modata.cull_curve

    #Create the initial model parameters ===============================
    msd_2 = msd + Month(3)
    msd_3 = msd - Month(6)
    msd_4 = msd - Month(3)

    current_stock = 0
    current_lactating = 0
    current_dry = 0
    current_heifers = 0
    current_dh = 0
    current_weaned = 0
    current_calves = 0
    optimal_dry = 0
    optimal_heifers = 0
    optimal_dh = 0
    optimal_weaned = 0
    optimal_calves = 0
    tradeable_stock = 0
    sending = receiving = Vector{Int}()
    rng = Xoshiro(seed)
    pop_p = 0
    pop_r = 0
    pop_s = 0
    pop_d = 0
    id_counter = 0
    positions = Array{Array{Int}}[]
    processed = false
    current_spring = current_autumn = optimal_spring = optimal_autumn = 0
    current_b1 = 0
    current_b2 = 0
    current_b3 = 0
    current_b4 = 0

    N = optimal_stock 
    optimal_b1 = optimal_b2 = optimal_b3 = optimal_b4 = ceil(0.25*N)
    #Set up the model ====================================================

    pen_counter = 0
    calf_pen = 8

    pen_decon = ifelse(pen_decon < 0.5, false, true)


    animalModel = AnimalModel(farmno, animals, timestep, date, rng, system, msd, msd_2, msd_3, msd_4, seed,  optimal_stock, treatment_prob, treatment_length, carrier_prob, current_stock, current_lactating, optimal_lactating, current_heifers, optimal_heifers, current_calves, optimal_calves, current_weaned, optimal_weaned, current_dh, optimal_dh, current_dry, optimal_dry, tradeable_stock, sending, receiving, density_lactating, density_calves, density_dry, positions, pop_r, pop_s, pop_p, pop_d, id_counter, vacc_rate, fpt_rate, prev_r, prev_p, prev_cr, prev_cp,   current_autumn, optimal_autumn, current_spring, optimal_spring, current_b1, current_b2, current_b3, current_b4, optimal_b1, optimal_b2, optimal_b3, optimal_b4, sim, environment, pen_counter, calf_pen, pen_decon, transmission, infections, treat_calf, treat_dry, treat_lac, vacc_protection, vacc_shedding, vacc_duration, p_retreatment, cull_curve, longevity, allData, pens, spare_bact)
    
    # Set the initial stock parameters
    animalModel.optimal_heifers = animalModel.optimal_weaned = animalModel.optimal_calves = animalModel.optimal_dh = animalModel.optimal_heifers = animalModel.optimal_dry = ceil(0.3*animalModel.optimal_lactating)
    
    # Add the dry cows ---------------------------------------------

    animalModel.id_counter = 0
    animalModel.id_counter = 0

    for b1dry in 1:ceil(optimal_b1*0.7)
      animalModel.id_counter += 1
      id = Int(animalModel.id_counter)
      stage = 6
      range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
      while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
      push!(animalModel.positions, pos)
      status = 0
      days_infected = 0
      days_exposed = Int(0)
      days_carrier = 0
      days_recovered = Int(0)
      days_treated = Int(0)
      treatment = false
      pop_p = Float64(0.0)
      bacteriaSubmodel = nothing
      pop_p = 0
      pop_r = 0
      dic =  Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(262), 199, 283))))   
      dim = 0
      stress = false
      sex = 1#Female
      calving_season = 1#Spring
      age = Int(ceil(rand(animalModel.rng, truncated(Rayleigh(animalModel.longevity*365*0.75),(3*365), (10*365))))) 
      lactation= 0
      pregstat = 1#Pregnant
      trade_status = 0#false
      neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
      carryover = false 
      fpt = false
      vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
      susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
      clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
  end

  # B1 dry heifers
  for b1heifers in 1:ceil(optimal_b1*0.3)
    animalModel.id_counter += 1
    id = Int(animalModel.id_counter)
    stage = 4
    range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry))))
    pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
    while pos in animalModel.positions == true
        pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
    end 
    push!(animalModel.positions, pos)
    status = 0
    days_infected = 0
    days_exposed = Int(0)
    days_carrier = 0
    days_recovered = Int(0)
    days_treated = Int(0)
    treatment = false
    pop_p = Float64(0.0)
    bacteriaSubmodel = nothing
    pop_p = 0
    pop_r = 0
    dic =  Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(272), 199, 283))))   
    dim = 0
    stress = false
    sex = 1#Female
    calving_season = 1#Spring
    age =  Int(ceil(rand(animalModel.rng, truncated(Rayleigh(2*365),(22*30), (25*30))))) 
    lactation= 0
    pregstat = 1#Pregnant
    trade_status = 0#false
    neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
    carryover = false 
    fpt = false
    vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
    susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
    clinical = false
    pen = 0
    times_treated = 0
    animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
    push!(animalModel.animals, animal)
end


for b1weaned in 1:ceil(optimal_b1*0.5)
  animalModel.id_counter += 1
  id = Int(animalModel.id_counter)
  stage = 2
  range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
      while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
  push!(animalModel.positions, pos)
  status = 0
  days_infected = 0
  days_exposed = Int(0)
  days_carrier = 0
  days_recovered = Int(0)
  days_treated = Int(0)
  treatment = false
  pop_p = Float64(0.0)
  bacteriaSubmodel = nothing
  pop_p = 0
  pop_r = 0
  dic =  0   
  dim = 0
  stress = false
  sex = 1#Female
  calving_season = 1#Spring
  age =  Int(ceil(rand(animalModel.rng, truncated(Rayleigh(315),(281), (365))))) 
  lactation= 0
  pregstat = 0#Pregnant
  trade_status = 0#false
  neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
  carryover = false 
  fpt = false
  vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
  susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
  clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
end

# B2 ------------------------------------------
# Lactating
for b2lac in 1:ceil(optimal_b2)
  animalModel.id_counter += 1
  id = Int(animalModel.id_counter)
  stage = 5
  range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_lactating)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_lactating))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
      while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
  push!(animalModel.positions, pos)
  status = 0
  days_infected = 0
  days_exposed = Int(0)
  days_carrier = 0
  days_recovered = Int(0)
  days_treated = Int(0)
  treatment = false
  pop_p = Float64(0.0)
  bacteriaSubmodel = nothing
  pop_p = 0
  pop_r = 0
  dim = Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(237), 189, 273))))
  stress = false
  sex = 1#Female
  calving_season = 2#Spring
  age = Int(ceil(rand(animalModel.rng, truncated(Rayleigh(animalModel.longevity*365*0.75),(2*365), (10*365))))) 
  lactation= 0
  pregstat = rand(animalModel.rng) < 0.85 ? 1 : 0#Pregnant
  dic = pregstat == 1 ? Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(153), 85, 188)))) : 0
  trade_status = 0#false
  neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
  carryover = false 
  fpt = false
  vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
  susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
  clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
end

for b2dheifers in 1:ceil(optimal_b2*0.3)
  animalModel.id_counter += 1
  id = Int(animalModel.id_counter)
  stage = 4
  range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
      while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
  push!(animalModel.positions, pos)
  status = 0
  days_infected = 0
  days_exposed = Int(0)
  days_carrier = 0
  days_recovered = Int(0)
  days_treated = Int(0)
  treatment = false
  pop_p = Float64(0.0)
  bacteriaSubmodel = nothing
  pop_p = 0
  pop_r = 0
  dim = 0
  stress = false
  sex = 1#Female
  calving_season = 2#Spring
  age =  Int(ceil(rand(animalModel.rng, truncated(Rayleigh(603),553, 638)))) 
  lactation= 0
  pregstat = 1#Pregnant
  dic = Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(174), 126, 209))))
  trade_status = 0#false
  neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
  carryover = false 
  fpt = false
  vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
  susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
  clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
end

for b2weaned in 1:ceil(optimal_b2*0.5)
  animalModel.id_counter += 1
  id = Int(animalModel.id_counter)
  stage = 2
  range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
      while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
  push!(animalModel.positions, pos)
  status = 0
  days_infected = 0
  days_exposed = Int(0)
  days_carrier = 0
  days_recovered = Int(0)
  days_treated = Int(0)
  treatment = false
  pop_p = Float64(0.0)
  bacteriaSubmodel = nothing
  pop_p = 0
  pop_r = 0
  dim = 0
  stress = false
  sex = 1#Female
  calving_season = 2#Spring
  age =   Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(237), 189, 273))))
  lactation= 0
  pregstat = 0#Pregnant
  dic = 0
  trade_status = 0#false
  neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
  carryover = false 
  fpt = false
  vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
  susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
  clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
end

# B3 -----------------------------------------------------------------------------
for b3lac in 1:ceil(optimal_b3)
  animalModel.id_counter += 1
  id = Int(animalModel.id_counter)
  stage = 5
  range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_lactating)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_lactating))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
      while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
  push!(animalModel.positions, pos)
  status = 0
  days_infected = status == 1 || status == 2 ? 1 : 0
  days_exposed = Int(0)
  days_carrier = status == 5 || status == 6 ? 1 : 0
  days_recovered = Int(0)
  days_treated = Int(0)
  treatment = false
  pop_p = Float64(0.0)
  bacteriaSubmodel = nothing
  pop_p = 0
  pop_r = 0
  dim =  Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(145), 97, 180))))
  stress = false
  sex = 1#Female
  calving_season = 3#Spring
  age = Int(ceil(rand(animalModel.rng, truncated(Rayleigh(animalModel.longevity*365*0.75),(2*365), (10*365))))) 
  lactation= 0
  pregstat = rand(animalModel.rng) < 0.85 ? 1 : 0#Pregnant
  dic = pregstat == 1 ? Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(153), 85, 188)))) : 0
  trade_status = 0#false
  neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
  carryover = false 
  fpt = false
  vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
  susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
  clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
end

for b3heifers in 1:ceil(N*0.3)
  animalModel.id_counter += 1
  id = Int(animalModel.id_counter)
  stage = 4
  range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
      while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
  push!(animalModel.positions, pos)
  status = 0
  days_infected = 0
  days_exposed = Int(0)
  days_carrier = 0
  days_recovered = Int(0)
  days_treated = Int(0)
  treatment = false
  pop_p = Float64(0.0)
  bacteriaSubmodel = nothing
  pop_p = 0
  pop_r = 0
  dim =  Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(145), 97, 180))))
  stress = false
  sex = 1#Female
  calving_season = 3#Spring
  age =  Int(ceil(rand(animalModel.rng, truncated(Rayleigh(511),463, 546))))
  lactation= 0
  pregstat = 1 #Pregnant
  dic = Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(82), 33, 117))))
  trade_status = 0#false
  neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
  carryover = false 
  fpt = false
  vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
  susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
  clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
end

for b3weaned in 1:ceil(optimal_b3*0.5)
  animalModel.id_counter += 1
  id = Int(animalModel.id_counter)
  stage = 2
  range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
      while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
  push!(animalModel.positions, pos)
  status = 0
  days_infected = 0
  days_exposed = Int(0)
  days_carrier = 0
  days_recovered = Int(0)
  days_treated = Int(0)
  treatment = false
  pop_p = Float64(0.0)
  bacteriaSubmodel = nothing
  pop_p = 0
  pop_r = 0
  dim =  0
  stress = false
  sex = 1#Female
  calving_season = 3#Spring
  age =  Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(145), 97, 180))))
  lactation= 0
  pregstat = 0 #Pregnant
  dic = 0
  trade_status = 0#false
  neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
  carryover = false 
  fpt = false
  vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
  susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
  clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
end

# B4 --------------------------------
for b4lac in 1:ceil(optimal_b4)
  animalModel.id_counter += 1
  id = Int(animalModel.id_counter)
  stage = 5
  range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_lactating)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_lactating))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
      while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
  push!(animalModel.positions, pos)
  status = 0
  days_infected = 0
  days_exposed = Int(0)
  days_carrier = 0
  days_recovered = Int(0)
  days_treated = Int(0)
  treatment = false
  pop_p = Float64(0.0)
  bacteriaSubmodel = nothing
  pop_p = 0
  pop_r = 0
  dim =  Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(55), 7, 90))))
  stress = false
  sex = 1#Female
  calving_season = 4#Spring
  age = Int(ceil(rand(animalModel.rng, truncated(Rayleigh(animalModel.longevity*365*0.75),(2*365), (10*365))))) 
  lactation= 0
  pregstat = 0 #Pregnant
  dic = 0
  trade_status = 0#false
  neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
  carryover = false 
  fpt = false
  vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
  susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
  clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
end

for b4dheifers in 1:ceil(optimal_b4*0.3)
  animalModel.id_counter += 1
  id = Int(animalModel.id_counter)
  stage = 4
  range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
      while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
  push!(animalModel.positions, pos)
  status =0
  days_infected = 0
  days_exposed = Int(0)
  days_carrier = 0
  days_recovered = Int(0)
  days_treated = Int(0)
  treatment = false
  pop_p = Float64(0.0)
  bacteriaSubmodel = nothing
  pop_p = 0
  pop_r = 0
  dim =  0
  stress = false
  sex = 1#Female
  calving_season = 4#Spring
  age =   Int(ceil(rand(animalModel.rng, truncated(Rayleigh(420),372, 455))))
  lactation= 0
  pregstat = 0 #Pregnant
  dic = 0
  trade_status = 0#false
  neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
  carryover = false 
  fpt = false
  vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
  susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
  clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
end


for b4calves in 1:ceil(optimal_b4*0.5)
  animalModel.id_counter += 1
  id = Int(animalModel.id_counter)
  stage = 1
  range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_calves)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_calves))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
      while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
  push!(animalModel.positions, pos)
  status = 0
  days_infected = 0
  days_exposed = Int(0)
  days_carrier = 0
  days_recovered = Int(0)
  days_treated = Int(0)
  treatment = false
  pop_p = Float64(0.0)
  bacteriaSubmodel = nothing
  pop_p = 0
  pop_r = 0
  dim =  0
  stress = false
  sex = 1#Female
  calving_season = 4#Spring
  age =  Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(55), 7, 90))))
  lactation= 0
  pregstat = 0 #Pregnant
  dic = 0
  trade_status = 0#false
  neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
  carryover = false 
  fpt = false
  vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
  susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5   : rand(animalModel.rng, 0.45:0.01:0.55)
  clinical = false
      pen = 8
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
end

#Distribute the initial infections
num_r = Int(ceil(animalModel.prev_r*length(animalModel.animals)))
num_p = Int(ceil(animalModel.prev_p*length(animalModel.animals)))
num_cr = Int(ceil(animalModel.prev_cr*length(animalModel.animals)))
num_cp = Int(ceil(animalModel.prev_cp*length(animalModel.animals)))

#Set resistant 
uninfected = findall(x-> x.status == 0, animalModel.animals)
resistant = sample(animalModel.rng, uninfected, num_r)
for resist in resistant
  animalModel.animals[resist].status = 4
  animalModel.animals[resist].bacteriaSubmodel = initialiseBacteria()
  animalModel.animals[resist].bacteriaSubmodel.rng = Xoshiro(resist)
  animalModel.animals[resist].bacteriaSubmodel.total_status = 4
  animalModel.animals[resist].days_exposed = 1
  animalModel.animals[resist].clinical = true
  bact_step!(animalModel.animals[resist].bacteriaSubmodel)
end

#Set pathogenic 
uninfected = findall(x-> x.status == 0, animalModel.animals)
pathogenic = sample(animalModel.rng, uninfected, num_p)
for pathogen in pathogenic
  animalModel.animals[pathogen].bacteriaSubmodel = initialiseBacteria( total_status = Int(0), days_treated = Int(0), days_exposed = Int(0), days_recovered = Int(0), stress = false, seed = Int(pathogen))
  animalModel.animals[pathogen].bacteriaSubmodel.rng = Xoshiro(pathogen)
  animalModel.animals[pathogen].status = 3
  animalModel.animals[pathogen].bacteriaSubmodel.total_status = 3
  animalModel.animals[pathogen].days_exposed = 1
  animalModel.animals[pathogen].clinical = true
  bact_step!(animalModel.animals[pathogen].bacteriaSubmodel)
end

#Set carrier resistant
uninfected = findall(x-> x.status == 0, animalModel.animals)
carrier_resistant = sample(animalModel.rng, uninfected, num_cr)
for resist in carrier_resistant
  animalModel.animals[resist].bacteriaSubmodel = initialiseBacteria( total_status = Int(0), days_treated = Int(0), days_exposed = Int(0), days_recovered = Int(0), stress = false, seed = Int(resist))
  animalModel.animals[resist].bacteriaSubmodel.rng = Xoshiro(resist)
  animalModel.animals[resist].status = 6
  animalModel.animals[resist].bacteriaSubmodel.total_status = 6
  animalModel.animals[resist].days_carrier = 1
  bact_step!(animalModel.animals[resist].bacteriaSubmodel)
end

#Set carrier pathogenic 
uninfected = findall(x-> x.status == 0, animalModel.animals)
carrier_pathogenic = sample(animalModel.rng, uninfected, num_cp)
for pathogen in carrier_pathogenic
  animalModel.animals[pathogen].bacteriaSubmodel = initialiseBacteria()
  animalModel.animals[pathogen].bacteriaSubmodel.rng = Xoshiro(pathogen)
  animalModel.animals[pathogen].status = 5
  animalModel.animals[pathogen].bacteriaSubmodel.total_status = 5
  animalModel.animals[pathogen].days_carrier = 1
  bact_step!(animalModel.animals[pathogen].bacteriaSubmodel)
end

for animal in 1:length(animalModel.animals)
  animalModel.animals[animal].bacteriaSubmodel === nothing && continue
  animalModel.animals[animal].bacteriaSubmodel.rng = Xoshiro(animalModel.animals[animal].id)
end



  count_animals!(animalModel)

  optimal_stock = length(animalModel.animals)
  animalModel.timestep = 1


  return animalModel


end

"""
initialiseSplit!(kwargs)
Initialise a split calving system (CS2), as above
"""
function initialiseSplit(;
  farmno::Int = FarmAgent.id,
  system::Int,
  msd::Date,
  seed::Int,
  optimal_stock::Int,
  optimal_lactating::Int,
  treatment_prob::Float64,
  treatment_length::Int,
  carrier_prob::Float64,
  timestep::Int,
  density_lactating::Int,
  density_dry::Int,
  density_calves::Float64,
  date::Date,
  vacc_rate::Float64,
  fpt_rate::Float64,
  prev_r::Float64,
  prev_p::Float64,
  prev_cr::Float64,
  prev_cp::Float64, 
  pen_decon::Float64,
  treat_calf::Float64,
  treat_dry::Float64,
  treat_lac::Float64,
  vacc_protection::Float64,
  vacc_shedding::Float64,
  vacc_duration::Float64, 
  p_retreatment::Float64,
  longevity::Float64
  )

  # Initialise the objects that will be used by the model 
  spare_bact = initialiseBacteria()
  allData = AllData()
  pens = Pens()
  animals = Array{AnimalAgent}[]
  sim = AnimalData()
  environment = Environment()
  transmission = Transmissions()
  infections = Infections()
  cull_curve = modata.cull_curve
  
  #Create the initial model parameters ===============================
  msd_2 = msd - Month(4)
  msd_3 = msd_4 = Date(0)
  current_stock = 0
  current_lactating = 0
  current_dry = 0
  current_heifers = 0
  current_dh = 0
  current_weaned = 0
  current_calves = 0
  optimal_dry = optimal_heifers = optimal_dh = optimal_weaned = optimal_calves = 0
  tradeable_stock = 0
  sending = receiving = Vector{Int}()
  rng = Xoshiro(seed)
  pop_p = 0
  pop_r = 0
  pop_s = 0
  pop_d = 0
  id_counter = 0
  positions = Array{Array{Int}}[]
  processed = false

  N = optimal_stock
  optimal_spring = optimal_autumn = Int(ceil(N*0.5))
  current_spring = Int(0)
  current_autumn = Int(0) 

  current_b1 = current_b2 = current_b3 = current_b4 = 0
  optimal_b1 = optimal_b2 = optimal_b3 = optimal_b4 = 0
  #Set up the model ====================================================
    
    pen_counter = 0
    calf_pen = 8
  
    pen_decon = ifelse(pen_decon < 0.5, false, true)


    animalModel = AnimalModel(farmno, animals, timestep, date, rng, system, msd, msd_2, msd_3, msd_4, seed,  optimal_stock, treatment_prob, treatment_length, carrier_prob, current_stock, current_lactating, optimal_lactating, current_heifers, optimal_heifers, current_calves, optimal_calves, current_weaned, optimal_weaned, current_dh, optimal_dh, current_dry, optimal_dry, tradeable_stock, sending, receiving, density_lactating, density_calves, density_dry, positions, pop_r, pop_s, pop_p, pop_d, id_counter, vacc_rate, fpt_rate, prev_r, prev_p, prev_cr, prev_cp,   current_autumn, optimal_autumn, current_spring, optimal_spring, current_b1, current_b2, current_b3, current_b4, optimal_b1, optimal_b2, optimal_b3, optimal_b4, sim, environment, pen_counter, calf_pen, pen_decon, transmission, infections, treat_calf, treat_dry, treat_lac, vacc_protection, vacc_shedding, vacc_duration, p_retreatment, cull_curve, longevity, allData, pens, spare_bact)
   
  #Set up the model ====================================================

  
  # Set the initial stock parameters
  animalModel.optimal_heifers = animalModel.optimal_weaned = animalModel.optimal_calves = animalModel.optimal_dh = animalModel.optimal_heifers = animalModel.optimal_dry = ceil(0.3*animalModel.optimal_lactating)
  
  
  #Set the animal number generator to 0
  animalModel.id_counter = 0
   # Add the dry cows ---------------------------------------------
  #Dry stage is 6, Dry plane is 6. Model opens on day before psc
   for cow in 1:Int(ceil(animalModel.optimal_spring*0.70))
      animalModel.id_counter += 1
      id = Int(animalModel.id_counter)
      stage = 6
      range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_dry))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
      while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
      push!(animalModel.positions, pos)
      status = 0
      days_infected = 0
      days_exposed = Int(0)
      days_carrier = 0
      days_recovered = Int(0)
      days_treated = Int(0)
      treatment = false
      pop_d = Float64(0.0)
      bacteriaSubmodel = nothing
      dic =  Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(240), 199, 280))))   
      dim = Int(0)
      pop_p = 0
      pop_r = 0
      stress = false
      sex = 1#Female
      calving_season = 1#Split1
      age = Int(ceil(rand(animalModel.rng, truncated(Rayleigh(animalModel.longevity*365*0.75),(3*365), (10*365))))) 
      lactation= round(age/365) - 1
      pregstat = 1#Pregnant
      trade_status = 0#false
      neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
      carryover = false
      fpt = false
      vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
      susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
      clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)

    end

# Add the heifers ---------------------------------------------
#Heifers about to calve, heifer stage 4
  for heifer in 1:ceil(animalModel.optimal_spring*0.3)
      animalModel.id_counter += 1
      id = Int(animalModel.id_counter)
      stage = 4
      range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_heifers)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_heifers))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
       while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
      push!(animalModel.positions, pos)
      status = 0
      days_infected = 0
      days_exposed = Int(0)
      days_carrier = 0
      days_recovered = Int(0)
      days_treated = Int(0)
      treatment = false
      pop_p = Float64(0.0)
      bacteriaSubmodel = nothing
      pop_p = 0
      pop_r = 0
      dic =  Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(240), 199, 280))))   
      dim = 0
      stress = false
      sex = 1#Female
      calving_season = 1#Split1
      age = Int(ceil(rand(animalModel.rng, truncated(Rayleigh(2*365),(22*30), (25*30))))) 
      lactation= 0
      pregstat = 1#Pregnant
      trade_status = 0#false
      neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
      carryover = false
      fpt = false
      vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
      susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
      clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
  end

   #Add weaned animals

  for weaned in 1:ceil(animalModel.optimal_spring*0.25)
      animalModel.id_counter += 1
      id = Int(animalModel.id_counter)
      stage = 2
      range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_weaned)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_weaned))))
      pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
       while pos in animalModel.positions == true
          pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
      end 
      push!(animalModel.positions, pos)
      status = 0
      days_infected = 0
      days_exposed = Int(0)
      days_carrier = 0
      days_recovered = Int(0)
      days_treated = Int(0)
      treatment = false
      pop_p = Float64(0.0)
      bacteriaSubmodel = nothing
      pop_p = 0
      pop_r = 0
      dic =  Int(0)   
      dim = 0
      stress = false
      sex = 1#Female
      calving_season = 1#Spring
      age = Int(ceil(rand(animalModel.rng, truncated(Rayleigh(365),(295), (350))))) 
      lactation= 0
      pregstat = 0#Pregnant
      trade_status = 0#false
      neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
      carryover = false 
      fpt = false
      vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
      susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5   : rand(animalModel.rng, 0.45:0.01:0.55)
      
      clinical = false
      pen = 0
      times_treated = 0
      animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
      push!(animalModel.animals, animal)
  end


#Calving period 2  ----------------

 #Lactating autumn cows
 for cow in 1:ceil(optimal_autumn)
    animalModel.id_counter += 1
    id = Int(animalModel.id_counter)
    stage = Int(5)
    range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_lactating)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_lactating))))
    pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
     while pos in animalModel.positions == true
        pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
    end 
    push!(animalModel.positions, pos)
    status = 0
    days_infected = 0
    days_exposed = Int(0)
    days_carrier = 0
    days_recovered = Int(0)
    days_treated = Int(0)
    treatment = false
    pop_d = Float64(0.0)
    bacteriaSubmodel = nothing
    dic =  ifelse(rand(animalModel.rng) < 0.85, Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(33), 31, 123)))), 0)   
    dim = Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(100), 37, 121))))
    pop_p = 0
    pop_r = 0
    stress = false
    sex = 1#Female
    calving_season = 2#Split2
    age = Int(ceil(rand(animalModel.rng, truncated(Rayleigh(animalModel.longevity*365*0.75),(2*365), (10*365))))) 
    lactation= round(age/365) - 1
    pregstat = ifelse(dic == 0, 0, 1)#Pregnant
    trade_status = 0#false
    neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
    carryover = false
    fpt = false
    vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
    susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
    clinical = false
    pen = 0
    times_treated = 0
    animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
    push!(animalModel.animals, animal)
end

#Split 2 heifers ------------------

for cow in 1:ceil(animalModel.optimal_autumn*0.25)
  animalModel.id_counter += 1
  id = Int(animalModel.id_counter)
  stage = Int(4)
  range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_heifers)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_heifers))))
  pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
  while pos in animalModel.positions == true
      pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
  end 
  push!(animalModel.positions, pos)
  status = 0
  days_infected = 0
  days_exposed = Int(0)
  days_carrier = 0
  days_recovered = Int(0)
  days_treated = Int(0)
  treatment = false
  pop_d = Float64(0.0)
  bacteriaSubmodel = nothing
  dic =  Int(ceil(rand(animalModel.rng, truncated(Rayleigh(42),(1), (60)))))   
  dim = 0
  pop_p = 0
  pop_r = 0
  stress = false
  sex = 1#Female
  calving_season = 2#Split2
  age = Int(ceil(rand(animalModel.rng, truncated(Rayleigh(2*365 - 4*30),(22*30 - 4*30), (25*30 - 4*30)))))
  lactation= round(age/365) - 1
  pregstat = 1#Pregnant
  trade_status = 0#false
  neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
  carryover = false
  fpt = false
  vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
  susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)
  clinical = false
  pen = 0
  times_treated = 0
  animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
  push!(animalModel.animals, animal)
end

# Split2 weaned

for cow in 1:ceil(animalModel.optimal_autumn*0.25)
  animalModel.id_counter += 1
  id = Int(animalModel.id_counter)
  stage = Int(2)
  range = Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_weaned)))) > 250 ? 250 : Int(ceil(√(abs(animalModel.density_dry*animalModel.optimal_weaned))))
  pos =Tuple( [rand(animalModel.rng,  1:range, 2)..., stage])
   while pos in animalModel.positions == true
      pos = Tuple([rand(animalModel.rng, 1:range, 2)..., stage])
  end 
  push!(animalModel.positions, pos)
  status = 0
  days_infected = 0
  days_exposed = Int(0)
  days_carrier = 0
  days_recovered = Int(0)
  days_treated = Int(0)
  treatment = false
  pop_d = Float64(0.0)
  bacteriaSubmodel = nothing
  dic =  0   
  dim = 0
  pop_p = 0
  pop_r = 0
  stress = false
  sex = 1#Female
  calving_season = 2#Split2
  age = Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(100), 37, 121))))
  lactation= round(age/365) - 1
  pregstat = 0#Pregnant
  trade_status = 0#false
  neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
  carryover = false
  fpt = false
  vaccinated = rand(animalModel.rng) < animalModel.vacc_rate ? true : false
  susceptibility = vaccinated == true ?  animalModel.vacc_protection*0.5  : rand(animalModel.rng, 0.45:0.01:0.55)

  clinical = false
  pen = 0
  times_treated = 0
  animal = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
  push!(animalModel.animals, animal)
end
 #Distribute the initial infections
 num_r = Int(ceil(animalModel.prev_r*length(animalModel.animals)))
 num_p = Int(ceil(animalModel.prev_p*length(animalModel.animals)))
 num_cr = Int(ceil(animalModel.prev_cr*length(animalModel.animals)))
 num_cp = Int(ceil(animalModel.prev_cp*length(animalModel.animals)))
 
 #Set resistant 
 uninfected = findall(x-> x.status == 0, animalModel.animals)
 resistant = sample(animalModel.rng, uninfected, num_r)
 for resist in resistant
   animalModel.animals[resist].status = 4
   animalModel.animals[resist].bacteriaSubmodel = initialiseBacteria()
   animalModel.animals[resist].bacteriaSubmodel.rng = Xoshiro(resist)
   animalModel.animals[resist].bacteriaSubmodel.total_status = 4
   animalModel.animals[resist].days_exposed = 1
   animalModel.animals[resist].clinical = true
   bact_step!(animalModel.animals[resist].bacteriaSubmodel)
 end

 #Set pathogenic 
 uninfected = findall(x-> x.status == 0, animalModel.animals)
 pathogenic = sample(animalModel.rng, uninfected, num_p)
 for pathogen in pathogenic
   animalModel.animals[pathogen].bacteriaSubmodel = initialiseBacteria()
   animalModel.animals[pathogen].bacteriaSubmodel.rng = Xoshiro(pathogen)
   animalModel.animals[pathogen].status = 3
   animalModel.animals[pathogen].bacteriaSubmodel.total_status = 3
   animalModel.animals[pathogen].days_exposed = 1
   animalModel.animals[pathogen].clinical = true
   bact_step!(animalModel.animals[pathogen].bacteriaSubmodel)
 end

 #Set carrier resistant
 uninfected = findall(x-> x.status == 0, animalModel.animals)
 carrier_resistant = sample(animalModel.rng, uninfected, num_cr)
 for resist in carrier_resistant
   animalModel.animals[resist].bacteriaSubmodel = initialiseBacteria()
   animalModel.animals[resist].bacteriaSubmodel.rng = Xoshiro(resist)
   animalModel.animals[resist].status = 6
   animalModel.animals[resist].bacteriaSubmodel.total_status = 6
   animalModel.animals[resist].days_carrier = 1
   bact_step!(animalModel.animals[resist].bacteriaSubmodel)
 end

 #Set carrier pathogenic 
 uninfected = findall(x-> x.status == 0, animalModel.animals)
 carrier_pathogenic = sample(animalModel.rng, uninfected, num_cp)
 for pathogen in carrier_pathogenic
   animalModel.animals[pathogen].bacteriaSubmodel = initialiseBacteria( total_status = Int(0), days_treated = Int(0), days_exposed = Int(0), days_recovered = Int(0), stress = false, seed = Int(pathogen))
   animalModel.animals[pathogen].bacteriaSubmodel.rng = Xoshiro(pathogen)
   animalModel.animals[pathogen].status = 5
   animalModel.animals[pathogen].bacteriaSubmodel.total_status = 5
   animalModel.animals[pathogen].days_carrier = 1
   bact_step!(animalModel.animals[pathogen].bacteriaSubmodel)
 end
 
 for animal in 1:length(animalModel.animals)
   animalModel.animals[animal].bacteriaSubmodel === nothing && continue
   animalModel.animals[animal].bacteriaSubmodel.rng = Xoshiro(animalModel.animals[animal].id)
  end
 count_animals!(animalModel)

 optimal_stock = length(animalModel.animals)
 animalModel.timestep = 1

 return animalModel
end



"""
update_animal!(animalModel)
Increment animal parameters
"""

  function update_animal!(animal::AnimalAgent)
    if animal.dim > 0 
        animal.dim += 1 
    end 

    if animal.dic > 0 
        animal.dic += 1
    end

    if animal.days_infected > 0
        animal.days_infected += 1
    end 
    
    if animal.days_treated > 0
        animal.days_treated += 1
    end

    if animal.days_carrier > 0
        animal.days_carrier += 1
    end

    if animal.days_exposed > 0
        animal.days_exposed += 1
    end

    if animal.days_recovered > 0
      animal.days_recovered += 1
    end

    #Advance age
    animal.age += 1

end


"""
run_submodel!(animalModel)
Run the bacterial submodel for each animalModel
"""
function run_submodel!(animal::AnimalAgent)
    if animal.days_recovered ≥ 15
      animal.bacteriaSubmodel = nothing
    end
    animal.bacteriaSubmodel === nothing && return
    animal.status == 0 && return
    bacteriaSubmodel = animal.bacteriaSubmodel
    #Update the submodel parameters
    bacteriaSubmodel.total_status = animal.status
    bacteriaSubmodel.days_treated = animal.days_treated
    bacteriaSubmodel.days_exposed = animal.days_exposed
    bacteriaSubmodel.days_recovered = animal.days_recovered

    if (animal.status in [5,6] && animal.days_exposed == 0) || (animal.status in [7,8] && animal.days_recovered < 15) || (animal.status in [1,2, 3, 4])  || (animal.status in [5,6] && animal.stress == true)
      bact_step!(animal.bacteriaSubmodel)
    end 

        
        if animal.bacteriaSubmodel !== nothing
          animal.pop_r = bacteriaSubmodel.pop_r
          animal.pop_p = bacteriaSubmodel.pop_p
        end
end


"""
animal_mortality!(animalModel. position)
Determine animal mortality if infected
"""
  function animal_mortality!(animalModel::AnimalModel, animal::AnimalAgent)
    animal.status ∉ [1,2] && return
    if animal.stage == 1 
      if animal.clinical == true
        if animal.days_treated == 0
          bernoulli_more(rand(animalModel.rng,  0.10:0.001:0.12), animalModel.rng) && return
          cull!(animal, animalModel)
        else
          bernoulli_more(rand(animalModel.rng,  0.6:0.001:0.72), animalModel.rng) && return
          cull!(animal, animalModel)
        end
      end
    elseif animal.stage != 1
      if animal.clinical == true
        if animal.days_treated == 0 
          bernoulli_more(rand(animalModel.rng,  0.04:0.001:0.05), animalModel.rng) && return
          cull!(animal, animalModel)
        else
          bernoulli_more(rand(animalModel.rng,  0.024:0.001:0.03), animalModel.rng) && return
          cull!(animal, animalModel)
        end
      end
    end
end

"""
animal_processed!(animalModel, position)
Reset the animal processed flag
"""
function animal_processed!(animal::AnimalAgent)
    animal.processed = false
end



"""
animal_recovery!(animal)
Animals recover from infection
"""
function animal_recovery!(animal::AnimalAgent, animalModel::AnimalModel)
  animal.status != 1 && animal.status != 2 && return
    animal.days_infected == 0 && return
    rec_test = animal.vaccinated == false ? Int(round(rand(animalModel.rng,  truncated(Rayleigh(5),(3), (10))))) : (round(rand(animalModel.rng,  truncated(Rayleigh(5),(3), (10))))*animalModel.vacc_duration)
    if animal.days_infected >= rec_test
      if rand(animalModel.rng) > animalModel.carrier_prob
            if animal.status == 1
              animal.days_infected = 0
              animal.days_recovered = 1
              animal.bacteriaSubmodel.days_recovered = 1
              animal.status = animal.bacteriaSubmodel.total_status =  7
              bact_step!(animal.bacteriaSubmodel)
            elseif animal.status == 2
              animal.days_infected = 0
              animal.days_recovered = 1
              animal.bacteriaSubmodel.days_recovered = 1
              animal.status = animal.bacteriaSubmodel.total_status = 8
              bact_step!(animal.bacteriaSubmodel)
            end
      elseif rand(animalModel.rng) <= animalModel.carrier_prob
          if animal.status == 2
            animal.days_infected = 0
            animal.days_carrier = 1
            animal.status = animal.bacteriaSubmodel.total_status = 6
            bact_step!(animal.bacteriaSubmodel)
          elseif animal.status == 1
            animal.days_infected = 0
            animal.days_carrier = 1
            animal.status = animal.bacteriaSubmodel.total_status =  5
            bact_step!(animal.bacteriaSubmodel)
          end
      end
    end
    end


"""
contamination!
Environmental contamination by infected agents.
"""
function contamination!(animal::AnimalAgent, animalModel::AnimalModel)
(animal.status == 0 || animal.status == 3 || animal.status == 4) && return
    existing_contam = animalModel.environment.contamination[animal.pos[3]][animal.pos[1], animal.pos[2]]
    new_contam = ifelse(animal.vaccinated == false, exp(animal.pop_r + animal.pop_p)/2, exp(animal.pop_r + animal.pop_p)/2*animalModel.vacc_shedding)

    if new_contam > existing_contam 
        animalModel.environment.contamination[animal.pos[3]][animal.pos[1], animal.pos[2]] = new_contam
        animalModel.environment.contam_time[animal.pos[3]][animal.pos[1], animal.pos[2]] = 1
        animalModel.environment.contam_type[animal.pos[3]][animal.pos[1], animal.pos[2]] = ifelse((animal.status != 0 && animal.status % 2) == 0, 2, 1)
    end 

end

"""
warm!
Determine warm months (this influences environmental contamination)
"""
function warm(month::Int, animalModel::AnimalModel)
  ((month == 1 || month == 2 || month == 10 || month == 11 || month == 12) && bernoulli_more(0.1, animalModel.rng))
end

"""
chilly!
Determine colder months (this influences environmental contamination)
"""
function chilly(month::Int, animalModel::AnimalModel)
  ((month == 3 || month == 4 || month == 5 || month == 9 || month == 6 || month == 7 || month == 8) && bernoulli_more(0.8, animalModel.rng))
end

"""
environmental_transmission!
Indirect transmission via the environment
"""
function environmental_transmission!(animal::AnimalAgent, animalModel::AnimalModel)
  animal.status != 0 && animal.status != 7 && animal.status != 8 && return
  animalModel.environment.contamination[animal.pos[3]][animal.pos[1], animal.pos[2]] == 0.0 && return
  bernoulli_more(animal.susceptibility, animalModel.rng) && return
  bernoulli_more(animalModel.environment.contamination[animal.pos[3]][animal.pos[1], animal.pos[2]], animalModel.rng) && return
  month = Dates.month(animalModel.date)
  warm(month, animalModel) && return
  chilly(month, animalModel) && return
  animal.bacteriaSubmodel = initialiseBacteria()
  animal.bacteriaSubmodel.rng = Xoshiro(animal.id)
  contam_type = animalModel.environment.contam_type[animal.pos[3]][animal.pos[1], animal.pos[2]]
  
  # Contaminate by host infection type
  if (contam_type != 0 && contam_type % 2 == 0)
    if animal.stage == 1
          bernoulli_less(rand(animalModel.rng, 0.6:0.01:1.0), animalModel.rng) ? animal.clinical = true : animal.clinical = false
    else
          bernoulli_less(rand(animalModel.rng, 0.01:0.001:0.05), animalModel.rng) ? animal.clinical = true : animal.clinical = false
    end
    animal.days_recovered = 0
    animal.status = 4
    animal.days_exposed = 1
    animal.bacteriaSubmodel.clinical = animal.clinical
    animal.bacteriaSubmodel.days_exposed = 1
    animal.bacteriaSubmodel.total_status = 4
    bact_step!(animal.bacteriaSubmodel)
  elseif (contam_type != 0 && contam_type % 2 != 0)
    if animal.stage == 1
      bernoulli_less(rand(animalModel.rng, 0.6:0.01:1.0), animalModel.rng) ? animal.clinical = true : animal.clinical = false
    else
      bernoulli_less(rand(animalModel.rng, 0.01:0.001:0.05), animalModel.rng) ? animal.clinical = true : animal.clinical = false
    end
    animal.days_recovered = 0
    animal.status = 3
    animal.days_exposed = 1
    animal.bacteriaSubmodel.days_exposed = 1
    animal.bacteriaSubmodel.total_status = 3
    bact_step!(animal.bacteriaSubmodel)
  end
end


"""
calfeteria!
Infection at calf feeding.
"""
function calfeteria!(animalModel::AnimalModel, sex::Int)
  for pen in 7:30
        calves = findall(x-> x.stage == 1 && x.sex == sex && x.pen == pen, animalModel.animals)
        length(calves) == 0 && continue
    for calf in shuffle(animalModel.rng, calves)
      transmitter = animalModel.animals[calf]
      (transmitter.status == 3 || transmitter.status == 4 || transmitter.status == 0) && continue
      transmitter.clinical == false && rand(animalModel.rng,Bool) && continue
      transmitter.vaccinated == true && bernoulli_less(animalModel.vacc_shedding, animalModel.rng) && continue
      buddies = unique(sample(animalModel.rng, calves,rand(animalModel.rng, 0:5)))
      while calf in buddies
        buddies = unique(sample(animalModel.rng, calves,rand(animalModel.rng, 0:5)))
      end
      
      for bud in buddies
        atrisk = animalModel.animals[bud]
        atrisk.status != 0 && atrisk.status != 7 && atrisk.status != 8 && continue
        bernoulli_more(atrisk.susceptibility, animalModel.rng) && continue
        (transmitter.status != 0 && transmitter.status % 2 == 0) && bernoulli_more(transmitter.pop_r, animalModel.rng) && continue
        (transmitter.status != 0 && transmitter.status % 2 != 0) && bernoulli_more(transmitter.pop_p, animalModel.rng) && continue
        (transmitter.status != 0 && transmitter.status % 2 == 0) ? atrisk.status = 4 : atrisk.status = 3
        bernoulli_less(rand(animalModel.rng, 0.6:0.01:1.0), animalModel.rng) ? atrisk.clinical = true : atrisk.clinical = false
        atrisk.bacteriaSubmodel = initialiseBacteria()
        atrisk.bacteriaSubmodel.rng = Xoshiro(atrisk.id)
        atrisk.days_exposed = 1
        atrisk.days_recovered = 0
        atrisk.bacteriaSubmodel.total_status = atrisk.status
        atrisk.bacteriaSubmodel.clinical = atrisk.clinical
        bact_step!(atrisk.bacteriaSubmodel)
      end
    end
  end
end



"""
animal_transmission!(animal)
Transmit infection between animals.
Only infected, recovering or carrier animals can transmit to their neighbours
"""
function animal_transmission!(animal::AnimalAgent, animalModel::AnimalModel)

  animal.bacteriaSubmodel === nothing && return
   animal.status == 0 && return
   animal.status == 4 && return
   animal.status == 3 && return
   
   # Short circuit evaluation based on infection type and dominant bacterial population
    (animal.status == 1 && (rand(animalModel.rng) > animalModel.pop_p)) && return
    (animal.status == 2 && (rand(animalModel.rng) > animalModel.pop_r)) && return
    (animal.status == 5 && (rand(animalModel.rng) > animalModel.pop_p)) && return
    (animal.status == 6 && (rand(animalModel.rng) > animalModel.pop_r)) && return
    (animal.status == 7 && (rand(animalModel.rng) > animalModel.pop_p)) && return
    (animal.status == 8 && (rand(animalModel.rng) > animalModel.pop_r)) && return 


    pos = animal.pos
    animal.neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
    animal.processed == true && return

    # Infection based on competition between neighbours
        competing_neighbours = animalModel.animals[findall(x->x.pos in animal.neighbours, animalModel.animals)]

        for competing_neighbour in shuffle(animalModel.rng, competing_neighbours)
        competing_neighbour === nothing && continue
        competing_neighbour.status != 0 && competing_neighbour.status != 7 && competing_neighbour.status != 8 && continue

        if bernoulli_less(competing_neighbour.susceptibility, animalModel.rng) 
              animal.clinical == false && rand(animalModel.rng,Bool) && continue
              animal.vaccinated == true && bernoulli_less(animalModel.vacc_shedding, animalModel.rng) && continue
              if animal.stage == 1
                bernoulli_less(rand(animalModel.rng, 0.3:0.001:0.6), animalModel.rng) ? competing_neighbour.clinical = true : competing_neighbour.clinical = false
              else
                bernoulli_less(rand(animalModel.rng, 0.01:0.001:0.05), animalModel.rng) ? competing_neighbour.clinical = true : competing_neighbour.clinical = false
              end
              competing_neighbour.bacteriaSubmodel = initialiseBacteria()
              competing_neighbour.bacteriaSubmodel.rng = Xoshiro(competing_neighbour.id)
              (animal.status != 0 && animal.status % 2 == 0) ? competing_neighbour.status = 4 : competing_neighbour.status = 3
              (animal.status != 0 && animal.status % 2 == 0) ? competing_neighbour.bacteriaSubmodel.total_status = 4 : competing_neighbour.bacteriaSubmodel.total_status = 3
              competing_neighbour.days_exposed = 1
              competing_neighbour.days_recovered = 0
              competing_neighbour.bacteriaSubmodel.days_exposed = 1
              competing_neighbour.bacteriaSubmodel.clinical = competing_neighbour.clinical
              bact_step!(competing_neighbour.bacteriaSubmodel)
              competing_neighbour.processed = true
            end
          end
    end


"""
record_transmission!
Record transmissions into the Transmissions struct
"""
 function record_transmission!(animalModel, from_id, to_id, stage, from, to, type, effective, clinical)
  push!(animalModel.transmissions.step, animalModel.timestep)
  push!(animalModel.transmissions.from_id, from_id)
  push!(animalModel.transmissions.to_id, to_id)
  push!(animalModel.transmissions.stage, stage)
  push!(animalModel.transmissions.from, from)
  push!(animalModel.transmissions.to, to)
  push!(animalModel.transmissions.type, type)
  push!(animalModel.transmissions.effective, effective)
  push!(animalModel.transmissions.clinical, clinical)
end 

"""
record_infections!
Record infections (individual animal) as they occur
"""
function record_infections!(animalModel::AnimalModel, animal::AnimalAgent, dead::Bool, cull_reason::Symbol)
  push!(animalModel.infections.id, animal.id)
  push!(animalModel.infections.status, animal.status)
  push!(animalModel.infections.step, animalModel.timestep)
  push!(animalModel.infections.stage, animal.stage)
  push!(animalModel.infections.clin, animal.clinical)
  push!(animalModel.infections.death, dead)
  push!(animalModel.infections.days_inf, animal.days_infected)
  push!(animalModel.infections.days_exposed, animal.days_exposed)
  push!(animalModel.infections.vaccinated, animal.vaccinated)
  push!(animalModel.infections.fpt, animal.fpt)
  push!(animalModel.infections.age, animal.age)
  push!(animalModel.infections.cull_reason, cull_reason)

end

"""
animal_shedding!(animal)
Recrudescent infection from carrier animals
"""
function animal_shedding!(animal::AnimalAgent)
  animal.status != 5 && animal.status != 6 && return
    if animal.stress == true
      if animal.status == 5
        animal.bacteriaSubmodel.days_exposed = 1
        animal.bacteriaSubmodel.total_status = 3
        animal.clinical = true
      elseif animal.status == 6
        animal.bacteriaSubmodel.days_exposed = 1
        animal.bacteriaSubmodel.total_status = 4
        animal.clinical = true
      end
    else
      if animal.status == 5
        animal.bacteriaSubmodel.days_exposed = 0
        animal.bacteriaSubmodel.total_status = 5
        animal.bacteriaSubmodel.days_carrier = 1
        animal.clinical = false
      elseif animal.status == 6
        animal.bacteriaSubmodel.days_exposed = 0
        animal.bacteriaSubmodel.total_status = 6
        animal.bacteriaSubmodel.days_carrier = 1
        animal.clinical = false
      end
    end
  end

"""
animal_susceptiblility(animal, animalModel)
Animals return to susceptibility at a variable interval after recovery, simulates waning immunity
"""
  function animal_susceptiblility!(animal::AnimalAgent, animalModel::AnimalModel)
      animal.status != 7 && animal.status != 8 && return
      animal.susceptibility = ℯ^((-animalModel.longevity*365)/animal.days_recovered)/2   

end

"""
check_retreat!
Check if animal has been previously treated, if so, treat according to the simulation's retreatment policy.
"""
check_retreat(animal::AnimalAgent, animalModel::AnimalModel) = animal.times_treated > 0 && bernoulli_more(animalModel.p_retreatment, animalModel.rng) 

"""
treat_sick!
Treat sick animals for Salmonella infections
"""
function treat_sick!(animal::AnimalAgent, animalModel::AnimalModel)
  animal.status != 1 && animal.status != 2 && return
  animal.treatment == true && return
  check_retreat(animal, animalModel) && return  
    if animal.clinical == true 
      bernoulli_more(animalModel.treatment_prob, animalModel.rng) && return
      animal.days_treated = 1
      animal.treatment = true
      animal.bacteriaSubmodel.days_treated = 1
    elseif animal.clinical == false
      animal.stage == 1 && bernoulli_more(animalModel.treat_calf, animalModel.rng) && return
      animal.stage == 5 && bernoulli_more(animalModel.treat_lac, animalModel.rng) && return
      (animal.stage != 5 && animal.stage != 1) && bernoulli_more(animalModel.treat_dry, animalModel.rng) && return
      animal.days_treated = 1
      animal.treatment = true 
      animal.bacteriaSubmodel === nothing && return
      animal.bacteriaSubmodel.days_treated = 1
    end
end

"""
treat_calf!
Incidental treatment of calves.
"""
function treat_calf!(animal::AnimalAgent, animalModel::AnimalModel)
  animal.status != 1 && animal.status != 2 && return
  animal.stage != 1 && return
  animal.times_treated > 0 && bernoulli_more(animalModel.p_retreatment, animalModel.rng) && return
  bernoulli_more(animalModel.treat_calf, animalModel.rng) && return
    animal.days_treated = 1
    animal.treatment = true
    animal.bacteriaSubmodel === nothing && return
    animal.bacteriaSubmodel.days_treated = 1 
end


"""
treat_lac!
Incidental treatment of lactating cattle.
"""
function treat_lac!(animal::AnimalAgent, animalModel::AnimalModel)
  animal.status != 1 && animal.status != 2 && return
  animal.stage != 5 && return
  check_retreat(animal, animalModel) && return  
  bernoulli_more(animalModel.treat_lac, animalModel.rng) && return
    animal.days_treated = 1
    animal.treatment = true
    animal.bacteriaSubmodel === nothing & return
    animal.bacteriaSubmodel.days_treated =1 
  end


"""
animal_treatment!(animal, animalModel)
Decide to treat animals
"""
function treat_drystock!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.status != 1 && animal.status != 2 && return
    check_retreat(animal, animalModel) && return  
    (animal.stage == 5 || animal.stage == 1) && return
    bernoulli_more(animalModel.treat_dry, animalModel.rng) && return
      animal.days_treated = 1
      animal.treatment = true
      animal.bacteriaSubmodel === nothing & return
      animal.bacteriaSubmodel.days_treated =1 
end

"""
animal_fpt_vacc(animal)
Adapt animal susceptibility based on vaccination status
"""
function animal_fpt_vacc!(animal::AnimalAgent, animalModel::AnimalModel)
  animal.stage != 1 && return
  if (animal.fpt == true && animal.age <= 10)
    if animal.age >1 
      animal.susceptibility = (animal.susceptibility*(1-0.05))^animal.age
    end
  end

  if (animal.fpt == false && (animal.age <= rand(animalModel.rng,  14:28)))
    animal.susceptibility = (animal.susceptibility*(1+0.05))^animal.age
  elseif (animal.fpt == true && animal.age > 28)
    animal.vaccinated == true && return
    animal.susceptibility = rand(animalModel.rng, 0.45:0.01:0.55)
  end

  if (animal.age == 54 && (rand(animalModel.rng) < animalModel.vacc_rate))
    animal.vaccinated = true
    animal.susceptibility = animalModel.vacc_protection*0.5
  end

end

"""
end_treatment!(animal, animalModel)
End treatment after course duration.
"""
function end_treatment!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.treatment == false && return
    (animal.days_treated < animalModel.treatment_length) && return
    animal.treatment = false
    animal.days_treated = 0
    animal.bacteriaSubmodel === nothing && return
    animal.bacteriaSubmodel.days_treated = 0
    animal.times_treated += 1
end

"""
getloc(animalModel, pos)
Find the animal's current position in the model
"""
function getloc(animalModel,pos)
  for i in eachindex(animalModel.positions)
    if animalModel.positions[i] == pos
      return i
    end
  end
end

"""move_animal(animal, animalModel)
Shuffle animals at each step with ranges determined by their stocking density.
"""
  function move_animal!(animal::AnimalAgent, animalModel::AnimalModel, stage::Int, density::Union{Int, Float64},  stock_in_class::Int)
    stock_in_class <= 0 ? range = 10 : range = Int(ceil(√(abs(density*stock_in_class))))

    range > 250 ? range = 250 : range = range
    
    oldpos = animal.pos
    newpos = Tuple([rand(animalModel.rng,  1:range, 2)...,stage])

    while newpos in animalModel.positions == true
        newpos = Tuple([rand(animalModel.rng,  1:range, 2)...,stage])
    end
    
    ind = getloc(animalModel, oldpos)
    ind === nothing && return
    deleteat!(animalModel.positions,ind)

    animal.pos = newpos

    push!(animalModel.positions, newpos)
end


"""
move_calf!(animal, animalModel)
Move Calves
"""
function move_calf!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.stage != 1 && return
    if (animal.sex == 1 && animal.age == 1)
      if animalModel.pen_counter > 10
        animalModel.pen_counter = 0
        if animalModel.pen_decon == true 
          animalModel.calf_pen < 25 ? animalModel.calf_pen += 1 : animalModel.calf_pen = 8
        else
          animalModel.calf_pen < 18 ? animalModel.calf_pen += 1 : animalModel.calf_pen = 8
        end
        animal.pen = animalModel.calf_pen
        calves = count(x-> x.stage == 1 & x.sex == 1 & x.pen == animalModel.calf_pen, animalModel.animals)
        move_animal!(animal, animalModel, animalModel.calf_pen, animalModel.density_calves, calves)
      else 
        animal.pen = animalModel.calf_pen
        calves = count(x-> x.stage == 1 & x.sex == 1 & x.pen == animalModel.calf_pen, animalModel.animals)

        move_animal!(animal, animalModel, animalModel.calf_pen, animalModel.density_calves, calves)
      end
    elseif (animal.sex == 0 && animal.age == 1)
      animal.pen = 7
      bobbies = count(x-> x.stage == 1 && x.sex == 0, animalModel.animals)
      move_animal!(animal, animalModel, 7, animalModel.density_calves, bobbies)
    end

    if animal.age > 1
      calves = count(x-> x.stage == 1 & x.sex == 1 & x.pen == animalModel.calf_pen, animalModel.animals)

      move_animal!(animal, animalModel, animal.pen, animalModel.density_calves, calves)
    end
end 

"""
move_weaned!(animal, animalModel)
Move weaned
"""
function move_weaned!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.stage != 2 && return
    move_animal!(animal, animalModel, 2, animalModel.density_dry, animalModel.current_weaned)
end

"""
move_dh!(animal, animalModel)
Move dh
"""
function move_dh!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.stage != 3 && return
    move_animal!(animal, animalModel, 4, animalModel.density_dry, animalModel.current_dh)
end

"""
move_heifer!(animal, animalModel)
Move heifer
"""
function move_heifer!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.stage != 4 && return
    move_animal!(animal, animalModel, 3, animalModel.density_dry, animalModel.current_heifers)
end

"""
move_lactating!(animal, animalModel)
Move lactating
"""
function move_lactating!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.stage != 5 && return
    move_animal!(animal, animalModel, 5, animalModel.density_lactating, animalModel.current_lactating)
end

"""
move_dry!(animal, animalModel)
Move dry
"""
function move_dry!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.stage != 6 && return
    move_animal!(animal, animalModel, 6, Int(animalModel.density_dry), Int(animalModel.current_dry))
end

"""
animal_shuffle!(animal, animalModel)
Randomly move animals.
"""
function animal_shuffle!(animal::AnimalAgent, animalModel::AnimalModel)
        move_calf!(animal, animalModel)
        move_weaned!(animal, animalModel)
        move_dh!(animal, animalModel)
        move_heifer!(animal, animalModel)
        move_lactating!(animal, animalModel)
        move_dry!(animal, animalModel)
end

"""
cull!(animal, animalModel)
Move an animal to level 10, culled
"""
function cull!(animal::AnimalAgent, animalModel::AnimalModel)

    if animal.stage == 1
        animalModel.current_calves -= 1
    elseif animal.stage == 2
        animalModel.current_weaned -=1
    elseif animal.stage == 3
        animalModel.current_heifers -=1 
    elseif animal.stage == 4
        animalModel.current_dh -= 1 
    elseif animal.stage == 6
        animalModel.current_dry -= 1
    end

    if animal.stage == 5
      
      if animalModel.system == 1
        animalModel.current_lactating -= 1
      end
      
      if animalModel.system == 2
        if animal.calving_season == 1
          animalModel.current_spring -= 1
        else
          animalModel.current_autumn -= 1
        end
      end
    
      if animalModel.system == 3
        if animal.calving_season == 1
          animalModel.current_b1 -= 1
        elseif animal.calving_season == 2
          animalModel.current_b2 -= 1
        elseif animal.calving_season == 3
          animalModel.current_b3 -= 1  
        elseif animal.calving_season == 4
          animalModel.current_b4 -= 1
      end
    end
  end
  

    deleteat!(animalModel.animals, findall(x -> x == animal, animalModel.animals))
    deleteat!(animalModel.positions, findall(x -> x == animal.pos, animalModel.positions))
end

"""
cull_empty_dry!(animal, animalModel)
"""
function cull_empty_dry!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.stage != 6 && return
    animal.pregstat != 0 && return
    cull!(animal, animalModel)
    animalModel.current_dry -= 1
end

"""
background_mortality
Background death rates in a herd.
"""
function background_mortality!(animal::AnimalAgent, animalModel::AnimalModel)
  (animal.status == 1 || animal.status == 2) && return
  bernoulli_more(animalModel.cull_curve[animal.age], animalModel.rng) && return
  cull!(animal, animalModel)
end

"""
cull_slipped!(animal, animalModel)
cull animals more than 320 dic that have not calved
"""
function cull_slipped!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.dic < 320 && return
    cull!(animal, animalModel)
end

"""
age_cull!(animal)
Age-based culling
"""
function age_cull!(animal::AnimalAgent, animalModel::AnimalModel)
    bernoulli_more(cdf(truncated(Rayleigh(animalModel.longevity*365), 2.5*365, 10*365), Int(animal.age)), animalModel.rng) && return
    cull!(animal, animalModel)
end

"""
fertility_cull!(animal, animalModel)
cull for fertility
"""
function fertility_cull!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.dim < 150 && return
    animal.dic ≥ 200 && return
    cull!(animal, animalModel)
end

"""
do_culls!(animal, animalModel, system)
Perform both cull types
"""
function do_culls!(animal::AnimalAgent, animalModel::AnimalModel)
    if animal.age > 5000
      cull!(animal, animalModel)
    end
    age_cull!(animal, animalModel)
    fertility_cull!(animal, animalModel)
    cull_empty_carryover!(animal, animalModel)
end

"""
cull_seasonal(animal, animalModel)
Cull for seasonal systems (system = 1)
"""
function cull_seasonal!(animal, animalModel)
    animalModel.timestep % 30 != 0 && return
    animalModel.system != 1 && return
    animal.stage != 5 && return
    (animalModel.current_lactating <= animalModel.optimal_lactating*rand(animalModel.rng,  0.9:0.01:1.1)) && return
    do_culls!(animal, animalModel)
end

"""
milking!
Infection of animals at milking
"""
function milking!(animalModel::AnimalModel)
    milkers = findall(x-> x.stage == 5, animalModel.animals)
    for i in 1:2

      
 # For each milker
 for milker in shuffle(animalModel.rng, milkers)
      transmitter = animalModel.animals[milker]
      transmitter.processed == true && continue   
      (transmitter.status == 3 || transmitter.status == 4 || transmitter.status == 0) && continue#######FLAG######
      transmitter.clinical == false && (rand(animalModel.rng,Bool)) && continue
      transmitter.vaccinated == true && bernoulli_less(animalModel.vacc_shedding, animalModel.rng) && continue

      # Pick a buddy
      buddies = sample(animalModel.rng, milkers,rand(animalModel.rng, 0:2))

      # Transmission between buddies at milking
      for bud in buddies
        bud == milker && continue
        atrisk = animalModel.animals[bud]
        atrisk.status != 0 && atrisk.status != 7 && atrisk.status != 8 && continue
        bernoulli_more(atrisk.susceptibility, animalModel.rng) && continue
        (transmitter.status != 0 && transmitter.status % 2 == 0) && bernoulli_more(transmitter.pop_r, animalModel.rng) && continue
        (transmitter.status != 0 && transmitter.status % 2 != 0) && bernoulli_more(transmitter.pop_p, animalModel.rng) && continue
        bernoulli_less(rand(animalModel.rng, 0.01:0.01:0.05), animalModel.rng) ? atrisk.clinical = true : atrisk.clinical = false
        atrisk.bacteriaSubmodel = initialiseBacteria()
        atrisk.bacteriaSubmodel = initialiseBacteria()
        atrisk.bacteriaSubmodel.rng = Xoshiro(atrisk.id)
        atrisk.days_recovered = 0
        (transmitter.status != 0 && transmitter.status % 2 == 0) ? atrisk.status = 4 : atrisk.status = 3
        atrisk.days_exposed = 1
        atrisk.bacteriaSubmodel.total_status = atrisk.status
        atrisk.bacteriaSubmodel.clinical = atrisk.clinical
        bact_step!(atrisk.bacteriaSubmodel)
      end
    end
  end


end

"""
cull_split!(animal, animalModel)
Cull for split systems (system 1)
"""
function cull_split!(animal, animalModel)
    animalModel.timestep % 30 != 0 && return

    animalModel.system != 2 && return
    animal.stage!= 5 && return
    if animalModel.current_spring > animalModel.optimal_spring
        animal.calving_season != 1 && return
        do_culls!(animal, animalModel)
    end

    if animalModel.current_autumn > animalModel.optimal_autumn
        animal.calving_season != 2 && return
        do_culls!(animal, animalModel)
    end
end

"""
cull_batch!
Apply culling functions in batch calving herds.
"""
function cull_batch!(animal::AnimalAgent, animalModel::AnimalModel)
  animalModel.timestep % 30 != 0 && return

  animalModel.system != 3 && return
  animal.stage!= 5 && return
  animalModel.current_lactating < animalModel.optimal_lactating && return
  
  # Trigger to determine cull season

  trigger = rand(animalModel.rng, 1:4)

  if trigger == 1

    if animalModel.current_b1 > animalModel.optimal_b1
        animal.calving_season != 1 && return
        do_culls!(animal, animalModel)
    end 
  
  elseif trigger == 2 

    if animalModel.current_b2 > animalModel.optimal_b2
        animal.calving_season != 2 && return
        do_culls!(animal, animalModel)
    end
  
  elseif trigger == 3
  
    if animalModel.current_b3 > animalModel.optimal_b3
        animal.calving_season != 3 && return
        do_culls!(animal, animalModel)
    end

  elseif trigger == 4

      if animalModel.current_b4 > animalModel.optimal_b4
          animal.calving_season != 4 && return
          do_culls!(animal, animalModel)
      end
  end    


end
"""
calving!(animal, animalModel)
Calve cows, create calf.
"""
function calving!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.dic != 283  && return
        animal.pregstat = 0
        animal.dic = 0
        animal.stage = 5
        animal.dim = 1
        animal.lactation += 1
        animal.carryover = false
        animal_birth!(animal, animalModel)
        animal.times_treated = 0
        move_animal!(animal, animalModel, 5, animalModel.density_lactating, animalModel.current_lactating)
end

"""
animal_birth!(animal,animalModel)
Create a calf
"""
function animal_birth!(animal::AnimalAgent, animalModel::AnimalModel)
        animalModel.id_counter += 1
        animalModel.pen_counter += 1
        id = Int(animalModel.id_counter)
        stage = 1
        pos = animal.pos
        push!(animalModel.positions, pos)
        status = 0
        days_infected = 0
        days_exposed = 0
        days_carrier = Int(0)
        days_recovered = Int(0)
        days_treated = Int(0)
        treatment = false

        seed = animalModel.seed

          bacteriaSubmodel = nothing
          pop_p = 0
          pop_d = 0
          pop_r = 0
          dic = 0
        dim = 0
        stress = false
        sex = (rand(animalModel.rng) > 0.5) ? 1 : 0
        calving_season = animal.calving_season
        age = Int(1)
        lactation = 0
        pregstat = 0
        trade_status = 0#false
        neighbours = modata.possible_neighbours[CartesianIndex(Tuple(pos))]
        processed = true
        carryover = false
        fpt = (rand(animalModel.rng) < animalModel.fpt_rate) ? true : false
        vaccinated = false
        if fpt == true
          susceptibility = rand(animalModel.rng, 0.9:0.01:1.0)
        elseif animal.vaccinated == true
          susceptibility = animalModel.vacc_protection*0.5  
        else 
          susceptibility = rand(animalModel.rng, 0.65:0.01:0.75)
        end
        clinical = false
        pen = 0        
        times_treated = 0
        calf = AnimalAgent(id, pos, status, stage, days_infected, days_exposed, days_carrier, days_recovered, days_treated, treatment, pop_p, pop_d, pop_r, bacteriaSubmodel, dic, dim, stress, sex, calving_season, age, lactation, pregstat, trade_status, neighbours, processed, carryover, fpt, vaccinated, susceptibility, clinical, pen,times_treated)    
        environmental_transmission!(calf, animalModel)
        move_calf!(calf, animalModel)
        push!(animalModel.animals, calf)
end


"""
bobby_cull!(animal, animalModel)
Cull bobby calves
"""
function bobby_cull!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.stage != 1 && return
    animal.sex != 0 && return
    animal.age < 4 && return
    cull!(animal, animalModel)
    animalModel.current_calves -= 1
end


"""
cull_empty_carryover!
Cull carryover cows that have not conceived.
"""
function cull_empty_carryover!(animal::AnimalAgent, animalModel::AnimalModel)
  animalModel.timestep % 30 != 0 && return

  if animal.carryover == true & animal.pregstat == 0
    animal.dim > 150 && return
    cull!(animal, animalModel)
  end
end 



"""
join_seasonal!(animal, animalModel)
Join animals in seasonal systems (CS1)
"""
function join_seasonal!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.stage != 5 && return
    animal.pregstat != 0 && return
    if (animalModel.date == (animalModel.msd + Month(3)))
        if rand(animalModel.rng) <= 0.85
            animal.pregstat = 1
            animal.dic = Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(63), 1, 84))))
        end

    end
end

"""
join_split!(animal, animalModel)
Join animals in split calving herds (CS2)
"""
function join_split!(animal::AnimalAgent, animalModel::AnimalModel)
  animal.pregstat != 0 && return
  animal.stage != 5 && return
  (rand(animalModel.rng) ≤ 0.15) && return

  if (animal.calving_season == 1 && (animalModel.date == (animalModel.msd + Month(3))))
    animal.dic = Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(63), 1, 84)))) 
    animal.pregstat = 1
  elseif (animal.calving_season == 2 && (animalModel.date == (animalModel.msd_2 + Month(3))))
    animal.dic =  Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(63), 1, 84))))
    animal.pregstat = 1
  end
end 

"""
join_batch!(animal, animalModel)
Join animals in batch calving systems (CS3)
"""
function join_batch!(animal::AnimalAgent, animalModel::AnimalModel)

  animal.pregstat != 0 && return
  animal.stage != 5 && return
  (rand(animalModel.rng)) ≤ 0.15 && return

  if (animal.calving_season == 1 && (animalModel.date == animalModel.msd + Month(3)))
    animal.dic = Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(63), 1, 84))))
    animal.pregstat = 1
  elseif (animal.calving_season == 2 && (animalModel.date == animalModel.msd_2 - Year(1) + Month(3)))
    animal.dic = Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(63), 1, 84))))
    animal.pregstat = 1
  elseif (animal.calving_season == 3 && (animalModel.date == animalModel.msd_3 + Month(3)))
    animal.dic = Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(63), 1, 84))))
    animal.pregstat = 1
  elseif (animal.calving_season == 4 && (animalModel.date == animalModel.msd_4 + Month(3)))
    animal.dic = Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(63), 1, 84))))
    animal.pregstat = 1
  end


end

"""
animal_joining!(animal, animalModel)
Wrapper function for assigning pregnancy status.
"""
function animal_joining!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.pregstat != 0 && return
    animal.stage != 5 && return
    if animalModel.system == 1
        join_seasonal!(animal, animalModel)
    elseif animalModel.system == 2
        join_split!(animal, animalModel)
    elseif animalModel.system == 3
        join_batch!(animal, animalModel)
    end

end

"""
animal_status!(animal)
Update the status of each animal depending on its bacterial population.
"""
function animal_status!(animal::AnimalAgent)
    animal.status != 3 && animal.status != 4 && return
      if animal.pop_r ≥ 0.5
          animal.status = 2 
          animal.days_infected = 1
          animal.days_exposed = 0
      elseif animal.pop_p ≥ 0.5
          animal.status = 1
          animal.days_infected = 1
          animal.days_exposed = 0
      end
end

"""
animal_wean!(animal, animalModel)
Wean calves to next lifestage
"""
function animal_wean!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.stage != 1 && return
    (animal.age <= rand(animalModel.rng,  56:90)) && return
    animal.stage = 2
    animal.times_treated = 0
    move_animal!(animal, animalModel, 2, animalModel.density_dry, animalModel.current_weaned)

end

"""
animal_heifer!(animal, animalModel)
Transition to the heifer lifestage
"""
function animal_heifer!(animal::AnimalAgent, animalModel::AnimalModel)  
    animal.stage != 2 && return
    (animal.age <= 13*30*rand(animalModel.rng,  0.9:0.01:1.1))  && return
    animal.stage = 3
    animal.times_treated = 0
    move_animal!(animal, animalModel, 3, animalModel.density_dry, animalModel.current_heifers)
end

"""
heifer_pregnancy!(animal, animalModel)
Set heifer pregnancy status
"""
function heifer_pregnancy!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.pregstat = 1
    animal.stage = 4
    animal.dic = Int(ceil(rand(animalModel.rng,  truncated(Rayleigh(42), 1, 63))))
    animal.times_treated = 0
    move_animal!(animal, animalModel, 4, animalModel.density_dry, animalModel.current_dh)
end

"""
join_heifer_seasonal!(animal, animalModel)
Join heifers for seasonal systems 
"""
function join_heifer_seasonal!(animal::AnimalAgent, animalModel::AnimalModel)
    animalModel.system != 1 && return
    if animalModel.date == (animalModel.msd + Day(42))
        heifer_pregnancy!(animal, animalModel)
    end
end

"""
join_heifer_split!(animal, animalModel)
Join heifers in split systems
"""
function join_heifer_split!(animal::AnimalAgent, animalModel::AnimalModel)
    animalModel.system != 2 && return
    if animal.calving_season == 1 && (animalModel.date == (animalModel.msd + Day(42)))
        heifer_pregnancy!(animal, animalModel)
    elseif animal.calving_season == 2 && (animalModel.date == (animalModel.msd_2 + Day(42)))
        heifer_pregnancy!(animal, animalModel)
    end
end

"""
join_heifer_batch!(animal, animalModel)
Join heifers in batch systems
"""
function join_heifer_batch!(animal::AnimalAgent, animalModel::AnimalModel)
    animalModel.system != 3 && return
    if animal.calving_season == 1 && animalModel.date == (animalModel.msd + Day(42))
        heifer_pregnancy!(animal, animalModel)
    elseif animal.calving_season == 2 && animalModel.date == (animalModel.msd_2  - Year(1) + Day(42))
        heifer_pregnancy!(animal, animalModel)
    elseif animal.calving_season == 3 && animalModel.date == (animalModel.msd_3 + Day(42))
        heifer_pregnancy!(animal, animalModel)
    elseif animal.calving_season == 4 && animalModel.date == (animalModel.msd_4 + Day(42))
        heifer_pregnancy!(animal, animalModel)
    end
end

"""
join_heifers!(animal, animalModel)
Join heifers for all systems
"""
function join_heifers!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.pregstat == 1 && return
    animal.stage != 3 && return
    join_heifer_seasonal!(animal, animalModel)
    join_heifer_split!(animal, animalModel)
    join_heifer_batch!(animal, animalModel)
end

"""
cull_empty_heifers!(animal, animalModel)
Cull empty heifers
"""
function cull_empty_heifers!(animal, animalModel)
    animalModel.timestep % 30 != 0 && return
    animal.stage != 3 && return
    animal.pregstat == 1 && return
    (animal.age <= rand(animalModel.rng,  550:650)) && return
    cull!(animal, animalModel)
end

"""
cull_surplus_heifers_spring!
Remove heifers in excess of replacement rate based on pregnancy status (i.e. keep animals that are more heavily in calf.)
"""
function cull_surplus_heifers_spring!(animal::AnimalAgent, animalModel::AnimalModel)
  animalModel.timestep % 30 != 0 && return
  animalModel.system != 1 && return
  animal.stage != 4 && return
  animal.pregstat == 0 && return
  (animalModel.date != animalModel.msd + Day(84)) && return
  animal.dic > 63 && return
  cull!(animal, animalModel)
end

"""
cull_surplus_heifers_split!
Remove heifers in excess of replacement rate based on pregnancy status (i.e. keep animals that are more heavily in calf.)
"""
function cull_surplus_heifers_split!(animal::AnimalAgent, animalModel)
  animalModel.timestep % 30 != 0 && return
  animalModel.system != 2 && return
  animal.stage != 4 && return
  animal.pregstat == 0 && return
  if (animal.calving_season == 1 && (animalModel.date == (animalModel.msd + Day(84))))
    animal.dic > 63 && return
    cull!(animal, animalModel)
  elseif (animal.calving_season == 2 && (animalModel.date == (animalModel.msd_2 + Day(84))))
    animal.dic > 63 && return
    cull!(animal, animalModel)
  end
end

"""
cull_surplus_heifers_batch!
Remove heifers in excess of replacement rate based on pregnancy status (i.e. keep animals that are more heavily in calf.)
"""
function cull_surplus_heifers_batch!(animal::AnimalAgent, animalModel::AnimalModel)
  animalModel.timestep % 30 != 0 && return
  animalModel.system != 3 && return
  animal.stage != 4 && return
  animal.pregstat == 0 && return
  if (animal.calving_season == 1 && animalModel.date == (animalModel.msd + Day(84)))
    animal.dic > 63 && return
    cull!(animal, animalModel)
  end

if (animal.calving_season == 2 && animalModel.date == (animalModel.msd_2  - Year(1) + Day(84)))
  animal.dic > 63 && return
  cull!(animal, animalModel)
end

if (animal.calving_season == 3 && animalModel.date == (animalModel.msd_3 + Day(84)))
  animal.dic > 63 && return
  cull!(animal, animalModel)
end

if (animal.calving_season == 4 && animalModel.date == (animalModel.msd_4 + Day(84)))
  animal.dic > 63 && return
  cull!(animal, animalModel)
end

end

"""
set_dry!(animal, animalModel)
Set an animal's status to dry
"""
function set_dry!(animal::AnimalAgent, animalModel::AnimalModel)
    
    if rand(animalModel.rng) < animalModel.vacc_rate
      animal.vaccinated = true
      animal.susceptibility = animalModel.vacc_protection*0.5
    end
  
    animal.stage = 6
    animal.dim = 0
    move_animal!(animal, animalModel, 6, animalModel.density_dry, animalModel.current_dry)
end


"""
dryoff_seasonal!(animal, animalModel)
Dry off lactating cows in seasonal systems
"""
function dryoff_seasonal!(animal::AnimalAgent, animalModel::AnimalModel)
    animalModel.system != 1 && return
    if (animal.pregstat == 0 && animal.dim < 330) && (rand(animalModel.rng) < 0.4 && animal.carryover == false)
      animal.carryover = true
  else
      set_dry!(animal, animalModel)
  end
end

"""
dryoff_split!(animal, animalModel)
Dry off lactating cows in split systems
"""
function dryoff_split!(animal::AnimalAgent, animalModel::AnimalModel)
    animalModel.system != 2 && return
    if (animal.pregstat == 0 && animal.dim < 330) && (rand(animalModel.rng) < 0.4 && animal.carryover == false)
        animal.calving_season == 1 ? animal.calving_season = 2 : animal.calving_season = 1
        animal.carryover = true
    else
        set_dry!(animal, animalModel)
    end
end

"""
dryoff_batch!(animal, animalModel)
Dryoff animals in batch systems
"""
function dryoff_batch!(animal, animalModel)
    animalModel.system != 3 && return

    if (animal.pregstat == 0 && animal.dim < 330) && (rand(animalModel.rng) < 0.4 && animal.carryover == false)
        animal.carryover = true
        
        if animal.calving_season < 4
            animal.calving_season += 1
        else
            animal.calving_season = 1
        end 

    else
        set_dry!(animal, animalModel)
    end
end

"""
export_alldata!
Export animal level data into the allData constructor if required.
"""
function export_alldata!(animal::AnimalAgent, animalModel::AnimalModel)
  push!(animalModel.allData.id, animal.id)
  push!(animalModel.allData.step, animalModel.timestep)
  push!(animalModel.allData.stage, animal.stage)
  push!(animalModel.allData.status, animal.status)

  push!(animalModel.allData.pregstat, animal.pregstat)
  push!(animalModel.allData.dic, animal.dic)
  push!(animalModel.allData.dim, animal.dim)
  push!(animalModel.allData.age, animal.age)
  push!(animalModel.allData.days_exposed, animal.days_exposed)
  push!(animalModel.allData.days_infected, animal.days_infected)
  push!(animalModel.allData.days_recovered, animal.days_recovered)
  push!(animalModel.allData.days_carrier, animal.days_carrier)
  push!(animalModel.allData.days_treated, animal.days_treated)

  push!(animalModel.allData.treatment, animal.treatment)
  push!(animalModel.allData.pop_p, animal.pop_p)
  push!(animalModel.allData.pop_r, animal.pop_r)
  push!(animalModel.allData.susceptibility, animal.susceptibility)
end


"""
animal_dryoff!(animal, animalModel)
Dryoff function for all calving systems
"""
function animal_dryoff!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.stage  != 5 && return
    animal.dim % 305 != 0 && return
    animal.times_treated = 0
    dryoff_seasonal!(animal, animalModel)
    dryoff_batch!(animal, animalModel)
    dryoff_split!(animal, animalModel)
end

"""
flag_trades!(animal, animalModel)
Flag animals eligible for trading
"""
function flag_trades!(animal::AnimalAgent, animalModel::AnimalModel)
    animal.stage != 3 && animal.stage != 4 && animal.stage != 5 && animal.stage != 6 && return
    if (animal.stage == 3 || animal.stage == 4) && (animal.age >= 13*30 && animal.age <= 18*30)
        animal.trade_status = true
    elseif (animal.stage == 4 || animal.stage ==  6) && animal.dic <= 241
        animal.trade_status = true
    elseif animal.stage == 5 && animal.dim >= 120
        animal.trade_status = true
    else
        animal.trade_status = false
    end
end

"""
update_msd!(animalModel)
Update the msd for each year
"""
function update_msd!(animalModel::AnimalModel)
    if Year(animalModel.date) > Year(animalModel.msd)
        animalModel.msd += Year(1)
    end
      
    if (animalModel.system == 2 || animalModel.system == 3)
        if Year(animalModel.date) > Year(animalModel.msd_2)
            animalModel.msd_2 += Year(1)
        end
    end
    
  if animalModel.system == 3
        if Year(animalModel.date) > Year(animalModel.msd_3)
            animalModel.msd_3 += Year(1)
        elseif Year(animalModel.date) > Year(animalModel.msd_4)
              animalModel.msd_4 += Year(1)
        end
  end
end

"""
animal_mstep!(animal, animalModel)
Update some parameters once per day
"""
function animal_mstep!(animalModel::AnimalModel)
    animalModel.timestep += 1
    animalModel.date += Day(1)
    update_msd!(animalModel)
    count_animals!(animalModel)
    animalData = animalModel.sim
    animal_export!(animalModel,animalData)
  end

"""
animal_stress!(animal,animalModel)
Apply a shedding-inducing stress event depending on age and stage.
"""
function animal_stress!(animal::AnimalAgent, animalModel::AnimalModel)
    if animal.dic >= 223 || animal.dim < 21 || (animal.age <= 2.5*365 && animal.stage == 5)
      animal.stress = true
      animal.status in [7,8] && return
      animal.susceptibility = rand(animalModel.rng,  0.65:0.01:0.75)
    else
      animal.stress = false
      animal.susceptibility = animal.vaccinated == true ? animalModel.vacc_protection*0.5  : animal.susceptibility
    end
end


"""
traing_need!
Determine if the herd needs to do a little cattle trading at time t.
"""
function trading_need!(animalModel::AnimalModel)
  animalModel.sending = Vector{Int}()
    # Tradeable stock
    tradeable_heifers = shuffle(findall(x-> x.trade_status == true && x.stage == 3, animalModel.animals))
    tradeable_dh = shuffle(findall(x-> x.trade_status == true && x.stage == 4, animalModel.animals))
    tradeable_lactating = shuffle(findall(x-> x.trade_status == true && x.stage == 5, animalModel.animals))
    tradeable_dry = shuffle(findall(x-> x.trade_status == true && x.stage == 6, animalModel.animals))
    
    # The truck
    truckers = Vector{Int}()

    # Draft them out for the truck

    if animalModel.current_heifers > animalModel.optimal_heifers
      length(tradeable_heifers) == 0 && return
      maxsample = (animalModel.current_heifers - animalModel.optimal_heifers) > length(tradeable_heifers) ? length(tradeable_heifers) : (animalModel.current_heifers - animalModel.optimal_heifers)
      for i in 1:maxsample
        push!(truckers, tradeable_heifers[i])
      end

    end

    if animalModel.current_dh > animalModel.optimal_dh
      length(tradeable_dh) == 0 && return
      maxsample = (animalModel.current_dh - animalModel.optimal_dh) > length(tradeable_dh) ? length(tradeable_dh) : (animalModel.current_dh - animalModel.optimal_dh)
      for i in 1:maxsample
        push!(truckers, tradeable_dh[i])
      end
    end

    if animalModel.current_lactating > animalModel.optimal_lactating
      length(tradeable_lactating) == 0 && return
      maxsample = (animalModel.current_lactating - animalModel.optimal_lactating) > length(tradeable_lactating) ? length(tradeable_lactating) : (animalModel.current_lactating - animalModel.optimal_lactating)
      
      for i in 1:maxsample
        push!(truckers, tradeable_lactating[i])
      end
    end

    if animalModel.current_dry > animalModel.optimal_dry
      length(tradeable_dry) == 0 && return
      maxsample = (animalModel.current_dry - animalModel.optimal_dry) > length(tradeable_dry) ? length(tradeable_dry) : (animalModel.current_dry - animalModel.optimal_dry)
      for i in 1:maxsample
        push!(truckers, tradeable_dry[i])
      end
    end

    min = round(0.0025*animalModel.current_stock)
    most_likely = round(0.025*animalModel.current_stock)
    max = round(0.04*animalModel.current_stock)

    # Now we're trucking
    animalModel.sending = truckers
end

"""
pendata!
Export data on the animals in each pen if required
"""
function pendata!(animalModel::AnimalModel)
  for pen in 7:33
    num = length(findall(x-> x.stage == 1 && x.pen == pen, animalModel.animals))
    push!(animalModel.pens.step, animalModel.timestep)
    push!(animalModel.pens.pen, pen)
    push!(animalModel.pens.num, num)
  end
end


"""
step_contam!
Step contamination on each square at each timestep.
"""
function step_contam!(animalModel::AnimalModel)
  #Threads.@threads 
  for i in 1:length(animalModel.environment.contamination)
        nnz(animalModel.environment.contamination[i]) == 0 && continue
        @inbounds    animalModel.environment.contamination[i]    .= ifelse.(animalModel.environment.contamination[i] .!=0, animalModel.environment.contamination[i].*(1 ./ animalModel.environment.contam_time[i]), 0)
        @inbounds    animalModel.environment.contam_time[i]      .= ifelse.(animalModel.environment.contam_time[i] .>1, animalModel.environment.contam_time[i] .+=1, 0)
   end
end

"""
animal_step!
Container function for stepping the animal's in the model at time t
"""
function animal_step!(animalModel::AnimalModel)

    # Run the submodel first. You can thread this if you want to run a singe sim.
   #Threads.@threads 
   for animal in findall(x->x.bacteriaSubmodel !== nothing, animalModel.animals)
        @inbounds run_submodel!(animalModel.animals[animal])
    end 

    # Run the functions on the shuffled critters.
     for x in eachindex(shuffle(animalModel.rng, animalModel.animals))
      checkbounds(Bool, animalModel.animals, x) == false && continue   
      @inbounds animal = animalModel.animals[x]
        animal.processed = false
        
        # Trading function
        flag_trades!(animal, animalModel) 
        #  vaccination and transmission
        animal_fpt_vacc!(animal, animalModel) 
        environmental_transmission!(animal, animalModel)  
        animal_transmission!(animal, animalModel)  #
        animal_recovery!(animal, animalModel) 
        #Treatment
        treat_sick!(animal, animalModel) 
        treat_lac!(animal, animalModel)   
        treat_calf!(animal, animalModel)
        treat_drystock!(animal, animalModel)
        # Shedding, susceptibility
        animal_stress!(animal, animalModel) #
        animal_shedding!(animal)  
        animal_susceptiblility!(animal, animalModel)  #
        end_treatment!(animal, animalModel) 
        contamination!(animal, animalModel) 

        # Calf housekeeping functions
         bobby_cull!(animal, animalModel)
         animal_wean!(animal, animalModel)
        # Reproductive functions
          animal_joining!(animal, animalModel)
          join_heifers!(animal, animalModel)
          animal_heifer!(animal, animalModel)
          calving!(animal, animalModel)  
          background_mortality!(animal, animalModel)
          animal_mortality!(animalModel, animal)     
          animal_dryoff!(animal, animalModel)
      
        # Culling
          cull_empty_dry!(animal, animalModel)
          cull_seasonal!(animal, animalModel) 
          cull_split!(animal, animalModel)
          cull_batch!(animal, animalModel)
          cull_slipped!(animal, animalModel)
          cull_empty_heifers!(animal, animalModel)
          cull_surplus_heifers_spring!(animal, animalModel)
          cull_surplus_heifers_split!(animal, animalModel)
          cull_surplus_heifers_batch!(animal, animalModel)
          # Housekeeping functions
          animal_shuffle!(animal, animalModel)
          animal_status!(animal)
          update_animal!(animal)
        #  export_alldata!(animal, animalModel)
    end 

    # Progress contamination
    step_contam!(animalModel) 
 
  # Run milking and calf infections
  milking!(animalModel) 
  calfeteria!(animalModel,1)
  calfeteria!(animalModel,0)
  # Update time-sensitive steps
  animal_mstep!(animalModel)
  # Determine the trading need of each farm 
  trading_need!(animalModel)

end


"""
animal_export!(animalModel, animalData)data from the model at each timestep.
Export the data at each model tep.
"""
function animal_export!(animalModel::AnimalModel,animalData::AnimalData)

    
    alive = animalModel.animals

    lifestages = countmap([animal.stage for animal in alive])
    statuses = countmap([animal.status for animal in alive])
    
    data = animalData
    step = animalModel.timestep
    data.id[step] = animalModel.farmno
    data.timestep[step] = step
    data.pop_r[step] = animalModel.pop_r
    data.pop_s[step] = animalModel.pop_s
    data.pop_d[step] = animalModel.pop_d
    data.pop_p[step] = get(statuses, 1, 0)

    # This could be functionalised/vectorised, but it is (even) more confusing

    lifestages = countmap([animal.stage for animal in alive])

    data.num_calves[step] = get(lifestages, 1, 0)
    data.num_weaned[step] = get(lifestages, 2, 0)
    data.num_heifers[step] = get(lifestages, 3, 0)
    data.num_dh[step] = get(lifestages, 4, 0)
    data.num_lactating[step] = get(lifestages, 5, 0)
    data.num_dry[step] = get(lifestages, 6, 0)

    data.pop_rec_r[step] = get(statuses, 8, 0)
    data.pop_rec_p[step] = get(statuses, 7, 0)
    data.pop_car_r[step] = get(statuses, 6, 0)
    data.pop_car_p[step] = get(statuses, 5, 0)

    data.pop_er[step] = get(statuses, 4, 0)
    data.pop_ep[step] = get(statuses, 3, 0)

    status_stage = countmap([(i.status, i.stage) for i in alive])

    data.inf_calves[step] = get(status_stage, (1,1), 0)
    data.ri_calves[step] = get(status_stage, (2,1), 0)

    data.car_calves[step] =  get(status_stage, (5,1), 0)
    data.rc_calves[step] = get(status_stage, (6,1), 0)

    data.inf_weaned[step] = get(status_stage, (1,2), 0)
    data.ri_weaned[step] = get(status_stage, (2,2), 0)
    data.car_weaned[step] = get(status_stage, (5,2), 0)
    data.rc_weaned[step] = get(status_stage, (6,2), 0)

    data.inf_heifers[step] = get(status_stage, (1,3), 0)
    data.ri_heifers[step] = get(status_stage, (2,3), 0)
    data.car_heifers[step] = get(status_stage, (5,3), 0)
    data.rc_heifers[step] = get(status_stage, (6,3), 0)

    data.inf_dh[step] = get(status_stage, (1,4), 0)
    data.ri_dh[step] = get(status_stage, (2,4), 0)
    data.car_dh[step] = get(status_stage, (5,4), 0)
    data.rc_dh[step] = get(status_stage, (6,4), 0)

    data.inf_lac[step] = get(status_stage, (1,5), 0)
    data.ri_lac[step] = get(status_stage, (2,5), 0)
    data.car_lac[step] = get(status_stage, (5,5), 0)
    data.rc_lac[step] = get(status_stage, (6,5), 0)

    data.inf_dry[step] = get(status_stage, (1,6), 0)
    data.ri_dry[step] = get(status_stage, (2,6), 0)
    data.car_dry[step] = get(status_stage, (5,6), 0)
    data.rc_dry[step] = get(status_stage, (6,6), 0)


    data.clinical[step] = count(i->((i.status == 1 || i.status == 2) && i.clinical == true), alive)
    data.subclinical[step] = count(i->((i.status == 1 || i.status == 2) && i.clinical != true), alive)

    current_batch = countmap([(i.stage, i.calving_season) for i in alive])
    data.current_b1[step] = get(current_batch, (5,1), 0)
    data.current_b2[step] = get(current_batch, (5,2), 0)
    data.current_b3[step] = get(current_batch, (5,3), 0)
    data.current_b4[step] = get(current_batch, (5,4), 0)
  end

  """
Run a single farm with a selected calving system, with n cows

"""



function run_herd!(; optimal_stock, calving_system, pen_decon, treatment_prob, vacc_rate, fpt_rate, treat_calf, treat_lac, simdays, init_r, init_p, init_cr, init_cp, treatment_length, density_lactating, density_dry, density_calves,
  treat_dry,
  vacc_protection,
  vacc_shedding,
  vacc_duration,
  p_retreatment,
  longevity,
  run
  )
    
  
        prev_r = Float64(init_r)
        prev_p = Float64(init_p)
        prev_cr = Float64(init_cr)
        prev_cp = Float64(init_cp)
   

    # Create empty vectors for trades
    trades_from = []
    trades_to = []

    # Set the farm's attributes based on the parameters for that ID in svars (the farm file)
    ncows = optimal_stock
    system = calving_system
    pen_decon = pen_decon
    optimal_stock = Int(optimal_stock)
    treatment_prob = Float64(treatment_prob)
    vacc_rate = Float64(vacc_rate)
    fpt_rate = Float64(fpt_rate)
    treat_calf = treat_calf
    treat_lac = treat_lac
    saleyards = 500
    treat_dry = treat_dry
    vacc_protection = vacc_protection
    vacc_shedding = vacc_shedding
    vacc_duration = vacc_duration
    p_retreatment = p_retreatment
    longevity = longevity

    # Choose an initialisation function based on whether the farm is CS1 (spring), CS2 (split) or CS3 (batch)
    if calving_system == 1
        herd = initialiseSpring(
                farmno = Int(1),
                system = Int(1),
                msd = Date(2021,9,24),
                seed = Int(42),
                optimal_stock = optimal_stock,
                optimal_lactating = optimal_stock,
                treatment_prob = treatment_prob,
                treatment_length = Int(treatment_length),
                carrier_prob = Float64(0.1),
                timestep = Int(0),
                density_lactating = Int(density_lactating),
                density_dry = Int(density_dry),
                density_calves = Float64(density_calves),
                date = Date(2021,7,2),
                vacc_rate = vacc_rate,
                fpt_rate = fpt_rate,
                prev_r = Float64(prev_r),
                prev_p = Float64(prev_p),
                prev_cr = Float64(prev_cr),
                prev_cp = Float64(prev_cp),
                pen_decon = Float64(pen_decon),
                treat_calf = treat_calf/70,
                treat_dry = treat_dry/365,
                treat_lac = treat_lac/365,
                vacc_protection = vacc_protection,
                vacc_shedding = vacc_shedding,
                vacc_duration = vacc_duration,
                p_retreatment = p_retreatment,
                longevity = Float64(longevity)
            )
    elseif calving_system == 2
        herd = initialiseSplit(
            farmno = Int(1),
            system = Int(2),
            msd = Date(2021,9,24),
            seed = Int(42),
            optimal_stock = optimal_stock,
                optimal_lactating = optimal_stock,
                treatment_prob = treatment_prob,
                treatment_length = Int(treatment_length),
                carrier_prob = Float64(0.1),
                timestep = Int(0),
                density_lactating = Int(density_lactating),
                density_dry = Int(density_dry),
                density_calves = Float64(density_calves),
                date = Date(2021,7,2),
                vacc_rate = vacc_rate,
                fpt_rate = fpt_rate,
                prev_r = Float64(prev_r),
                prev_p = Float64(prev_p),
                prev_cr = Float64(prev_cr),
                prev_cp = Float64(prev_cp),
                pen_decon = Float64(pen_decon),
                treat_calf = treat_calf/70,
                treat_dry = treat_dry/365,
                treat_lac = treat_lac/365,
                vacc_protection = vacc_protection,
                vacc_shedding = vacc_shedding,
                vacc_duration = vacc_duration,
                p_retreatment = p_retreatment,
                longevity = Float64(longevity)
          )
    else 
        herd = initialiseBatch(
            farmno = Int(1),
            system = Int(3),
            msd = Date(2021,9,24),
            seed = Int(42),
            optimal_stock = optimal_stock,
                optimal_lactating = optimal_stock,
                treatment_prob = treatment_prob,
                treatment_length = Int(treatment_length),
                carrier_prob = Float64(0.1),
                timestep = Int(0),
                density_lactating = Int(density_lactating),
                density_dry = Int(density_dry),
                density_calves = Float64(density_calves),
                date = Date(2021,7,2),
                vacc_rate = vacc_rate,
                fpt_rate = fpt_rate,
                prev_r = Float64(prev_r),
                prev_p = Float64(prev_p),
                prev_cr = Float64(prev_cr),
                prev_cp = Float64(prev_cp),
                pen_decon = Float64(pen_decon),
                treat_calf = treat_calf/70,
                treat_dry = treat_dry/365,
                treat_lac = treat_lac/365,
                vacc_protection = vacc_protection,
                vacc_shedding = vacc_shedding,
                vacc_duration = vacc_duration,
                p_retreatment = p_retreatment,
                longevity = Float64(longevity)
          )
    end

   

    @time [animal_step!(herd) for _ in 1:simdays]

     animalModel = herd

    data = DataFrame(

    id = run,
    timestep= animalModel.sim.timestep,

    pop_r = animalModel.sim.pop_r,
    pop_s= animalModel.sim.pop_s,
    pop_p= animalModel.sim.pop_p,
    pop_d= animalModel.sim.pop_d,
    pop_rec_r= animalModel.sim.pop_rec_r,
    pop_rec_p= animalModel.sim.pop_rec_p,
    pop_car_p= animalModel.sim.pop_car_p,
    pop_car_r= animalModel.sim.pop_car_r,

    num_calves= animalModel.sim.num_calves,
    num_dh= animalModel.sim.num_dh,
    num_heifers= animalModel.sim.num_heifers,
    num_lactating= animalModel.sim.num_lactating,
    num_dry= animalModel.sim.num_dry,
    num_weaned = animalModel.sim.num_weaned,

    pop_er= animalModel.sim.pop_er,
    pop_ep= animalModel.sim.pop_ep,

    inf_calves= animalModel.sim.inf_calves,
    car_calves = animalModel.sim.car_calves,
    ri_calves = animalModel.sim.ri_calves,
    rc_calves = animalModel.sim.rc_calves,

    inf_weaned= animalModel.sim.inf_weaned,
    car_weaned = animalModel.sim.car_weaned,
    ri_weaned = animalModel.sim.ri_weaned,
    rc_weaned = animalModel.sim.rc_weaned,

    inf_heifers= animalModel.sim.inf_heifers,
    car_heifers = animalModel.sim.car_heifers,
    ri_heifers = animalModel.sim.ri_heifers,
    rc_heifers = animalModel.sim.rc_heifers,

    inf_dh= animalModel.sim.inf_dh,
    car_dh = animalModel.sim.car_dh,
    ri_dh = animalModel.sim.ri_dh,
    rc_dh = animalModel.sim.rc_dh,

    inf_lactating= animalModel.sim.inf_lac,
    car_lactating = animalModel.sim.car_lac,
    ri_lactating = animalModel.sim.ri_lac,
    rc_lactating = animalModel.sim.rc_lac,

    inf_dry= animalModel.sim.inf_dry,
    car_dry = animalModel.sim.car_dry,
    ri_dry = animalModel.sim.ri_dry,
    rc_dry = animalModel.sim.rc_dry,

   
    clinical= animalModel.sim.clinical,
    subclinical= animalModel.sim.subclinical,

    current_b1= animalModel.sim.current_b1,
    current_b2= animalModel.sim.current_b2,
    current_b3= animalModel.sim.current_b3,
    current_b4= animalModel.sim.current_b4, 

    current_spring = animalModel.current_spring,
    current_autumn = animalModel.current_autumn,

    
    optimal_stock = optimal_stock,
    treatment_prob = treatment_prob,
    vacc_rate = vacc_rate,
    fpt_rate = fpt_rate,
    pen_decon = pen_decon,
    treat_calf = treat_calf,
    treat_lac = treat_lac,

    iteration = run,
    system = system,
    simid = run


  )
  


  return data  

end




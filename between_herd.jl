
# The purpose of this script is to define the between-farm spread model
# and the functions required to mediate trades between farms.
# Import the packages required, installing if not already present ------------------------


# Create some structs that will be used to define model-specific objects -------------------------------

"""
Farm(id, neighbours_to, neighbours_from, trades_from, trades_to, ncows, system, herd, traded_to, traded_from, status)

Constructor, define a farm in the between-herd model along with a summary of farm-level attributes required for farm to farm trading.

# Fields
id = Farm ID, corresponds to farms in the farm_file.csv dataset \n
neighbours_to = Farms that the farm with this ID sends stock to \n
neighbours_from = Farms that the farm with this ID receives stock from \n
trades_from = Holding array of agents of type AnimalAgent that the farm with this ID can send to other farms. \n
trades_to = Holding array of agents of type AnimalAgent that the farm with this ID is receiving from another farm in a model step. \n
ncows = Optimal number of milking cattle this farm has. \n
system = Farm's calving system (1 = CS1, 2 = CS2, 3 = CS3) \n
herd = AnimalModel, within-herd model for this farm, stepped as each model day advances. \n
traded_to = Whether the farm has received animals in this trading step \n
traded_from = Whether the farm has sent animals in this trading step.\n
status = Infection status of the farm (0 = No Salmonella infection, 1 = Salmonella infection sensitive to antibiotics, 2 = AMR Salmonella present.)\n
"""
mutable struct Farm
    id::Int
    neighbours_to::Array{Int}
    neighbours_from::Array{Int}
    trades_from::Array{AnimalAgent}
    trades_to::Array{AnimalAgent}
    ncows::Int
    system::Int
    herd::AnimalModel
    traded_to::Bool 
    traded_from::Bool
    status::Int
    node::Int
    saleyards::Vector{Int}
end
"""
Saleyard(id, animals, neighbours_to, neighbours_from, trades_from, trades_to, heifers, dh, lactating, dry)
Constructor for making a Saleyard object, that will hold animals indirectly traded between farms
    # Fields
    * id = Identification number of the saleyard
    * animals = An array of AnimalAgents sold from farms that can be purchased by farms.
    * neighbours_to = farms that cattle can be sent to
    * neighbours_from = farms that cattle can be received from
    * trades_from = vector of sent cattle
    * trades_to = vector of received cattle
    * heifers...dry = vectors of cattle in each PLS.
"""
mutable struct Saleyard
    id::Int
    animals::Array{AnimalAgent}
    neighbours_to::Array{Int}
    neighbours_from::Array{Int}
    trades_from::Array{AnimalAgent}
    trades_to::Array{AnimalAgent}
    heifers::Array{Int}
    dh::Array{Int}
    lactating::Array{Int}
    dry::Array{Int}
end

"""
FarmModel(timestep, farms, rng, space, svars, movements, statuses, run)

A constructor for making an object that contains the between-herd model. A FarmModel object includes information on the progress of the model through time, holds ABMs for each farm, and data on the movements between farms and farm status.

# Fields
* timestep = Model timestep in days
* farms = An array of objects initialised with the type constructor Farm
* space = A DataFrame containing the edgelist of directional movements in the pre-specified graph loaded in network_structure.csv
* svars = A DataFrame containing pre-specified information on the parameters for each farm.
* movements = A DataFrame of movements between farms during each model run.
* statuses = A DataFrame of the infection status of each farm over the course of the model run.
* run = An integer tracking the stage of the simulation in days.
"""
mutable struct FarmModel
    timestep::Int
    farms::Array{Farm}
    saleyards::Array{Saleyard}
    rng::Xoshiro
    space::DataFrame
    svars::DataFrame
    movements::DataFrame
    statuses::DataFrame
    run::Int
end

"""
    initiliaseFarms(svars, run)

    A function to initialise the model at the beginning of each run. The initialiseFarms function:
        * Creates each farm according to the parameters in `farm_file.csv`
        * Determines the position of each farm in the network graph specified in `network_structure.csv`
        * Randomly assigns infections to farms in the network.
    
    # Arguments:
    * svars: A DataFrame of information on the parameters for each farm, containing the following fields.
        - vacc_rate: Probability that farm vaccinates an animal
        - fpt_rate: Probability that a calf born has failure of passive transfer
        - treatment_prob: Probability that a clinical animal with Salmonella is treated with antibiotics
        - optimal_stock: The optimal number of lactating cattle on this farm.
        - pen_decon: Probability that farm decontaminates calving pens between uses
        - treat_calf: Probability that a calf is treated with antibiotics other than for Salmonella (per day)
        - treat_lac: Probability that a lactating cow is treated with antibiotics other than for Salmonella (per day)
        - farm: ID of the farm
        - calving_system: Herd calving pattern (1 = CS1, 2 = CS2, 3 = CS3)
    * run: Integer specifying the model run

"""
function initialiseFarms(;
    svars::DataFrame,
    run::Int,
    space::DataFrame,
    farm_si::Int,
    farm_amrsi::Int
    )

    #Set up the network ------------------------------
    
    # The number of farms is specified in the farm file (loaded as svars)
    numfarms = length(svars.optimal_stock)
    
    # Define where farms send animals to and where they receive animals from
    
    # to_edges are the farms that a farm sends animals to, from edges are the farms that a farm receives animals from 
    to_edges = Array{Vector{Int}}(undef, numfarms)

  
    for src in 1:numfarms
       dat = filter(x -> x.src == src, space)
       to_edges[src] = dat.dst
    end


    from_edges = Array{Vector{Int}}(undef, numfarms)

  
    for dst in 1:numfarms
       dat = filter(x -> x.dst == dst, space)
       from_edges[dst] = dat.src
    end

    # Create DataFrames that will hold information on each model run. These will be held in the Farm object

    """
        movements

        A DataFrame of animal movements specifying the date of movement, the source farm, the destination farm and the status of each animal moved.
    """
    movements = DataFrame(
        date = Date(0),
        from = 0,
        to = 0,
        from_status = 0,
        to_status = 0,
        infected_movements = false,
        typeofinf = 999
    )

    """
        statuses

        A DataFrame of farm statuses at each model step, specifying the date the status was observed on, the farm id, the percentage of the population with different infection types and the number of stock in each class for that farm.
    """
    statuses = DataFrame(
        date = Date(0),
        farm = 0,
        status = 0,
        pop_p = 0.0,
        pop_r = 0.0,
        pop_cp = 0.0,
        pop_cr = 0.0,
        pop_s = 0.0,
        num_calves = 0,
        num_weaned = 0,
        num_dh = 0,
        num_heifers = 0,
        num_lac = 0,
        num_dry = 0
    )

    # Shuffle farm positions in the network at each model iteration

    nodes = shuffle!([1:1:numfarms;])

     
    svars.node_occupied = nodes

    # Export the attributes and positions of each farm for each model step

    # Initialise the farms

    # Create an array for farms to go into of length numfarms
    farms = Array{Farm}(undef,numfarms)

    saleyards = Array{Saleyard}(undef, 2)
    yards = [500, 501]
    for i in eachindex(yards)
        animals = Array{AnimalAgent}[]
        heifers = Array{AnimalAgent}[]
        dh = Array{AnimalAgent}[]
        lactating = Array{AnimalAgent}[]
        dry = Array{AnimalAgent}[]
        id = yards[i]
        trades_from = []
        trades_to = []
        neighbours_to = neighbours_from = svars[(svars.saleyard_1 .== yards[i]) .| (svars.saleyard_2 .== yards[i]), :].farm
        saleyards[i] = Saleyard(id,animals, neighbours_to, neighbours_from, trades_from, trades_to, heifers, dh, lactating, dry)
    end


    # Create the model, an object of type farmModel, with the farms, space, parameters and RNG specified. Supply empty DataFrames for farm information.
    farmModel = FarmModel(0,farms, saleyards, Xoshiro(run), space, svars, movements, statuses, run)

    numres = Int(ceil(farm_amrsi/100*numfarms))
    numsi = Int(round(farm_si/100*numfarms))

    @info numres
    @info numsi

    # Seed infection randomly into two farms, one sensitive farm (inf_farm), one resistant farm (res_farm)
    res_farm = sample(farmModel.rng, svars.farm, numres)
    inf_farms = sample(farmModel.rng, findall(x-> x != res_farm, svars.farm), numsi, replace = false)

    saleyard_connections = Vector{Vector{Int}}(undef, length(svars.farm))
    for i in eachindex(saleyard_connections)
        saleyard_connections[i] = filter!(x-> x != 999, [svars.saleyard_1[i], svars.saleyard_2[i]])
    end
    shuffle!(saleyard_connections)

    #Set up the initial population (MT)
   Threads.@threads for i in eachindex(svars.farm)

        id = svars.farm[i]

        node = nodes[id]
        
        if id in res_farm
            status = 2
        elseif id in inf_farms
            status = 1
        else 
            status = 0
        end 

        if status == 2 @info "The resistant farm is $id" end
        # Set initial prevalence in infected farms
        
        if status == 1
            prev_r = Float64(0.00)
            prev_p = Float64(svars.init_p[i])
            prev_cr = Float64(0.00)
            prev_cp = Float64(svars.init_cp[i])
        elseif status == 2
            prev_r = Float64(svars.init_r[i])
            prev_p = Float64(0.00)
            prev_cr = Float64(svars.init_cr[i])
            prev_cp = Float64(0)
        elseif status == 0
            prev_r = Float64(0.0)
            prev_p = Float64(0.00)
            prev_cr = Float64(0.0)
            prev_cp = Float64(0)
        end

        # Specify the farm's outgoing and incoming trading partners
        neighbours_to = to_edges[node]
        neighbours_from = from_edges[node]

        # Create empty vectors for trades
        trades_from = []
        trades_to = []

        # Set the farm's attributes based on the parameters for that ID in svars (the farm file)
        ncows = svars.optimal_stock[i]
        system = svars.calving_system[i]
        pen_decon = svars.pen_decon[i]
        optimal_stock = Int(svars.optimal_stock[i])
        treatment_prob = Float64(svars.treatment_prob[i])
        vacc_rate = Float64(svars.vacc_rate[i])
        fpt_rate = Float64(svars.fpt_rate[i])
        treat_calf = svars.treat_calf[i]
        treat_lac = svars.treat_lac[i]
        saleyards = saleyard_connections[i]
        treat_dry = svars.treat_dry[i]
        p_retreatment = svars.p_retreatment[i]
        treatment_length = svars.treatment_length[i]
        vacc_duration = svars.vacc_duration[i]
        vacc_shedding = svars.vacc_shedding[i]
        vacc_protection = svars.vacc_protection[i]
        density_lactating = svars.density_lactating[i]
        density_dry = svars.density_dry[i]
        density_calf = svars.density_calves[i]
        longevity = svars.longevity[i]

        # Choose an initialisation function based on whether the farm is CS1 (spring), CS2 (split) or CS3 (batch)
        if svars.calving_system[i] == 1
            herd = initialiseSpring(
                    farmno = Int(i),
                    system = Int(1),
                    msd = Date(2021,9,24),
                    seed = Int(id + run),
                    optimal_stock = optimal_stock,
                    optimal_lactating = optimal_stock,
                    treatment_prob = treatment_prob,
                    treatment_length = Int(treatment_length),
                    carrier_prob = Float64(0.1),
                    timestep = Int(0),
                    density_lactating = Int(density_lactating),
                    density_dry = Int(density_dry),
                    density_calves = Float64(density_calf),
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
        elseif svars.calving_system[i] == 2
            herd = initialiseSplit(
                farmno = Int(i),
                system = Int(2),
                msd = Date(2021,9,24),
                    seed = Int(id + run),
                    optimal_stock = optimal_stock,
                    optimal_lactating = optimal_stock,
                    treatment_prob = treatment_prob,
                    treatment_length = Int(treatment_length),
                    carrier_prob = Float64(0.1),
                    timestep = Int(0),
                    density_lactating = Int(density_lactating),
                    density_dry = Int(density_dry),
                    density_calves = Float64(density_calf),
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
                farmno = Int(i),
                system = Int(3),
                msd = Date(2021,9,24),
                seed = Int(id + run),
                optimal_stock = optimal_stock,
                optimal_lactating = optimal_stock,
                treatment_prob = treatment_prob,
                treatment_length = Int(treatment_length),
                carrier_prob = Float64(0.1),
                timestep = Int(0),
                density_lactating = Int(density_lactating),
                density_dry = Int(density_dry),
                density_calves = Float64(density_calf),
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

        # Farms have not traded in this step
        traded_to = false
        traded_from = false
        # Take the attributes specified in this step and create an object ot type Farm
        thisfarm = Farm(id, neighbours_to, neighbours_from, trades_from, trades_to, ncows, system, herd, traded_to, traded_from, status, node, saleyards)

        # Place that object into the array of farms in farmModel.
        farmModel.farms[i] = thisfarm
       

    end

    # Record the neighbours of each farm at each iteration and the indegrees and outdegrees of each farm at the node it currently occupies
    all_to = []
    all_from = []
    in_deg = []
    out_deg = []
    for farm in farmModel.farms 
        push!(all_to, farm.neighbours_to)
        push!(out_deg, length(farm.neighbours_to))
        push!(all_from, farm.neighbours_from)
        push!(in_deg, length(farm.neighbours_from))
    end

    svars.neighbours_to = all_to
    svars.neighbours_from = all_from
    svars.in_deg = in_deg
    svars.out_deg = out_deg
    svars.saleyard_connections = saleyard_connections

    CSV.write("./export/farm_attributes_run_$run.csv",DataFrame(svars))

    farmModel.timestep = 1

    return farmModel

    end


"""
sell_to_saleyard!
Farms selling cattle to saleyard.
"""
function sell_to_saleyard!(farmModel::FarmModel, sending_farm::Farm)
        sending_farm.traded_from == true && return
        sending_farm.trades_from = Vector{AnimalAgent}()
        length(sending_farm.herd.sending) == 0 && return
        length(sending_farm.saleyards) == 0 && return
        purchasing_yard = length(sending_farm.saleyards) == 2 ? sample(farmModel.rng, sending_farm.saleyards)[1] : sending_farm.saleyards[1]
        purchasing_yard == 999 && return
        
        herd = sending_farm.herd
        trading_abilities = Vector{Int}(undef, 4)
        trading_abilities[1] = herd.current_heifers - herd.optimal_heifers
        trading_abilities[2] = herd.current_dh - herd.optimal_dh
        trading_abilities[3] = herd.current_lactating - herd.optimal_lactating
        trading_abilities[4] = herd.current_dry - herd.optimal_dry

        disposal_power = Int(round(rand( Truncated(Rayleigh(round(0.02*herd.current_stock)), round(0.0025*herd.current_stock), round(0.04*herd.current_stock)))))
        
        trading_abilities = max.(trading_abilities, 0)

        
        fromeach = disposal_power > sum(trading_abilities) ? 1 : disposal_power/sum(trading_abilities)

        trading_proportions = round.(fromeach*trading_abilities)

        avaialable_stages = Vector{Int}()
        for animal in herd.sending
            isassigned(herd.animals, animal) == false && continue
            push!(avaialable_stages, herd.animals[animal].stage) 
        end 
        for_sale = Vector{Vector{Int}}(undef, 4)
        for_sale[1] = herd.sending[findall(x-> x == 3, avaialable_stages)]
        for_sale[2] = herd.sending[findall(x-> x == 4, avaialable_stages)]
        for_sale[3] = herd.sending[findall(x-> x == 5, avaialable_stages)]
        for_sale[4] = herd.sending[findall(x-> x == 6, avaialable_stages)]
        
        sum([length(i) for i in for_sale]) < 5 && return

        for i in eachindex(for_sale)
            length(for_sale[i]) == 0 && continue
            trading_proportions[i] <= 0 && continue
            to_trade = trading_proportions[i] < length(for_sale[i]) ? trading_proportions[i] : length(for_sale[i])
            to_trade == 0 && continue
            trading = unique(sample(farmModel.rng, for_sale[i], Int(to_trade)))
            # @info trading
            for trade in trading
                push!(sending_farm.trades_from, herd.animals[trade])
            end
        end
        
        length(sending_farm.trades_from) == 0 && return

        receiving_yard = farmModel.saleyards[findfirst(x-> x.id == purchasing_yard, farmModel.saleyards)]
        receiving_yard.trades_to = sending_farm.trades_from

        id_counter = length(receiving_yard.animals) == 0 ? 1 : maximum(x.id for x in receiving_yard.animals) 

        for bought in 1:length(receiving_yard.trades_to)
            isassigned(receiving_yard.trades_to, bought) == false && continue
            #@info "sold to saleyard"
            id_counter += 1
            receiving_yard.trades_to[bought].id = id_counter
            push!(receiving_yard.animals, receiving_yard.trades_to[bought])
        end

        for sold in 1:length(sending_farm.trades_from)
                isassigned(sending_farm.trades_from, sold) == false && continue
                findfirst(isequal(sending_farm.trades_from[sold]), sending_farm.herd.animals) === nothing && continue
                deleteat!(sending_farm.herd.animals, findfirst(isequal(sending_farm.trades_from[sold]), sending_farm.herd.animals))
                push!(farmModel.movements.from, sending_farm.id)
                push!(farmModel.movements.to, receiving_yard.id)
                push!(farmModel.movements.from_status, sending_farm.status)
                push!(farmModel.movements.to_status, 0)
                push!(farmModel.movements.infected_movements, ifelse(sending_farm.trades_from[sold].status ∉ [0,7,8], true, false))
                push!(farmModel.movements.typeofinf, sending_farm.trades_from[sold].status)
                push!(farmModel.movements.date, sending_farm.herd.date)
            end

        receiving_yard.trades_to = Vector{AnimalAgent}()
        draft_saleyard!(farmModel)# draft out the cattle in the yard by PLS.
        sending_farm.traded_from = true
            
end
    
"""
draft_saleyard!
Draft out cattle in the saleyard by PLS
"""
function draft_saleyard!(farmModel::FarmModel)
    for yard in farmModel.saleyards
        yard.heifers = findall(x-> x.stage == 3, yard.animals )
        yard.dh = findall(x -> x.stage == 4, yard.animals)
        yard.lactating = findall(x -> x.stage == 5, yard.animals)
        yard.dry = findall(x -> x.stage == 6, yard.animals)
    end
end
    
"""
buy_from_saleyard!
Farm purchasing of cattle from saleyard. Cattle shed when moved.
"""
function buy_from_saleyard!(farmModel::FarmModel, receiving_farm::Farm)
        receiving_farm.traded_to == true && return
        length(receiving_farm.saleyards) == 0 && return
        herd = receiving_farm.herd
        # Determine trading needs
        trading_needs = Vector{Int}(undef, 4)
        trading_needs[1] = herd.optimal_heifers - herd.current_heifers
        trading_needs[2] = herd.optimal_dh - herd.current_dh
        trading_needs[3] = herd.optimal_lactating - herd.current_lactating
        trading_needs[4] = herd.optimal_dry - herd.current_dry

        trading_needs = max.(trading_needs, 0)

        purchasing_power = Int(round(rand(farmModel.rng, Truncated(Rayleigh(round(0.02*herd.current_stock)), round(0.0025*herd.current_stock), round(0.04*herd.current_stock)))))
        
        toeach = purchasing_power > sum(trading_needs) ? 1 : purchasing_power/sum(trading_needs)

        trading_proportions = round.(toeach*trading_needs)

        
        
        selling_yard = sample(farmModel.rng, receiving_farm.saleyards, 1)[1]
        selling_yard == 999 && return
        selling_yard = farmModel.saleyards[findfirst(x-> x.id == selling_yard, farmModel.saleyards)]

        for_sale = Vector{Vector{Int}}(undef, 4)
        for_sale[1] = selling_yard.heifers
        for_sale[2] = selling_yard.dh
        for_sale[3] = selling_yard.lactating
        for_sale[4] = selling_yard.dry

        sum([length(i) for i in for_sale]) < 5 && return

        purchased = Vector{Int}()

        for i in eachindex(trading_needs)
            trading_needs[i] <= 0 && continue
            length(for_sale[i]) == 0 && continue
            to_purchase = trading_proportions[i] > length(for_sale[i]) ? length(for_sale[i]) : trading_proportions[i]
            to_purchase == 0 && continue
            sold = unique(sample(farmModel.rng, for_sale[i], Int(to_purchase)))
            for i in sold push!(purchased, i) end
        end
        
        # No single cow trades
        length(purchased) < 5 && return

        didtrades = 0
        yard = selling_yard
        for sent in purchased
            # @info "Purchased from saleyard"
            isassigned(yard.animals, sent) == false && continue
            push!(receiving_farm.trades_to, yard.animals[sent])
            deleteat!(yard.animals, sent)
        end

        for received in receiving_farm.trades_to
            herd.id_counter += 1
            received.id = herd.id_counter
            received.stress = true
            animal_shedding!(received)
            push!(herd.animals, received)
            push!(farmModel.movements.from, yard.id)
            push!(farmModel.movements.to, receiving_farm.id)
            push!(farmModel.movements.from_status, yard.id)
            push!(farmModel.movements.to_status, receiving_farm.status)
            push!(farmModel.movements.infected_movements, ifelse(received.status ∉ [0,7,8], true, false))
            push!(farmModel.movements.typeofinf, received.status)
            push!(farmModel.movements.date, herd.date)
            didtrades += 1
        end

        receiving_farm.trades_to = Vector{AnimalAgent}[]
        
        didtrades == 0 && return

        receiving_farm.traded_to = true
end

"""
selling
Find all sending animals from the selling farm.
"""
function selling(for_sale, index, selling_herd, avaialable_stages)
    for_sale[index] = selling_herd.sending[findall(x-> x == (index + 2), avaialable_stages)]
end


"""
between_farm_trade!
Broker needs-based trading between farms.
"""
function between_farm_trade!(farmModel::FarmModel, receiving_farm::Farm)
    receiving_farm.traded_to == true && return
    length(receiving_farm.neighbours_from) == 0 && return

    sending_trader = farmModel.farms[receiving_farm.neighbours_from[rand(farmModel.rng, eachindex(receiving_farm.neighbours_from))]]
    length(sending_trader.herd.sending) < 5 && return
    sending_trader.trades_from = Vector{AnimalAgent}()
    receiving_herd = receiving_farm.herd
    
    # Determine the acquisition needs of the receiving farm
    trading_needs = Vector{Int}(undef, 4)
    trading_needs[1] = receiving_herd.optimal_heifers - receiving_herd.current_heifers
    trading_needs[2] = receiving_herd.optimal_dh - receiving_herd.current_dh
    trading_needs[3] = receiving_herd.optimal_lactating - receiving_herd.current_lactating
    trading_needs[4] = receiving_herd.optimal_dry - receiving_herd.current_dry

    trading_needs = max.(trading_needs, 0)

    purchasing_power = Int(round(rand(farmModel.rng, Truncated(Rayleigh(round(0.025*receiving_herd.current_stock)), round(0.025*receiving_herd.current_stock), round(0.05*receiving_herd.current_stock)))))
    
    toeach = purchasing_power > sum(trading_needs) ? 1 : purchasing_power/sum(trading_needs)

    receiver_needs = round.(toeach*trading_needs)


    # Determine the selling abilities of the sending farm
    selling_herd = sending_trader.herd

    trading_abilities = Vector{Int}(undef, 4)
    trading_abilities[1] = selling_herd.current_heifers - selling_herd.optimal_heifers
    trading_abilities[2] = selling_herd.current_dh - selling_herd.optimal_dh
    trading_abilities[3] = selling_herd.current_lactating - selling_herd.optimal_lactating
    trading_abilities[4] = selling_herd.current_dry - selling_herd.optimal_dry

    trading_abilities = max.(trading_abilities, 0)

    disposal_power = Int(round(rand(farmModel.rng, Truncated(Rayleigh(round(0.025*selling_herd.current_stock)), round(0.025*selling_herd.current_stock), round(0.05*selling_herd.current_stock)))))
    
    fromeach = disposal_power > sum(trading_abilities) ? 1 : disposal_power/sum(trading_abilities)

    seller_abilities = round.(fromeach*trading_abilities)

    needs = sum(receiver_needs)
    abilities = sum(seller_abilities)
    
    tradeno =  needs > abilities ? abilities : needs

    tradeno = Int(tradeno > length(sending_trader.herd.sending) ? length(sending_trader.herd.sending) : tradeno)

    tradeno < 4 && return

    for_sale = shuffle(sending_trader.herd.sending)[1:tradeno]

    selling_herd = sending_trader.herd
    for i in eachindex(for_sale)
            isassigned(selling_herd.animals, for_sale[i]) == false && continue
            push!(sending_trader.trades_from, selling_herd.animals[for_sale[i]])
        # end 
    end

    length(sending_trader.trades_from) == 0 && return

    receiving_farm.trades_to = sending_trader.trades_from
        didtrades = 0
        
        # Move the purchased animals onto the receiving farm
        for bought in 1:length(receiving_farm.trades_to)
            isassigned(receiving_farm.trades_to, bought) == false && continue
            receiving_farm.herd.id_counter += 1
            receiving_farm.trades_to[bought].id = receiving_farm.herd.id_counter
            receiving_farm.trades_to[bought].stress = true
            animal_shedding!(receiving_farm.trades_to[bought])
            push!(receiving_farm.herd.animals, receiving_farm.trades_to[bought])
            didtrades += 1
        end

        didtrades == 0 && return

        # Set the status of each farm to prevent further trading in this step
        receiving_farm.traded_to = true
        sending_trader.traded_from = true
        for sold in 1:length(sending_trader.trades_from)
            isassigned(sending_trader.trades_from, sold) == false && continue
            findfirst(isequal(sending_trader.trades_from[sold]), sending_trader.herd.animals) === nothing && continue
            deleteat!(sending_trader.herd.animals, findfirst(isequal(sending_trader.trades_from[sold]), sending_trader.herd.animals))
            push!(farmModel.movements.from, sending_trader.id)
            push!(farmModel.movements.to, receiving_farm.id)
            push!(farmModel.movements.from_status, sending_trader.status)
            push!(farmModel.movements.to_status, receiving_farm.status)
            push!(farmModel.movements.infected_movements, ifelse(sending_trader.trades_from[sold].status ∉ [0,7,8], true, false))
            push!(farmModel.movements.typeofinf, sending_trader.trades_from[sold].status)
            push!(farmModel.movements.date, receiving_farm.herd.date)
        end

end
    

"""
trades!
Seasonally mediate trading probabilities.
"""
function trades!(farmModel::FarmModel)
    current_month = month(farmModel.farms[1].herd.date)
    summer = [12,1,2]
    autumn = [3,4,5]
    winter = [6,7,8]
    spring = [9,10,11]
    
        for farm in farmModel.farms
        current_month in summer && rand(farmModel.rng) > 0.03/2 && continue
        current_month in autumn && rand(farmModel.rng) > 0.04/2 && continue
        current_month in winter && rand(farmModel.rng) > 0.05/2 && continue
        current_month in spring && rand(farmModel.rng) > 0.06/2 && continue
        between_farm_trade!(farmModel, farm)
        if rand(farmModel.rng, Bool)
            buy_from_saleyard!(farmModel, farm)
        else 
            sell_to_saleyard!(farmModel, farm)
        end

    end
end

"""
export_statuses!(farmModel)

A function that iterates over farms in the farmModel object, exporting data to DataFrames recording farm status at each model step

#Args 
* farmModel: An object of type FarmModel.
"""
function export_statuses!(farmModel::FarmModel)
    # Iterate over each farm in farmModel. A farm is an object of type AnimalModel
    for farm in farmModel.farms
    statuses = farmModel.statuses
        #push!(statuses, farmModel.statuses)
        # Append date, farm_id, status to the statuses DataFrame
        push!(statuses.date, farm.herd.date)
        push!(statuses.farm, farm.id)
        push!(statuses.status, farm.status)
    

        critters = length(farm.herd.animals)
        # Append information on the within-herd prevalence of resistant and sensitive Salmonella infections at this model step.
        # Each prevalence is the number of animals with that infection type over the total in the herd

        push!(statuses.pop_p, 100*farm.herd.pop_p/critters)
        push!(statuses.pop_r, 100*farm.herd.pop_r/critters)
        push!(statuses.pop_cp, 100*farm.herd.sim.pop_car_p[farmModel.timestep]/critters)
        push!(statuses.pop_cr, 100*farm.herd.sim.pop_car_r[farmModel.timestep]/critters)
        push!(statuses.pop_s, 100*farm.herd.pop_s/critters)

        lifestage = countmap([animal.stage for animal in farm.herd.animals])

        # Append information on the number of animals in each stock class
        push!(statuses.num_calves, get(lifestage, 1, 0))
        push!(statuses.num_weaned, get(lifestage, 2, 0))
        push!(statuses.num_heifers, get(lifestage, 3, 0))
        push!(statuses.num_dh, get(lifestage, 4, 0))
        push!(statuses.num_lac, get(lifestage, 5, 0))
        push!(statuses.num_dry, get(lifestage, 6, 0))

    end
end

"""
farm_step!(farmModel)

The farm_step! function advances the farmModel one timestep, stepping the each farm's AnimalModel, adding trading information and determining the status of the each farm at each step.
"""
    function farm_step!(farmModel::FarmModel)
    

    # Clear previous trading status at the beginning of each step, empty the trading from and to vectors, set the status.
    Threads.@threads for farm in farmModel.farms
        farm.traded_from = false
        farm.traded_to = false
        farm.trades_from = []
        farm.trades_to = []
        
        if (farm.herd.pop_r > 0) | (farm.herd.sim.pop_car_r[farmModel.timestep] > 0)
            farm.status = 2
        elseif (farm.herd.pop_r == 0) & (farm.herd.sim.pop_car_r[farmModel.timestep] == 0) & 
                ((farm.herd.pop_p > 0) | (farm.herd.sim.pop_car_p[farmModel.timestep] > 0))
            farm.status = 1
        elseif (farm.herd.pop_r == 0) & (farm.herd.sim.pop_car_r[farmModel.timestep] == 0) & 
                (farm.herd.pop_p == 0) & (farm.herd.sim.pop_car_p[farmModel.timestep] == 0)
            farm.status = 0
        end

    end


    if farmModel.timestep > 365*2 
    trades!(farmModel) 
    export_statuses!(farmModel)
    end
    
    # Advance the timestep by one day

    farmstepper = Vector{Task}(undef, length(farmModel.farms))

    # Spawn tasks for each farm in farmModel.farms
    for i in eachindex(farmModel.farms)
        farmstepper[i] = Threads.@spawn animal_step!(farmModel.farms[i].herd)
    end

    
    # Wait until the spawned tasks have finished.
    for i in eachindex(farmstepper)
        wait(farmstepper[i])
    end

    farmModel.timestep += 1
end

function run_farm_model!(farm_file::DataFrame, numruns::Int, numdays::Int, farm_si::Int, farm_amrsi::Int, space::DataFrame)

    statuses = DataFrame[]
    movements = DataFrame[]
   # spaces = Vector{DataFrame}[]


    # Iterate running the model over the number of runs specified 
  @time  for run in 1:numruns
        run = 0000 + run#For large scale runs this is managed using a shell script.
         @time farmModel = initialiseFarms(svars = farm_file, run = run, space = space, farm_si = farm_si, farm_amrsi = farm_amrsi);
         @time  for _ in 1:numdays
                   farm_step!(farmModel);
                end
        
        @info "Run $run complete."

        @info farmModel.statuses
        @info farmModel.movements
        
        @info "before"
        push!(statuses, farmModel.statuses)
        push!(movements, farmModel.movements)
                
        @info "after"
       
      #  push!(spaces, farmModel.spaces)
       
        # At the end of each model run, export information on the run to disk.
       #=  CSV.write("./export/farm_run_statuses_$run.csv", farmModel.statuses)
        CSV.write("./export/farm_run_movements_$run.csv", farmModel.movements)
        CSV.write("./export/farm_run_space_$run.csv",farmModel.space) =#
        farmModel = nothing
        GC.gc() # Probably unecessary, but avoids excesive RAM usage with slightly more aggressive GC.
    end

    return statuses, movements
end 






    
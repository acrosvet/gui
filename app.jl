using GenieFramework, DataFrames, CSV, PlotlyBase, LatinHypercubeSampling
@genietools


include("model.jl")
include("between_herd.jl")
include("make_farm_network.jl")
@app begin


    #Tab selection 
    @in selected_component = "parmtab"
    @in left_drawer_open = true
    @in ministate = true


  
    #Default button and conditional states
    @in run_model::Bool = false
    @in lhs_draw::Bool = false
    @in run_ensemble::Bool = false
    @in gen_ensemble_plots::Bool = false
    @in gen_spread_plots::Bool = false
    @out plotting_completed::Bool = false
    @out single_sim_completed::Bool = false
    @out between_herd_plotting::Bool = false
    @in gen_farm_network::Bool = false
    @in gen_farm_parms::Bool = false
    @in run_between_herd::Bool = false
    @in gen_spread_networks::Bool = false
    @in generated_networks::Bool = false



    #Prefill values
    #Single herd simulation
    @in system = "Spring"
    @in optimal_stock = 273
    @in vacc_rate = 0.5
    @in fpt_rate = 0.6
    @in treatment_probability = 0.5
    @in pen_decon = "FALSE"
    @in simdays = 5
    @in init_r = 1.5
    @in init_p = 1.5
    @in init_cr = 2.5
    @in init_cp = 2.5
    @in treatment_length = 5
    @in density_lactating = 50
    @in density_dry = 250
    @in density_calf = 3
    @in treat_dry = 0.01
    @in vacc_protection = 0.5
    @in vacc_shedding = 0.5
    @in vacc_duration = 0.5
    @in p_retreatment = 0.2
    @in longevity = 8
    @in density_dry = 250

    #Ensemble simulations
    @in optimal_stock_min = 80
    @in optimal_stock_max = 1000

    @in pen_decon_min = 0
    @in pen_decon_max = 1
    @in fpt_rate_min = 0
    @in fpt_rate_max = 1
    @in treat_prob_min = 0
    @in treat_prob_max = 1
    @in treat_calf_min = 0
    @in treat_calf_max = 1
    @in treat_lac_min = 0
    @in treat_lac_max = 1
    @in treat_dry_min = 0
    @in treat_dry_max = 1
    @in treat_length_min = 1
    @in treat_length_max = 5
    @in density_calf_min = 2
    @in density_calf_max = 5
    @in density_lactating_min = 45
    @in density_lactating_max = 55
    @in density_dry_min = 200
    @in density_dry_max = 300
    @in vacc_protection_min = 0.2
    @in vacc_protection_max = 0.9
    @in vacc_shedding_min = 0.2
    @in vacc_shedding_max = 0.9
    @in vacc_duration_min = 0.2
    @in vacc_duration_max = 0.9
    @in p_retreatment_min = 0
    @in p_retreatment_max = 0.8
    @in longevity_min = 7
    @in longevity_max = 10
    @in simno = 10
    @in prob_spring = 0.33
    @in prob_split = 0.33
    @in prob_batch = 0.33
    @in vacc_rate_min = 0
    @in vacc_rate_max = 1
    @in ensemble_length = 5
    @in init_r_min = 1
    @in init_r_max = 1.5
    @in init_cr_min = 2.5
    @in init_cr_max = 3
    @in init_p_min = 1
    @in init_p_max = 1.5
    @in init_cp_min = 2.5
    @in init_cp_max = 3

    # Between herd simulations
    @in numsims = 1
    @in numfarms = 50
    @in numedges = 103
    @in farm_si = 70
    @in farm_amrsi = 1
    @in in_deg = 10.8
    @in out_deg = 11.3
    @in smallest_farm = 80
    @in mean_farm = 220
    @in largest_farm = 500
    @in vacc_rate_min_f = 0
    @in vacc_rate_max_f = 1
    @in pen_decon_f = 0.5
    @in fpt_rate_min_f = 0.5
    @in fpt_rate_max_f = 0.9
    @in treat_prob_min_f = 0.1
    @in treat_prob_max_f = 0.9
    @in treat_lac_min_f = 0.03
    @in treat_lac_max_f = 0.06
    @in treat_calf_min_f = 0.05
    @in treat_calf_max_f = 0.25
    @in treat_dry_min_f = 0.01
    @in treat_dry_max_f = 0.02
    @in treat_length_f = 5
    @in density_calf_f = 3
    @in density_lactating_f = 50
    @in density_dry_f = 250
    @in vacc_protection_f = 0.5
    @in vacc_shedding_f = 0.5
    @in vacc_duration_f = 0.5
    @in p_retreatment_min_f = 0.1
    @in p_retreatment_max_f = 0.5
    @in longevity_min_f = 7
    @in longevity_max_f = 9
    @in init_r_min_f = 1.5
    @in init_r_max_f = 3
    @in init_cr_min_f = 3
    @in init_cr_max_f  = 5
    @in init_p_min_f = 1.5
    @in init_p_max_f = 3
    @in init_cp_min_f = 3
    @in init_cp_max_f = 5 
    @in prob_spring_f = 0.7
    @in prob_split_f = 0.25
    @in prob_batch_f = 0.05
    
   # @load "iconuri.jld2"

    
    # Output vals for within-herd (these are static)
    @out statuses = ["Infected (SI)", "Infected (AMRSI)", "Uninfected", "Infected (AMRSI and SI)"]
    @out systems = ["Spring", "Split", "Batch"]
    @out decon_options = ["TRUE", "FALSE"]
    @in treat_lac = 0.05
    @in treat_calf = 0.1

    # Timer functions
    @out runstatus::String = "Click to 'Run simulation' (once) to run."
    @out ensemblestatus::String = "Click 'Run ensemble' (once) to run."
    @out betweenstatus::String = "Click 'Run between-herd model' (once) to run."


    #Tables 
    @out simtable_pagination = DataTablePagination(rows_per_page=100)
    @out simtable = DataTable()
  
        
    # Plot initialisers
    @out sim_trace = []
    @out demo_trace = []
    @out rec_trace = []
    @out calf_trace = []
    @out weaned_trace = []
    @out heifer_trace = []
    @out dh_trace = []
    @out lactating_trace = []
    @out dry_trace = []
    @out ensemble_amrsi_trace = []
    @out ensemble_amrsi_car_trace = []
    @out ensemble_si_car_trace = []
    @out ensemble_si_trace = []

    


    function pls_layout(pls)
        PlotlyBase.Layout(
        title = pls,
        xaxis_title = "Simulation day",
        yaxis_title = "Number of cattle",
        height = 350,
        plot_bgcolor = "white",
        paper_bgcolor = "white"
    )
    end

    @out post_sna_layout =  PlotlyBase.Layout(
        paper_bgcolor="rgba(255,255,255,1)",
        plot_bgcolor="rgba(255,255,255,1)",
        hovermode="closest",
        #title="Farm trading network structure",
        titlefont_size=16,
        showlegend=false,
        showarrow=false,
        xaxis=attr(showgrid=false, zeroline=false, showticklabels=false),
        yaxis=attr(showgrid=false, zeroline=false, showticklabels=false),
        height =  700
    )

    @out post_sna = []
    @out pre_sna = []

    @out calf_layout = pls_layout("Calves")
    @out weaned_layout = pls_layout("Weaned")
    @out heifer_layout = pls_layout("Heifers")
    @out dh_layout = pls_layout("Pregnant heifers")
    @out lactating_layout = pls_layout("Lactating")
    @out dry_layout = pls_layout("Dry")

    @out sim_layout = PlotlyBase.Layout(
        title = "Infected",
        xaxis_title = "Simulation day",
        yaxis_title = "Point prevalence (%)",
        height = 350,
        plot_bgcolor = "white",
        paper_bgcolor = "white"
    )

    @out rec_layout = PlotlyBase.Layout(
        title = "Susceptible and recovered",
        xaxis_title = "Simulation day",
        yaxis_title = "Percentage of cattle",
        height = 350,
        plot_bgcolor = "white",
        paper_bgcolor = "white",
        yaxis_range = [0,100]
    )

    @out demo_layout = PlotlyBase.Layout(
        title = "Population demographics",
        xaxis_title = "Simulation day",
        yaxis_title = "Number of cattle",
        height = 350,
        plot_bgcolor = "white",
        paper_bgcolor = "white"
    )

   
    function ensemble_plot(type)
        PlotlyBase.Layout(
            title = type,
            xaxis_title = "Simulation day",
            yaxis_title = "Point prevalence (%)",
            height = 350,
            plot_bgcolor = "white",
            paper_bgcolor = "white"#,
           # yaxis = attr(range=[0,100])
        )
    
    end

    @out between_prev_amrsi_layout = ensemble_plot("Between-herd prevalence: AMRSI")
    @out between_prev_amrsi = []
    @out between_prev_si = []
    @out between_prev_si_layout = ensemble_plot("Between-herd prevalence: SI")
    @out between_prev_uninf_layout = ensemble_plot("Between-herd prevalence: Uninfected")
    @out between_prev_uninf = []
    @out within_prev_amrsi = []
    @out within_prev_amrsi_layout = ensemble_plot("Within-herd prevalence: Active AMRSI")
    @out within_prev_si = []
    @out within_prev_si_layout = ensemble_plot("Within-herd prevalence: Active SI")
    @out within_prev_uninf = []
    @out within_prev_uninf_layout = ensemble_plot("Within-herd prevalence: Uninfected")

    @out ensemble_amrsi_layout = ensemble_plot("AMRSI (active)")
    @out ensemble_amrsi_car_layout = ensemble_plot("AMRSI (carrier)")
    @out ensemble_si_layout = ensemble_plot("SI (active)")
    @out ensemble_si_car_layout = ensemble_plot("SI (carrier)")
    @out ensemble_si_layout = ensemble_plot("SI (active)")
    @out ensemble_recovered_layout = ensemble_plot("Recovered")
    @out ensemble_susceptible_layout = ensemble_plot("Susceptible") 
    @out ensemble_recovered_trace = []
    @out ensemble_susceptible_trace = []

    #Bignumbers 
    @out max_prev_si::String = "0"
    @out max_prev_amrsi::String = "0"
    @out max_prev_active_si::String = "0"
    @out max_prev_active_amrsi::String = "0"
    @out max_prev_car_si::String = "0"
    @out max_prev_car_amrsi::String = "0"

    function bind_ensemble()
        ensemble_res = DataFrame()

        path = "./private/"
        files = readdir(path)  
        ensemble_results = filter(x -> occursin("raw", x), files)


        for file in ensemble_results
            @info file
            filepath = joinpath(path, file)
            res_df = DataFrame(CSV.File(filepath))
            run_number = parse(Int, match(r"\d+", file).match)
            res_df[!, :run] = fill(run_number, nrow(res_df))
                ensemble_res = vcat(ensemble_res, res_df, cols = :union)
        end

        return(ensemble_res)
    end 



    
    @onbutton run_model begin 
        running = true
        @info "System is $system"
        @info "Herd size is $optimal_stock"
        @info "Vacc rate is $vacc_rate"
        @info "FPT rate is $fpt_rate"

       
        
        decons = Dict("TRUE"=>1, "FALSE"=>0)
        sys_select = Dict("Spring"=>1, "Split"=>2, "Batch"=>3)
       
        @info sys_select[system]
        @info 365*(simdays+2)
        
        function start_timer(task_to_monitor::Task)
            start_time = Base.time()
            runstatus = "Simulation running"
            try
                while !istaskdone(task_to_monitor)
                    elapsed_time = Base.time() - start_time
                    runstatus = string("Running simulation, time elapsed: ", round(elapsed_time; digits=2), " sec")
                    sleep(0.1)
                end
            finally
                elapsed_time = Base.time() - start_time
                runstatus = string("Run completed in: ", round(elapsed_time; digits=2), " sec")
            end
            return runstatus  
        end
        
        function run_sim()
            simdata = run_herd!(optimal_stock = optimal_stock, calving_system = sys_select[system], pen_decon = decons[pen_decon], treatment_prob = treatment_probability, vacc_rate = vacc_rate, fpt_rate = fpt_rate, treat_calf = treat_calf, treat_lac = treat_lac, simdays = 365*(simdays+2), init_r = init_r/100, init_p = init_p/100, init_cr = init_cr/100, init_cp = init_cp/100, treatment_length = treatment_length, density_lactating = density_lactating, density_dry = density_dry, density_calves = density_calf, treat_dry = treat_dry, vacc_protection = vacc_protection, vacc_shedding = vacc_shedding, vacc_duration = vacc_duration, p_retreatment = p_retreatment, longevity = longevity, run = 1)
            
            filter!(row->row.timestep <= 365*(simdays+2) && row.timestep !=0, simdata)
            simdata.total_moos = simdata.num_calves .+ simdata.num_weaned .+ simdata.num_dh .+ simdata.num_heifers .+ simdata.num_lactating .+ simdata.num_dry
            simdata.prev_r = round.(100 .* (simdata.pop_r .+ simdata.pop_car_r) ./ simdata.total_moos, digits =2)
            simdata.prev_active_amrsi = round.(100 .* (simdata.pop_r) ./ simdata.total_moos, digits =2)
            simdata.prev_carrier_amrsi = round.(100 .* (simdata.pop_car_r) ./ simdata.total_moos, digits =2)
            simdata.prev_p = round.(100 .* (simdata.pop_p .+ simdata.pop_car_p)./ simdata.total_moos, digits =2)
            simdata.prev_active_si = round.(100 .* (simdata.pop_p) ./ simdata.total_moos, digits =2)
            simdata.prev_carrier_si = round.(100 .* (simdata.pop_car_p) ./ simdata.total_moos, digits =2)
            simdata.prev_s = round.(100 .* (simdata.pop_s) ./ simdata.total_moos, digits =2)
            simdata.prev_rec = round.(100 .* (simdata.pop_rec_p + simdata.pop_rec_r) ./ simdata.total_moos, digits =2)
            @save "private/cachedrun_sim.jld2" simdata
            @info "Model ran!"
        end


       simtask = Threads.@spawn run_sim()

        start_timer(simtask)


       wait(simtask)

       println("Simulation has completed.")
     
      
        @load "private/cachedrun_sim.jld2"
        prev_data = DataFrames.select(simdata, :timestep, :prev_active_amrsi, :prev_active_si, :prev_carrier_si, :prev_carrier_amrsi, :prev_rec, :prev_s)
        demographic_data = DataFrames.select(simdata, :timestep, :num_calves, :num_weaned, :num_heifers, :num_dh, :num_lactating, :num_dry)
        @info "Loaded the new simdata"
        sim_trace = [
            scattergl(
                x = prev_data[:, "timestep"],
                y = prev_data[:, "prev_active_amrsi"],
                mode = "markers+lines",
                marker = attr(size = 2, color = "rgba(31, 119, 180, 1.0))"),
                name = "AMRSI"
            ),
            scattergl(
                x = prev_data[:, "timestep"],
                y = prev_data[:, "prev_active_si"],
                mode = "markers+lines",
                marker = attr(size = 2, color = "rgba(44, 160, 44, 1.0)"),
                name = "SI"
            ),
            scattergl(
                x = prev_data[:, "timestep"],
                y = prev_data[:, "prev_carrier_si"],
                mode = "markers+lines",
                marker = attr(size = 2, color = "rgba(214, 39, 40, 1.0)"),
                name = "Carrier (SI)"
            ),
            scattergl(
                x = prev_data[:, "timestep"],
                y = prev_data[:, "prev_carrier_amrsi"],
                mode = "markers+lines",
                marker = attr(size = 2, color = "rgba(148, 103, 189, 1.0)"),
                name = "Carrier (AMRSI)"
            )
        ] 

        rec_trace = [
            scattergl(
                x = prev_data[:, "timestep"],
                y = prev_data[:, "prev_s"],
                mode = "markers+lines",
                marker = attr(size = 2),
                name = "Susceptible"
            ),
            scattergl(
                x = prev_data[:, "timestep"],
                y = prev_data[:, "prev_rec"],
                mode = "markers+lines",
                marker = attr(size = 2),
                name = "Recovered"
            )
        ]

        demo_trace = [
            scattergl(
                x = demographic_data[:, "timestep"],
                y = demographic_data[:, "num_calves"],
                mode = "markers+lines",
                marker = attr(size = 2),
                name = "Calves"
            ),
            scattergl(
                x = demographic_data[:, "timestep"],
                y = demographic_data[:, "num_weaned"],
                mode = "markers+lines",
                marker = attr(size = 2),
                name = "Weaned"
            ),
            scattergl(
                x = demographic_data[:, "timestep"],
                y = demographic_data[:, "num_heifers"],
                mode = "markers+lines",
                marker = attr(size = 2),
                name = "Heifers"
            ),
            scattergl(
                x = demographic_data[:, "timestep"],
                y = demographic_data[:, "num_dh"],
                mode = "markers+lines",
                marker = attr(size = 2),
                name = "Pregnant heifers"
            ),
            scattergl(
                x = demographic_data[:, "timestep"],
                y = demographic_data[:, "num_lactating"],
                mode = "markers+lines",
                marker = attr(size = 2, color = "rgba(31, 119, 180, 1.0)"),
                name = "Lactating"
            ),
            scattergl(
                x = demographic_data[:, "timestep"],
                y = demographic_data[:, "num_dry"],
                mode = "markers+lines",
                marker = attr(size = 2),
                name = "Dry"
            )
        ] 

       
            max_prev_si = string(maximum(simdata.prev_p[simdata.timestep .> 731]),"%")
            max_prev_amrsi = string(maximum(simdata.prev_r[simdata.timestep .> 731]), "%")
            max_prev_active_amrsi = string(maximum(simdata.prev_active_amrsi[simdata.timestep .> 731]), "%")
            max_prev_active_si = string(maximum(simdata.prev_active_si[simdata.timestep .> 731]), "%")
            max_prev_car_amrsi = string(maximum(simdata.prev_carrier_amrsi[simdata.timestep .> 731]), "%")
            max_prev_car_si = string(maximum(simdata.prev_carrier_si[simdata.timestep .> 731]), "%")

   

            simtable = DataTable(DataFrames.select(simdata, :timestep => :Timestep, :pop_r=>"Active AMRSI", :pop_s=>"Susceptible", :pop_p=>"Active SI", :pop_car_p=>"Carrier SI", :pop_car_r=>"Carrier AMRSI", :num_calves=>"Calves", :num_weaned=>"Weaned", :num_heifers=>"Heifers", :num_dh=>"Pregnant heifers", :num_lactating=>"Lactating", :num_dry=>"Dry", :total_moos=>"Total stock", :prev_r=>"Point prev, AMRSI (inc. carriers, %)", :prev_p=>"Point prev. SI (inc. carriers, %)", :prev_s=>"Susceptible (%)", :prev_rec=>"Recovered (%)" ))

            function ls_plot(simdata, lifestage, label)
                scattergl(
                    x = simdata[:, "timestep"],
                    y = simdata[:, string(lifestage)],
                    mode = "markers",
                    marker = attr(size = 4),
                    name = label
                )
            end

                calf_trace = [
                    ls_plot(simdata, "inf_calves", "Active SI"),
                    ls_plot(simdata, "ri_calves", "Active AMRSI"),
                    ls_plot(simdata, "car_calves", "Carrier SI"),
                    ls_plot(simdata, "rc_calves", "Carrier AMRSI")
                ]

                weaned_trace = [
                    ls_plot(simdata, "inf_heifers", "Active SI"),
                    ls_plot(simdata, "ri_heifers", "Active AMRSI"),
                    ls_plot(simdata, "car_heifers", "Carrier SI"),
                    ls_plot(simdata, "rc_heifers", "Carrier AMRSI")
                ]

                heifer_trace = [
                    ls_plot(simdata, "inf_weaned", "Active SI"),
                    ls_plot(simdata, "ri_weaned", "Active AMRSI"),
                    ls_plot(simdata, "car_weaned", "Carrier SI"),
                    ls_plot(simdata, "rc_weaned", "Carrier AMRSI")
                ]
                dh_trace = [
                    ls_plot(simdata, "inf_dh", "Active SI"),
                    ls_plot(simdata, "ri_dh", "Active AMRSI"),
                    ls_plot(simdata, "car_dh", "Carrier SI"),
                    ls_plot(simdata, "rc_dh", "Carrier AMRSI")
                ]

                lactating_trace = [
                    ls_plot(simdata, "inf_lactating", "Active SI"),
                    ls_plot(simdata, "ri_lactating", "Active AMRSI"),
                    ls_plot(simdata, "car_lactating", "Carrier SI"),
                    ls_plot(simdata, "rc_lactating", "Carrier AMRSI")
                ]

                dry_trace = [
                    ls_plot(simdata, "inf_dry", "Active SI"),
                    ls_plot(simdata, "ri_dry", "Active AMRSI"),
                    ls_plot(simdata, "car_dry", "Carrier SI"),
                    ls_plot(simdata, "rc_dry", "Carrier AMRSI")
                ]

            single_sim_completed = true

            @info single_sim_completed

        end

        @onbutton lhs_draw begin
            n_samples = simno
            dims = 21
            plan = randomLHC(n_samples, dims)
            scale_range = [
                (optimal_stock_min, optimal_stock_max),
                (vacc_rate_min, vacc_rate_max),
                (pen_decon_min, pen_decon_max),
                (fpt_rate_min, fpt_rate_max),
                (treat_prob_min, treat_prob_max),
                (treat_calf_min, treat_calf_max),
                (treat_lac_min, treat_lac_max),
                (treat_dry_min, treat_dry_max),
                (treat_length_min, treat_length_max),
                (density_calf_min, density_calf_max),
                (density_lactating_min, density_lactating_max),
                (density_dry_min, density_dry_max),
                (vacc_protection_min, vacc_protection_max),
                (vacc_shedding_min, vacc_shedding_max),
                (vacc_duration_min, vacc_duration_max),
                (p_retreatment_min, p_retreatment_max),
                (longevity_min, longevity_max),
                (init_r_min/100, init_r_max/100),
                (init_cr_min/100, init_cr_max/100),
                (init_p_min/100, init_p_max/100),
                (init_cp_min/100, init_cp_max/100)
            ]
        

        sampler = scaleLHC(plan, scale_range)
        @info scale_range
        varnames = [:optimal_stock, :vacc_rate, :pen_decon, :fpt_rate, :treatment_prob, :treat_calf, :treat_lac, :treat_dry, :treatment_length, :density_calves, :density_lactating, :density_dry, :vacc_protection, :vacc_shedding, :vacc_duration, :p_retreatment, :longevity, :init_r, :init_cr, :init_p, :init_cp]

            lhs_data = DataFrame(sampler, varnames)
            lhs_data.optimal_stock .= round.(lhs_data.optimal_stock, digits = 0)
            lhs_data.pen_decon .= round.(lhs_data.pen_decon, digits = 0)
            lhs_data.simdays .= 365*(ensemble_length + 2)
            calving_systems = [1,2,3]
            prob_batch_used = 1.0 - prob_split - prob_batch
            probs = [prob_spring, prob_split, prob_batch_used]
            assigned_cs = sample(calving_systems, Weights(probs), n_samples; replace = true)
            lhs_data.calving_system = assigned_cs
            lhs_data.run = 1:nrow(lhs_data)
            lhs_data.treatment_length = round.(lhs_data.treatment_length)
            lhs_data.density_lactating = round.(lhs_data.density_lactating)
            lhs_data.density_dry = round.(lhs_data.density_dry)
            lhs_data.density_calves = round.(lhs_data.density_calves)
            

            @info lhs_data

            CSV.write("./private/lhs.csv", lhs_data)

    end


    function ensemble_traces(results, var)
        traces = GenericTrace[]
    
       for (index, df) in enumerate(results)
            trace = scattergl(
                x = df[!, :timestep],  
                y = df[!, var],
                mode = "lines+markers",
                marker = attr(size = 1, color = "gray"),
                showlegend=false
            )
            push!(traces, trace)
        end
    
        
        return traces
    end
    
    @onbutton run_ensemble begin

        lhs_data = DataFrame(CSV.File("./private/lhs.csv"))
        ensemble_tasks = Vector{Task}(undef, nrow(lhs_data))

        exec_order = [
        :optimal_stock, 
        :calving_system, 
        :pen_decon, 
        :treatment_prob, 
        :vacc_rate, 
        :fpt_rate, 
        :treat_calf, 
        :treat_lac, 
        :simdays, 
        :init_r, 
        :init_p, 
        :init_cr, 
        :init_cp, 
        :treatment_length, 
        :density_lactating, 
        :density_dry, 
        :density_calves,
        :treat_dry,
        :vacc_protection,
        :vacc_shedding,
        :vacc_duration,
        :p_retreatment,
        :longevity,
        :run]

        function exec_ensemble(lhs_data)
            @sync for i in 1:nrow(lhs_data)
                kwargs = Dict{Symbol, Any}(col => lhs_data[i, col] for col in exec_order)
            
                ensemble_tasks[i] = Threads.@spawn begin
                    @info "Running run_herd! for row $i on thread $(Threads.threadid()) with parameters:", kwargs
            
                    result = run_herd!(; kwargs...)
                            @info "Completed run_herd! for row $i on thread $(Threads.threadid())"
                    return result
                end
            end
            
            results = [fetch(task) for task in ensemble_tasks if istaskdone(task) && !istaskfailed(task)]

            return results
        
        end

        function ensemble_timer(task_to_monitor::Task)
            start_time = Base.time()
            ensemblestatus = "Ensemble running"
            try
                while !istaskdone(task_to_monitor)
                    elapsed_time = Base.time() - start_time
                    ensemblestatus = string("Running ensemble simulation, time elapsed: ", round(elapsed_time; digits=2), " sec")
                    sleep(1)
                end
            finally
                elapsed_time = Base.time() - start_time
                ensemblestatus = string("Ensemble completed in: ", round(elapsed_time; digits=2), " sec")
            end
            return ensemblestatus 
        end

        ensemble_task = Threads.@spawn results = exec_ensemble(lhs_data)
        ensemble_timer(ensemble_task)
        #wait(ensemble_task)
        fetch(ensemble_task)

        @info "Arrived here"

        for (index, df) in enumerate(results)
            df.total_moos = df.num_calves .+ df.num_weaned .+ df.num_dh .+ df.num_heifers .+ df.num_lactating .+ df.num_dry
            df.prev_r = round.(100 .* (df.pop_r .+ df.pop_car_r) ./ df.total_moos, digits =2)
            df.prev_active_amrsi = round.(100 .* (df.pop_r) ./ df.total_moos, digits =2)
            df.prev_carrier_amrsi = round.(100 .* (df.pop_car_r) ./ df.total_moos, digits =2)
            df.prev_p = round.(100 .* (df.pop_p .+ df.pop_car_p)./ df.total_moos, digits =2)
            df.prev_active_si = round.(100 .* (df.pop_p) ./ df.total_moos, digits =2)
            df.prev_carrier_si = round.(100 .* (df.pop_car_p) ./ df.total_moos, digits =2)
            df.prev_s = round.(100 .* (df.pop_s) ./ df.total_moos, digits =2)
            df.prev_rec = round.(100 .* (df.pop_rec_p + df.pop_rec_r) ./ df.total_moos, digits =2)
            select!(df, [:timestep, :prev_r, :total_moos, :prev_active_amrsi, :prev_carrier_amrsi, :prev_p, :prev_active_si, :prev_carrier_si, :prev_s, :prev_rec])
        end

        @save "private/results.jld2" results

        @info "Has now saved"


    end



    @onbutton gen_ensemble_plots begin

        @load "private/results.jld2"
        
        ensemble_amrsi_trace = ensemble_traces(results, :prev_active_amrsi)
        ensemble_amrsi_car_trace = ensemble_traces(results, :prev_carrier_amrsi)
        ensemble_si_trace = ensemble_traces(results, :prev_active_si)
        ensemble_si_car_trace = ensemble_traces(results, :prev_carrier_si)
        ensemble_susceptible_trace = ensemble_traces(results, :prev_s)
        ensemble_recovered_trace = ensemble_traces(results, :prev_rec)

         plotting_completed = true

    end

  

    @onbutton gen_farm_network begin
        @info in_deg
        make_static_graph(numfarms, numedges, in_deg, out_deg)

        @info "Network generated"
 
    end

    @onbutton gen_farm_parms begin
        n_samples = numfarms
        dims = 12
        plan = randomLHC(n_samples, dims)
        scale_range = [
            (vacc_rate_min_f, vacc_rate_max_f),
            (fpt_rate_min_f, fpt_rate_max_f),
            (treat_prob_min_f, treat_prob_max_f),
            (treat_calf_min, treat_calf_max),
            (treat_lac_min_f, treat_lac_max_f),
            (treat_dry_min_f, treat_dry_max_f),
            (p_retreatment_min_f, p_retreatment_max_f),
            (longevity_min_f, longevity_max_f),
            (init_r_min_f/100, init_r_max_f/100),
            (init_cr_min_f/100, init_cr_max_f/100),
            (init_p_min_f/100, init_p_max_f/100),
            (init_cp_min_f/100, init_cp_max_f/100)
        ]
   
    sampler = scaleLHC(plan, scale_range)
    varnames = [:vacc_rate, :fpt_rate, :treatment_prob, :treat_calf, :treat_lac, :treat_dry, :p_retreatment, :longevity, :init_r, :init_cr, :init_p, :init_cp]
    @info length(varnames)    
        lhs_data = DataFrame(sampler, varnames)
 
       

        herd_size_dist = Truncated(Normal(mean_farm, round(mean_farm/4)), smallest_farm, largest_farm)
        lhs_data.optimal_stock = round.(rand(herd_size_dist, numfarms))
        
       
        lhs_data.pen_decon = round.(rand(numfarms))
        @info "Got here"
        lhs_data.density_calf .= density_calf_f
        lhs_data.density_lactating .= density_lactating_f
        lhs_data.density_dry .= density_dry_f
        lhs_data.treatment_length .= treat_length_f
        lhs_data.vacc_duration .= vacc_duration_f
        lhs_data.vacc_protection .= vacc_protection_f
        lhs_data.vacc_shedding .= vacc_shedding_f
        calving_systems = [1,2,3]
        prob_batch_used = 1.0 - prob_split_f - prob_batch_f
        probs = [prob_spring, prob_split, prob_batch_used]
        assigned_cs = sample(calving_systems, Weights(probs), n_samples; replace = true)
        lhs_data.calving_system = assigned_cs
        function weighted_sample(choices, weights)
            return choices[sample(1:length(choices), Weights(weights))]
        end
        
        lhs_data.saleyard_1 = [weighted_sample([500, 501, 999], [0.45, 0.45, 0.1]) for _ in 1:nrow(lhs_data)]
        lhs_data.saleyard_2 = [weighted_sample([500, 501, 999], [0.05, 0.05, 0.85]) for _ in 1:nrow(lhs_data)]
        for i in 1:nrow(lhs_data)
            if lhs_data.saleyard_1[i] == 999 && lhs_data.saleyard_2[i] == 999
                lhs_data.saleyard_1[i] = weighted_sample([500, 501], [0.5, 0.5])  
            end
        end

      
        lhs_data.optimal_stock .= round.(lhs_data.optimal_stock, digits = 0)
        lhs_data.pen_decon .= round.(lhs_data.pen_decon, digits = 0)
        lhs_data.run = 1:nrow(lhs_data)
        lhs_data.treatment_length = round.(lhs_data.treatment_length)
        lhs_data.density_lactating = round.(lhs_data.density_lactating)
        lhs_data.density_dry = round.(lhs_data.density_dry)
        lhs_data.density_calves = round.(lhs_data.density_calf)
        lhs_data.farm = 1:nrow(lhs_data)
        farm_file = lhs_data
        

        @save "private/farm_file.jld2" farm_file

    end



    

    @onbutton run_between_herd begin    
        @info "Triggered run"
        function run_between()
            @load "private/farm_file.jld2"
            @load "data/static_network.jld2"
    
            farm_statuses, farm_movements = run_farm_model!(farm_file, numsims, 2555, farm_si, farm_amrsi, space);
    
            @info statuses
    
            @save "private/statuses.jld2" farm_statuses
            @save "private/movements.jld2" farm_movements
    
            @info "rundone"
        end

        function between_timer(task_to_monitor::Task)
            start_time = Base.time()
            betweenstatus = "Simulation running"
            try
                while !istaskdone(task_to_monitor)
                    elapsed_time = Base.time() - start_time
                    betweenstatus = string("Running simulation, time elapsed: ", round(elapsed_time; digits=2), " sec")
                    sleep(0.1)
                end
            finally
                elapsed_time = Base.time() - start_time
                betweenstatus = string("Run completed in: ", round(elapsed_time; digits=2), " sec")
            end
            return betweenstatus  
        end

        
        between_task = Threads.@spawn run_between()
        between_timer(between_task)
        wait(between_task)
        @info "Between herd simulation completed!"
    end

    function process_between_herd(df::DataFrame)
        filter!(row -> row.farm != 0, df)
        transform!(df, :date => (x -> x .- minimum(x)) => :timestep)
        df.timestep =  getfield.(df.timestep, :value)

       
        nfarms = length(unique(df.farm))
      
       
        result = combine(groupby(df, [:timestep, :status]), nrow => :count)
        transform!(result, :status => ByRow(status -> if status == 2 "amrsi" elseif status == 1 "si" else "uninf" end) => :status)
        result = unstack(result, :status, :count)
        result.nfarms .= nfarms
        transform!(result, [:amrsi, :nfarms] => ByRow((amrsi, nfarms) -> 100 * amrsi / nfarms) => :prev_amrsi)
        transform!(result, [:si, :nfarms] => ByRow((si, nfarms) -> 100 * si / nfarms) => :prev_si)
        transform!(result, [:uninf, :nfarms] => ByRow((uninf, nfarms) -> 100 * uninf / nfarms) => :prev_uninf)

        

        return result
    end

    function process_movements(df::DataFrame)
        filter!(row -> row.from != 0, df)
        transform!(df, :date => (x -> x .- minimum(x)) => :timestep)
        df.timestep =  getfield.(df.timestep, :value)
        movement_summary = groupby(df, [:from, :to])
        movement_summary = combine(movement_summary, nrow => :nmoved)

        return movement_summary

    end

    function get_node_attributes()
        @load "private/farm_file.jld2"
        @load "private/statuses.jld2"

        simend_status = filter(row -> row.date == maximum(farm_statuses[1].date), farm_statuses[1])
        simend_status = simend_status[:, [:farm, :status]]
        farm_info = farm_file[:, [:farm, :vacc_rate, :fpt_rate, :treatment_prob, :optimal_stock]]
        farm_info =  innerjoin(farm_info, simend_status, on = :farm)
        return farm_info
    end

    function create_post_network(df::DataFrame)
        nfarms = length(unique(df.from))
        g = Graphs.SimpleDiGraph(nfarms)

        node_info = get_node_attributes()
        node_info.node = node_info.farm

        saleyard_info = DataFrame(farm = [500,501], vacc_rate = 0, fpt_rate = 0, treatment_prob = 0, optimal_stock = 0, status = 999, node = [500,501])
        node_info = vcat(node_info, saleyard_info)

        function gen_label(row)
            return """
                Farm $(row.farm) 
                Farm size: $(Int(row.optimal_stock))
            """
        end
        info_map = Dict(row.node => [gen_label(row)] for row in eachrow(node_info))
        hover_texts = [info_map[node] for node in node_info.node]  # Extract hover texts using node identifiers


        unique_nodes = unique(vcat(df.from, df.to))
        node_index_map = Dict(node => idx for (idx, node) in enumerate(unique_nodes))
        node_nmoved = Dict(node => 0 for node in unique_nodes)
        for row in eachrow(df)
            from_index = node_index_map[row.from]
            to_index = node_index_map[row.to]
            Graphs.add_edge!(g, from_index, to_index)
            # Increment nmoved sums for both 'from' and 'to' nodes
            node_nmoved[row.from] += row.nmoved
            node_nmoved[row.to] += row.nmoved
        end

        for i in 1:nrow(df)
            from_index = node_index_map[df.from[i]]
            to_index = node_index_map[df.to[i]]
            Graphs.add_edge!(g, from_index, to_index)
        end
        pos = GraphPlot.spring_layout(g,
                    C=18.6,
                    MAXITER=100,
                    INITTEMP=2.0)


        node_sizes = [node_nmoved[node] for node in unique_nodes]
        max_size = maximum(node_sizes)
        sizes = map(s -> 5 + 15 * (s / max_size), node_sizes)*2  # Adjust size for visibility


        x_positions, y_positions = pos

        # Initialize arrays for edge coordinates
        edge_x, edge_y = Float64[], Float64[]

        color_map = [size(Graphs.neighbors(g, node))[1] for node in 1:nfarms]
        
        # Populate edge coordinates
        for edge in Graphs.edges(g)
            push!(edge_x, x_positions[Graphs.src(edge)])
            push!(edge_y, y_positions[Graphs.src(edge)])
            push!(edge_x, x_positions[Graphs.dst(edge)])
            push!(edge_y, y_positions[Graphs.dst(edge)])
            push!(edge_x, NaN)  # This ensures edges are broken between segments
            push!(edge_y, NaN)
        end
        
        # Create scatter plot for edges
        edge_trace = scattergl(
            x=edge_x, 
            y=edge_y, 
            mode="lines",
            #line=attr(color="black", width=1),
            line=attr(
                width=1,
             color="#888"
                     )
        )
        
        # Create scatter plot for nodes
        node_trace = scatter(
            x=x_positions, 
            y=y_positions, 
            mode="markers",
            marker=attr(
                showscale=true,
                colorscale=colors.imola,
                color=color_map,
                size=sizes,
                text = hover_texts,
                hoverinfo = "text",
                colorbar=attr(
                    thickness=15,
                    title="Number of trading connections",
                    xanchor="left",
                    titleside="right"
                            )
                            
                ))
        
              #  layout = PlotlyBase.layout(title="Dynamic network (post simulation)")


        plot =  [edge_trace, node_trace]
                    
            return plot


    end

    function create_pre_network(df::DataFrame)
        nfarms = length(unique(df.src))
        g = Graphs.SimpleDiGraph(nfarms)
        df = transform(groupby(df, :src), :dst => (x -> length(unique(x))) => :nmoved)

        unique_nodes = unique(vcat(df.src, df.dst))
        node_index_map = Dict(node => idx for (idx, node) in enumerate(unique_nodes))
        node_nmoved = Dict(node => 0 for node in unique_nodes)
        for row in eachrow(df)
            from_index = node_index_map[row.src]
            to_index = node_index_map[row.dst]
            Graphs.add_edge!(g, from_index, to_index)
            # Increment nmoved sums for both 'from' and 'to' nodes
            node_nmoved[row.src] += row.nmoved
            node_nmoved[row.dst] += row.nmoved
        end

        for i in 1:nrow(df)
            from_index = node_index_map[df.src[i]]
            to_index = node_index_map[df.dst[i]]
            Graphs.add_edge!(g, from_index, to_index)
        end

        pos = GraphPlot.spring_layout(g,
                    C=18.6,
                    MAXITER=100,
                    INITTEMP=2.0)


        node_sizes = [node_nmoved[node] for node in unique_nodes]
        max_size = maximum(node_sizes)
        sizes = map(s -> 5 + 15 * (s / max_size), node_sizes)*2  # Adjust size for visibility


        for i in 1:nrow(df)
            from_index = node_index_map[df.src[i]]
            to_index = node_index_map[df.dst[i]]
            Graphs.add_edge!(g, from_index, to_index)
        end
        pos = GraphPlot.spring_layout(g,
                    C=18.6,
                    MAXITER=100,
                    INITTEMP=2.0)


        x_positions, y_positions = pos

        # Initialize arrays for edge coordinates
        edge_x, edge_y = Float64[], Float64[]

        color_map = [size(Graphs.neighbors(g, node))[1] for node in 1:nfarms]
        
        # Populate edge coordinates
        for edge in Graphs.edges(g)
            push!(edge_x, x_positions[Graphs.src(edge)])
            push!(edge_y, y_positions[Graphs.src(edge)])
            push!(edge_x, x_positions[Graphs.dst(edge)])
            push!(edge_y, y_positions[Graphs.dst(edge)])
            push!(edge_x, NaN)  # This ensures edges are broken between segments
            push!(edge_y, NaN)
        end
        
        # Create scatter plot for edges
        edge_trace = scattergl(
            x=edge_x, 
            y=edge_y, 
            mode="lines",
            #line=attr(color="black", width=1),
            line=attr(
                width=1,
             color="#888"
                     )
        )
        
        # Create scatter plot for nodes
        node_trace = scatter(
            x=x_positions, 
            y=y_positions, 
            mode="markers",
            marker=attr(
                showscale=true,
                colorscale=colors.imola,
                color=color_map,
                size=sizes,
               # text = hover_texts,
             #   hoverinfo = "text",
                colorbar=attr(
                    thickness=15,
                    title="Number of trading connections",
                    xanchor="left",
                    titleside="right"
                            )
                            
                ))
        
        #layout = PlotlyBase.layout(title="Static network (pre simulation)")

        plot =  [edge_trace, node_trace]
                    
            return plot
    end

    @onbutton gen_spread_plots begin
        @info "Triggered event"
        @load "private/statuses.jld2"
        nits = length(farm_statuses)
        prevalences = Vector{DataFrame}(undef, nits)

        for status in eachindex(farm_statuses)
            prevalence = process_between_herd(farm_statuses[status])
            prevalences[status] = prevalence
        end
    
    between_prev_amrsi = ensemble_traces(prevalences, :prev_amrsi)
    between_prev_si = ensemble_traces(prevalences, :prev_si)
    between_prev_uninf = ensemble_traces(prevalences, :prev_uninf)
    
    within_prev_amrsi = ensemble_traces(farm_statuses, :pop_r)
    within_prev_si = ensemble_traces(farm_statuses, :pop_p)
    within_prev_uninf = ensemble_traces(farm_statuses, :pop_s)

    between_herd_plotting = true
    
    end

    @onbutton gen_spread_networks begin
        
        @load "private/movements.jld2"

        movements = Vector{DataFrame}()

        for movement in farm_movements
            #@show movement
            movement_df = process_movements(movement)
            push!(movements, movement_df)
        end

        post_sna = create_post_network(movements[1])

        @load "data/static_network.jld2"

        pre_sna = create_pre_network(space)

        generated_networks = true

    end


end

@page("/", "ui.jl")
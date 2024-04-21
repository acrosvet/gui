StippleUI.layout(view="hHh Lpr lff", 


"""
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.title = "SARMS-PDH"; // Sets the title

  var link = document.querySelector("link[rel~='icon']"); // Finds an existing favicon
  if (!link) {
    link = document.createElement('link');
    link.rel = 'icon';
    document.getElementsByTagName('head')[0].appendChild(link);
  }
  link.href = './flipped.png'; // Sets the favicon path
});
</script>
                    <div style="position: fixed; top: 0; left: 0; height: 100vh; width: 100vw; background-image: url('/sarms_pdh.png'); background-size: cover; background-repeat: no-repeat; background-position: center; z-index: 0;">
                      </div>
                      """,
                     [ 
        quasar(:header, style="background:primary;",  toolbar(
            [
            btn(; dense=true, flat=true, round=true, icon="menu", @click("left_drawer_open = !left_drawer_open")),
            toolbartitle(
                #icon = :sarms_pdh_icon,
               # "SARMS-PDH | <small>A Simulation of AMR <i>Salmonella enterica</i> for pastoral dairy herds</small>")
               """
               <div style="display: flex; align-items: center;">
                <b>SARMS-PDH</b> &nbsp;| <small>&nbsp;A simulation of antimicrobial resistant <i>Salmonella enterica</i> within and between pastoral dairy herds</small>
                <img src="flipped.png" style="height: 48px; margin-left: auto; margin-right: 10px;">
              </div>

               """)
       
        ],
        ),
        ),
                      drawer(bordered="", 
                      fieldname="left_drawer_open", 
                      side="left", 
                      #var":mini"="ministate", 
                      #var"@mouseover"="ministate = false", 
                      #var"@mouseout"="ministate = true", 
                      #var"mini-to-overlay"=true, width="170", 
                      #breakpoint=200,
                      width="200",
                      permanent="true",
                      overlay="true",
                             list(bordered=true, separator=true, overlap=false,
                                  [
                                   item(clickable="", vripple="", @click("selected_component = 'parmtab'"),
                                   
                                        [
                                         #itemsection(avatar=true, icon("tune")),
                                         itemsection("<b>SET PARAMETERS</b>  Single simulation<i>Within-herd model</i>")
                                        ]),
                                   item(clickable="", vripple="", @click("selected_component = 'restab'"), 
                                   #@iif("single_sim_completed == true"),
                                        [
                                         #itemsection(avatar=true, icon("show_chart")),
                                         itemsection("<b>VIEW RESULTS</b> Single simulation")
                                        ]),
                                    item(clickable="", vripple="", @click("selected_component = 'tabtab'"),  @iif("single_sim_completed == true"),
                                    [
                                        #itemsection(avatar=true, icon("table_restaurant")),
                                        itemsection("<b>RESULTS TABLE</b> Single simulation <i>Within-herd model</i>")
                                    ]),
                                    item(clickable="", vripple="", @click("selected_component = 'lifestage'"),  @iif("single_sim_completed == true"),
                                    [
                                       # itemsection(avatar=true, icon("grass")),
                                        itemsection("<b>LIFESTAGE PLOTS</b> Single simulation")
                                    ]),
                                    item(clickable="", vripple="", @click("selected_component = 'ensemble'"),
                                    [
                                        #itemsection(avatar=true, icon("tune")),
                                        itemsection("<b>SET PARAMETERS</b> Simulation ensemble <i>Within-herd model</i>")
                                    ]),
                                    item(clickable="", vripple="", @click("selected_component = 'ensembleresults'"),
                                    [
                                        #itemsection(avatar=true, icon("show_chart")),
                                        itemsection("<b>RESULTS</b> Simulation ensemble")
                                    ]),
                                    item(clickable="", vripple="", @click("selected_component = 'spreadparms'"),
                                    [
                                        #itemsection(avatar=true, icon("show_chart")),
                                        itemsection("<b>SET PARAMETERS</b>Between-herd spread model")
                                    ]),
                                    item(clickable="", vripple="", @click("selected_component = 'spreadresults'"),
                                    [
                                        #itemsection(avatar=true, icon("show_chart")),
                                        itemsection("<b>RESULTS</b>Between-herd spread model")
                                    ]),
                                    item(clickable="", vripple="", @click("selected_component = 'spreadnetworks'"),
                                    [
                                        #itemsection(avatar=true, icon("show_chart")),
                                        itemsection("<b>RESULTS</b>Movement networks")
                                    ])
                                  ]
                                 )),

                      page_container(  
                                     [
                                      Html.div(class="", 
                                     
                                      style="margin-left: 200px; background-color: red",
                                      @iif("selected_component == 'parmtab'"), 
                                               [
                                                row([
                                                    cell([
                                                    h6("Set model parameters")
                                                    ]),
                                                    cell([
                                                    Html.div(class = "", style="transform: scale(0.5);transform-origin: top left;",
                                                    [bignumber(number=:runstatus, title = "", style="font-weight: 600; text-transform: none; width: 100%;")]
                                                    )
                                                    ])
                                                    ]),
                                                
                                                    btn("Run simulation", color="primary", style="font-weight: 600; text-transform: none; width: 100%;", @click("run_model = true; alert('Simulation has started');")),
                                                row([
                                                    cell([
                                                        Stipple.select(:system; options=:systems, label = "Choose calving system")
                                                    ])
                                                ]),
                                                row([
                                                    cell([
                                                        numberfield(class = "q-my-md", "Set herd size:", :optimal_stock, hint = "Optimal number of lactating cows", min = 80, max = 1500)
                                                    ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Simulation length:", :simdays, hint = "Sim length in years", min = 1, max = 5)
                                                    ])
                                                ]),
                                                row([
                                                    cell([
                                                        Stipple.select( :pen_decon, label = "Pen decontamination:", options=:decon_options, hint = "Herd decontaminates calf pens")
                                                    ])
                                                ]),
                                                row([
                                                    cell([
                                                        numberfield(class = "q-my-md", "FPT probability:", :fpt_rate, hint = "Failure of passive transfer in calves", min = 0, max = 1, step = 0.01)
                                                    ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Vaccination probability:", :vacc_rate, hint = "Probability animals vaccinated", min = 0, max = 1, step = 0.01)
                                                    ])
                                                ]),
                                                row([
                                                    cell([
                                                        numberfield(class = "q-my-md", "Treatment probability:", :treatment_probability, hint = "Probability clinical cases receive antibiotics", min = 0, max = 1, step = 0.01)
                                                    ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Incidental treatment (lactating):", :treat_lac, hint = "Probability of incidental antibiotic treatment (per 305 days at risk)", min = 0, max = 1, step = 0.01)
                                                    ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Incidental treatment (calves):", :treat_calf, hint = "Probability of incidental antibiotic treatment (per 90 days at risk)", min = 0, max = 1, step = 0.01)
                                                    ])
                                                ]),
                                                 row([
                                                    cell([
                                                        numberfield(class = "q-my-md", "Initial prev: AMRSI (active)", :init_r, hint = "Active AMRSI at start", min = 0, max = 1, step = 0.01)
                                                    ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Initial prev: SI (active):", :init_p, hint = "Active SI at start", min = 0, max = 1, step = 0.01)
                                                    ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Initial prev: AMRSI (carrier):", :init_cr, hint = "Carrier AMRSI at start", min = 0, max = 1, step = 0.01)
                                                    ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Initial prev: SI (carrier)", :init_cp, hint = "Carrier SI", min = 0, max = 1, step = 0.01)
                                                    ])
                                                ]), 
                                                row([
                                                    cell([
                                                        numberfield(class = "q-my-md", "Treatment length", :treatment_length, hint = "Treatment length (in days)", min = 1, max = 5, step = 1)
                                                    ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Stocking density (calves)", :density_calf, hint = "Space per animal (square metres)", min = 1, max = 5, step = 1)
                                                    ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Stocking density (lactating)", :density_lactating, hint = "Space per animal (square meters)", min = 40, max = 80, step = 1)
                                                    ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Stocking density (lactating)", :density_dry, hint = "Space per animal (square meters)", min = 200, max = 300, step = 1)
                                                ])
                                                ]),
                                                row([
                                                    cell([
                                                        numberfield(class = "q-my-md", "Vacc protection", :vacc_protection, hint = "Protection from infection", min = 0, max = 1, step = 0.01)
                                                    ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Vacc shedding", :vacc_shedding, hint = "Reduction in shedding", min = 0, max = 1, step = 0.01)
                                                    ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Vacc duration", :vacc_duration, hint = "Reduction in duration of infectiousness", min = 0, max = 1, step = 0.01)
                                                    ])
                                                ]),
                                                row([
                                                    cell([
                                                        numberfield(class = "q-my-md", "Retreatment probability", :p_retreatment, hint = "Probability of receiving a second course of antibiotics", min = 0, max = 1, step = 0.01)
                                                ]),
                                                    cell([
                                                        numberfield(class = "q-my-md", "Cow longevity", :longevity, hint = "Target productive life (in years)", min = 0, max = 1, step = 0.01)

                                                    ])
                                                ])
                                                
                                               ]),
                                      Html.div(class="", @iif("selected_component == 'restab'"), 
                                      style="margin-left: 200px;",

                                               [ 
                                                h6("Simulation results"),
                                                h6("Run an example simulation to show results", @iif("single_sim_completed == false"), style="transform: scale(0.75);transform-origin: top left;"),
                                                row([ h6("Maximum prevalence", style="text-align: center; position: relative;", @iif("single_sim_completed == true"))]),
                                                row([
                                                    cell([bignumber(number=:max_prev_si, title = "All SI", @iif("single_sim_completed == true"),style="position: relative;")]),
                                                    cell([bignumber(number=:max_prev_amrsi, title = "All AMRSI", @iif("single_sim_completed == true"),style="position: relative;")]),     
                                                    cell([bignumber(number=:max_prev_active_si, title = "Active SI", @iif("single_sim_completed == true"),style="position: relative;")]),
                                                    cell([bignumber(number=:max_prev_active_amrsi, title = "Active AMRSI", @iif("single_sim_completed == true"),style="position: relative;")]),
                                                    cell([bignumber(number=:max_prev_car_si, title = "Carrier SI", @iif("single_sim_completed == true"),style="position: relative;")]),
                                                    cell([bignumber(number=:max_prev_car_amrsi, title = "Carrier AMRSI", @iif("single_sim_completed == true"),style="position: relative;")])
                                                ]),
                                                plot(:sim_trace, layout=:sim_layout, @iif("single_sim_completed == true")),
                                                plot(:rec_trace, layout=:rec_layout, @iif("single_sim_completed == true")),
                                                plot(:demo_trace, layout=:demo_layout, @iif("single_sim_completed == true"))
                                                

                                               ]),
                                        Html.div(class="", @iif("selected_component == 'tabtab'"), 
                                        style="margin-left: 200px;",

                                        [ 
                                        h6("Daily status"),
                                        h6("Run an example simulation to show results", @iif("single_sim_completed == false"), style="transform: scale(0.75);transform-origin: top left;"),

                                        row([
                                            GenieFramework.table(:simtable; dense=true, flat=true, style="height: 600px;", pagination=:simtable_pagination, @iif("single_sim_completed == true"))
                     
                                        ])
                                        ]),
                                        Html.div(class="", @iif("selected_component == 'lifestage'"), 
                                        style="margin-left: 200px;",

                                        [ 
                                        h6("Results by production lifestage"),
                                        h6("Run an example simulation to show results", @iif("single_sim_completed == false"), style="transform: scale(0.75);transform-origin: top left;"),

                                        row([
                                            cell([plot(:calf_trace, layout=:calf_layout, @iif("single_sim_completed == true"))]),
                                            cell([plot(:weaned_trace, layout=:weaned_layout, @iif("single_sim_completed == true"))])
                                        ]),
                                        row([
                                            cell([plot(:heifer_trace, layout=:heifer_layout, @iif("single_sim_completed == true"))]),
                                            cell([plot(:dh_trace, layout=:dh_layout, @iif("single_sim_completed == true"))])
                                        ]),
                                        row([
                                            cell([plot(:lactating_trace, layout=:lactating_layout, @iif("single_sim_completed == true"))]),
                                            cell([plot(:dry_trace, layout=:dry_layout, @iif("single_sim_completed == true"))])
                                        ])
                                              
                                        
                                        ]),
                                        Html.div(class="", @iif("selected_component == 'ensemble'"),
                                        style="margin-left: 200px;",

                                        [
                                            h6("Parameterise a simulation ensemble"),
                                            row([
                                            btn("Draw parameters", color="teal", style="font-weight: 600; text-transform: none; width: 100%;", @click("lhs_draw = true; alert('Parameters drawn and saved');")),

                                            ]),
                                            row([
                                                cell([], class="q-pa-xs q-ma-xs q-mb-md")
                                            ]),
                                            row([

                                            cell([btn("Run ensemble", color="red", style="font-weight: 600; text-transform: none; width: 100%;", @click("run_ensemble = true; alert('The ensemble simulation is starting...');"))]),
                                            cell([
                                                Html.div(class = "", style="transform: scale(0.5);transform-origin: top left;",
                                                [bignumber(number=:ensemblestatus, title = "", style="font-weight: 600; text-transform: none; width: 100%;")]
                                                )
                                                ])
                                            
                                            ]),
                                            row([
                                                cell([numberfield(class = "q-ma-xs", "Number of simulations", :simno, hint = "Expected duration: 10 sims = Cup of tea, 50 sims = Pot of tea, 100 sims = Samovar, 1000 sims = Plant, harvest and ferment Camellia sinensis (recommend using an HPC for large scale sims).", min = 10, max = 1500, step = 1)]),
                                                cell([numberfield(class = "q-ma-xs", "Simulation length", :ensemble_length, min = 1, max = 5)])

                                            ]),
                                            row([
                                                cell([numberfield(class = "q-ma-xs", "Herd size (min)", :prob_spring, hint = "Probablity drawn herd is spring", min = 0, max = 1, step = 0.01)]),
                                                cell([numberfield(class = "q-ma-xs", "Herd size (min)", :prob_split, hint = "Probablity drawn herd is split", min = 0, max = 1, step = 0.01)]),
                                                cell([numberfield(class = "q-ma-xs", "Herd size (min)", :prob_batch, hint = "Probablity drawn herd is batch", min = 0, max = 1, step = 0.01)])
                                            ]),

                                            row([
                                                        cell([numberfield(class = "q-ma-xs", "Herd size (min)", :optimal_stock_min, hint = " ", min = 80, max = 1500, step = 1)]),
                                                        cell([numberfield(class = "q-ma-xs", "Herd size (max)", :optimal_stock_max, hint = "", min = 80, max = 1500, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Vacc rate (min)", :vacc_rate_min, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Vacc rate (max)", :vacc_rate_max, hint = "", min = 0, max = 1, step = 0.01)])
                                                            ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Pen decon prob (min)", :pen_decon_min, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Pen decon prob (max)", :pen_decon_max, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "FPT prob (min)", :fpt_rate_min, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "FPT prob (max)", :fpt_rate_max, hint = "", min = 0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Treat prob (min)", :treat_prob_min, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Treat prob (max)", :treat_prob_max, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Incidental treatment lactating (min)", :treat_lac_min, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Incidental treatment lactating (max)", :treat_lac_max, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Incidental treatment calf (min)", :treat_calf_min, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Incidental treatment calf (max)", :treat_calf_max, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Incidental treatment dry (min)", :treat_dry_min, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Incidental treatment dry (max)", :treat_dry_max, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Treatment length (min)", :treat_length_min, hint = " ", min = 1, max = 5, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Treatment length (max)", :treat_length_max, hint = "", min =  1, max = 5, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Density calves (min)", :density_calf_min, hint = " ", min = 1, max = 10, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Density calves (max)", :density_calf_max, hint = "", min =  1, max = 10, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Density lactating (min)", :density_lactating_min, hint = " ", min = 10, max = 100, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Density lactating (max)", :density_lactating_max, hint = "", min =  10, max = 100, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Density dry (min)", :density_dry_min, hint = " ", min = 100, max = 500, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Density dry (max)", :density_dry_max, hint = "", min =  100, max = 500, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Vacc protection (min)", :vacc_protection_min, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Vacc protection (max)", :vacc_protection_max, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Vacc shedding (min)", :vacc_shedding_min, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Vacc shedding (max)", :vacc_shedding_max, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Vacc duration (min)", :vacc_duration_min, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Vacc duration (max)", :vacc_duration_max, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Retreatment probability (min)", :p_retreatment_min, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Retreatment probability", :p_retreatment_max, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Cow longevity (min)", :longevity_min, hint = " ", min = 5, max = 15, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Cow longevity (max)", :longevity_max, hint = "", min =  5, max = 15, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev active AMRSI (min)", :init_r_min, hint = " ", min = 0, max = 100, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev active AMRSI (max)", :init_r_max, hint = "", min =  0, max = 100, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev carrier AMRSI (min)", :init_cr_min, hint = " ", min = 5, max = 15, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev carrier AMRSI (max)", :init_cr_max, hint = "", min =  5, max = 15, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev active SI (min)", :init_p_min, hint = " ", min = 0, max = 100, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev active SI (max)", :init_p_max, hint = "", min =  0, max = 100, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev carrier SI (min)", :init_cp_min, hint = " ", min = 5, max = 15, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev carrier SI (max)", :init_cp_max, hint = "", min =  5, max = 15, step = 1)])
                                                        ])


                                                            ]
                                        ),
                                        Html.div(class="", @iif("selected_component == 'ensembleresults'"), 
                                        style="margin-left: 200px;",

                                        [ 
                                         h6("Ensemble results"),
                                         btn("Plot results", color="teal", style="font-weight: 600; text-transform: none; width: 100%;", @click(:gen_ensemble_plots)),
                                         h6("Click to refresh plots. These will take some time to generate, depending on how many simulations were run.", @iif("plotting_completed == false"), style="transform: scale(0.75);transform-origin: top left;"),
                                         #row([ h6("Maximum prevalence", style="text-align: center;")]),
                                         plot(:ensemble_amrsi_trace, layout=:ensemble_amrsi_layout, @iif("plotting_completed == true")),
                                         plot(:ensemble_amrsi_car_trace, layout=:ensemble_amrsi_car_layout, @iif("plotting_completed == true")),
                                         plot(:ensemble_si_trace, layout=:ensemble_si_layout, @iif("plotting_completed == true")),
                                         plot(:ensemble_si_car_trace, layout=:ensemble_si_car_layout, @iif("plotting_completed == true")),
                                         plot(:ensemble_recovered_trace, layout=:ensemble_recovered_layout, @iif("plotting_completed == true")),
                                         plot(:ensemble_susceptible_trace, layout=:ensemble_susceptible_layout, @iif("plotting_completed == true"))
                                         

                                        ]),
                                        Html.div(class="", @iif("selected_component == 'spreadparms'"), 
                                        style="margin-left: 200px;",
                                        [ 
                                            h6("Set parameters for a between-herd spread simulation"),
                                            row([h6("First specify a static trading network structure, next generate parameters for individual farms and then run the model.", style="transform: scale(0.85);transform-origin: top left;")]),
                                            
                                            row([
                                            cell([btn("Run between-herd model", color="red", style="font-weight: 600; text-transform: none; width: 100%;", @click("run_between_herd = true; alert('Running between herd model. This will take 15 to 30 min ...');"))]),
                                            cell([
                                                Html.div(class = "", style="transform: scale(0.5);transform-origin: top left;", 
                                                bignumber(number=:betweenstatus, title = "", style="font-weight: 600; text-transform: none; width: 100%;"))]),

                                           
                                            ]),
                                            row([h6("Farm network", style="transform: scale(0.85);transform-origin: top left;")]),
                                            btn("Generate farm network", color="orange", style="font-weight: 600; text-transform: none; width: 100%;", @click("gen_farm_network = true; alert('Generated trading network structure');")),
                                            row([
                                               
                                               # cell([numberfield(class = "q-ma-xs", "Number of simulations", :numsims, hint = "", min =  1, max = 15, step = 1)]),
                                                cell([numberfield(class = "q-ma-xs", "Number of farms", :numfarms, hint = "", min =  1, max = 15, step = 1)]),
                                                cell([numberfield(class = "q-ma-xs", "Initial prevalence SI", :farm_si, min = 1, max = 100, step = 1)]),
                                                cell([numberfield(class = "q-ma-xs", "Initial prevalence AMRSI", :farm_amrsi, min = 1, max = 100, step = 1)])
                                                ]),
                                            row([
                                                cell([numberfield(class = "q-ma-xs", "In-degree", :in_deg, hint = "", min =  1, max = 15, step = 1)]),
                                                cell([numberfield(class = "q-ma-xs", "Out-degree", :out_deg, hint = "", min =  1, max = 15, step = 1)]),
                                                cell([numberfield(class = "q-ma-xs", "Edges", :numedges, hint = "", min =  1, max = 15, step = 1)])
                                            ]),
                                            row([
                                                cell([numberfield(class = "q-ma-xs", "Probability spring calving", :prob_spring_f, hint = "Probablity drawn herd is spring", min = 0, max = 1, step = 0.01)]),
                                                cell([numberfield(class = "q-ma-xs", "Probability split calving", :prob_split_f, hint = "Probablity drawn herd is split", min = 0, max = 1, step = 0.01)]),
                                                cell([numberfield(class = "q-ma-xs", "Probability batch calving", :prob_batch_f, hint = "Probablity drawn herd is batch", min = 0, max = 1, step = 0.01)])
                                            ]),
                                            row([h6("Farm parameters", style="transform: scale(0.85);transform-origin: top left;")]),
                                            btn("Generate farm parameters", color="primary", style="font-weight: 600; text-transform: none; width: 100%;", @click("gen_farm_parms = true; alert('Farm parameters saved');")),

                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Minimum farm size", :smallest_farm, hint = "", min =  80, max = 100, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Average farm size", :mean_farm, hint = "", min =  250, max = 300, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Maximum farm size", :largest_farm, hint = "", min =  850, max = 1500, step = 1)])
                                                            ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Vacc rate (min)", :vacc_rate_min_f, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Vacc rate (max)", :vacc_rate_max_f, hint = "", min = 0, max = 1, step = 0.01)])
                                                            ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Pen decon prob", :pen_decon_f, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "FPT prob (min)", :fpt_rate_min_f, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "FPT prob (max)", :fpt_rate_max_f, hint = "", min = 0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Treat prob (min)", :treat_prob_min_f, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Treat prob (max)", :treat_prob_max_f, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Incidental treatment lactating (min)", :treat_lac_min_f, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Incidental treatment lactating (max)", :treat_lac_max_f, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Incidental treatment calf (min)", :treat_calf_min_f, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Incidental treatment calf (max)", :treat_calf_max_f, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Incidental treatment dry (min)", :treat_dry_min_f, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Incidental treatment dry (max)", :treat_dry_max_f, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Treatment length", :treat_length_f, hint = " ", min = 1, max = 5, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Density calves", :density_calf_f, hint = " ", min = 1, max = 10, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Density lactating", :density_lactating_f, hint = " ", min = 10, max = 100, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Density dry", :density_dry_f, hint = "", min =  100, max = 500, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Vacc protection", :vacc_protection_f, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Vacc shedding", :vacc_shedding_f, hint = "", min =  0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Vacc duration", :vacc_duration_f, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Retreatment probability (min)", :p_retreatment_min_f, hint = " ", min = 0, max = 1, step = 0.01)]),
                                                            cell([numberfield(class = "q-ma-xs", "Retreatment probability", :p_retreatment_max_f, hint = "", min =  0, max = 1, step = 0.01)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Cow longevity (min)", :longevity_min_f, hint = " ", min = 5, max = 15, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Cow longevity (max)", :longevity_max_f, hint = "", min =  5, max = 15, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev active AMRSI (min)", :init_r_min_f, hint = " ", min = 0, max = 100, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev active AMRSI (max)", :init_r_max_f, hint = "", min =  0, max = 100, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev carrier AMRSI (min)", :init_cr_min_f, hint = " ", min = 5, max = 15, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev carrier AMRSI (max)", :init_cr_max_f, hint = "", min =  5, max = 15, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev active SI (min)", :init_p_min_f, hint = " ", min = 0, max = 100, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev active SI (max)", :init_p_max_f, hint = "", min =  0, max = 100, step = 1)])
                                                        ]),
                                                        row([
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev carrier SI (min)", :init_cp_min_f, hint = " ", min = 5, max = 15, step = 1)]),
                                                            cell([numberfield(class = "q-ma-xs", "Initial prev carrier SI (max)", :init_cp_max_f, hint = "", min =  5, max = 15, step = 1)])
                                                        ])
                                         
                                        ]),
                                        Html.div(class="", @iif("selected_component == 'spreadresults'"), 
                                        style="margin-left: 200px;",

                                        [ 
                                         h6("Between herd spread results"),
                                         btn("Plot results", color="teal", style="font-weight: 600; text-transform: none; width: 100%;", @click(:gen_spread_plots)),
                                         h6("Click to refresh plots. These will take some time to generate, depending on how many simulations were run.", style="transform: scale(0.75);transform-origin: top left;"),
                                         h6("Between herd prevalence", @iif("between_herd_plotting == true")),
                                         plot(:between_prev_amrsi, layout=:between_prev_amrsi_layout, @iif("between_herd_plotting == true")),
                                         plot(:between_prev_si, layout=:between_prev_si_layout, @iif("between_herd_plotting == true")),
                                         plot(:between_prev_uninf, layout=:between_prev_uninf_layout, @iif("between_herd_plotting == true")),
                                         h6("Within-herd prevalence", @iif("between_herd_plotting == true")),
                                         plot(:within_prev_amrsi, layout=:within_prev_amrsi_layout, @iif("between_herd_plotting == true")),
                                         plot(:within_prev_si, layout=:within_prev_si_layout, @iif("between_herd_plotting == true")),
                                         plot(:within_prev_uninf, layout=:within_prev_uninf_layout, @iif("between_herd_plotting == true"))
                                        ]),
                                        Html.div(class="", @iif("selected_component == 'spreadnetworks'"), 
                                        style="margin-left: 200px;",

                                        [ 
                                         h6("Between herd spread networks", style="position: relative;"),
                                         btn("Generate networks", color="teal", style="font-weight: 600; text-transform: none; width: 100%;", @click(:gen_spread_networks)),
                                         h6("Click to generate networks. These will take some time to generate, depending on how many simulations were run.", style="transform: scale(0.75);transform-origin: top left;"),
                                         h6("Static network", @iif("generated_networks == true"),style="position: relative;"),
                                         plot(:pre_sna, layout = :post_sna_layout, @iif("generated_networks == true")),
                                         h6("Post-simulation network", @iif("generated_networks == true"), style="position: relative;"),
                                         plot(:post_sna, layout = :post_sna_layout, @iif("generated_networks == true"))
                                        ])
 
                                     ])
                     ]
                    )




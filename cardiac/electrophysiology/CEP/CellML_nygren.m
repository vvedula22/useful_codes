
%function [VOI, STATES, ALGEBRAIC, CONSTANTS] = mainFunction()
% This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
[VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
%end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".  
    algebraicVariableCount =49;
end
% There are a total of 29 entries in each of the rate and state variable arrays.
% There are a total of 51 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel()
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over 
    tspan = [0, 10];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

    % Solve model with ODE solver
    [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI);

    % Plot state variables against variable of integration
    [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends();
    figure();
    plot(VOI, STATES);
    xlabel(LEGEND_VOI);
    l = legend(LEGEND_STATES);
    set(l,'Interpreter','none');
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('time in component environment (second)');
    LEGEND_STATES(:,1) = strpad('V in component membrane (millivolt)');
    LEGEND_CONSTANTS(:,1) = strpad('R in component membrane (millijoule_per_mole_kelvin)');
    LEGEND_CONSTANTS(:,2) = strpad('T in component membrane (kelvin)');
    LEGEND_CONSTANTS(:,3) = strpad('F in component membrane (coulomb_per_mole)');
    LEGEND_CONSTANTS(:,4) = strpad('Cm in component membrane (nanoF)');
    LEGEND_ALGEBRAIC(:,27) = strpad('i_Na in component sodium_current (picoA)');
    LEGEND_ALGEBRAIC(:,29) = strpad('i_Ca_L in component L_type_Ca_channel (picoA)');
    LEGEND_ALGEBRAIC(:,31) = strpad('i_t in component Ca_independent_transient_outward_K_current (picoA)');
    LEGEND_ALGEBRAIC(:,32) = strpad('i_sus in component sustained_outward_K_current (picoA)');
    LEGEND_ALGEBRAIC(:,36) = strpad('i_K1 in component inward_rectifier (picoA)');
    LEGEND_ALGEBRAIC(:,35) = strpad('i_Kr in component delayed_rectifier_K_currents (picoA)');
    LEGEND_ALGEBRAIC(:,33) = strpad('i_Ks in component delayed_rectifier_K_currents (picoA)');
    LEGEND_ALGEBRAIC(:,37) = strpad('i_B_Na in component background_currents (picoA)');
    LEGEND_ALGEBRAIC(:,39) = strpad('i_B_Ca in component background_currents (picoA)');
    LEGEND_ALGEBRAIC(:,40) = strpad('i_NaK in component sodium_potassium_pump (picoA)');
    LEGEND_ALGEBRAIC(:,41) = strpad('i_CaP in component sarcolemmal_calcium_pump_current (picoA)');
    LEGEND_ALGEBRAIC(:,42) = strpad('i_NaCa in component Na_Ca_ion_exchanger_current (picoA)');
    LEGEND_ALGEBRAIC(:,1) = strpad('i_Stim in component membrane (picoA)');
    LEGEND_CONSTANTS(:,5) = strpad('stim_start in component membrane (second)');
    LEGEND_CONSTANTS(:,6) = strpad('stim_end in component membrane (second)');
    LEGEND_CONSTANTS(:,7) = strpad('stim_period in component membrane (second)');
    LEGEND_CONSTANTS(:,8) = strpad('stim_duration in component membrane (second)');
    LEGEND_CONSTANTS(:,9) = strpad('stim_amplitude in component membrane (picoA)');
    LEGEND_ALGEBRAIC(:,13) = strpad('E_Na in component sodium_current (millivolt)');
    LEGEND_CONSTANTS(:,10) = strpad('P_Na in component sodium_current (nanolitre_per_second)');
    LEGEND_STATES(:,2) = strpad('Na_c in component cleft_space_ion_concentrations (millimolar)');
    LEGEND_STATES(:,3) = strpad('Na_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES(:,4) = strpad('m in component sodium_current_m_gate (dimensionless)');
    LEGEND_STATES(:,5) = strpad('h1 in component sodium_current_h1_gate (dimensionless)');
    LEGEND_STATES(:,6) = strpad('h2 in component sodium_current_h2_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,2) = strpad('m_infinity in component sodium_current_m_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,14) = strpad('tau_m in component sodium_current_m_gate (second)');
    LEGEND_ALGEBRAIC(:,3) = strpad('h_infinity in component sodium_current_h1_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,15) = strpad('tau_h1 in component sodium_current_h1_gate (second)');
    LEGEND_ALGEBRAIC(:,16) = strpad('tau_h2 in component sodium_current_h2_gate (second)');
    LEGEND_CONSTANTS(:,11) = strpad('g_Ca_L in component L_type_Ca_channel (nanoS)');
    LEGEND_CONSTANTS(:,12) = strpad('E_Ca_app in component L_type_Ca_channel (millivolt)');
    LEGEND_ALGEBRAIC(:,28) = strpad('f_Ca in component L_type_Ca_channel (dimensionless)');
    LEGEND_CONSTANTS(:,13) = strpad('k_Ca in component L_type_Ca_channel (millimolar)');
    LEGEND_STATES(:,7) = strpad('Ca_d in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES(:,8) = strpad('d_L in component L_type_Ca_channel_d_L_gate (dimensionless)');
    LEGEND_STATES(:,9) = strpad('f_L_1 in component L_type_Ca_channel_f_L1_gate (dimensionless)');
    LEGEND_STATES(:,10) = strpad('f_L_2 in component L_type_Ca_channel_f_L2_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,4) = strpad('d_L_infinity in component L_type_Ca_channel_d_L_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,17) = strpad('tau_d_L in component L_type_Ca_channel_d_L_gate (second)');
    LEGEND_ALGEBRAIC(:,5) = strpad('f_L_infinity in component L_type_Ca_channel_f_L1_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,18) = strpad('tau_f_L1 in component L_type_Ca_channel_f_L1_gate (second)');
    LEGEND_ALGEBRAIC(:,19) = strpad('tau_f_L2 in component L_type_Ca_channel_f_L2_gate (second)');
    LEGEND_ALGEBRAIC(:,30) = strpad('E_K in component Ca_independent_transient_outward_K_current (millivolt)');
    LEGEND_CONSTANTS(:,14) = strpad('g_t in component Ca_independent_transient_outward_K_current (nanoS)');
    LEGEND_STATES(:,11) = strpad('K_c in component cleft_space_ion_concentrations (millimolar)');
    LEGEND_STATES(:,12) = strpad('K_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES(:,13) = strpad('r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless)');
    LEGEND_STATES(:,14) = strpad('s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,20) = strpad('tau_r in component Ca_independent_transient_outward_K_current_r_gate (second)');
    LEGEND_ALGEBRAIC(:,6) = strpad('r_infinity in component Ca_independent_transient_outward_K_current_r_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,21) = strpad('tau_s in component Ca_independent_transient_outward_K_current_s_gate (second)');
    LEGEND_ALGEBRAIC(:,7) = strpad('s_infinity in component Ca_independent_transient_outward_K_current_s_gate (dimensionless)');
    LEGEND_CONSTANTS(:,15) = strpad('g_sus in component sustained_outward_K_current (nanoS)');
    LEGEND_STATES(:,15) = strpad('r_sus in component sustained_outward_K_current_r_sus_gate (dimensionless)');
    LEGEND_STATES(:,16) = strpad('s_sus in component sustained_outward_K_current_s_sus_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,22) = strpad('tau_r_sus in component sustained_outward_K_current_r_sus_gate (second)');
    LEGEND_ALGEBRAIC(:,8) = strpad('r_sus_infinity in component sustained_outward_K_current_r_sus_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,23) = strpad('tau_s_sus in component sustained_outward_K_current_s_sus_gate (second)');
    LEGEND_ALGEBRAIC(:,9) = strpad('s_sus_infinity in component sustained_outward_K_current_s_sus_gate (dimensionless)');
    LEGEND_CONSTANTS(:,16) = strpad('g_Ks in component delayed_rectifier_K_currents (nanoS)');
    LEGEND_CONSTANTS(:,17) = strpad('g_Kr in component delayed_rectifier_K_currents (nanoS)');
    LEGEND_STATES(:,17) = strpad('n in component delayed_rectifier_K_currents_n_gate (dimensionless)');
    LEGEND_STATES(:,18) = strpad('p_a in component delayed_rectifier_K_currents_pa_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,34) = strpad('p_i in component delayed_rectifier_K_currents_pi_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,24) = strpad('tau_n in component delayed_rectifier_K_currents_n_gate (second)');
    LEGEND_ALGEBRAIC(:,10) = strpad('n_infinity in component delayed_rectifier_K_currents_n_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,25) = strpad('tau_p_a in component delayed_rectifier_K_currents_pa_gate (second)');
    LEGEND_ALGEBRAIC(:,11) = strpad('p_a_infinity in component delayed_rectifier_K_currents_pa_gate (dimensionless)');
    LEGEND_CONSTANTS(:,18) = strpad('g_K1 in component inward_rectifier (nanoS)');
    LEGEND_CONSTANTS(:,19) = strpad('g_B_Na in component background_currents (nanoS)');
    LEGEND_CONSTANTS(:,20) = strpad('g_B_Ca in component background_currents (nanoS)');
    LEGEND_ALGEBRAIC(:,38) = strpad('E_Ca in component background_currents (millivolt)');
    LEGEND_STATES(:,19) = strpad('Ca_c in component cleft_space_ion_concentrations (millimolar)');
    LEGEND_STATES(:,20) = strpad('Ca_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,21) = strpad('k_NaK_K in component sodium_potassium_pump (millimolar)');
    LEGEND_CONSTANTS(:,22) = strpad('k_NaK_Na in component sodium_potassium_pump (millimolar)');
    LEGEND_CONSTANTS(:,23) = strpad('i_NaK_max in component sodium_potassium_pump (picoA)');
    LEGEND_CONSTANTS(:,24) = strpad('i_CaP_max in component sarcolemmal_calcium_pump_current (picoA)');
    LEGEND_CONSTANTS(:,25) = strpad('k_CaP in component sarcolemmal_calcium_pump_current (millimolar)');
    LEGEND_CONSTANTS(:,26) = strpad('k_NaCa in component Na_Ca_ion_exchanger_current (picoA_per_millimolar_4)');
    LEGEND_CONSTANTS(:,27) = strpad('d_NaCa in component Na_Ca_ion_exchanger_current (per_millimolar_4)');
    LEGEND_CONSTANTS(:,28) = strpad('gamma in component Na_Ca_ion_exchanger_current (dimensionless)');
    LEGEND_CONSTANTS(:,29) = strpad('phi_Na_en in component intracellular_ion_concentrations (picoA)');
    LEGEND_CONSTANTS(:,30) = strpad('Vol_i in component intracellular_ion_concentrations (nanolitre)');
    LEGEND_CONSTANTS(:,50) = strpad('Vol_d in component intracellular_ion_concentrations (nanolitre)');
    LEGEND_ALGEBRAIC(:,43) = strpad('i_di in component intracellular_ion_concentrations (picoA)');
    LEGEND_CONSTANTS(:,31) = strpad('tau_di in component intracellular_ion_concentrations (second)');
    LEGEND_ALGEBRAIC(:,47) = strpad('i_up in component Ca_handling_by_the_SR (picoA)');
    LEGEND_ALGEBRAIC(:,49) = strpad('i_rel in component Ca_handling_by_the_SR (picoA)');
    LEGEND_ALGEBRAIC(:,44) = strpad('dOCdt in component intracellular_Ca_buffering (per_second)');
    LEGEND_ALGEBRAIC(:,45) = strpad('dOTCdt in component intracellular_Ca_buffering (per_second)');
    LEGEND_ALGEBRAIC(:,46) = strpad('dOTMgCdt in component intracellular_Ca_buffering (per_second)');
    LEGEND_STATES(:,21) = strpad('O_C in component intracellular_Ca_buffering (dimensionless)');
    LEGEND_STATES(:,22) = strpad('O_TC in component intracellular_Ca_buffering (dimensionless)');
    LEGEND_STATES(:,23) = strpad('O_TMgC in component intracellular_Ca_buffering (dimensionless)');
    LEGEND_STATES(:,24) = strpad('O_TMgMg in component intracellular_Ca_buffering (dimensionless)');
    LEGEND_CONSTANTS(:,32) = strpad('Mg_i in component intracellular_Ca_buffering (millimolar)');
    LEGEND_CONSTANTS(:,51) = strpad('Vol_c in component cleft_space_ion_concentrations (nanolitre)');
    LEGEND_CONSTANTS(:,33) = strpad('tau_Na in component cleft_space_ion_concentrations (second)');
    LEGEND_CONSTANTS(:,34) = strpad('tau_K in component cleft_space_ion_concentrations (second)');
    LEGEND_CONSTANTS(:,35) = strpad('tau_Ca in component cleft_space_ion_concentrations (second)');
    LEGEND_CONSTANTS(:,36) = strpad('Na_b in component cleft_space_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,37) = strpad('Ca_b in component cleft_space_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,38) = strpad('K_b in component cleft_space_ion_concentrations (millimolar)');
    LEGEND_ALGEBRAIC(:,48) = strpad('i_tr in component Ca_handling_by_the_SR (picoA)');
    LEGEND_CONSTANTS(:,39) = strpad('I_up_max in component Ca_handling_by_the_SR (picoA)');
    LEGEND_CONSTANTS(:,40) = strpad('k_cyca in component Ca_handling_by_the_SR (millimolar)');
    LEGEND_CONSTANTS(:,41) = strpad('k_srca in component Ca_handling_by_the_SR (millimolar)');
    LEGEND_CONSTANTS(:,42) = strpad('k_xcs in component Ca_handling_by_the_SR (dimensionless)');
    LEGEND_CONSTANTS(:,43) = strpad('alpha_rel in component Ca_handling_by_the_SR (picoA_per_millimolar)');
    LEGEND_STATES(:,25) = strpad('Ca_rel in component Ca_handling_by_the_SR (millimolar)');
    LEGEND_STATES(:,26) = strpad('Ca_up in component Ca_handling_by_the_SR (millimolar)');
    LEGEND_CONSTANTS(:,44) = strpad('Vol_up in component Ca_handling_by_the_SR (nanolitre)');
    LEGEND_CONSTANTS(:,45) = strpad('Vol_rel in component Ca_handling_by_the_SR (nanolitre)');
    LEGEND_ALGEBRAIC(:,12) = strpad('r_act in component Ca_handling_by_the_SR (per_second)');
    LEGEND_ALGEBRAIC(:,26) = strpad('r_inact in component Ca_handling_by_the_SR (per_second)');
    LEGEND_CONSTANTS(:,46) = strpad('r_recov in component Ca_handling_by_the_SR (per_second)');
    LEGEND_STATES(:,27) = strpad('O_Calse in component Ca_handling_by_the_SR (dimensionless)');
    LEGEND_STATES(:,28) = strpad('F1 in component Ca_handling_by_the_SR (dimensionless)');
    LEGEND_STATES(:,29) = strpad('F2 in component Ca_handling_by_the_SR (dimensionless)');
    LEGEND_CONSTANTS(:,47) = strpad('tau_tr in component Ca_handling_by_the_SR (second)');
    LEGEND_CONSTANTS(:,48) = strpad('k_rel_i in component Ca_handling_by_the_SR (millimolar)');
    LEGEND_CONSTANTS(:,49) = strpad('k_rel_d in component Ca_handling_by_the_SR (millimolar)');
    LEGEND_RATES(:,1) = strpad('d/dt V in component membrane (millivolt)');
    LEGEND_RATES(:,4) = strpad('d/dt m in component sodium_current_m_gate (dimensionless)');
    LEGEND_RATES(:,5) = strpad('d/dt h1 in component sodium_current_h1_gate (dimensionless)');
    LEGEND_RATES(:,6) = strpad('d/dt h2 in component sodium_current_h2_gate (dimensionless)');
    LEGEND_RATES(:,8) = strpad('d/dt d_L in component L_type_Ca_channel_d_L_gate (dimensionless)');
    LEGEND_RATES(:,9) = strpad('d/dt f_L_1 in component L_type_Ca_channel_f_L1_gate (dimensionless)');
    LEGEND_RATES(:,10) = strpad('d/dt f_L_2 in component L_type_Ca_channel_f_L2_gate (dimensionless)');
    LEGEND_RATES(:,13) = strpad('d/dt r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless)');
    LEGEND_RATES(:,14) = strpad('d/dt s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless)');
    LEGEND_RATES(:,15) = strpad('d/dt r_sus in component sustained_outward_K_current_r_sus_gate (dimensionless)');
    LEGEND_RATES(:,16) = strpad('d/dt s_sus in component sustained_outward_K_current_s_sus_gate (dimensionless)');
    LEGEND_RATES(:,17) = strpad('d/dt n in component delayed_rectifier_K_currents_n_gate (dimensionless)');
    LEGEND_RATES(:,18) = strpad('d/dt p_a in component delayed_rectifier_K_currents_pa_gate (dimensionless)');
    LEGEND_RATES(:,3) = strpad('d/dt Na_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,12) = strpad('d/dt K_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,20) = strpad('d/dt Ca_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,7) = strpad('d/dt Ca_d in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,21) = strpad('d/dt O_C in component intracellular_Ca_buffering (dimensionless)');
    LEGEND_RATES(:,22) = strpad('d/dt O_TC in component intracellular_Ca_buffering (dimensionless)');
    LEGEND_RATES(:,23) = strpad('d/dt O_TMgC in component intracellular_Ca_buffering (dimensionless)');
    LEGEND_RATES(:,24) = strpad('d/dt O_TMgMg in component intracellular_Ca_buffering (dimensionless)');
    LEGEND_RATES(:,2) = strpad('d/dt Na_c in component cleft_space_ion_concentrations (millimolar)');
    LEGEND_RATES(:,11) = strpad('d/dt K_c in component cleft_space_ion_concentrations (millimolar)');
    LEGEND_RATES(:,19) = strpad('d/dt Ca_c in component cleft_space_ion_concentrations (millimolar)');
    LEGEND_RATES(:,27) = strpad('d/dt O_Calse in component Ca_handling_by_the_SR (dimensionless)');
    LEGEND_RATES(:,25) = strpad('d/dt Ca_rel in component Ca_handling_by_the_SR (millimolar)');
    LEGEND_RATES(:,26) = strpad('d/dt Ca_up in component Ca_handling_by_the_SR (millimolar)');
    LEGEND_RATES(:,28) = strpad('d/dt F1 in component Ca_handling_by_the_SR (dimensionless)');
    LEGEND_RATES(:,29) = strpad('d/dt F2 in component Ca_handling_by_the_SR (dimensionless)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    STATES(:,1) = -74.2525;
    CONSTANTS(:,1) = 8314;
    CONSTANTS(:,2) = 306.15;
    CONSTANTS(:,3) = 96487;
    CONSTANTS(:,4) = 0.05;
    CONSTANTS(:,5) = 0.1;
    CONSTANTS(:,6) = 100000000;
    CONSTANTS(:,7) = 1;
    CONSTANTS(:,8) = 0.006;
    CONSTANTS(:,9) = -280;
    CONSTANTS(:,10) = 0.0016;
    STATES(:,2) = 130.011;
    STATES(:,3) = 8.5547;
    STATES(:,4) = 0.0032017;
    STATES(:,5) = 0.8814;
    STATES(:,6) = 0.8742;
    CONSTANTS(:,11) = 6.75;
    CONSTANTS(:,12) = 60;
    CONSTANTS(:,13) = 0.025;
    STATES(:,7) = 7.2495e-5;
    STATES(:,8) = 1.3005e-5;
    STATES(:,9) = 0.9986;
    STATES(:,10) = 0.9986;
    CONSTANTS(:,14) = 7.5;
    STATES(:,11) = 5.3581;
    STATES(:,12) = 129.435;
    STATES(:,13) = 0.0010678;
    STATES(:,14) = 0.949;
    CONSTANTS(:,15) = 2.75;
    STATES(:,15) = 0.00015949;
    STATES(:,16) = 0.9912;
    CONSTANTS(:,16) = 1;
    CONSTANTS(:,17) = 0.5;
    STATES(:,17) = 0.0048357;
    STATES(:,18) = 0.0001;
    CONSTANTS(:,18) = 3;
    CONSTANTS(:,19) = 0.060599;
    CONSTANTS(:,20) = 0.078681;
    STATES(:,19) = 1.8147;
    STATES(:,20) = 6.729e-5;
    CONSTANTS(:,21) = 1;
    CONSTANTS(:,22) = 11;
    CONSTANTS(:,23) = 70.8253;
    CONSTANTS(:,24) = 4;
    CONSTANTS(:,25) = 0.0002;
    CONSTANTS(:,26) = 0.0374842;
    CONSTANTS(:,27) = 0.0003;
    CONSTANTS(:,28) = 0.45;
    CONSTANTS(:,29) = -1.68;
    CONSTANTS(:,30) = 0.005884;
    CONSTANTS(:,31) = 0.01;
    STATES(:,21) = 0.0275;
    STATES(:,22) = 0.0133;
    STATES(:,23) = 0.1961;
    STATES(:,24) = 0.7094;
    CONSTANTS(:,32) = 2.5;
    CONSTANTS(:,33) = 14.3;
    CONSTANTS(:,34) = 10;
    CONSTANTS(:,35) = 24.7;
    CONSTANTS(:,36) = 130;
    CONSTANTS(:,37) = 1.8;
    CONSTANTS(:,38) = 5.4;
    CONSTANTS(:,39) = 2800;
    CONSTANTS(:,40) = 0.0003;
    CONSTANTS(:,41) = 0.5;
    CONSTANTS(:,42) = 0.4;
    CONSTANTS(:,43) = 200000;
    STATES(:,25) = 0.6465;
    STATES(:,26) = 0.6646;
    CONSTANTS(:,44) = 0.0003969;
    CONSTANTS(:,45) = 4.41e-5;
    CONSTANTS(:,46) = 0.815;
    STATES(:,27) = 0.4369;
    STATES(:,28) = 0.4284;
    STATES(:,29) = 0.0028;
    CONSTANTS(:,47) = 0.01;
    CONSTANTS(:,48) = 0.0003;
    CONSTANTS(:,49) = 0.003;
    CONSTANTS(:,50) =  0.0200000.*CONSTANTS(:,30);
    CONSTANTS(:,51) =  0.136000.*CONSTANTS(:,30);
    if (isempty(STATES)), warning('Initial values for states not set');, end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
        utilOnes = ones(statesRowCount, 1);
    end
    RATES(:,24) =  2000.00.*CONSTANTS(:,32).*((1.00000 - STATES(:,23)) - STATES(:,24)) -  666.000.*STATES(:,24);
    ALGEBRAIC(:,12) =  203.800.*(power(STATES(:,20)./(STATES(:,20)+CONSTANTS(:,48)), 4.00000)+power(STATES(:,7)./(STATES(:,7)+CONSTANTS(:,49)), 4.00000));
    RATES(:,28) =  CONSTANTS(:,46).*((1.00000 - STATES(:,28)) - STATES(:,29)) -  ALGEBRAIC(:,12).*STATES(:,28);
    ALGEBRAIC(:,2) = 1.00000./(1.00000+exp((STATES(:,1)+27.1200)./ - 8.21000));
    ALGEBRAIC(:,14) =  4.20000e-05.*exp( - power((STATES(:,1)+25.5700)./28.8000, 2.00000))+2.40000e-05;
    RATES(:,4) = (ALGEBRAIC(:,2) - STATES(:,4))./ALGEBRAIC(:,14);
    ALGEBRAIC(:,3) = 1.00000./(1.00000+exp((STATES(:,1)+63.6000)./5.30000));
    ALGEBRAIC(:,15) = 0.0300000./(1.00000+exp((STATES(:,1)+35.1000)./3.20000))+0.000300000;
    RATES(:,5) = (ALGEBRAIC(:,3) - STATES(:,5))./ALGEBRAIC(:,15);
    ALGEBRAIC(:,16) = 0.120000./(1.00000+exp((STATES(:,1)+35.1000)./3.20000))+0.00300000;
    RATES(:,6) = (ALGEBRAIC(:,3) - STATES(:,6))./ALGEBRAIC(:,16);
    ALGEBRAIC(:,4) = 1.00000./(1.00000+exp((STATES(:,1)+9.00000)./ - 5.80000));
    ALGEBRAIC(:,17) =  0.00270000.*exp( - power((STATES(:,1)+35.0000)./30.0000, 2.00000))+0.00200000;
    RATES(:,8) = (ALGEBRAIC(:,4) - STATES(:,8))./ALGEBRAIC(:,17);
    ALGEBRAIC(:,5) = 1.00000./(1.00000+exp((STATES(:,1)+27.4000)./7.10000));
    ALGEBRAIC(:,18) =  0.161000.*exp( - power((STATES(:,1)+40.0000)./14.4000, 2.00000))+0.0100000;
    RATES(:,9) = (ALGEBRAIC(:,5) - STATES(:,9))./ALGEBRAIC(:,18);
    ALGEBRAIC(:,19) =  1.33230.*exp( - power((STATES(:,1)+40.0000)./14.2000, 2.00000))+0.0626000;
    RATES(:,10) = (ALGEBRAIC(:,5) - STATES(:,10))./ALGEBRAIC(:,19);
    ALGEBRAIC(:,20) =  0.00350000.*exp( - power(STATES(:,1)./30.0000, 2.00000))+0.00150000;
    ALGEBRAIC(:,6) = 1.00000./(1.00000+exp((STATES(:,1) - 1.00000)./ - 11.0000));
    RATES(:,13) = (ALGEBRAIC(:,6) - STATES(:,13))./ALGEBRAIC(:,20);
    ALGEBRAIC(:,21) =  0.481200.*exp( - power((STATES(:,1)+52.4500)./14.9700, 2.00000))+0.0141400;
    ALGEBRAIC(:,7) = 1.00000./(1.00000+exp((STATES(:,1)+40.5000)./11.5000));
    RATES(:,14) = (ALGEBRAIC(:,7) - STATES(:,14))./ALGEBRAIC(:,21);
    ALGEBRAIC(:,22) = 0.00900000./(1.00000+exp((STATES(:,1)+5.00000)./12.0000))+0.000500000;
    ALGEBRAIC(:,8) = 1.00000./(1.00000+exp((STATES(:,1)+4.30000)./ - 8.00000));
    RATES(:,15) = (ALGEBRAIC(:,8) - STATES(:,15))./ALGEBRAIC(:,22);
    ALGEBRAIC(:,23) = 0.0470000./(1.00000+exp((STATES(:,1)+60.0000)./10.0000))+0.300000;
    ALGEBRAIC(:,9) = 0.400000./(1.00000+exp((STATES(:,1)+20.0000)./10.0000))+0.600000;
    RATES(:,16) = (ALGEBRAIC(:,9) - STATES(:,16))./ALGEBRAIC(:,23);
    ALGEBRAIC(:,24) = 0.700000+ 0.400000.*exp( - power((STATES(:,1) - 20.0000)./20.0000, 2.00000));
    ALGEBRAIC(:,10) = 1.00000./(1.00000+exp((STATES(:,1) - 19.9000)./ - 12.7000));
    RATES(:,17) = (ALGEBRAIC(:,10) - STATES(:,17))./ALGEBRAIC(:,24);
    ALGEBRAIC(:,25) = 0.0311800+ 0.217180.*exp( - power((STATES(:,1)+20.1376)./22.1996, 2.00000));
    ALGEBRAIC(:,11) = 1.00000./(1.00000+exp((STATES(:,1)+15.0000)./ - 6.00000));
    RATES(:,18) = (ALGEBRAIC(:,11) - STATES(:,18))./ALGEBRAIC(:,25);
    ALGEBRAIC(:,26) = 33.9600+ 339.600.*power(STATES(:,20)./(STATES(:,20)+CONSTANTS(:,48)), 4.00000);
    RATES(:,29) =  ALGEBRAIC(:,12).*STATES(:,28) -  ALGEBRAIC(:,26).*STATES(:,29);
    ALGEBRAIC(:,30) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*log(STATES(:,11)./STATES(:,12));
    ALGEBRAIC(:,31) =  CONSTANTS(:,14).*STATES(:,13).*STATES(:,14).*(STATES(:,1) - ALGEBRAIC(:,30));
    ALGEBRAIC(:,32) =  CONSTANTS(:,15).*STATES(:,15).*STATES(:,16).*(STATES(:,1) - ALGEBRAIC(:,30));
    ALGEBRAIC(:,36) = ( CONSTANTS(:,18).*power(STATES(:,11)./1.00000, 0.445700).*(STATES(:,1) - ALGEBRAIC(:,30)))./(1.00000+exp(( 1.50000.*((STATES(:,1) - ALGEBRAIC(:,30))+3.60000).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,34) = 1.00000./(1.00000+exp((STATES(:,1)+55.0000)./24.0000));
    ALGEBRAIC(:,35) =  CONSTANTS(:,17).*STATES(:,18).*ALGEBRAIC(:,34).*(STATES(:,1) - ALGEBRAIC(:,30));
    ALGEBRAIC(:,33) =  CONSTANTS(:,16).*STATES(:,17).*(STATES(:,1) - ALGEBRAIC(:,30));
    ALGEBRAIC(:,40) = ( (( (( CONSTANTS(:,23).*STATES(:,11))./(STATES(:,11)+CONSTANTS(:,21))).*power(STATES(:,3), 1.50000))./(power(STATES(:,3), 1.50000)+power(CONSTANTS(:,22), 1.50000))).*(STATES(:,1)+150.000))./(STATES(:,1)+200.000);
    RATES(:,12) =  - ((ALGEBRAIC(:,31)+ALGEBRAIC(:,32)+ALGEBRAIC(:,36)+ALGEBRAIC(:,35)+ALGEBRAIC(:,33)) -  2.00000.*ALGEBRAIC(:,40))./( CONSTANTS(:,30).*CONSTANTS(:,3));
    RATES(:,11) = (CONSTANTS(:,38) - STATES(:,11))./CONSTANTS(:,34)+((ALGEBRAIC(:,31)+ALGEBRAIC(:,32)+ALGEBRAIC(:,36)+ALGEBRAIC(:,35)+ALGEBRAIC(:,33)) -  2.00000.*ALGEBRAIC(:,40))./( CONSTANTS(:,51).*CONSTANTS(:,3));
    ALGEBRAIC(:,13) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*log(STATES(:,2)./STATES(:,3));
    ALGEBRAIC(:,27) = ( (( CONSTANTS(:,10).*power(STATES(:,4), 3.00000).*( 0.900000.*STATES(:,5)+ 0.100000.*STATES(:,6)).*STATES(:,2).*STATES(:,1).*power(CONSTANTS(:,3), 2.00000))./( CONSTANTS(:,1).*CONSTANTS(:,2))).*(exp(( (STATES(:,1) - ALGEBRAIC(:,13)).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2))) - 1.00000))./(exp(( STATES(:,1).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2))) - 1.00000);
    ALGEBRAIC(:,28) = STATES(:,7)./(STATES(:,7)+CONSTANTS(:,13));
    ALGEBRAIC(:,29) =  CONSTANTS(:,11).*STATES(:,8).*( ALGEBRAIC(:,28).*STATES(:,9)+ (1.00000 - ALGEBRAIC(:,28)).*STATES(:,10)).*(STATES(:,1) - CONSTANTS(:,12));
    ALGEBRAIC(:,37) =  CONSTANTS(:,19).*(STATES(:,1) - ALGEBRAIC(:,13));
    ALGEBRAIC(:,38) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./( 2.00000.*CONSTANTS(:,3))).*log(STATES(:,19)./STATES(:,20));
    ALGEBRAIC(:,39) =  CONSTANTS(:,20).*(STATES(:,1) - ALGEBRAIC(:,38));
    ALGEBRAIC(:,41) = ( CONSTANTS(:,24).*STATES(:,20))./(STATES(:,20)+CONSTANTS(:,25));
    ALGEBRAIC(:,42) = ( CONSTANTS(:,26).*( power(STATES(:,3), 3.00000).*STATES(:,19).*exp(( CONSTANTS(:,28).*CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2))) -  power(STATES(:,2), 3.00000).*STATES(:,20).*exp(( (CONSTANTS(:,28) - 1.00000).*STATES(:,1).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2)))))./(1.00000+ CONSTANTS(:,27).*( power(STATES(:,2), 3.00000).*STATES(:,20)+ power(STATES(:,3), 3.00000).*STATES(:,19)));
    ALGEBRAIC(:,1) = piecewise({VOI>=CONSTANTS(:,5)&VOI<=CONSTANTS(:,6)&(VOI - CONSTANTS(:,5)) -  floor((VOI - CONSTANTS(:,5))./CONSTANTS(:,7)).*CONSTANTS(:,7)<=CONSTANTS(:,8), CONSTANTS(:,9) }, 0.00000);
    RATES(:,1) =  ( - 1.00000./CONSTANTS(:,4)).*(ALGEBRAIC(:,1)+ALGEBRAIC(:,27)+ALGEBRAIC(:,29)+ALGEBRAIC(:,31)+ALGEBRAIC(:,32)+ALGEBRAIC(:,36)+ALGEBRAIC(:,35)+ALGEBRAIC(:,33)+ALGEBRAIC(:,37)+ALGEBRAIC(:,39)+ALGEBRAIC(:,40)+ALGEBRAIC(:,41)+ALGEBRAIC(:,42));
    RATES(:,3) =  - (ALGEBRAIC(:,27)+ALGEBRAIC(:,37)+ 3.00000.*ALGEBRAIC(:,40)+ 3.00000.*ALGEBRAIC(:,42)+CONSTANTS(:,29))./( CONSTANTS(:,30).*CONSTANTS(:,3));
    RATES(:,2) = (CONSTANTS(:,36) - STATES(:,2))./CONSTANTS(:,33)+(ALGEBRAIC(:,27)+ALGEBRAIC(:,37)+ 3.00000.*ALGEBRAIC(:,40)+ 3.00000.*ALGEBRAIC(:,42)+CONSTANTS(:,29))./( CONSTANTS(:,51).*CONSTANTS(:,3));
    RATES(:,19) = (CONSTANTS(:,37) - STATES(:,19))./CONSTANTS(:,35)+((ALGEBRAIC(:,29)+ALGEBRAIC(:,39)+ALGEBRAIC(:,41)) -  2.00000.*ALGEBRAIC(:,42))./( 2.00000.*CONSTANTS(:,51).*CONSTANTS(:,3));
    RATES(:,21) =  200000..*STATES(:,20).*(1.00000 - STATES(:,21)) -  476.000.*STATES(:,21);
    ALGEBRAIC(:,43) = ( (STATES(:,7) - STATES(:,20)).*2.00000.*CONSTANTS(:,3).*CONSTANTS(:,50))./CONSTANTS(:,31);
    RATES(:,7) =  - (ALGEBRAIC(:,29)+ALGEBRAIC(:,43))./( 2.00000.*CONSTANTS(:,50).*CONSTANTS(:,3));
    RATES(:,22) =  78400.0.*STATES(:,20).*(1.00000 - STATES(:,22)) -  392.000.*STATES(:,22);
    RATES(:,23) =  200000..*STATES(:,20).*((1.00000 - STATES(:,23)) - STATES(:,24)) -  6.60000.*STATES(:,23);
    ALGEBRAIC(:,47) = ( CONSTANTS(:,39).*(STATES(:,20)./CONSTANTS(:,40) - ( power(CONSTANTS(:,42), 2.00000).*STATES(:,26))./CONSTANTS(:,41)))./((STATES(:,20)+CONSTANTS(:,40))./CONSTANTS(:,40)+( CONSTANTS(:,42).*(STATES(:,26)+CONSTANTS(:,41)))./CONSTANTS(:,41));
    ALGEBRAIC(:,48) = ( (STATES(:,26) - STATES(:,25)).*2.00000.*CONSTANTS(:,3).*CONSTANTS(:,45))./CONSTANTS(:,47);
    RATES(:,26) = (ALGEBRAIC(:,47) - ALGEBRAIC(:,48))./( 2.00000.*CONSTANTS(:,44).*CONSTANTS(:,3));
    ALGEBRAIC(:,49) =  CONSTANTS(:,43).*power(STATES(:,29)./(STATES(:,29)+0.250000), 2.00000).*(STATES(:,25) - STATES(:,20));
    ALGEBRAIC(:,44) = RATES(:,21);
    ALGEBRAIC(:,45) = RATES(:,22);
    ALGEBRAIC(:,46) = RATES(:,23);
    RATES(:,20) =  - (((( - ALGEBRAIC(:,43)+ALGEBRAIC(:,39)+ALGEBRAIC(:,41)) -  2.00000.*ALGEBRAIC(:,42))+ALGEBRAIC(:,47)) - ALGEBRAIC(:,49))./( 2.00000.*CONSTANTS(:,30).*CONSTANTS(:,3)) - ( 0.0800000.*ALGEBRAIC(:,45)+ 0.160000.*ALGEBRAIC(:,46)+ 0.0450000.*ALGEBRAIC(:,44));
    RATES(:,27) =  480.000.*STATES(:,25).*(1.00000 - STATES(:,27)) -  400.000.*STATES(:,27);
    RATES(:,25) = (ALGEBRAIC(:,48) - ALGEBRAIC(:,49))./( 2.00000.*CONSTANTS(:,45).*CONSTANTS(:,3)) -  31.0000.*RATES(:,27);
   RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI)
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        utilOnes = ones(statesRowCount, 1);
    end
    ALGEBRAIC(:,12) =  203.800.*(power(STATES(:,20)./(STATES(:,20)+CONSTANTS(:,48)), 4.00000)+power(STATES(:,7)./(STATES(:,7)+CONSTANTS(:,49)), 4.00000));
    ALGEBRAIC(:,2) = 1.00000./(1.00000+exp((STATES(:,1)+27.1200)./ - 8.21000));
    ALGEBRAIC(:,14) =  4.20000e-05.*exp( - power((STATES(:,1)+25.5700)./28.8000, 2.00000))+2.40000e-05;
    ALGEBRAIC(:,3) = 1.00000./(1.00000+exp((STATES(:,1)+63.6000)./5.30000));
    ALGEBRAIC(:,15) = 0.0300000./(1.00000+exp((STATES(:,1)+35.1000)./3.20000))+0.000300000;
    ALGEBRAIC(:,16) = 0.120000./(1.00000+exp((STATES(:,1)+35.1000)./3.20000))+0.00300000;
    ALGEBRAIC(:,4) = 1.00000./(1.00000+exp((STATES(:,1)+9.00000)./ - 5.80000));
    ALGEBRAIC(:,17) =  0.00270000.*exp( - power((STATES(:,1)+35.0000)./30.0000, 2.00000))+0.00200000;
    ALGEBRAIC(:,5) = 1.00000./(1.00000+exp((STATES(:,1)+27.4000)./7.10000));
    ALGEBRAIC(:,18) =  0.161000.*exp( - power((STATES(:,1)+40.0000)./14.4000, 2.00000))+0.0100000;
    ALGEBRAIC(:,19) =  1.33230.*exp( - power((STATES(:,1)+40.0000)./14.2000, 2.00000))+0.0626000;
    ALGEBRAIC(:,20) =  0.00350000.*exp( - power(STATES(:,1)./30.0000, 2.00000))+0.00150000;
    ALGEBRAIC(:,6) = 1.00000./(1.00000+exp((STATES(:,1) - 1.00000)./ - 11.0000));
    ALGEBRAIC(:,21) =  0.481200.*exp( - power((STATES(:,1)+52.4500)./14.9700, 2.00000))+0.0141400;
    ALGEBRAIC(:,7) = 1.00000./(1.00000+exp((STATES(:,1)+40.5000)./11.5000));
    ALGEBRAIC(:,22) = 0.00900000./(1.00000+exp((STATES(:,1)+5.00000)./12.0000))+0.000500000;
    ALGEBRAIC(:,8) = 1.00000./(1.00000+exp((STATES(:,1)+4.30000)./ - 8.00000));
    ALGEBRAIC(:,23) = 0.0470000./(1.00000+exp((STATES(:,1)+60.0000)./10.0000))+0.300000;
    ALGEBRAIC(:,9) = 0.400000./(1.00000+exp((STATES(:,1)+20.0000)./10.0000))+0.600000;
    ALGEBRAIC(:,24) = 0.700000+ 0.400000.*exp( - power((STATES(:,1) - 20.0000)./20.0000, 2.00000));
    ALGEBRAIC(:,10) = 1.00000./(1.00000+exp((STATES(:,1) - 19.9000)./ - 12.7000));
    ALGEBRAIC(:,25) = 0.0311800+ 0.217180.*exp( - power((STATES(:,1)+20.1376)./22.1996, 2.00000));
    ALGEBRAIC(:,11) = 1.00000./(1.00000+exp((STATES(:,1)+15.0000)./ - 6.00000));
    ALGEBRAIC(:,26) = 33.9600+ 339.600.*power(STATES(:,20)./(STATES(:,20)+CONSTANTS(:,48)), 4.00000);
    ALGEBRAIC(:,30) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*log(STATES(:,11)./STATES(:,12));
    ALGEBRAIC(:,31) =  CONSTANTS(:,14).*STATES(:,13).*STATES(:,14).*(STATES(:,1) - ALGEBRAIC(:,30));
    ALGEBRAIC(:,32) =  CONSTANTS(:,15).*STATES(:,15).*STATES(:,16).*(STATES(:,1) - ALGEBRAIC(:,30));
    ALGEBRAIC(:,36) = ( CONSTANTS(:,18).*power(STATES(:,11)./1.00000, 0.445700).*(STATES(:,1) - ALGEBRAIC(:,30)))./(1.00000+exp(( 1.50000.*((STATES(:,1) - ALGEBRAIC(:,30))+3.60000).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,34) = 1.00000./(1.00000+exp((STATES(:,1)+55.0000)./24.0000));
    ALGEBRAIC(:,35) =  CONSTANTS(:,17).*STATES(:,18).*ALGEBRAIC(:,34).*(STATES(:,1) - ALGEBRAIC(:,30));
    ALGEBRAIC(:,33) =  CONSTANTS(:,16).*STATES(:,17).*(STATES(:,1) - ALGEBRAIC(:,30));
    ALGEBRAIC(:,40) = ( (( (( CONSTANTS(:,23).*STATES(:,11))./(STATES(:,11)+CONSTANTS(:,21))).*power(STATES(:,3), 1.50000))./(power(STATES(:,3), 1.50000)+power(CONSTANTS(:,22), 1.50000))).*(STATES(:,1)+150.000))./(STATES(:,1)+200.000);
    ALGEBRAIC(:,13) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*log(STATES(:,2)./STATES(:,3));
    ALGEBRAIC(:,27) = ( (( CONSTANTS(:,10).*power(STATES(:,4), 3.00000).*( 0.900000.*STATES(:,5)+ 0.100000.*STATES(:,6)).*STATES(:,2).*STATES(:,1).*power(CONSTANTS(:,3), 2.00000))./( CONSTANTS(:,1).*CONSTANTS(:,2))).*(exp(( (STATES(:,1) - ALGEBRAIC(:,13)).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2))) - 1.00000))./(exp(( STATES(:,1).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2))) - 1.00000);
    ALGEBRAIC(:,28) = STATES(:,7)./(STATES(:,7)+CONSTANTS(:,13));
    ALGEBRAIC(:,29) =  CONSTANTS(:,11).*STATES(:,8).*( ALGEBRAIC(:,28).*STATES(:,9)+ (1.00000 - ALGEBRAIC(:,28)).*STATES(:,10)).*(STATES(:,1) - CONSTANTS(:,12));
    ALGEBRAIC(:,37) =  CONSTANTS(:,19).*(STATES(:,1) - ALGEBRAIC(:,13));
    ALGEBRAIC(:,38) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./( 2.00000.*CONSTANTS(:,3))).*log(STATES(:,19)./STATES(:,20));
    ALGEBRAIC(:,39) =  CONSTANTS(:,20).*(STATES(:,1) - ALGEBRAIC(:,38));
    ALGEBRAIC(:,41) = ( CONSTANTS(:,24).*STATES(:,20))./(STATES(:,20)+CONSTANTS(:,25));
    ALGEBRAIC(:,42) = ( CONSTANTS(:,26).*( power(STATES(:,3), 3.00000).*STATES(:,19).*exp(( CONSTANTS(:,28).*CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2))) -  power(STATES(:,2), 3.00000).*STATES(:,20).*exp(( (CONSTANTS(:,28) - 1.00000).*STATES(:,1).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2)))))./(1.00000+ CONSTANTS(:,27).*( power(STATES(:,2), 3.00000).*STATES(:,20)+ power(STATES(:,3), 3.00000).*STATES(:,19)));
    ALGEBRAIC(:,1) = piecewise({VOI>=CONSTANTS(:,5)&VOI<=CONSTANTS(:,6)&(VOI - CONSTANTS(:,5)) -  floor((VOI - CONSTANTS(:,5))./CONSTANTS(:,7)).*CONSTANTS(:,7)<=CONSTANTS(:,8), CONSTANTS(:,9) }, 0.00000);
    ALGEBRAIC(:,43) = ( (STATES(:,7) - STATES(:,20)).*2.00000.*CONSTANTS(:,3).*CONSTANTS(:,50))./CONSTANTS(:,31);
    ALGEBRAIC(:,47) = ( CONSTANTS(:,39).*(STATES(:,20)./CONSTANTS(:,40) - ( power(CONSTANTS(:,42), 2.00000).*STATES(:,26))./CONSTANTS(:,41)))./((STATES(:,20)+CONSTANTS(:,40))./CONSTANTS(:,40)+( CONSTANTS(:,42).*(STATES(:,26)+CONSTANTS(:,41)))./CONSTANTS(:,41));
    ALGEBRAIC(:,48) = ( (STATES(:,26) - STATES(:,25)).*2.00000.*CONSTANTS(:,3).*CONSTANTS(:,45))./CONSTANTS(:,47);
    ALGEBRAIC(:,49) =  CONSTANTS(:,43).*power(STATES(:,29)./(STATES(:,29)+0.250000), 2.00000).*(STATES(:,25) - STATES(:,20));
    ALGEBRAIC(:,44) = RATES(:,21);
    ALGEBRAIC(:,45) = RATES(:,22);
    ALGEBRAIC(:,46) = RATES(:,23);
end

% Compute result of a piecewise function
function x = piecewise(cases, default)
    set = [0];
    for i = 1:2:length(cases)
        if (length(cases{i+1}) == 1)
            x(cases{i} & ~set,:) = cases{i+1};
        else
            x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
        end
        set = set | cases{i};
        if(set), break, end
    end
    if (length(default) == 1)
        x(~set,:) = default;
    else
        x(~set,:) = default(~set);
    end
end

% Pad out or shorten strings to a set length
function strout = strpad(strin)
    req_length = 160;
    insize = size(strin,2);
    if insize > req_length
        strout = strin(1:req_length);
    else
        strout = [strin, blanks(req_length - insize)];
    end
end

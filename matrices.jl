#Reactions
# "V1: ATP + Citrulline + Aspartate -> AMP + PPi + 2-(Nomega-L-arginine)succinate"
# "V2: 2-(Nomega-L-arginine)succinate -> Fumarate + Arginine"
# "V3: Arginine + H2O -> Ornithine + Urea"
# "V4: Carbamoyl phosphate + Ornithine -> Orthophosphate + Citrulline"
# "V5: 2 Arginine + 4 O2 + 3 NADPH + 3 H+ <-> 2 NO + 2 Citruilline + 3 NADP+ + 4 H2O"
# "B1: Carbamoyl phosphate (b) -> Carbamoyl phosphate"
# "B2: Aspartate (b) -> Aspartate"
# "B3: Fumarate -> Fumarate (b)"
# "B4: Urea -> Urea (b)"
# "B5: AMP -> AMP (b)"
# "B6: PPi -> PPi (b)"
# "B7: ATP (b) -> ATP"
# "B8: H2O (b) -> H2O"       # from V3
# "B9: H2O -> H2O (b)"       # from V5
# "B10: NO -> NO (b)"
# "B11: O2 (b) -> O2"
# "B12: NADPH (b) -> NADPH"
# "B13: H+ (b) -> H+"
# "B14: NADP+ -> NADP+ (b)"
# "B15: Pi -> Pi (b)"

include("VarnerFlux.jl")

Stoich=[
                -1.0    0     0      0     0      0    0    0     0    0    0   1.0  0    0     0   0     0     0    0     0     0
                -1.0    0     0     1.0   2.0   -2.0   0    0     0    0    0    0   0    0     0   0     0     0    0     0     0
                -1.0    0     0      0     0      0    0   1.0    0    0    0    0   0    0     0   0     0     0    0     0     0
                 1.0    0     0      0     0      0    0    0     0    0  -1.0   0   0    0     0   0     0     0    0     0     0
                 1.0    0     0      0     0      0    0    0     0    0    0  -1.0  0    0     0   0     0     0    0     0     0
                 1.0  -1.0    0      0     0      0    0    0     0    0    0    0   0    0     0   0     0     0    0     0     0
                  0    1.0    0      0     0      0    0    0   -1.0   0    0    0   0    0     0   0     0     0    0     0     0
                  0    1.0  -1.0     0   -2.0    2.0   0    0     0    0    0    0   0    0     0   0     0     0    0     0     0
                  0     0   -1.0     0    4.0   -4.0   0    0     0    0    0    0   0   1.0  -1.0  0     0     0    0     0     0
                  0     0    1.0   -1.0    0      0    0    0     0    0    0    0   0    0    0    0     0     0    0     0     0
                  0     0    1.0     0     0      0    0    0     0  -1.0   0    0   0    0    0    0     0     0    0     0     0
                  0     0     0      0   -3.0    3.0   0    0     0    0    0    0   0    0    0    0     0    1.0   0     0     0
                  0     0     0      0   -3.0    3.0   0    0     0    0    0    0   0    0    0    0     0     0   1.0    0     0
                  0     0     0      0   -4.0    4.0   0    0     0    0    0    0   0    0    0    0    1.0    0    0     0     0
                  0     0     0      0    2.0   -2.0   0    0     0    0    0    0   0    0    0  -1.0    0     0    0     0     0
                  0     0     0      0    3.0   -3.0   0    0     0    0    0    0   0    0    0    0     0     0    0   -1.0    0
                  0     0     0    -1.0    0      0   1.0   0     0    0    0    0   0    0    0    0     0     0    0     0     0
                  0     0     0     1.0    0      0    0    0     0    0    0    0   0    0    0    0     0     0    0     0   -1.0
];

balance=[
    10  6   4  10  0   10  4  6   0   5  1  21 0 0 0 21 1 0
    16  13  7  14  4   18  4  14  2   12 4  30 1 0 0 29 4 3
    5   3   1  5   0   4   0  4   0   2  2  7  0 0 1 7  1 0
    13  3   4  7   7   6   4  2   1   2  1  17 0 2 1 17 5 4
    3   0   0  1   2   0   0  0   0   0  0  3  0 0 0 3  1 1
    0   0   0  0   0   0   0  0   0   0  0  0  0 0 0 0  0 0
];

atom_blc = balance * Stoich

# Species Balance
species_bounds_array = [
         0.0    0.0	;       # ATP
         0.0    0.0	;       # Citrulline
         0.0    0.0	;       # Aspartate
         0.0    0.0	;       # AMP
         0.0    0.0	;       # PPi
         0.0    0.0	;       # Succinate
         0.0    0.0	;       # Fumarate
         0.0    0.0	;       # Arginine
         0.0    0.0	;       # H2O
         0.0    0.0	;       # Ornithine
         0.0    0.0	;       # Urea
         0.0    0.0	;       # NADPH
         0.0    0.0	;       # H+
         0.0    0.0	;       # O2
         0.0    0.0	;       # NO
         0.0    0.0	;       # NADP+
         0.0    0.0	;       # Carbamoyl phosphate
         0.0    0.0	;       # Pi
    ];

# Objective coefficient array
objective_coefficient_array = [
    	0.0    ;    # V1
        0.0    ;    # V2
        0.0    ;    # V3
        0.0    ;    # V4
        0.0    ;    # V5
        0.0    ;    # B1
        0.0    ;    # B2
        0.0    ;    # B3
        1.0    ;    # B4: maximizing Urea production
        0.0    ;    # B5
        0.0    ;    # B6
        0.0    ;    # B7
        0.0    ;    # B8
        0.0    ;    # B9
        0.0    ;    # B10
        0.0    ;    # B11
        0.0    ;    # B12
        0.0    ;    # B13
        0.0    ;    # B14
        0.0    ;    # B15
    ];

# Defining V_max for default bounds array
V_max1 = 7.308    ;  # mmol/g-DW*hr
V_max2 = 1.242    ;  # mmol/g-DW*hr
V_max3 = 8.964    ;  # mmol/g-DW*hr
V_max4 = 3.1716   ;  # mmol/g-DW*hr
V_max5 = 0.49319  ;  # mmol/g-DW*hr

# Defining Control  for default bounds array
Cntrl1 = 0.9225 * 0.9898  ;
Cntrl2 = 1                ;
Cntrl3 = 0.1418           ;
Cntrl4 = 0.7373           ;
Cntrl5 = 0.9865           ;

# Defining Default bounds array
default_bounds_array = [
       0	V_max1*Cntrl1	;	# V1
       0	V_max2*Cntrl2	;	# V2
       0	V_max3*Cntrl3	;	# V3
       0	V_max4*Cntrl4	;	# V4
       -V_max5*Cntrl5	V_max5*Cntrl5	;	# V5
       0	10     ;       # B1
       0	10     ;       # B2
       0	10     ;       # B3
       0	10     ;       # B4
       0	10     ;       # B5
       0	10     ;       # B6
       0	10     ;       # B7
       0	10     ;       # B8
     -10	10     ;       # B9
     -10	10     ;       # B10
     -10	10     ;       # B11
     -10	10     ;       # B12
     -10	10     ;       # B13
     -10	10     ;       # B14
       0	10     ;       # B15
   ];

(objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculte_optimal_flux_distribution(S,default_bounds_array,species_bounds_array,objective_coefficient_array,min_flag=false)

# [
#  [ MPIStateArray Name, Field Name, Maximum, Minimum, Mean, Standard Deviation ],
#  [         :                :          :        :      :          :           ],
# ]

parr = [
    ["state", :ρ, 12, 12, 12, 12],
    ["state", "ρu[1]", 12, 12, 12, 12],
    ["state", "ρu[2]", 12, 12, 12, 12],
    ["state", :ρθ, 12, 12, 10, 12],
]

#! format: off
polyorder1_rusanov = (
[
 [ "state",     "ρ",   9.93537454151033561089e-01,  1.00757504372666661929e+00,  9.99999999999999777955e-01,  1.02263433378261691005e-03 ],
 [ "state", "ρu[1]",  -2.52149323222216004137e-01,  5.89035169767154198617e-01,  1.59153776752221198798e-01,  1.08271170994963500389e-01 ],
 [ "state", "ρu[2]",  -2.83904021692772789986e-01,  3.49961171703760953733e-01, -7.54387871615413985182e-15,  8.42032266471293422416e-02 ],
 [ "state", "ρu[3]",  -1.63932566753800940118e-01,  1.23965597965165222338e-01, -1.13011813449026554679e-14,  3.89889444446707769676e-02 ],
 [ "state",    "ρθ",  -2.76247620490824052908e+00,  3.10549331157095087619e+00,  3.81639164714897560771e-17,  3.92101634706877422154e-01 ],
],
    parr,
)

polyorder1_roeflux = (
[
 [ "state",     "ρ",   9.97585895647369769179e-01,  1.00236035990324023892e+00,  9.99999999999999666933e-01,  5.00528832279226957441e-04 ],
 [ "state", "ρu[1]",  -3.84580135492388203167e-02,  3.77035244545484926615e-01,  1.59153776752221198798e-01,  5.73715254456575937669e-02 ],
 [ "state", "ρu[2]",  -2.43621178091514400954e-01,  2.47293992442084220595e-01, -7.68525867944624963002e-15,  6.82648495080709211136e-02 ],
 [ "state", "ρu[3]",  -2.09420608163400190360e-01,  2.17831173702847658014e-01, -1.13068191961995800909e-14,  5.54495074411652175139e-02 ],
 [ "state",    "ρθ",  -1.26980084672239956767e+00,  1.08666605998342746808e+00,  9.54097911787243901927e-17,  1.55448762892842345940e-01 ],
],
    parr,
)

polyorder1_rusanov_OI = (
[
 [ "state",     "ρ",   9.95958037734723578005e-01,  1.00495855045460968924e+00,  1.00000000000000466294e+00,  5.57944398922684704073e-04 ],
 [ "state", "ρu[1]",  -4.34993659111327646283e-02,  4.21153217882477615142e-01,  1.59153832831525815505e-01,  5.75287520097671686847e-02 ],
 [ "state", "ρu[2]",  -2.62346444903813313942e-01,  1.89066914031029154053e-01,  1.73563549078204747485e-13,  5.57705873789293080089e-02 ],
 [ "state", "ρu[3]",  -2.07290593236049203174e-01,  2.19447710086280806108e-01,  1.70645133825944280891e-13,  5.13975251802195609585e-02 ],
 [ "state",    "ρθ",  -1.12053362422573310475e+00,  1.22788292508450758156e+00, -5.60765425863910049406e-16,  1.84284312273684497407e-01 ],
],
    parr,
)

polyorder1_roeflux_OI = (
[
 [ "state",     "ρ",   9.96392461713754284958e-01,  1.00273194462648151948e+00,  1.00000000000000421885e+00,  4.75134340693243468340e-04 ],
 [ "state", "ρu[1]",  -8.54600194713670824331e-02,  4.02538964605315685574e-01,  1.59153832831525593461e-01,  5.89090417737306304424e-02 ],
 [ "state", "ρu[2]",  -3.21940872783143539060e-01,  2.64398003261648806284e-01,  1.73801174069904768761e-13,  6.08279803152186812620e-02 ],
 [ "state", "ρu[3]",  -2.23230857311768421392e-01,  2.39411196621541449980e-01,  1.70854328627343568357e-13,  5.50776248221017880602e-02 ],
 [ "state",    "ρθ",  -5.04681716215467313091e-01,  6.30768576983098072652e-01, -5.30182893472615246298e-16,  7.40683592050967315457e-02 ],
],
    parr,
)

polyorder4_rusanov_OI = (
[
 [ "state",     "ρ",   9.81081119687905900406e-01,  1.02546863200421922713e+00,  9.99981537278744170294e-01,  9.66849142176961053388e-04 ],
 [ "state", "ρu[1]",  -4.31180237030100166340e-01,  6.83271932198187825769e-01,  1.59197324162233116995e-01,  6.93659410491637984375e-02 ],
 [ "state", "ρu[2]",  -4.64680386632805531022e-01,  4.07526712445438032972e-01,  5.05220081165664268441e-04,  7.02023022289678111374e-02 ],
 [ "state", "ρu[3]",  -4.17302889932255671734e-01,  4.22237897112829496660e-01, -2.23936105152157209890e-04,  5.88107393397426930770e-02 ],
 [ "state",    "ρθ",  -5.30642712966413237154e+05,  3.05809315733645984437e+05, -4.26409605695937869996e+01,  1.38833821676920106256e+04 ],
],
    parr,
)

polyorder4_roeflux_OI = (
[
 [ "state",     "ρ",   9.94742609290946155909e-01,  1.00383107477793220852e+00,  9.99994986973244737172e-01,  3.83131689751475921091e-04 ],
 [ "state", "ρu[1]",  -1.43632325541406702385e-01,  5.51335819502245194634e-01,  1.59250404849606036484e-01,  5.28425027374051434204e-02 ],
 [ "state", "ρu[2]",  -2.53378560432881150266e-01,  3.01289778488659565348e-01,  2.01684398256973883116e-04,  6.56919957209414762112e-02 ],
 [ "state", "ρu[3]",  -2.50004518195408198533e-01,  2.78664912597020086871e-01, -1.46870695308643243668e-04,  4.39167259046288754876e-02 ],
 [ "state",    "ρθ",  -1.01085368604178116314e+00,  8.88723796619741435165e-01, -2.72714617833441435538e-04,  5.56687696608747506488e-02 ],
],
    parr,
)

#! format: on

refVals = (;
    polyorder1_rusanov,
    polyorder1_roeflux,
    polyorder1_rusanov_OI,
    polyorder1_roeflux_OI,
    polyorder4_rusanov_OI,
    polyorder4_roeflux_OI,
)

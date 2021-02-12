# [
#  [ MPIStateArray Name, Field Name, Maximum, Minimum, Mean, Standard Deviation ],
#  [         :                :          :        :      :          :           ],
# ]

parr = [
    ["state", :ρ, 12, 12, 12, 12],
    ["state", "ρu[1]", 12, 12, 12, 12],
    ["state", "ρu[2]", 12, 12, 0, 12],
    ["state", "ρu[3]", 12, 12, 0, 12],
    ["state", :ρθ, 12, 12, 0, 12],
]

#! format: off
first_order = (
[
 [ "state",     "ρ",   9.96392461713754284958e-01,  1.00273194462648151948e+00,  1.00000000000000421885e+00,  4.75134340693243468340e-04 ],
 [ "state", "ρu[1]",  -8.54600194713670824331e-02,  4.02538964605315685574e-01,  1.59153832831525593461e-01,  5.89090417737306304424e-02 ],
 [ "state", "ρu[2]",  -3.21940872783143539060e-01,  2.64398003261648806284e-01,  1.73801174069904768761e-13,  6.08279803152186812620e-02 ],
 [ "state", "ρu[3]",  -2.23230857311768421392e-01,  2.39411196621541449980e-01,  1.70854328627343568357e-13,  5.50776248221017880602e-02 ],
 [ "state",    "ρθ",  -5.04681716215467313091e-01,  6.30768576983098072652e-01, -5.30182893472615246298e-16,  7.40683592050967315457e-02 ],
],
    parr,
)

fourth_order = (
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

refVals = (; first_order, fourth_order)

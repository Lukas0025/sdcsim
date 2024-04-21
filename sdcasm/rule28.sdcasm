#
# RULE 110 cellular automaton implementation in DNA|SIMD
# @autor Lukáš Plevač <xpleva07@vutbr.cz>
# @date 11.21.2023
#

define:
    0 [AB][CDE]
    1 [ABCD]{E}

data:
    1001111010000

instructions: # O(6)
    {E*A*B*F*}      # mark 10
    {D*E*A*B*C*G*}  # mark 11
    {DEABCG}        # remove mark 11
    {A*B*} {C*D*E*} # write 0
    {EABF}          # remove mark 10
    {A*B*C*D*}      # write 1
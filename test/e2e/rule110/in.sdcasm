#
# RULE 110 cellular automaton implementation in DNA|SIMD
# @autor Lukáš Plevač <xpleva07@vutbr.cz>
# @date 11.21.2023
#

define:
    0 [ABC][DE]
    1 {A}[BCDE]

data: # O(6)
    1001111010

instructions:
    {D*E*A*F*}      # mark 01
    {D*E*A*B*C*G*}  # mark 11
    {DEABCG}        # remove mark 11
    {A*B*C*} {D*E*} # write 0
    {DEAF}          # remove mark 01
    {B*C*D*E*}      # write 1
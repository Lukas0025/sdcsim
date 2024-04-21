#
# Shift left register in DNA|SIMD
# @autor Lukáš Plevač <xpleva07@vutbr.cz>
# @date 11.28.2023
#

define:
    #   BIT     NOT selector
    0 [ABC][DE]{NSL}
    1 [ABCD]{ENSL}

data:
    #    Implicit zero
    00110[ABC][DE][NSL]
    01100[ABC][DE][NSL]
    11000[ABC][DE][NSL]
    10000[ABC][DE][NSL]
    10101[ABC][DE][NSL]
    01010[ABC][DE][NSL]
    11010[ABC][DE][NSL]
    01101[ABC][DE][NSL]
    00000[ABC][DE][NSL]
    11111[ABC][DE][NSL]

instructions: # O(33)
    {E*N*S*L*A*}      # unvrap DE from zero and bind here on tother site onvrap first base
    {D*E*N*S*L*G*}    # replace unvraper if on left site is 0
    {ENSLA}           # remove all existing unvrapers free when left 1
    {E*N*S*L*A*B*}    # unvrap 0 when is on right is 0 on free space and fully binde here if is 1 right 
    {ENSLAB}          # unvrap unwraper if is 1 on right
    {E*N*S*L*I*}      # bind not selector for 11 on free space (only free when 1 right and 1 left) (prevent E binding to)
    {N*S*L*A*B*C*F*}  # REMOVE {E*N*S*L*A*B*}
    {DENSLG}          # REMOVE {D*E*N*S*L*G*}
    {NSLABCF}         # REMOVE {N*S*L*A*B*C*F*} 
    {A*B*C*}          # write zero fisrt part back

    {E*N*S*L*J*}      # try fit in free space between cells (when 1 on right thare is no free base when 0 there is one free base)
    {D*E*N*S*L*G*}    # replace {E*N*S*L*J*} when 0 on rigth
    {DENSLG}          # REMOVE {D*E*N*S*L*G*} (now free only when have 0 on rigth)
    {D*E*N*S*L*A*B*}  # unvrap 0 when is on left is 0 on free space and fully binde here if is 0 right 
    {DENSLAB}         # unvrap {D*E*N*S*L*A*B*} when 1 on rigth

    {E*N*S*L*Y*}         # pad NOT selector to prevent bind not selector here (01)
    {B*C*D*G*} {DENSLAB} # unvrap {D*E*N*S*L*A*B*} using {B*C*D*G*} and remove 

    {ENSLI}          # remove temp not selector for 11
    {N*S*L*U*}       # bind not selector for 00

    {ENSLY} {ENSLJ}  # remove PAD {E*N*S*L*Y*} and PAD {E*N*S*L*J*}
    
    # here is problem with 000
    # {BCDG} remove unvraper {B*C*D*G*} not needed becasuse is replaced by 0

    {A*B*C*} {D*E*}    # write 0
    
    # selected not subprogram

    {G*D*E*N*}         # mark NOT 0 and NOT 1
    {ABCD}             # remove unwraped 1
    {GDEN}             # remove mark
    {C*D*E*N*}         # mark write 0
    {CDEN}             # remowe mark write 0 (is only posible when is unvraped for second part of zero)
    {A*B*C*D*}         # write 1
    {ABCD}             # remove not writed 1 (is unwraped by mark write 0)
    {N*S*L*G*}         # unwrap write 0 mark
    {CDEN}             # remove write 0 mark
    {NSLG}             # remove unwraper
    {A*B*C*} {D*E*}    # write 0
    {NSLU}             # remove all not selectors

    
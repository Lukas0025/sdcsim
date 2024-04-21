#
# Selected NOT in DNA|SIMD
# @autor Lukáš Plevač <xpleva07@vutbr.cz>
# @date 11.27.2023
#

define:
    #   BIT     NOT selector
    0 [ABC][DE][NSE].<U*>
    1 [ABCD]{E}[NSE].<U*>
    NO [ABC][DE]{NSE}
    Nl [ABCD]{E}{NSE}

data:
    0 NO 1 Nl

instructions: # O(12)
    {G*D*E*N*}         # mark NOT 0 and NOT 1
    {ABCD}             # remove unwraped 1
    {GDEN}             # remove mark
    {C*D*E*N*}         # mark write 0
    {CDEN}             # remowe mark write 0 (is only posible when is unvraped for second part of zero)
    {A*B*C*D*}         # write 1
    {ABCD}             # remove not writed 1 (is unwraped by mark write 0)
    {N*S*E*G*}         # unwrap write 0 mark
    {CDEN}             # remove write 0 mark
    {NSEG}             # remove unwraper
    {A*B*C*} {D*E*}    # write 0
    {N*S*E*U*}         # lock all NOT selectors
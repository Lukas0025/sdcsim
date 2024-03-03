#
# Increment instruction implementation in DNA|SIMD
# @autor Lukáš Plevač <xpleva07@vutbr.cz>
# @date 11.21.2023
#

define:
    0 [ABC][DE]
    1 [AB][CDE]

data: # O(6) + O(N)
    # {F} is TUE BINDING hold
    1011{F}

instructions:
    # mark last if is 1 if exist and replace this and open base C
    #
    #                                               /
    # - > - - > - - > - > - > - - > - >   - - - - /  
    # | | | | | | | | | | | | | | | | |   | | | | |  
    # - - - - - - - - - - - - - - - - - - - - - - -  
    # A B C D E A B C D E A B C D E A B C D E A B C  
    #
    # if next bit is 1 it again replace this and open base C this is chain reaction
    #
    {D*E*F*G*}   {D*E*A*B*C*H*}

    # remove all markers
    {DEFG}        {DEABCH}

    # set 1 end this unvrap last 0
    {C*D*E*}

    # shift 1 end to center of register cell if possible. This allow unvrap by 0
    {B*C*D*}

    # set 0
    {A*B*C*} {D*E*}

    # remove last 0
    {ABC}

    # set last 0 to 1
    {A*B*}
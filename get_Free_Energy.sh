#!/bin/bash
#Get F vs i_s; from final configurations m vs z in fort.1*
for fl in fort.1*; do awk 'NR==1{z1=$1;m1=$2;z_prev=$1;m_prev=$2;F=0}
    NR>1{
        dm=$2-m_prev;dz=$1-z_prev;F+=dz*( -.5*$2^2 + .25*$2^4 + .5*(dm/dz)^2);
        z_prev=$1;m_prev=$2}
    END{dm=m_prev-m1;F+=dz*( -.5*m1^2 + .25*m1^4 + .5*(dm/dz)^2);
    print FILENAME,F
        }' $fl;done | awk -F"fort." '{print $2}' | awk '{print $1-100,$2}'

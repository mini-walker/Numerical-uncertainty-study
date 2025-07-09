%****************************************************************************
%
%  PROGRAM: Data_Process_For_Plate
%
%  PURPOSE:  Process the data of flow over plate
%
%  Programer: Shanqin Jin
%
%  Time: Mar.25.2017
%
%  Location: MUN
%****************************************************************************

%Input variables: 1 r21----h2/h1
%                            2 r32----h3/h2
%                            3 eplsion21----Fai_2-Fai_1
%                            4 eplsion32----Fai_3-Fai_2
%                            5 P----the observed order of grid convergence

function F = P_Calculater_non_linear(p,r21,r32,eplsion21,eplsion32)

    s=sign(eplsion32/eplsion21);
    F=p-abs(log(abs(eplsion32/eplsion21))+log((r21^p-s)/(r32^p-s)))/log(r21);

end




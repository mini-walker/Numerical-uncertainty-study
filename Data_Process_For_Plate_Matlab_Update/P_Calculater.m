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

function [P] = P_Calculater(r21,r32,eplsion21,eplsion32)

     s=sign(eplsion32/eplsion21);

    %write(*,*) r21,r32,eplsion21,eplsion32
    p_ini=2.0;
    times=0;
    P=0;
    
    while(times<1000)

        p_mid=abs(log(abs(eplsion32/eplsion21))+log((r21^p_ini-s)/(r32^p_ini-s)))/log(r21);

        times=times+1;
   
        if(abs(p_mid-p_ini)<0.0001)
            P=p_mid;
            break;
        end
        
         p_ini=p_mid;
         
         continue;
   
    end 

    if(times==1000) 
         fprintf('%s\n', '10000 is not enough for convergence in Simple P Calculater!');  
     end

end


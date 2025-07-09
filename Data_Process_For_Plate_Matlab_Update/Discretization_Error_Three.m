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


%Input variables: 1 Fai_1----desired variable
%                            2 h_1----typical cell size
%                            3 Fai_0(x(1))----the estimated exact solution
%                            4 alpha(x(2))----coefficient
%                            5 P(x(3))----------the observed order of grid convergence

function F = Discretization_Error_Three(x,Fai_1,h_1,Fai_2,h_2,Fai_3,h_3)

    F(1)=Fai_1-x(1)-x(2)*h_1^x(3);
    F(2)=Fai_2-x(1)-x(2)*h_2^x(3);
    F(3)=Fai_3-x(1)-x(2)*h_3^x(3);

end


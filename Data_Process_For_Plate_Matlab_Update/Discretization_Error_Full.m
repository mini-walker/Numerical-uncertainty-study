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

function F = Discretization_Error_Full(x,n_g,w,Fai,h)

    w_i_Fai_i=0.0;
    w_i_h_i_P=0.0;
    w_i_h_i_2P=0.0;
    w_i_Fai_i_h_i_P=0.0;
    w_i_Fai_i_h_i_P_log_h_i=0.0;
    w_i_h_i_P_log_h_i=0.0;
    w_i_h_i_2P_log_h_i=0.0;
    for i=1:n_g
        w_i_Fai_i=w_i_Fai_i+w(i)*Fai(i);
        w_i_h_i_P=w_i_h_i_P+w(i)*h(i)^x(3);
        w_i_h_i_2P=w_i_h_i_2P+w(i)*h(i)^(2*x(3));
        w_i_Fai_i_h_i_P=w_i_Fai_i_h_i_P+w(i)*Fai(i)*h(i)^x(3);
        w_i_Fai_i_h_i_P_log_h_i=w_i_Fai_i_h_i_P_log_h_i+w(i)*Fai(i)*h(i)^x(3)*log(h(i));
        w_i_h_i_P_log_h_i=w_i_h_i_P_log_h_i+w(i)*h(i)^x(3)*log(h(i));
        w_i_h_i_2P_log_h_i=w_i_h_i_2P_log_h_i+w(i)*h(i)^(2*x(3))*log(h(i));
    end


    F(1)=w_i_Fai_i-x(2)*w_i_h_i_P-x(1);
    F(2)=(w_i_Fai_i_h_i_P-w_i_Fai_i*w_i_h_i_P)/(w_i_h_i_2P-w_i_h_i_P*w_i_h_i_P)-x(2);
    F(3)=w_i_Fai_i_h_i_P_log_h_i-x(1)*w_i_h_i_P_log_h_i-x(2)*w_i_h_i_2P_log_h_i;

end


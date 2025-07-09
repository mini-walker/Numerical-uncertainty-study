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

% Clear all the before data
clear all;
close all;
clc;

% save the running record
if exist('.\running record.txt') 
    delete('.\running record.txt');
end 

diary('.\running record.txt');
diary on;    % Start the recording

 % Body of Data_Process_For_Plate
 % Record the CPU time  
 fprintf('%s\n', '************************Calculation Starting***********************');  
start_time=clock;
    
% Assign the mesh ratio
 mesh_ratio=[1.000,1.231,1.455,1.600,2.000,2.462,2.909,3.200,4.000,4.923,5.818,6.400,8.000];
 h=mesh_ratio;
 % CaseI
 mesh_number=[294912,194688,139392,115200,73728,48672,34848,28800,18432,12168,8712,7200,4608];
 % CaseII
 %mesh_number=[393216,259584,185856,153600,98304,64896,46464,38400,24576,16224,11616,9600,6144];
 % CaseIII
 %mesh_number=[491520,324480,232320,192000,12880,81120,58080,48000,30720,20280,14520,12000,7680];
 mesh_ratio_char=sprintf('%0.3f',mesh_ratio);           %.n n is the number after point
        
% Check the time the default value is 200000
 Change_desired_time_or_not = input('Do you want to change the desired time ? (y/n):','s');
 disp(['You input is:',Change_desired_time_or_not]); 
  
 desired_time=100000;
 if(Change_desired_time_or_not=='y') 
      desired_time = input('The time you want to check is:');
      fprintf('%s %d\n', 'The current check time is :',desired_time);   
 else
      fprintf('%s %d\n', 'The current check time is :',desired_time);  
 end 
    
 % Transfer the desired time to character
 time_file_name=num2str(desired_time);
    
% Open the save file 
exist_or_not=exist('conv_data_cd.dat','file');
if(exist_or_not>0) 
    fid_conv_data_cd=fopen('conv_data_cd.dat','r');
    if(fid_conv_data_cd>0)
         fprintf('%s\n', 'conv_data_cd.dat is opened successfull !'); 
    end
else
     fprintf('%s\n', 'conv_data_cd.dat is not exist and you need to add it !');  
end

exist_or_not=exist('conv_surface.dat','file');
if(exist_or_not>0) 
    fid_conv_surface=fopen('conv_surface.dat','r');
    if(fid_conv_surface>0)
         fprintf('%s\n', 'conv_surface.dat is opened successfull !'); 
    end
else
     fprintf('%s\n', 'conv_surface.dat is not exist and you need to add it !');  
end
    
% Read the data in the file and just save the number in the temp file
fid_temp_data_file=fopen('Temp_Data.dat','w');
while ~feof(fid_conv_data_cd)
    
    tline=fgetl(fid_conv_data_cd);     % Read the file line by line
     
    %The first character of the line is number or not
    if double(tline(1))>=48&&double(tline(1))<=57  
       
        % if this is not the character line, just saving the data in the
        % temp data file
        fprintf(fid_temp_data_file,'%s\n\n', tline);  
        continue
    end
end

fclose(fid_conv_data_cd);
double Drag_Coefficient;
Drag_Coefficient=importdata('Temp_Data.dat');      %input the data to the workspace
%format long;

%*************************************************************************************************
% Calculate the P base on every vicnity three sets
Fai=Drag_Coefficient(:,2);
n_g=13;                                                    % Number of the grid
oscillatory_convergence_number=0;
for i=1:n_g
    if((i-1)>=1&&(i+1)<=n_g)
        N_1(i)=mesh_number(i-1);
        N_2(i)=mesh_number(i);
        N_3(i)=mesh_number(i+1);
        r_21(i)=h(i)/h(i-1);
        r_32(i)=h(i+1)/h(i);
        Fai_1(i)=Fai(i-1);
        Fai_2(i)=Fai(i);
        Fai_3(i)=Fai(i+1);
        Eplsion_32(i)=Fai_3(i)-Fai_2(i);
        Eplsion_21(i)=Fai_2(i)-Fai_1(i);
        
         if(abs(r_21(i)-r_32(i))<1.0e-8)
            P_Simple(i)=abs(log(abs(Eplsion_32(i)/Eplsion_21(i))))/log(r_21(i));
         else     
            %***********************Solve the non-linear equations*****************
            options=optimset;
            options.MaxFunEvals=20000;
            options.MaxIter=10000;
            options.Display='off';
            
            [p,fval]=fzero(@(p) P_Calculater_non_linear(p,r_21(i),r_32(i),Eplsion_21(i),Eplsion_32(i)),2.000,options);
            %[p,fval]=fsolve(@(p) P_Calculater_non_linear(p,r_21(i),r_32(i),Eplsion_21(i),Eplsion_32(i)),0.100,options);
            [P_Simple(i)]=p;
            fval;
        end 

        if((Eplsion_32/Eplsion_21)<0.0d0)
            oscillatory_convergence_number=oscillatory_convergence_number+1;
            fprintf('%s %8.3f %s %15.8f\n', 'Oscillatory convergence happened at :',h(i),'and p is',P_Simple(i));  
        end 
        
        Fai_ext_21(i)=(r_21(i)^P_Simple(i)*Fai_1(i)-Fai_2(i))/(r_21(i)^P_Simple(i)-1);
        Error_a_21(i)=abs((Fai_1(i)-Fai_2(i))/Fai_1(i));
        Error_ext_21(i)=abs((Fai_ext_21(i)-Fai_1(i))/Fai_ext_21(i));
        GCI_fine_21(i)=1.25*Error_ext_21(i)/(r_21(i)^P_Simple(i)-1);
    end 
end

% Save the simple P and relative data
fid_Simple_P=fopen('Simple_P.dat','w');
Title={'N_1';'N_2';'N_3';'r_21';'r_32';'Fai_1';'Fai_2';'Fai_3';'Eplsion_32';'Eplsion_21';'P_Simple';'Fai_ext_21';'Error_a_21';'Error_ext_21';'GCI_fine_21'};
for i=2:12
      Transfer_Data(1,i)=N_1(i);
      Transfer_Data(2,i)=N_2(i);
      Transfer_Data(3,i)=N_3(i);
      Transfer_Data(4,i)=r_21(i);
      Transfer_Data(5,i)=r_32(i);
      Transfer_Data(6,i)=Fai_1(i);
      Transfer_Data(7,i)=Fai_2(i);
      Transfer_Data(8,i)=Fai_3(i);
      Transfer_Data(9,i)=Eplsion_32(i);
      Transfer_Data(10,i)=Eplsion_21(i);
      Transfer_Data(11,i)=P_Simple(i);
      Transfer_Data(12,i)=Fai_ext_21(i);
      Transfer_Data(13,i)=Error_a_21(i);
      Transfer_Data(14,i)=Error_ext_21(i);
      Transfer_Data(15,i)=GCI_fine_21(i);
end
% size(Title);

% Output the data to the file
for i=1:15
      fprintf(fid_Simple_P,'%12s',Title{i});  % The cell should also be also in {}
      for j=2:12
          if(i<=3)
              if(j==12) 
                  fprintf(fid_Simple_P,'%15d\n',floor(Transfer_Data(i,j))); 
              else
                   fprintf(fid_Simple_P,'%15d',floor(Transfer_Data(i,j))); 
              end 
          else
              if(j==12) 
                  fprintf(fid_Simple_P,'%15.8f\n',Transfer_Data(i,j)); 
              else
                  fprintf(fid_Simple_P,'%15.8f',Transfer_Data(i,j)); 
              end
          end
      end 
end

% Output the data to the oscillatory convergence number
fprintf(fid_Simple_P,'%s %d\n\n', 'Note: the oscillatory convergence number is',oscillatory_convergence_number);
fprintf('%s\n', 'Simple_P.dat is :'); 
type('Simple_P.dat'); % display the data file
fclose(fid_Simple_P);
%%*******************************************************************************************************
   
% %*************************************************************************************************
low_P=-25;
upper_P=25;
% Solve the observed convergence order p based on complex algorithm      
 for i=1:13
    if((i-1)>=1&&(i+1)<=13)
        N_1(i)=mesh_number(i-1);
        N_2(i)=mesh_number(i);
        N_3(i)=mesh_number(i+1);
        r_21(i)=h(i)/h(i-1);
        r_32(i)=h(i+1)/h(i);
        Fai_1(i)=Fai(i-1);
        Fai_2(i)=Fai(i);
        Fai_3(i)=Fai(i+1);
        Eplsion_32(i)=Fai_3(i)-Fai_2(i);
        Eplsion_21(i)=Fai_2(i)-Fai_1(i);
        
        %Solve the non-linear equations
        options=optimset;
        options.MaxFunEvals=20000;
        options.MaxIter=10000;
        options.Display='off';
        
        % Loop the different initial value of P(-2.5---2.5)
        fval=1.0;
        for j=low_P:upper_P
            p_initial=j/10.0;
            before_sum_fval=sum(fval);
            x0=[(Fai_1(i)+Fai_2(i)+Fai_3(i))/3.0;(Fai_1(i)-(Fai_2(i)+Fai_2(i)+Fai_3(i))/3.0)/h(i)^p_initial;p_initial];

            [x,fval,exitflag]=fsolve(@(x) Discretization_Error_Three(x,Fai_1(i),h(i-1),Fai_2(i),h(i),Fai_3(i),h(i+1)),x0,options);

            if(j==low_P)          %initial x
                ture_x=x;
            end
            
            after_sum_fval=sum(fval);
            if(after_sum_fval<before_sum_fval)
                ture_x=x;
            end
        end

        Fai_0(i)=ture_x(1);
        alpha(i)=ture_x(2);
        P_Complex_Three(i)=ture_x(3);
        Delta_RE(i)=alpha(i)*h(i)^P_Complex_Three(i);
        Error_ext_21(i)=abs((Fai_0(i)-Fai_1(i))/Fai_0(i));
    end 
 end  
 
% Save the data obtained from the three non-linear equations
fid_Non_Linear_Three_P=fopen('P_Three_non_linear_equs.dat','w');
Title={'N_1';'N_2';'N_3';'r_21';'r_32';'Fai_1';'Fai_2';'Fai_3';'Fai_0';'Alpha';'P';'Delta_RE';'Error_a_21';'Error_ext_21'};
for i=2:12
      Transfer_Data(1,i)=N_1(i);
      Transfer_Data(2,i)=N_2(i);
      Transfer_Data(3,i)=N_3(i);
      Transfer_Data(4,i)=r_21(i);
      Transfer_Data(5,i)=r_32(i);
      Transfer_Data(6,i)=Fai_1(i);
      Transfer_Data(7,i)=Fai_2(i);
      Transfer_Data(8,i)=Fai_3(i);
      Transfer_Data(9,i)=Fai_0(i);
      Transfer_Data(10,i)=alpha(i);
      Transfer_Data(11,i)=P_Complex_Three(i);
      Transfer_Data(12,i)=Delta_RE(i);
      Transfer_Data(13,i)=Error_a_21(i);
      Transfer_Data(14,i)=Error_ext_21(i);
end
% size(Title);

% Output the data to the file
for i=1:14
      fprintf(fid_Non_Linear_Three_P,'%12s',Title{i});  % The cell should also be also in {}
      for j=2:12
          if(i<=3)
              if(j==12) 
                  fprintf(fid_Non_Linear_Three_P,'%15d\n',floor(Transfer_Data(i,j))); 
              else
                   fprintf(fid_Non_Linear_Three_P,'%15d',floor(Transfer_Data(i,j))); 
              end 
          else
              if(j==12) 
                  fprintf(fid_Non_Linear_Three_P,'%15.8f\n',Transfer_Data(i,j)); 
              else
                  fprintf(fid_Non_Linear_Three_P,'%15.8f',Transfer_Data(i,j)); 
              end
          end
      end 
end

% Output the data to the oscillatory convergence number
%fprintf(fid_Non_Linear_Three_P,'%s %d\n\n', 'Note: the oscillatory convergence number is',oscillatory_convergence_number);
fprintf('%s\n', 'P_Three_non_linear_equs.dat is :');
type('P_Three_non_linear_equs.dat') % display the data file
fclose(fid_Non_Linear_Three_P);

%*************************************************************************************************

% %*************************************************************************************************
% % Slove the P based on all the grids results in same set
 P_Calculater_Model={'P=1','P=2','P1=1 and P2=2','P_Calaulater'};
 P_Model=P_Calculater_Model{4};
%***********************************Wighted Procedure(P_1)**************************************
%Solve the non-linear equations
options=optimset;
options.MaxFunEvals=20000;
options.MaxIter=10000;
options.Display='off';

% Loop the different initial value of P(-2.5---2.5)
fval=1.0;

for i=1:2
    
    % Wighted Procedure(P(1))
    if(i==1)
        % Calculate the sum of one over h
        sum_one_over_h=0.0;
        for k=1:n_g
            sum_one_over_h=1.0/h(k)+sum_one_over_h;
        end

        % Calculate the weighted number
        for k=1:n_g
            w(k)=1.0/h(k)/sum_one_over_h;
            nw(k)=n_g*w(k);  
        end
        
    % No-Wighted Procedure(P(2))
    elseif(i==2)
        
         % Calculate the weighted number
         for k=1:n_g
               w(k)=1.0;
               nw(k)=1.0;  
         end
    end
    
    for j=low_P:upper_P
        p_initial=j/10.0;
        before_standard_deviation(i)=0.0;
        x0=[sum(Fai)/n_g;(Fai(1)-sum(Fai)/n_g)/h(1)^p_initial;p_initial];

        if(j==low_P)           % Initial standard deviation
            for k=1:n_g
                 before_standard_deviation(i)=before_standard_deviation(i)+nw(k)*(Fai(k)-(x0(1)+x0(2)*h(k)^x0(3)))^2;
            end
        else     
            for k=1:n_g
                 before_standard_deviation(i)=before_standard_deviation(i)+nw(k)*(Fai(k)-(x(1)+x(2)*h(k)^x(3)))^2;
            end
        end

        before_standard_deviation(i)=sqrt(before_standard_deviation(i)/(n_g-3.0));

        [x,fval,exitflag]=fsolve(@(x) Discretization_Error_Full(x,n_g,w,Fai,h),x0,options);

        if(j==low_P)          %initial x
            ture_x=x;
        end

        after_standard_deviation(i)=0.0;
        for k=1:n_g
              after_standard_deviation(i)=after_standard_deviation(i)+nw(k)*(Fai(k)-(x(1)+x(2)*h(k)^x(3)))^2;
        end

        after_standard_deviation(i)=sqrt(after_standard_deviation(i)/(n_g-3.0));
        if(after_standard_deviation(i)<before_standard_deviation(i))
            ture_x=x;
        end
    end

    Fai_0_Full(i)=ture_x(1);
    alpha_Full(i)=ture_x(2);
    P_Complex_Full(i)=ture_x(3);

end
%*************************************************************************************************

%***********************Define P Base on the non-linear equations************************
% Both P larger than 0, choose the minimun standard deviation
 if(P_Complex_Full(1)>0.0 && P_Complex_Full(2)>0.0)
     [min_standard_deviation,i]=min(after_standard_deviation);  % Find the minimum standard deviation
     Delta=min_standard_deviation; 
     Final_Fai_0=Fai_0_Full(i);
     Final_alpha=alpha_Full(i);
     Final_P=P_Complex_Full(i);
     
 % Both P less than 0, report this is an anomalou data
 elseif(P_Complex_Full(1)<0.0 && P_Complex_Full(2)<0.0)
     fprintf('%s\n', 'This is an anomalou data !');  
     Final_P=-2.0;
 % One P less than 0 and One P larger than 0, Choose the positive one
 elseif(P_Complex_Full(1)<0.0 && P_Complex_Full(2)>0.0)
     Delta=after_standard_deviation(2); 
     Final_Fai_0=Fai_0_Full(2);
     Final_alpha=alpha_Full(2);
     Final_P=P_Complex_Full(2);
 elseif(P_Complex_Full(1)>0.0 && P_Complex_Full(2)<0.0)
     Delta=after_standard_deviation(1); 
     Final_Fai_0=Fai_0_Full(1);
     Final_alpha=alpha_Full(1);
     Final_P=P_Complex_Full(1); 
 end
 %*************************************************************************************************
 
 %*************************Choose the alternative model or not*****************************
 Alternative_P=[1,2,3];
 % Does the P located in the reliabe zone
 %Final_P=0.1;
 Temp_Final_P=Final_P;                                    %Save the Temp Final when 2.0<P<2.1 is useful
 n_g=13;
 
 if(Final_P>2)

    for i=1:2  % Weighted or non-weighted

         % Wighted Procedure(P(1))
         if(i==1)
            % Calculate the sum of one over h
            sum_one_over_h=0.0;
            for k=1:n_g
                sum_one_over_h=1.0/h(k)+sum_one_over_h;
            end

            % Calculate the weighted number
            for k=1:n_g
                w(k)=1.0/h(k)/sum_one_over_h;
                nw(k)=n_g*w(k);  
            end

        % No-Wighted Procedure(P(2))
        elseif(i==2)

             % Calculate the weighted number
             for k=1:n_g
                   w(k)=1.0;
                   nw(k)=1.0;  
             end
         end
        
         for j=1:2              % Alternative mode 1 or mode 2
               
               w_i_Fai_i=0.0; 
               w_i_h_i_P=0.0;
               w_i_h_i_2P=0.0;
               w_i_Fai_i_h_i_P=0.0;
               for k=1:n_g
                     w_i_Fai_i=w_i_Fai_i+w(k)*Fai(k);
                     w_i_Fai_i_h_i_P=w_i_Fai_i_h_i_P+w(k)*Fai(k)*h(k)^Alternative_P(j);
                     w_i_h_i_P=w_i_h_i_P+w(k)*h(k)^Alternative_P(j);
                     w_i_h_i_2P=w_i_h_i_2P+w(k)*h(k)^(2*Alternative_P(j));
               end
               
               % Set the linear system
               coefficient_1=[1,w_i_h_i_P;w_i_h_i_P,w_i_h_i_2P];
               right_hand_1=[w_i_Fai_i;w_i_Fai_i_h_i_P];
               result=coefficient_1\right_hand_1;
               result_1(i,j,1)=result(1);
               result_1(i,j,2)=result(2);
                
               % Calculate the standard deviation
               alternative_standard_deviation(i,j)=0.0;
                for k=1:n_g
                      alternative_standard_deviation(i,j)=alternative_standard_deviation(i,j)+nw(k)*(Fai(k)-(result(1)+result(2)*h(k)^Alternative_P(j)))^2;
                end
                alternative_standard_deviation(i,j)=sqrt(alternative_standard_deviation(i,j)/(n_g-2.0));
                
         end
          
    end
    
    % Search the minimum of alternative standard deviation
     Delta=min(min(alternative_standard_deviation)); 
     [i,j]= find(alternative_standard_deviation==min(min(alternative_standard_deviation)));    
     Final_Fai_0=result_1(i,j,1);
     Final_alpha=result_1(i,j,2);
     Final_P=Alternative_P(j);
     P_Model=P_Calculater_Model{j};
         
 end
 
 % P is too small the value is not reliable
 if(Final_P<0.5&&Final_P>0.0)
     
     for i=1:2
     
         % Wighted Procedure(P(1))
         if(i==1)
            % Calculate the sum of one over h
            sum_one_over_h=0.0;
            for k=1:n_g
                sum_one_over_h=1.0/h(k)+sum_one_over_h;
            end

            % Calculate the weighted number
            for k=1:n_g
                w(k)=1.0/h(k)/sum_one_over_h;
                nw(k)=n_g*w(k);  
            end

        % No-Wighted Procedure(P(2))
        elseif(i==2)

             % Calculate the weighted number
             for k=1:n_g
                   w(k)=1.0;
                   nw(k)=1.0;  
             end
         end
        
         for j=1:3              % Alternative mode 1, mode 2 and mode 3
             
              if(j==3)          % Alternative mode 3
                   
                   w_i_h_i_1=0.0;
                   w_i_h_i_2=0.0;
                   w_i_h_i_3=0.0;
                   w_i_h_i_4=0.0;
                   w_i_Fai_i=0.0; 
                   w_i_Fai_i_h_i_1=0.0;
                   w_i_Fai_i_h_i_2=0.0;
                   for k=1:n_g
                       w_i_h_i_1=w_i_h_i_1+w(k)*h(k);
                       w_i_h_i_2=w_i_h_i_2+w(k)*h(k)^2;
                       w_i_h_i_3=w_i_h_i_3+w(k)*h(k)^3;
                       w_i_h_i_4=w_i_h_i_4+w(k)*h(k)^4;
                       w_i_Fai_i=w_i_Fai_i+w(k)*Fai(k);
                       w_i_Fai_i_h_i_1=w_i_Fai_i_h_i_1+w(k)*Fai(k)*h(k);
                       w_i_Fai_i_h_i_2=w_i_Fai_i_h_i_2+w(k)*Fai(k)*h(k)^2;
                   end

                   % Set the linear system
                   coefficient_2=[1,w_i_h_i_1,w_i_h_i_2;w_i_h_i_1,w_i_h_i_2,w_i_h_i_3;w_i_h_i_2,w_i_h_i_3,w_i_h_i_4];
                   right_hand_2=[w_i_Fai_i;w_i_Fai_i_h_i_1;w_i_Fai_i_h_i_2];
                   result_prime=coefficient_2\right_hand_2;
                   result_2(i,j,1)=result_prime(1);
                   result_2(i,j,2)=result_prime(2);
                   result_2(i,j,3)=result_prime(3);

                   % Calculate the standard deviation
                   alternative_standard_deviation(i,j)=0.0;
                    for k=1:n_g
                          alternative_standard_deviation(i,j)=alternative_standard_deviation(i,j)+nw(k)*(Fai(k)-(result_prime(1)+result_prime(2)*h(k)+result_prime(3)*h(k)^2))^2;
                    end
                    alternative_standard_deviation(i,j)=sqrt(alternative_standard_deviation(i,j)/(n_g-3.0));
             
              else      % Alternative mode 1, mode 2
                 
                   w_i_Fai_i=0.0; 
                   w_i_h_i_P=0.0;
                   w_i_h_i_2P=0.0;
                   w_i_Fai_i_h_i_P=0.0;
                   for k=1:n_g
                         w_i_Fai_i=w_i_Fai_i+w(k)*Fai(k);
                         w_i_Fai_i_h_i_P=w_i_Fai_i_h_i_P+w(k)*Fai(k)*h(k)^Alternative_P(j);
                         w_i_h_i_P=w_i_h_i_P+w(k)*h(k)^Alternative_P(j);
                         w_i_h_i_2P=w_i_h_i_2P+w(k)*h(k)^(2*Alternative_P(j));
                   end

                   % Set the linear system
                   coefficient_1=[1,w_i_h_i_P;w_i_h_i_P,w_i_h_i_2P];
                   right_hand_1=[w_i_Fai_i;w_i_Fai_i_h_i_P];
                   result=coefficient_1\right_hand_1;
                   result_1(i,j,1)=result(1);
                   result_1(i,j,2)=result(2);

                   % Calculate the standard deviation
                   alternative_standard_deviation(i,j)=0.0;
                    for k=1:n_g
                          alternative_standard_deviation(i,j)=alternative_standard_deviation(i,j)+nw(k)*(Fai(k)-(result(1)+result(2)*h(k)^Alternative_P(j)))^2;
                    end
                    alternative_standard_deviation(i,j)=sqrt(alternative_standard_deviation(i,j)/(n_g-2.0));

              end
         end
     end

    % Search the minimum of alternative standard deviation
     Delta=min(min(alternative_standard_deviation)); 
     [i,j]= find(alternative_standard_deviation==min(min(alternative_standard_deviation))); 
     if(j==3)
         Final_Fai_0=result_2(i,j,1);
         Final_alpha_1=result_2(i,j,2);
         Final_alpha_2=result_2(i,j,3);
         Final_P=12;
     else    
         Final_Fai_0=result_1(i,j,1);
         Final_alpha=result_1(i,j,2);
         Final_P=Alternative_P(j);
     end
     
     P_Model=P_Calculater_Model{j};
     
 end 
%*********************Finished Choose the alternative model or not**********************

%*************************************************************************************************
 % Output the P relative vaule
 fprintf('%s %s\n', 'The P Calculater model is :',P_Model);
 if strcmp(P_Model,P_Calculater_Model{3})
     fprintf('%s %15.8f\n', 'The Fai_0 is :',Final_Fai_0);
     fprintf('%s %15.8f\n', 'The alpha_1 is :',Final_alpha_1);
     fprintf('%s %15.8f\n', 'The alpha_2 is :',Final_alpha_2);
     fprintf('%s %15.8f\n', 'The Delta_12 of the finest grid is :',Final_alpha_1*h(1)+Final_alpha_2*h(2)^2);
 elseif strcmp(P_Model,P_Calculater_Model{1})
     fprintf('%s %15.8f\n', 'The Fai_0 is :',Final_Fai_0);
     fprintf('%s %15.8f\n', 'The alpha is :',Final_alpha);
     fprintf('%s %15.8f\n', 'The Delta_1 of the finest grid is :',Final_alpha*h(1));
 elseif strcmp(P_Model,P_Calculater_Model{2}) 
     fprintf('%s %15.8f\n', 'The Fai_0 is :',Final_Fai_0);
     fprintf('%s %15.8f\n', 'The alpha is :',Final_alpha);
     fprintf('%s %15.8f\n', 'The Delta_2 of the finest grid is :',Final_alpha*h(1)^2);
 elseif strcmp(P_Model,P_Calculater_Model{4}) 
     fprintf('%s %15.8f\n', 'The P is :',Final_P);
     fprintf('%s %15.8f\n', 'The Fai_0 is :',Final_Fai_0);
     fprintf('%s %15.8f\n', 'The alpha is :',Final_alpha);
     fprintf('%s %15.8f\n', 'The Delta_2 of the finest grid is :',Final_alpha*h(1)^Final_P);
 end
%*************************************************************************************************

%*************************************************************************************************
 % Calculate the Extrapolated values  at each point
 fid_Extrapolated_value_Full=fopen('Extrapolated_value_Full.dat','w');
 fprintf('%s %s\n', 'The P Calculater model is :',P_Model);
 if strcmp(P_Model,P_Calculater_Model{3})
     low_mesh_ratio=min(h);
     current_mesh_ratio=low_mesh_ratio;
     upper_mesh_ratio=max(h);
     i=0;
     while(current_mesh_ratio<=upper_mesh_ratio)
         i=i+1;
         save_mesh_ratio(i)=current_mesh_ratio;
         Extrapolated_value_Full(i)=Final_Fai_0+Final_alpha_1*current_mesh_ratio+Final_alpha_2*current_mesh_ratio^2;
         fprintf(fid_Extrapolated_value_Full,'%15.8f %15.8f\n',current_mesh_ratio,Extrapolated_value_Full(i));
         current_mesh_ratio=current_mesh_ratio+0.1;
         continue
     end
    
 elseif strcmp(P_Model,P_Calculater_Model{1})
     low_mesh_ratio=min(h);
     current_mesh_ratio=low_mesh_ratio;
     upper_mesh_ratio=max(h);
     i=0;
     while(current_mesh_ratio<=upper_mesh_ratio)
         i=i+1;
         save_mesh_ratio(i)=current_mesh_ratio;
         Extrapolated_value_Full(i)=Final_Fai_0+Final_alpha*current_mesh_ratio;
         fprintf(fid_Extrapolated_value_Full,'%15.8f %15.8f\n',current_mesh_ratio,Extrapolated_value_Full(i));
         current_mesh_ratio=current_mesh_ratio+0.1;
         continue
     end
 elseif strcmp(P_Model,P_Calculater_Model{2}) 
     low_mesh_ratio=min(h);
     current_mesh_ratio=low_mesh_ratio;
     upper_mesh_ratio=max(h);
     i=0;
     while(current_mesh_ratio<=upper_mesh_ratio)
         i=i+1;
         save_mesh_ratio(i)=current_mesh_ratio;
         Extrapolated_value_Full(i)=Final_Fai_0+Final_alpha*current_mesh_ratio^2;
         fprintf(fid_Extrapolated_value_Full,'%15.8f %15.8f\n',current_mesh_ratio,Extrapolated_value_Full(i));
         current_mesh_ratio=current_mesh_ratio+0.1;
         continue
     end
 elseif strcmp(P_Model,P_Calculater_Model{4}) 
     low_mesh_ratio=min(h);
     current_mesh_ratio=low_mesh_ratio;
     upper_mesh_ratio=max(h);
     i=0;
     while(current_mesh_ratio<=upper_mesh_ratio)
         i=i+1;
         save_mesh_ratio(i)=current_mesh_ratio;
         Extrapolated_value_Full(i)=Final_Fai_0+Final_alpha*current_mesh_ratio^Final_P;
         fprintf(fid_Extrapolated_value_Full,'%15.8f %15.8f\n',current_mesh_ratio,Extrapolated_value_Full(i));
         current_mesh_ratio=current_mesh_ratio+0.1;
         continue
     end

 end
%*************************************************************************************************

%***********************Determine a data range parameter********************************
Delta_Fai=(max(Fai)-min(Fai))/(n_g-1.0);

if(Temp_Final_P>=0.5&&Temp_Final_P<2.1 && Delta<Delta_Fai)
    Fs=1.25;
else
    Fs=3;
end

for k=1:n_g

    % Calculate Eplsion_Fai and Fai_fit
     if strcmp(P_Model,P_Calculater_Model{1})
         Eplsion_Fai=Final_alpha*h(k);
         Fai_fit(k)=Final_Fai_0+Final_alpha*h(k);
     elseif strcmp(P_Model,P_Calculater_Model{2}) 
         Eplsion_Fai=Final_alpha*h(k)^2;
         Fai_fit(k)=Final_Fai_0+Final_alpha*h(k)^2;
     elseif strcmp(P_Model,P_Calculater_Model{3})
         Eplsion_Fai=Final_alpha_1*h(k)+Final_alpha_2*h(k)^2;
         Fai_fit(k)=Final_Fai_0+Final_alpha_1*h(k)+Final_alpha_2*h(k)^2;
     elseif strcmp(P_Model,P_Calculater_Model{4})
         Eplsion_Fai=Final_alpha*h(k)^Final_P;
         Fai_fit(k)=Final_Fai_0+Final_alpha*h(k)^Final_P;
     end

    % Calculate the Uncertainty
    if(Delta<Delta_Fai)
        U_Fai(k)=Fs*Eplsion_Fai+Delta+abs(Fai(k)-Fai_fit(k));
    else
        U_Fai(k)=3.0*Delta/Delta_Fai*(Eplsion_Fai+Delta+abs(Fai(k)-Fai_fit(k)));
    end
end

%************************************************************************************************
%*******************************Output the Uncertainty Data******************************
fid_Complex_Full_Uncertainty_Data=fopen('Complex_Full_Uncertainty_Data.dat','w');
for i=1:n_g
      fprintf(fid_Complex_Full_Uncertainty_Data,'%15.8f %15.8f\n',h(i),abs(U_Fai(i)));
end
fclose(fid_Complex_Full_Uncertainty_Data);

fid_Complex_Full=fopen('Complex_Full.dat','w');
fprintf(fid_Complex_Full,'%s %15.8f\n', 'The_Final_P is :',Final_P);
fprintf(fid_Complex_Full,'%s %15.8f\n', 'The_Fai_0 is :',Final_Fai_0);
fclose(fid_Complex_Full);


%*************************************************************************************************
%*******************************************Output the figure*********************************
% P_Simple(2:12)=1.0;
% P_Complex_Three(2:12)=2.0;
% Final_P=3.0;
plotx=h(2:12);
ploty1=P_Simple(2:12);
ploty2=P_Complex_Three(2:12);
Final_P_Complex(2:12)=Final_P;
ploty3=Final_P_Complex(2:12);

%make the Observed Order of Grid Convergence (P) plot
figure;
plot(plotx,ploty1,'r^',plotx,ploty2,'ko',plotx,ploty3,'-m');grid on;
% Set the coordinate font size
set(gca,'FontName','Times New Roman','FontSize',12);
xlabel('Grid Ratio','FontName','Times New Roman','FontSize',13);
ylabel('P','FontName','Times New Roman','FontSize',13);
title('Observed Order of Grid Convergence (P)','FontName','Times New Roman','FontSize',13);       %plot title
legend('Procedure 1','Procedure 2','Procedure 3');
print(gcf,'-dpng','Observed_convergence_order.png')       %output the figure

%make the Extrapolated values plot
% Fai_ext_21(2:12)=1;
% Fai_0(2:12)=2;
% Final_Fai_0=3;
ploty1=Fai_ext_21(2:12);
ploty2=Fai_0(2:12);
Final_P_Complex(2:12)=Final_Fai_0;
ploty3=Final_P_Complex(2:12);

figure;
plot(plotx,ploty1,'r^',plotx,ploty2,'ko',plotx,ploty3,'-m');grid on;
% Set the coordinate font size
set(gca,'FontName','Times New Roman','FontSize',12);
xlabel('Grid Ratio','FontName','Times New Roman','FontSize',13);
ylabel('\phi_0(\phi_{ext})','FontName','Times New Roman','FontSize',13);
title('Extrapolated Value (\phi)','FontName','Times New Roman','FontSize',13);       %plot title
legend('Procedure 1','Procedure 2','Procedure 3');
print(gcf,'-dpng','Extrapolated_Value.png')       %output the figure

%U_Fai(1:13)=1.0;
plotx=h;
ploty1=U_Fai;

figure;
plot(plotx,ploty1,'-m*');grid on;
% Set the coordinate font size
set(gca,'FontName','Times New Roman','FontSize',12);
xlabel('Grid Ratio','FontName','Times New Roman','FontSize',13);
ylabel('Uncertainty (U_{\phi})','FontName','Times New Roman','FontSize',13);
title('Uncertainty (U_{\phi})','FontName','Times New Roman','FontSize',13);       %plot title
legend('Uncertainty');
print(gcf,'-dpng','Uncertainty.png')       %output the figure
 
% fprintf('%s %f\n', 'The time of the program cost is :',etime(end_time,start_time));

fprintf('%s\n', '************************Calculation Finished***********************');  
% Record the calculation start time
end_time=clock;

fprintf('%s\n', 'The program have run to the end, Please press "enter" to finish it!');
fprintf('%s %f\n', 'The time of the program cost is :',etime(end_time,start_time));

diary off;   % End the recording




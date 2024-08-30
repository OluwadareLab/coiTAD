
Mean_Intra = [];
Mean_Inter = [];
Mean_Pcc = [];
Mean_diff = [];
D = 0;
current_TADQ_Index = 1;

    for  radius = min_radius : max_radius
            tdname = [out_path,'HDBSCAN_',num2str(radius),'_TAD_BinID.txt'];             
        if (TAD_Q(current_TADQ_Index, 1) ~= 0)
            TDname = strcat(tdname); % path to recognized domain
            TD = dlmread(TDname);
            % The average intra-TD and inter-TD contact frequencies
            len = length(TD(:,1)); % find the number of TD
            INTER  = [];
            INTRA = [];
            PCC = [];
            
            for i = 1:len
                %=============================================
                % For each Domain find intra-TD and inter-TD
                % Intra(i) denote the average of contact frequencies between bins within the same TD i
                % Inter(i, j) denote the% average of contact frequencies between a bin in TD i and a bin in adjacent TD j, where |i–j | = 1
                %---------------------------------------------
                % No 1: Intra    
                domain1 = TD(i,:);     
                Intra = intra(domain1,Chr_Data);
                %=============================================
                % IntraPCC 
                pcc = intraPCC(domain1,Chr_Data);  
                PCC = [PCC;pcc]; % Concatenates the intraPCC results
                %---------------------------------------------
        
                % No 2: Inter
                if (i==1)
                     % Adjacent td is just the next one            
                    if (len == i)
                        domain2 = TD(i,:); 
                    else
                        domain2 = TD(i+1,:); 
                    end
                    
                   [Sum,count]= inter_v2(domain1,domain2,Chr_Data);
                   Inter = Sum/count;
                elseif(i>1 && i < len)
                     % Adjacent td includes prev and next
                    domain2 = TD(i-1,:);
                    [Sum1,count1]= inter_v2(domain2,domain1,Chr_Data);
                    domain2 = TD(i+1,:);
                    [Sum2,count2]= inter_v2(domain1,domain2,Chr_Data);
                    Inter = (Sum1 + Sum2) / (count1 + count2);
                elseif(i==len)
                     % Adjacent td id just the  prev one     
                    domain2 = TD(i-1,:);
                    [Sum,count] = inter_v2(domain2,domain1,Chr_Data);
                    Inter = Sum/count;
                end
                INTER = [INTER;Inter];
                INTRA = [INTRA;Intra];
        
        
            end
        
            out = INTRA - INTER; % Concaenate the Intra anf Inter Results
            
        
            %outputa
            fprintf('The average PCC over all TD= %f\n', mean(PCC));
            fprintf('The average IF for over each TD= %f\n', mean(INTRA));
            fprintf('The average_Inter IF for over all TDs= %f\n', mean(INTER));
            fprintf('The average difference= %f\n', mean(out));
        
            Mean_Intra = [ Mean_Intra ; mean(INTRA)];
            Mean_Inter = [Mean_Inter;mean(INTER)];
            Mean_Pcc = [Mean_Pcc;mean(PCC)];
            Mean_diff = [ Mean_diff;mean(out)];
            
            if (D < mean(out))
                    D =  mean(out);
                    Str = ['KM_',num2str(radius),'_domain.txt'];
                    bestRadius = radius; 
             end
        else 
            fprintf("Skipping radius %d", radius)
        end
    
        intername = [out_path_1 ,'Overall_Inter_',Resolution,'.csv'];
        intraname = [out_path_1 ,'Overall_Intra_',Resolution,'.csv'];
        diffname = [out_path_1 ,'Overall_Diff_',Resolution,'.csv'];
        pccname = [out_path_1 ,'Overall_pcc_',Resolution,'.csv'];
        
        %dlmwrite(intername,Mean_Inter);
        %dlmwrite(intraname, Mean_Intra);  
        %dlmwrite(diffname, Mean_diff);
        %dlmwrite(pccname, Mean_Pcc);

        current_TADQ_Index = current_TADQ_Index + 1; 
    end
        %----------------------------------------------------------------------
        % Find the maximum difference
        %----------------------------------------------------------------------
        f= max(Mean_diff);
        readme=strcat(out_path_1,'Readme_',name,'.txt');
        msg = ['Recommended TAD = ', Str,' with value ',num2str(D)];
        fid = fopen(readme,'wt');
        fprintf(fid, msg');
        fclose(fid);
	    
	     if (D >MaxQuality)
           MaxQuality=  D;       
           Outputmsg = ['The name of Best TAD identified Can be found Here:' readme '\n' 'Recommended TAD = ', Str,' with value ',num2str(D) '\n'];
           
         end
    
    % and the difference between the average intra-TD and inter-TD contact frequencies
    
    % This find the intra-TDcontact frequency and returns average
    function [Intra_Average,Sum,count] = intra(domain,Chr_Data)
    Start = domain(:,1);
    End = domain(:,2);
    Sum = 0;
    count = 0 ;
    for i = Start:End
        for j = i+1 : End
           count =count  + 1;
           Sum = Sum + Chr_Data(i,j); 
        end
        
    end
    
    % Find the average
    
    if (Sum >  0 && count > 0)
        Intra_Average = Sum/count;
    else
        Intra_Average = 0;
    end
    
    end

function [ Intra_Average] = intraPCC( domain,Chr_Data )
%% Pearson’s correlation coefficient (PCC) between the contact profiles of bins within a TD as a quality measurement
Start = domain(:,1);
End = domain(:,2);
Mat = [];
for i = Start:End
    temp = [];
    for j = Start : End
        temp = [temp Chr_Data(i,j)];       
    end
    Mat = [Mat;temp];
    
end
% Now that Matrix has been formed, do a correlation
Matrix = corr(Mat); % Matrix contains the pearsons correlation of each bin element
Matrix(isnan(Matrix)) = -0.2;
%-----------------------------------------------------------------------------
% The average intra-TD Pearson’s correlation coefficient (PCC) and weighted PCC
%--------------------------------------------------------------------------
%1. Intra PCC
 num = size(Mat,2);
 Intra_Average = sum(sum(triu(Matrix)))/((num+1)*num/2); % Average

end

% This progrms can find the inter Inter(i, j) denote the
% average of contact frequencies between a bin in TD i and a
% bin in adjacent TD j, where |i–j | = 1
function [Sum,count] = inter_v2(domain_i, domain_j,Chr_Data)

%domain i
Start_i = domain_i(:,1);
End_i = domain_i(:,2);
%domain j
Start_j = domain_j(:,1);
End_j = domain_j(:,2);

Sum = 0;
count = 0 ;
incr = 0;
for i = Start_i:Start_j-1
%     fprintf('Value of i= %d\n',i);
    incr = incr + 1;
    c = 0;
    for j = End_i+1 :  End_j
%     fprintf('Value of i= %d\n',j);
    c = c + 1;
    count =count  + 1;
%     fprintf('%d and %d\n',i,j);
    Sum = Sum + Chr_Data(i,j); 
    if(c==incr)
        break;
    end     
      
    end
end


end

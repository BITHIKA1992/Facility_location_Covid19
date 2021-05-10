%%  Simulation program to maximize social distance in an area considering it as a graph with nodes and vertex
%%  by opening facility centres to cater the need of the people.

clearvars ; clc
tic

%% Load Input Data
load input_data ;
P_org = P;
final_social_distant = [];

nop_wsd = 4 ;             % nop_wsd = 4
A = 10;                   % A = 2000
b = 0.5;                  % b = 0.5

facility_ot = 7 ;
facility_ct = 12 ;
tof = 3;                  % Type of Facilities

mean_int_arrv = 1.00;
mean_int_service = 0.70;
nosim = 1;
social_distn = zeros(nosim,5) ;


%% Sorted Nodes based on Demand 
noN = 100;
dist_mat = dist_com(1:noN,1:noN);   % C
P = P_org(1:noN);
[B, I] = sort(P, 'descend');

%% Arrival Time Creation
arrival_int_time = zeros(noN, max(P));
for i=1:noN
    arrival_int_time(i,1:P(i)) = exprnd(mean_int_arrv,[1 P(i)]);
end
arrival_time = cumulative_function(arrival_int_time, facility_ot);

%% Nodes Assignment among Open Nodes
output_run = [];
for noF = 20:20:noN
%noF = 100;
asign_nodes = zeros(noF,10);
asign_demand = zeros(noF,10);
unassigned_nodes = 1:noN ;

sort_nodes_dem = I;

for i=1:noN
    if rem(i,noF)==0
        j=noF;
    else
        j= rem(i,noF) ;
    end
    k = ceil(i/noF);
    if rem(k,2)==1
    asign_nodes(j,k)= I(i) ;
    asign_demand(j,k)= B(i);
    pos = find(unassigned_nodes(1,:)== I(i)) ;
    unassigned_nodes(pos)=[];
    else
       asign_nodes(j,k)= I(min(k*noF,noN) - (i-(k-1)*noF -1)) ; 
       asign_demand(j,k)= B(min(k*noF,noN) - (i-(k-1)*noF -1));
       pos = find(unassigned_nodes(1,:)== I(min(k*noF,noN) - (i-(k-1)*noF -1))) ;
       unassigned_nodes(pos)=[];
    end    
end


%% Arrival Time Sorting accross Open Facilities

open_facility = I(1:noF);   

max_col = sum(P);
queue_memb = zeros(noN, max_col);
for i=1:noN
    if sum(i==open_facility)==1
       queue_memb(i,1:length(arrival_int_time(i,:)))= arrival_time(i,:);
    end
end

for i=1:noN
    node = I(i);
    [row,col] = find(asign_nodes==node) ;
    if sum(node==open_facility)==0
       %% For Grid Network
%        dist_temp = abs(dist_mat(i,i) - dist_mat(i,open_facility));
%        [~, fac_use_idx] = min(dist_temp) ;
       
       %% For Complete Network
       %[~, fac_use_idx] = min(dist_mat(i,open_facility)) ;
       
       facility_used = asign_nodes(row,1) ;
       
       for j=1:sum(arrival_time(node,:)>0)
           [~, ~, cur_mem] = find(queue_memb(facility_used,:)) ;
           [~, nCol, max_value] = find(max(arrival_time(node,j) - cur_mem,0));
           
           if isempty(nCol)
               queue_memb(facility_used,2:length(cur_mem)+1)= queue_memb(facility_used,1:length(cur_mem));
               queue_memb(facility_used,1)= arrival_time(node,j);        
           else
               %queue_memb(facility_used,1:max(nCol))= queue_memb(facility_used,1:max(nCol));
               queue_memb(facility_used,max(nCol)+2:length(cur_mem)+1)= queue_memb(facility_used,max(nCol)+1:length(cur_mem));
               queue_memb(facility_used,max(nCol)+1)= arrival_time(node,j);
           end           
       end
    end
end

people_cnt_fw = zeros(noN,1) ;

for i =1:noN
    people_cnt_fw(i,1) = sum(queue_memb(i,:)>0) ;
end

%% Service Time Generation

inter_service_time = zeros(noN, max_col);
service_time = zeros(noN, max_col);

for i=1:noN
    if sum(i==open_facility)==1
       inter_service_time(i,1:max_col) = exprnd(mean_int_service,[1 max_col]);
    end

    if sum(i==open_facility)==1
       for j=1:max_col
           if j==1
              intg = floor(queue_memb(i,j));
              frac = queue_memb(i,j) - intg ;
              frac = frac*100 + inter_service_time(i,j) ;
              hh = intg + floor(frac/60);
              mm = rem(frac,60);     
              service_time(i,j) = hh + mm/100 ;
           else
              intg = floor(max(queue_memb(i,j), service_time(i,j-1)));
              frac = max(queue_memb(i,j), service_time(i,j-1)) - intg ;
              frac = frac*100 + inter_service_time(i,j) ;
              hh = intg + floor(frac/60);
              mm = rem(frac,60); 
              service_time(i,j) = hh + mm/100 ; 
           end
       end        
    end
end
% service_time = cumulative_function(inter_service_time, facility_ot);

%% Number of People in Queue

people_queue = zeros(noN, max_col);

for i=1:noN
    if sum(i==open_facility)==1
        for j=1:people_cnt_fw(i,1) %sum(P)
            if j>1 && queue_memb(i,j)==0
               %people_queue(i,j) = 0 ;
               people_queue(i,j) = max(people_queue(i,j-1)-1,0) ;
            else
               people_queue(i,j) = sum(queue_memb(i,j)<service_time(i,1:j)) ;
            end
        end
    end
end

%% Number of People Served in a given time

people_served_fw = zeros(noN, 1);       % People served facility wise
people_queue_tb = zeros(noN, max_col);   % People in queue up to a given time

for i=1:noN
    if sum(i==open_facility)==1
        for j=1:max_col
            if service_time(i,j) < facility_ct
                if queue_memb(i,j)==0
                    people_served_fw(i,1) = people_served_fw(i,1) + 0 ;
                else
                    people_served_fw(i,1) = people_served_fw(i,1) + 1 ;
                end
                
                if j>1 && queue_memb(i,j)==0
                   people_queue_tb(i,j) = max(people_queue_tb(i,j-1)-1,0) ;
                else
                   people_queue_tb(i,j) = sum(queue_memb(i,j)<service_time(i,1:j)) ;
                end
            else
                break ;
            end
        end
    end
end


%% Social Distancing (Linear and Exponential Factor)
social_distn_le = zeros(noN, sum(P));
c = 2.00;
d = 0.25;

for i=1:noN
    if sum(i==open_facility)==1
        for j=1:sum(P)
            if people_queue(i,j) <= nop_wsd
               social_distn_le(i,j) = A ;
            elseif people_queue(i,j) <= c*nop_wsd
               social_distn_le(i,j) = A - b*(people_queue(i,j)- nop_wsd) ;
            else
               social_distn_le(i,j) = A - b*(c-1)*nop_wsd - b*exp((people_queue(i,j)- c*nop_wsd)^0.25) ;  
            end
        end
    end
end

%% Summation over Queue Length
r = 1;
social_distn(r,1) = noF ;
social_distn(r,2) = sum(sum(social_distn_le,2)) ;

%% Average over Queue Length
social_distn(r,3) = mean(mean(people_queue,2)) ;

%% People Served and Not Served
social_distn(r,4) = sum(people_served_fw);
social_distn(r,5) = sum(people_cnt_fw) - sum(people_served_fw) ;

%% Final Output
location = zeros(1, noN);
location(1, 1:noF) = open_facility;
output = [social_distn location];
output_run = [output_run; output];
end
toc


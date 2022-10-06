function [track_exposed, exposed_listSum, who_infect_list,  ...
    hospital_incidence, death_incidence, V2_list, vaccine_ever_list3, track_exposed_earlest] = ...
    Github_run_booster(endDays_Sim,  beta_others, beta_house, house_neighs, neighborsOthers, ...
    sigma, gamma, gammaA, prop_sym, omega_asym, omega_preasym, epsilon, VEDs, VEIs, VEHs, VESs, ...
    V1t, V2t, V3t, Rt, firstdose_rate, seconddose_rate, thirdose_rate, Ages1518, Ages5, end_day_dose, end_num_infections, daysBoosterTaken)

% x100 for % people initially infected and in exposed status
popImmunity = 0.01;
daysProtectV = [28, daysBoosterTaken];
% for infection
relativeVEI1 = 1.-VEIs(1,:);
relativeVEI2 = 1.-VEIs(2,:);
relativeVEI3 = 1.-VEIs(3,:);
relativeVEIR = 1.-VEIs(4,:);

% for symptom
relativeVES1 = 1.-VESs(1,:);
relativeVES2 = 1.-VESs(2,:);
relativeVES3 = 1.-VESs(3,:);
relativeVESR = 1.-VESs(4,:);

% for hospital
relativeVEH1 = (1.-VEHs(1,:))./(1.-VESs(1,:));
relativeVEH2 = (1.-VEHs(2,:))./(1.-VESs(2,:));
relativeVEH3 = (1.-VEHs(3,:))./(1.-VESs(3,:));
relativeVEHR = (1.-VEHs(4,:)+10^(-10))./(1.-VESs(4,:)+10^(-10));

% for death
relativeVED1 = (1.-VEDs(1,:))./(1.-VEHs(1,:));
relativeVED2 = (1.-VEDs(2,:))./(1.-VEHs(2,:));
relativeVED3 = (1.-VEDs(3,:))./(1.-VEHs(3,:));
relativeVEDR = (1.-VEDs(4,:)+10^(-10))./(1.-VEHs(4,:)+10^(-10));

%% symptomatic case hospitalization ratio
p_sym_hosp = [0, 0.025, 2.672, 9.334, 15.465]*0.01;
eta = 1/5.9;

deathRate = [0.0390, 0.1208, 0.0304, 0.1049, 0.2269];
death_mu = 0.12821;
gamma_H = 0.0912409;

N = length(V1t);
suscept_list = ones(N,1); 
exposed_list = zeros(N,1); 
presymp_list = zeros(N,1); 
symp_list = zeros(N,1); 
asym_list = zeros(N,1); 
recover_list = zeros(N,1); 
hospital_list = zeros(N,1); 
death_list = zeros(N,1); 
hospital_incidence = zeros(N,1);
symp_ever_list = zeros(N,1);
V1_time_list = V1t-end_day_dose;
V1_time_list( find(V1_time_list==-end_day_dose) ) = 10^5;
V2_time_list = V2t-end_day_dose;
V2_time_list( find(V2_time_list==-end_day_dose) ) = 10^5;
V3_time_list = V3t-end_day_dose;
V3_time_list( find(V3_time_list==-end_day_dose) ) = 10^5;

% timing of people recovered
recover_time_list = Rt-end_day_dose;
recover_time_list( find(recover_time_list==-end_day_dose) ) = 10^5;

% history of people vaccinated
vaccine_ever_list1 = zeros(N,1);
vaccine_ever_list2 = zeros(N,1);
vaccine_ever_list3 = zeros(N,1);
vaccine_ever_list1( find(V1_time_list<10^5) ) = 1;
vaccine_ever_list2( find(V2_time_list<10^5) ) = 1;
vaccine_ever_list3( find(V3_time_list<10^5) ) = 1;

% timing of vaccination vaccinated
V1_list = zeros(N,1);
V2_list = zeros(N,1);
V3_list = zeros(N,1); 
V1_list(find(V1t>0)) = 1;
V2_list(find(V2t>0)) = 1;
V3_list(find(V3t>0)) = 1;
recover_list(find(Rt>0)) = 1;
suscept_list(find(V1_list==1)) = 1;
suscept_list(find(V2_list==1)) = 1;
suscept_list(find(V3_list==1)) = 1;
suscept_list(find(recover_list==1)) = 1;

who_infect_list = 0;
who_infect_count = 0;

% Track who is exposed at time t
track_exposed = 10^5*ones(N,1); %tracks time at infection for each node.

% Initialization
p_zero = randperm(N, ceil(N*popImmunity));
track_time = 0;
suscept_list(p_zero) = 0;
exposed_list(p_zero) = 1;
track_exposed(p_zero) = 1;
V1_list(p_zero) = 0;
V2_list(p_zero) = 0;
V3_list(p_zero) = 0;
recover_list(p_zero) = 0;

% Track who is exposed at the first time
track_exposed_earlest = track_exposed;

% Daily number for first, second, and third dose
numVaccPerDay1 = round(N * firstdose_rate/100);
numVaccPerDay2 = round(N * seconddose_rate/100);
numVaccPerDay3 = round(N * thirdose_rate/100);

% Only can vaccinate
rankNodes = find( (Ages1518==1 | Ages5>=3) & (V1t == 0 & V2t ==0 & V3t ==0 & Rt==0) );

%%
exposed_listSum = zeros(endDays_Sim,1);

ToH = -1*ones(N,1); % if -1: not applied, if 0, to R, if 1 to H
ToHrate = zeros(N,1);
end_flag = 0;
% presymp_list_new = zeros(N,1);
for ic=1:endDays_Sim
    
    if end_flag==0
        track_time = track_time + 1;
        % S to V1
        % Can vaccinate S, E, A P R, which never in V1
        temp_numVaccPerDay = numVaccPerDay1(ic);
        temp = (suscept_list|exposed_list|asym_list|presymp_list|recover_list) ...
            & ~vaccine_ever_list1;
        temp = find(temp==1);
        for j=1:length(temp)
            i = temp(j);
            if temp_numVaccPerDay>0
                vaccine_ever_list1(i) = 1;
                temp_numVaccPerDay = temp_numVaccPerDay-1;
                if suscept_list(i) == 1
                    suscept_list(i) = 0;
                    V1_list(i) = 1;
                    V1_time_list(i) = track_time;
                end
            end
        end

        % initilize
        temp_numVaccPerDay = numVaccPerDay2(ic);
        temp1 = (track_time-V1_time_list)>=daysProtectV(1);
        temp = (V1_list|exposed_list|asym_list|presymp_list|recover_list) & temp1...
            & vaccine_ever_list1 & ~vaccine_ever_list2;
        temp = find(temp==1);
        for j=1:length(temp)
            i = temp(j);
            if temp_numVaccPerDay>0
                vaccine_ever_list2(i) = 1;
                temp_numVaccPerDay = temp_numVaccPerDay-1;
                if V1_list(i) == 1
                    V1_list(i) = 0;
                    V2_list(i) = 1;
                    V2_time_list(i) = track_time;
                end
            end
        end
        
        % from V2 to V3
        temp_numVaccPerDay = numVaccPerDay3(ic);
        temp1 = (track_time-V2_time_list)>=daysProtectV(2);
        temp = (V2_list|exposed_list|asym_list|presymp_list|recover_list) & temp1...
            & vaccine_ever_list2 & ~vaccine_ever_list3;
        temp = find(temp==1);
        for j=1:length(temp)
            i = temp(j);
            if temp_numVaccPerDay>0
                vaccine_ever_list3(i) = 1;
                temp_numVaccPerDay = temp_numVaccPerDay-1;
                if V2_list(i) == 1
                    V2_list(i) = 0;
                    V3_list(i) = 1;
                    V3_time_list(i) = track_time;
                end
            end
        end

        %% infect the susceptibles with the probability 1-exp(-bk) where k is the number of infeced neighbors.
        exposed_list_new = zeros(N,1);
        NodeList =  find(suscept_list==1);
        for i=1:length(NodeList)
            node = NodeList(i);
            pTemp = 0;

            % for household and non-househod transmission from those
            % non-isolated contacts.
            pTemp = pTemp+(1- exp(-beta_house*infected_neighborsSEIR(house_neighs, node, symp_list)...
                - omega_asym*beta_house*infected_neighborsSEIR(house_neighs, node, asym_list) ...
                - omega_preasym*beta_house*infected_neighborsSEIR(house_neighs, node, presymp_list)));
            pTemp = pTemp+(1- exp(-beta_others*infected_neighborsSEIR(neighborsOthers, node, symp_list)...
                - omega_asym*beta_others*infected_neighborsSEIR(neighborsOthers, node, asym_list) ...
                - omega_preasym*beta_others*infected_neighborsSEIR(neighborsOthers, node, presymp_list)));

            %% mainly for Re estimationw when without testing.
            %  So do not consdier the impact of isolated people no matter they
            %   have infectiousness or not.
            if rand<pTemp %&& rel_beta_isolated<=0
                exposed_list_new(node)=1;
                temp01 = house_neighs{node};
                temp11 = symp_list(temp01)+omega_asym*asym_list(temp01)+omega_preasym*presymp_list(temp01);

                temp02 = neighborsOthers{node};
                temp12 = symp_list(temp02)+omega_asym*asym_list(temp02)+omega_preasym*presymp_list(temp02);

                temp = [temp01, temp02];
                temp1 = [beta_house.*temp11; beta_others.*temp12];

                temp2 = temp(find( mnrnd(1, temp1./sum(temp1)) == 1));
            end
        end

        for ii=1:4
            if ii==1; NodeList = find(V1_list==1); end
            if ii==2; NodeList = find(V2_list==1); end
            if ii==3; NodeList = find(V3_list==1); end
            if ii==4; NodeList = find(recover_list==1); end

            for i=1:length(NodeList)
                node=NodeList(i);
                pTemp = 0;

                temp_beta = 1;
                if ii==1
                    if ic>V1_time_list(node)
                        temp_beta = relativeVEI1(ic-V1_time_list(node));
                    end
                end
                if ii==2
                    if ic>V2_time_list(node)
                        temp_beta = relativeVEI2(ic-V2_time_list(node));
                    end
                end
                if ii==3
                    if ic>V3_time_list(node)
                        temp_beta = relativeVEI3(ic-V3_time_list(node));
                    end
                end

                if ii==4
                    if ic>recover_time_list(node)
                        temp_beta = relativeVEIR(ic-recover_time_list(node));
                    end
                end

                % for household and non-househod transmission from those
                % non-isolated contacts.
                pTemp = pTemp+(1- exp(-temp_beta*beta_house*infected_neighborsSEIR(house_neighs, node, symp_list)...
                    - omega_asym*temp_beta*beta_house*infected_neighborsSEIR(house_neighs, node, asym_list) ...
                    - omega_preasym*temp_beta*beta_house*infected_neighborsSEIR(house_neighs, node, presymp_list)));

                pTemp = pTemp+(1- exp(-temp_beta*beta_others*infected_neighborsSEIR(neighborsOthers, node, symp_list)...
                    - omega_asym*temp_beta*beta_others*infected_neighborsSEIR(neighborsOthers, node, asym_list) ...
                    - omega_preasym*temp_beta*beta_others*infected_neighborsSEIR(neighborsOthers, node, presymp_list)));

                if rand<pTemp
                    exposed_list_new(node)=1;
                    temp01 = house_neighs{node};
                    temp11 = symp_list(temp01)+omega_asym*asym_list(temp01)+omega_preasym*presymp_list(temp01);
                
                    temp02 = neighborsOthers{node};
                    temp12 = symp_list(temp02)+omega_asym*asym_list(temp02)+omega_preasym*presymp_list(temp02);
                
                    temp = [temp01, temp02];
                    temp1 = [beta_house.*temp11; beta_others.*temp12];
                
                    temp2 = temp(find( mnrnd(1, temp1./sum(temp1)) == 1));
                    if end_num_infections~=-1
                        if who_infect_count <= end_num_infections
                            who_infect_count = who_infect_count+1;
                            %                         who_infect_list(who_infect_count) = temp2;
                            if (ismember(temp2, p_zero))
                                who_infect_list = who_infect_list+1;
                            end
                        else
                            end_flag = 1;
                        end
                    end

                end

            end
        end


        %% Store the time of infection for the newly infected nodes in track_exposed list
        track_exposed(find(exposed_list_new==1))= track_time;

        %All the infected nodes recover with a probailty
        asym_list_new = zeros(N,1);
        presymp_list_new = zeros(N,1);
        symp_list_new = zeros(N,1);
        recovered_list_new = zeros(N,1);
        asymrecovered_list_new = zeros(N,1);
        hospital_list_new = zeros(N,1);
        death_list_new = zeros(N,1);


        %%
        NodeList = find(exposed_list==1);
        for i=1:length(NodeList)
            node=NodeList(i);
            if rand<epsilon
                temp_prop_sym = prop_sym;
                if V3_time_list(node)<10^5;
                    temp_psi  = relativeVES3(ic-V3_time_list(node)+1);
                    temp_prop_sym = temp_psi*prop_sym;
                else
                    if V2_time_list(node)<10^5;
                        temp_psi  = relativeVES2(ic-V2_time_list(node)+1);
                        temp_prop_sym = temp_psi*prop_sym;
                    else
                        if V1_time_list(node)<10^5;
                            temp_psi  = relativeVES1(ic-V1_time_list(node)+1);
                            temp_prop_sym = temp_psi*prop_sym;
                        end
                    end
                end

                if recover_time_list(node)<10^5;
                    temp_psi  = relativeVESR(ic-recover_time_list(node)+1);
                    temp_prop_sym = temp_psi*prop_sym;
                end

                if rand < temp_prop_sym
                    presymp_list_new(node)=1;
                else
                    asym_list_new(node)=1;
                end
            end
        end

        %%
        NodeList = find(presymp_list==1);
        for i=1:length(NodeList)
            node=NodeList(i);
            if rand<sigma
                symp_list_new(node)=1;
            end
        end


        %%  ToH = -1*ones(N,1); % if -1: not applied, if 0, to R, if 1 to H
        NodeList = find(symp_list==1);
        for i=1:length(NodeList)
            node=NodeList(i);

            temp_phosp = 1;
                
            if V3_time_list(node)<10^5;
                temp_phosp= relativeVEH3(ic-V3_time_list(node)+1);
            else
                if V2_time_list(node)<10^5;
                    temp_phosp= relativeVEH2(ic-V2_time_list(node)+1);
                else
                    if V1_time_list(node)<10^5;
                        temp_phosp= relativeVEH1(ic-V1_time_list(node)+1);
                    end
                end
            end

            if recover_time_list(node)<10^5;
                temp_phosp= relativeVEHR(ic-recover_time_list(node)+1);
            end

            if ToH(node) == -1
                if rand< temp_phosp*p_sym_hosp(Ages5(node) )
                    ToH(node) = 1;
                    ToHrate(node) = eta;
                else
                    ToH(node) = 0;
                    ToHrate(node) = gamma;
                end
            end
            if ToH(node) == 1
                if rand<ToHrate(node)
                    hospital_list_new(node)=1;
                end
            end
            if ToH(node) == 0
                if rand<ToHrate(node)
                    recovered_list_new(node)=1;
                end
            end
        end

        %
        NodeList = find(hospital_list==1);
        for i=1:length(NodeList)
            node=NodeList(i);

            temp_pdeath = 1;
                
            if V3_time_list(node)<10^5;
                temp_pdeath= relativeVED3(ic-V3_time_list(node)+1);
            else
                if V2_time_list(node)<10^5;
                    temp_pdeath= relativeVED2(ic-V2_time_list(node)+1);
                else
                    if V1_time_list(node)<10^5;
                        temp_pdeath= relativeVED1(ic-V1_time_list(node)+1);
                    end
                end
            end

            if recover_time_list(node)<10^5;
                temp_pdeath= relativeVEDR(ic-recover_time_list(node)+1);
            end


            temprand = rand;
            tempRandt = (temp_pdeath*deathRate(Ages5(node))*death_mu ...
                + (1-temp_pdeath*deathRate(Ages5(node)))*gamma_H );
            if temprand< tempRandt
                if temprand< (temp_pdeath*deathRate(Ages5(node))*death_mu)
                    death_list_new(node)=1;
                else
                    recovered_list_new(node)=1;
                end
            end
        end

        %%
        NodeList = find(asym_list==1);
        for i=1:length(NodeList)
            node=NodeList(i);
            if rand<gammaA
                asymrecovered_list_new(node)=1;
            end
        end


        %% Store the time of recovery for the newly  nodes in  list
        suscept_list((exposed_list_new==1)) = 0;

        % update the original exposed_list
        exposed_list((exposed_list_new==1)) = 1;
        exposed_list((presymp_list_new==1)) = 0;
        exposed_list((symp_list_new==1)) = 0;
        exposed_list((asym_list_new==1)) = 0;
        exposed_list((recovered_list_new==1)) = 0;

        % update orignial presymp_list
        presymp_list((exposed_list_new==1)) = 0;
        presymp_list((presymp_list_new==1)) = 1;
        presymp_list((symp_list_new==1)) = 0;
        presymp_list((asym_list_new==1)) = 0;
        presymp_list((recovered_list_new==1)) = 0;


        % update orignial symp_list
        symp_list((exposed_list_new==1)) = 0;
        symp_list((presymp_list_new==1)) = 0;
        symp_list((symp_list_new==1)) = 1;
        symp_list((asym_list_new==1)) = 0;
        symp_list((recovered_list_new==1)) = 0;
        symp_list((hospital_list_new==1)) = 0;
        symp_ever_list(find(symp_list==1)) = 1;

        %
        asym_list((exposed_list_new==1)) = 0;
        asym_list((presymp_list_new==1)) = 0;
        asym_list((symp_list_new==1)) = 0;
        asym_list((asym_list_new==1)) = 1;
        asym_list((asymrecovered_list_new==1)) = 0;

        hospital_list((hospital_list_new==1)) = 1;
        hospital_list((death_list_new==1)) = 0;
        hospital_list((recovered_list_new==1)) = 0;

        death_list((death_list_new==1)) = 1;

        recover_list((recovered_list_new==1)) = 1;
        recover_list((asymrecovered_list_new==1)) = 1;
        recover_list((exposed_list==1)) = 0;
        recover_list((asym_list==1)) = 0;
        recover_list((presymp_list==1)) = 0;
        recover_list(exposed_list_new==1) = 0;

        % recover timing
        recover_time_list((recovered_list_new==1)) = track_time;
        recover_time_list((asymrecovered_list_new==1)) = track_time;

        
        V1_list((asym_list_new==1)) = 0;
        V1_list((presymp_list_new==1)) = 0;
        V1_list((exposed_list_new==1)) = 0;

        %
        V2_list((exposed_list==1)) = 0;
        V2_list((asym_list==1)) = 0;
        V2_list((presymp_list==1)) = 0;
        V2_list((exposed_list_new==1)) = 0;

        %
        V3_list((exposed_list==1)) = 0;
        V3_list((asym_list==1)) = 0;
        V3_list((presymp_list==1)) = 0;
        V3_list((exposed_list_new==1)) = 0;

        hospital_incidence(hospital_list_new==1) = hospital_incidence(hospital_list_new==1)+1;
        temp = find( track_exposed_earlest>track_exposed);
        track_exposed_earlest(temp) = track_exposed(temp);

        exposed_listSum(ic) = sum(exposed_list_new);
    end
end
death_incidence = death_list;
who_infect_list = who_infect_list/length(p_zero);

fid1 = fopen('th_real.txt', 'w');
fid2 = fopen('th_paper.txt', 'w');
fid3 = fopen('model.txt', 'w');
fclose(fid1);
fclose(fid2);
fclose(fid3);
users = 100;
preambs = 20;
grants = 8;
p_det = 1;
p_det_col = 1;

for j = 1:users
    sumResult = 0;
    for i = 1:min(preambs,j)
        result = Pr_k_v(j, preambs, i);  
        sumResult = sumResult + (result*min(i, grants)); %(result*i);
    end
    %sumResult = sumResult * ((1-1/preambs)^(users - 1));
    fid = fopen('th_real.txt', 'a');
    if fid ~= -1
        fprintf(fid, '%f', sumResult);
        fprintf(fid, '%s', ', ');
        fclose(fid);
    else
        error('Cannot open output file');
    end
    calculate_probability(10000,j, preambs, grants, p_det, p_det_col);
    papers_theory(j, preambs, grants, p_det, p_det_col);
end

%for i=0.1:0.1:1
%    calculate_probability(10000,users, preambs, grants, p_det,i);
%    papers_theory(users, preambs, grants, p_det, i);
%end
th_real = load('th_real.txt');
th_paper = load('th_paper.txt');
model = load('model.txt');
x = 1:length(model);
figure; % create a new figure
plot(x, model, "b",x,th_paper,"r",x,th_real,"g"); 
%,
xlabel('количсетво абонентов'); 
ylabel('Требуемое количсетво грантов'); 
title('Plot'); 
grid on; 


function [res] = Pr_k_v(n, V, k)
res = 0;
%if k<=n && k<=V
    for j=k:min(V,n)
        res = res + (-1)^j*(V-j)^(n-j)/(factorial(j-k)*factorial(V-j)*factorial(n-j));
    end
    res = res * (-1)^k*factorial(V)*factorial(n)/(V^n*factorial(k));
%end
end

function final_number = count_successes(number_of_users, number_of_preambles, number_of_grants, p_dec, p_dec_col)
    preambles = zeros(1,number_of_preambles);
    successful_preambles = 0;
    final_number = 0;
    for i = 1:number_of_users
        random_pa = randi([1, number_of_preambles]);
        preambles(random_pa) = preambles(random_pa) + 1;
    end
    for i = 1:number_of_preambles
        decode =rand(1,1);
        detect_col = rand(1,1);
        if preambles(i) == 1 && decode<p_dec
            successful_preambles = successful_preambles + 1;
        end
        if preambles(i)>1 && detect_col>p_dec_col
            successful_preambles = successful_preambles + 1;
            preambles(i) = -1;
        end
    end
    succ_preambles_after_errors = zeros(1,successful_preambles);
    counter = 1;
    for i = preambles
        if i == 1 || i == -1
            succ_preambles_after_errors(counter)= i;
            counter = counter + 1;
        end
    end
    grants = min(successful_preambles, number_of_grants);
    for i = 1:grants
        random_choice = randi([1, successful_preambles]);
        if succ_preambles_after_errors(random_choice)== 1
            final_number = final_number + 1;
        end
    end
    %number = successful_preambles;
end

function optimal_grants = calculate_probability(number_of_iterations, number_of_users, number_of_preambles, number_of_grants, p_dec, p_dec_col)
    %probability_old = 0;
    optimal_grants = 0;
    %for j = 1:number_of_preambles
        
        successes_every_iteration = zeros(1,number_of_iterations);
        for i = 1:number_of_iterations
            numbers = count_successes(number_of_users, number_of_preambles, number_of_grants, p_dec, p_dec_col);
            successes_every_iteration(i) = numbers;
        end
        probability = sum(successes_every_iteration)/(number_of_iterations);%*number_of_users);
        %if (probability) == number_of_users
            %optimal_grants = j;
        %end
    %probability_old = probability_new;
    %end
    efficiency = probability/number_of_grants;
    fid = fopen('model.txt', 'a');
    if fid ~= -1
        fprintf(fid, '%f', probability);
        fprintf(fid, '%s', ', ');
        fclose(fid);
    else
        error('Cannot open output file');
    end
    %    fid = fopen('th_paper.txt', 'a');
    %if fid ~= -1
    %    fprintf(fid, '%f', efficiency);
     %   fprintf(fid, '%s', ', ');
     %   fclose(fid);
    %else
    %    error('Cannot open output file');
    %end
end

function p_tags = papers_theory(number_of_users, number_of_preambles, number_of_grants, p_dec, p_dec_col)
    p1 = number_of_users * (1/number_of_preambles)*((1 - 1/number_of_preambles)^(number_of_users-1));
    p0 = (1-1/number_of_preambles)^number_of_users;
    p_c = 1 - p1 -p0;
    p_without_errors = ((1-1/number_of_preambles)^(number_of_users - 1)) * min(1, number_of_grants/(p1*number_of_preambles));
    p_with_det = p_dec*((1-1/number_of_preambles)^(number_of_users - 1)) * min(1, number_of_grants/(p1*number_of_preambles*p_dec));
    p_with_dec_col = p_dec*((1-1/number_of_preambles)^(number_of_users - 1)) * min(1, number_of_grants/(p1*number_of_preambles*p_dec + p_c*number_of_preambles*(1-p_dec_col)));
    p_tags = p_with_dec_col*number_of_users;%*((p1*p_dec)/(p1*p_dec + p_c*(1-p_dec_col)));
    
    fid = fopen('th_paper.txt', 'a');
    if fid ~= -1
        fprintf(fid, '%f', p_tags);
        fprintf(fid, '%s', ', ');
        fclose(fid);
    else
        error('Cannot open output file');
    end
end


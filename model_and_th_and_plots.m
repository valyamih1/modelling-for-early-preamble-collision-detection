fid1 = fopen('th_real.txt', 'w');
fid2 = fopen('th_paper.txt', 'w');
fid3 = fopen('model.txt', 'w');
fclose(fid1);
fclose(fid2);
fclose(fid3);
users = 20;
preambs = 20;
grants = 8;
p_det = 0.8;
p_det_col = 1;
experiments = 1000;
mods = zeros(1,users+1);
grants_total = zeros(1, preambs+1);
%for j = 1:users
    sumResult = 0;
    for i = 0:min(preambs,users)
        result = sum(Pr_d_v(users, preambs, i, p_det));  
        %res_plus_min = result*min(i, grants);
        %sumResult = sumResult + (result*min(i, grants)); %(result*i);
        grants_total = calculate_probability(experiments,users, preambs, i, p_det,p_det_col, mods);
    %sumResult = sumResult * ((1-1/preambs)^(users - 1));
        fid = fopen('th_real.txt', 'a');
        if fid ~= -1
            fprintf(fid, '%f', result);
            fprintf(fid, '%s', ', ');
            fclose(fid);
        else
            error('Cannot open output file');
        end
        %sum_model = sum(mods(2:users+1))/(1000*(users+1));
    end
    %calculate_probability(1000,j, preambs, grants, p_det, p_det_col);
    %papers_theory(j, preambs, grants, p_det, p_det_col);
%end

%for i=0.1:0.1:1
%    calculate_probability(10000,users, preambs, grants, p_det,i);
%    papers_theory(users, preambs, grants, p_det, i);
%end
for i = 1:preambs+1
    grants_total(i) = grants_total(i)/(experiments);
end

th_real = load('th_real.txt');
th_paper = load('th_paper.txt');
model = load('model.txt');
sum_model = sum(grants_total)
sum_th = sum(th_real)
x = 1:length(th_real);
figure; % create a new figure
plot(x,th_real,"g",x, grants_total, "b"); 
%,x, model, "b",x,th_paper,"r",
xlabel('Количество грантов'); 
ylabel('Вероятность D = k'); 
title('Plot'); 
grid on; 

function result = Pr_d_v(n, V, k, p_det)
result = zeros(1,n+1);
for i=k:min(V,n)
    b = nchoosek(i,k);
    b1 = b*(p_det^k);
    b2 = (1-p_det)^(i-k);
    pr_dk = b1*b2;
    res = Pr_k_v(n, V, i);
    result(i+1) =res*pr_dk;
end
end

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

function grants_counter = calculate_probability(number_of_iterations, number_of_users, number_of_preambles, number_of_grants, p_dec, p_dec_col, grants_counter)
    %probability_old = 0;
    optimal_grants = 0;
    pr_d = 0;
    grants_counter = zeros(1, number_of_preambles+1);
    
    %for j = 1:number_of_preambles
        
    successes_every_iteration = zeros(1,number_of_iterations);
    for i = 1:number_of_iterations
        numbers = count_successes(number_of_users, number_of_preambles, number_of_grants, p_dec, p_dec_col);
        successes_every_iteration(i) = numbers;
        grants_counter(numbers+1) = grants_counter(numbers+1) + 1;
    end
    pr_d = sum(grants_counter(2:number_of_grants+1))/(number_of_iterations*(number_of_preambles));
    probability = sum(successes_every_iteration)/((number_of_iterations)*number_of_users);
    required_grants = sum(successes_every_iteration)/(number_of_iterations);

    efficiency = probability/number_of_grants;
    fid = fopen('model.txt', 'w');
    if fid ~= -1
        
        fprintf(fid, '%f', grants_counter);
        fprintf(fid, '%s', ', ');
        fclose(fid);
    else
        error('Cannot open output file');
    end
end

function p_tags = papers_theory(number_of_users, number_of_preambles, number_of_grants, p_dec, p_dec_col)
    p1 = number_of_users * (1/number_of_preambles)*((1 - 1/number_of_preambles)^(number_of_users-1));
    p0 = (1-1/number_of_preambles)^number_of_users;
    p_c = 1 - p1 -p0;
    p_without_errors = ((1-1/number_of_preambles)^(number_of_users - 1)) * min(1, number_of_grants/(p1*number_of_preambles));
    p_with_det = p_dec*((1-1/number_of_preambles)^(number_of_users - 1)) * min(1, number_of_grants/(p1*number_of_preambles*p_dec));
    p_with_dec_col = p_dec*((1-1/number_of_preambles)^(number_of_users - 1)) * min(1, number_of_grants/(p1*number_of_preambles*p_dec + p_c*number_of_preambles*(1-p_dec_col)));
    p_tags = p_with_dec_col*number_of_users *((p1*p_dec)/(p1*p_dec + p_c*(1-p_dec_col)));
    
    fid = fopen('th_paper.txt', 'a');
    if fid ~= -1
        fprintf(fid, '%f', p_tags);
        fprintf(fid, '%s', ', ');
        fclose(fid);
    else
        error('Cannot open output file');
    end
end


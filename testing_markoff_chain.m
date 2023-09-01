users = 20;
preambs = 20;
grants = 3;
p_det = 0.7;
p_det_col = 0.8;


model_results = double_model(50000, preambs, users,  p_det, p_det_col, grants)
theory_results = markoff_chain(preambs, users, p_det, p_det_col, grants)

function theory_res = markoff_chain(N, M, p_det, p_det_col, grants)

    p_1 = M * (1 / N) * (1 - (1 / N))^(M - 1);
    p_0 = (1 - (1 / N))^M;
    p_2 = 1 - p_0 - p_1;
    D = N * p_1 * p_det;
    C = N * p_2 * (1 - p_det_col);

    num_of_nodes = round(D + 1) * round(C + 1);
    probab_matrix = zeros(num_of_nodes, num_of_nodes);
    b = zeros(1, num_of_nodes);
    b(num_of_nodes) = 1;
    counter = 1;
    rest_of_d = D;
    rest_of_c = C;

    for i = num_of_nodes:-1:1
        if rest_of_c <= 0
            rest_of_d = rest_of_d - 1;
            rest_of_c = C;
        end
        if rest_of_d >= 0 && rest_of_c > 0 && counter< num_of_nodes - round(C)
            probab_matrix(i, num_of_nodes - counter) = rest_of_c / (rest_of_c + rest_of_d);
            probab_matrix(i, num_of_nodes - counter - round(C)) = rest_of_d / (rest_of_c + rest_of_d);
            rest_of_c = rest_of_c - 1;
            counter = counter + 1;
        end
    end
    if N<grants
        grants = N;
    end
    if round(D) + round(C)<grants
        grants = round(D) + round(C);
    end
    result = probab_matrix^grants;
    total = b * result;
    total_matrix = zeros(round(D + 1), round(C + 1));
    line_counts = 0;
    rows_counts = 0;
    for i = 1:num_of_nodes
        if rows_counts >= min(round(D)+1, round(C)+1)
            line_counts = line_counts + 1;
            rows_counts = 0;
        end
        total_matrix(line_counts+1, rows_counts+1) = total(i);
        rows_counts = rows_counts + 1;
    end
    %total_matrix = reshape(total, round(D + 1), round(C + 1));

    %disp(result);
    %disp(total);
    %disp(total_matrix);
    theory_res = total_matrix;

end

function model_res = double_model(experiments, N, M, p_dec, p_dec_col, grants)
    success = 0;
    false_success = 0;
    collision = 0;
    detection_error = 0;
    idle = 0;
    for i =1:experiments
        preambls = zeros(1,N);
        %random_pa = randi([1, N], 1, M);
        for j = 1:M
            random_pa = randi([1, N]);
            preambls(random_pa) = preambls(random_pa) + 1;
        end

        for j = 1:N
            decode =rand(1,1);
            detect_col = rand(1,1);
            if preambls(j) == 1 && decode<p_dec
                success = success + 1;
            elseif preambls(j) == 1 && decode>=p_dec
                detection_error = detection_error + 1;
            elseif preambls(j) > 1 && detect_col>=p_dec_col
                false_success = false_success + 1;
            elseif preambls(j) > 1 && detect_col<p_dec_col
                collision = collision + 1;
            elseif preambls(j) == 0
                idle = idle + 1;
            end
        end
        %matrix(success+1, collision+1) = matrix(success+1, collision+1) + 1;
    end
    success = success/experiments;
    false_success = false_success/experiments;
    collision = collision/experiments;
    detection_error = detection_error/experiments;
    idle = idle/experiments;
    success = round(success);
    false_success = round(false_success);
    if N<grants
        grants = N;
    end
    if success + false_success<grants
        grants = success + false_success;
    end
    grants_given = zeros(success + 1, false_success + 1);
    for i =1:experiments
        new_succ = success;
        new_false_succ = false_success;
        for j = 1:grants
            random_grant = randi([1, new_succ + new_false_succ]);
            %grants_given(random_grant) = 1;
            if random_grant <= new_succ
                new_succ = new_succ - 1;
            else
                new_false_succ = new_false_succ -1;
            end
        end
        grants_given(new_succ+1,new_false_succ+1) = grants_given(new_succ+1,new_false_succ+1) + 1;
    end
    for i = 1:length(grants_given)
        for j = 1:size(grants_given, 2)
            grants_given(i,j) = grants_given(i,j)/experiments;
        end
    end
    %disp(grants_given);
    model_res = grants_given;
end

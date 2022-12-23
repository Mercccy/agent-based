clear
clc
middle={};
difficulty_matrix=[];
change_iteration_num=10;
% deadline for agent
t_max_l=100;
% deadline for agent h
t_max_h=100;
possiblity_threshold=0.5;
negotiation_times=[];
difficulty_error=[];
knowledge_error=[];
knowledge_coverage=[];

load final_each_recommandation_collection_question
load each_student_knowledege_matrix
% 每个问题包含的知识点和难度
load collect_knowledge_point
load collect_factor_difficulty_matrix
compute_max_value=collect_factor_difficulty_matrix(12,3);

%对题目的难度进行归一化
max_difficulty=max(cell2mat(collect_knowledge_point(:,6)));
for i=1:size(collect_knowledge_point,1)
    collect_knowledge_point{i,6}=collect_knowledge_point{i,6}/max_difficulty;
end

%教师设置知识点的权重
ans=collect_knowledge_point(:,2:4);
total_knowledge_item=string(zeros(814,3));
for i=1:814
    for j=1:3
        total_knowledge_item(i,j)=cell2mat(ans(i,j));
    end
end
total_knowledge_item=unique(total_knowledge_item);
total_knowledge_item = total_knowledge_item(~ismember(total_knowledge_item,'空'));
for i=1:size(total_knowledge_item,1)
    if strcmp (total_knowledge_item(i,1),'两位数减两位数')==1 || strcmp (total_knowledge_item(i,1),'两位数加两位数')==1
        total_knowledge_item(i,2)=0.8;
    elseif strcmp (total_knowledge_item(i,1),'两位数减一位数')==1 || strcmp (total_knowledge_item(i,1),'两位数加一位数')==1
        total_knowledge_item(i,2)=0.6;
    elseif strcmp (total_knowledge_item(i,1),'枚举填空')==1
        total_knowledge_item(i,2)=0.7;
    else
        total_knowledge_item(i,2)=0;
    end
end


for i=1:size(each_student_knowledege_matrix)
    middle_value=each_student_knowledege_matrix{i,1,1};
    count=1;
    middle_value_matrix={};
    for j=1:size(middle_value,1)
        if middle_value{j,4}~=0
            middle_value_matrix(count,:)=middle_value(j,:);
            count=count+1;
        else
            count=count;
        end
    end
    middle{i,1}=middle_value_matrix;
end
each_student_knowledege_matrix=middle;
[~,~,matirx]=xlsread('迭代所用.xlsx',3);%输入矩阵，只有包括题目做错还是做对
global x equation
ability_factor=[];
leaner_assignment_statement={};
%计算当前对于该类题目学生能力值的变化
for i=1: size(final_each_recommandation_collection_question,1)
    %获取每个学生的做题情况和难度矩阵
    learner_assignment=final_each_recommandation_collection_question{i,1,1};%调入生成的作业矩阵
    difficluty_matrix=[];
    collect_wrong_matrix=[];
    count=1;
    for j=1:size(learner_assignment,1)
        if strcmp(learner_assignment(j,1),'空')==1
            count=count;
        else
            difficluty_matrix(count,1)= learner_assignment{j,6}/6;
            collect_wrong_matrix(count,1)=matirx{i,count};
            count=count+1;
        end
    end
    %计算每个学习者的学习能力
    syms x
    equation=1;
    for j=1:size(collect_wrong_matrix,1)
        u=collect_wrong_matrix(j,1);
        difficulty=difficluty_matrix(j,1);
        middle_equation=exp(x-difficulty)/(1+exp(x-difficulty));
        final_equation=(middle_equation^u)*((1-middle_equation)^(1-u));
        equation=equation*(final_equation);
    end
    equation=diff(equation,x);
    fun=@fitness;
    z=lsqnonlin(fun,0,-3,3);
    z=round(z,2);
    ability_factor(i,1)=z;
    %统计每个学习者的每个知识点出现的次数
    if size(learner_assignment,1)==0
        leaner_assignment_statement{i,1}=[];
        compare_number_of_knowledge_point(i,1)=size(learner_assignment,1);
    else
        knowledge_point=learner_assignment(:,2:4);
        knowledge_point=unique(knowledge_point);
        count=1;
        final_knowledge_point={};
        for j=1:size(knowledge_point,1)
            if strcmp(knowledge_point(j,1),'空')==1
                count=count;
            else
                final_knowledge_point{count,1}=knowledge_point{j,1};
                count=count+1;
            end
        end
        knowledge_point=learner_assignment(:,2:4);
        for l=1:size(final_knowledge_point,1)
            count_total=0;
            for j=1:size(knowledge_point,1)
                for k=1:3
                    if strcmp(final_knowledge_point(l,1),knowledge_point(j,k))==1
                        count_total=count_total+1;
                    end
                end
            end
            final_knowledge_point{l,2}=count_total;
        end
        compare_number_of_knowledge_point(i,1)=size(final_knowledge_point,1);
        %根据matrix的情况建立每个的错误的知识点矩阵
        wrong_knowledge_point=learner_assignment(:,2:4);
        for j=1:size(collect_wrong_matrix,1)
            if collect_wrong_matrix(j,1)==1
                for k=1:3
                    wrong_knowledge_point{j,k}='空';
                end
            end
        end
        %计算每个知识点的错误率
        for l=1:size(final_knowledge_point,1)
            count_total=0;
            for j=1:size(knowledge_point,1)
                for k=1:3
                    if strcmp(final_knowledge_point(l,1),wrong_knowledge_point(j,k))==1
                        count_total=count_total+1;
                    end
                end
            end
            final_knowledge_point{l,3}=count_total;
            final_knowledge_point{l,4}=count_total/final_knowledge_point{l,2};
        end
        middle_value={};
        middle_value=each_student_knowledege_matrix{i,1,1};
        middle={};
        count=1;
        for j=1:size(middle_value,1)
            middle_middle=middle_value{j,1};
            for k=1:size(final_knowledge_point,1)
                if strcmp(middle_middle,final_knowledge_point(k,1))==1 && final_knowledge_point{k,4}~=0
                    middle(count,:)=final_knowledge_point(k,:);
                    count=count+1;
                end
            end
        end
        count=size(middle,1);
        for j=1:size(total_knowledge_item,1)
            ans=total_knowledge_item(j,1);
            for j_1=1:size(middle,1)
                if strcmp(middle(j_1,1),ans)==1 && str2num(total_knowledge_item(j,2))~=0
                    middle{j_1,5}=str2num(total_knowledge_item(j,2));
                elseif strcmp(middle(j_1,1),ans)==1 && str2num(total_knowledge_item(j,2))==0
                    middle{j_1,5}=0;
                end
            end
            if str2num(total_knowledge_item(j,2))~=0 && strcmp(middle(j_1,1),ans)==0
                count=count+1;
                middle{count,1}=total_knowledge_item(j,1);
                middle{count,4}=0;
                middle{count,5}=str2num(total_knowledge_item(j,2));
            end
        end
        
        c_max_h=[];
        c_min_l=[];
        for j=1:size(middle,1)
            value1= middle{j,5};
            value2= middle{j,4};
            c_max_h(j,1)=max(value1,value2);
            c_min_l(j,1)=min(value1,value2);
        end
        
        c_result=c_max_h;
        
        c_max_h=sum(c_max_h)/size(c_max_h,1);
        c_min_l=sum(c_min_l)/size(c_min_l,1);
        
        
        c_min_h=(c_max_h+c_min_l)/2-(c_max_h-c_min_l)/4;
        c_max_l=(c_max_h+c_min_l)/2+(c_max_h-c_min_l)/4;
        
        r_l=0.3+0.1*rand(1,1);%c_max_l
        r_h=0.3+0.1*rand(1,1); %c_min_h
        %beta_l=1.5+2.5*rand(1,1);
        beta_l=1.5+2.5*rand(1,1);
        beta_h=1.5+2.5*rand(1,1);
        
        data=[];
        data_value=[];
        C_t=c_min_l+(c_min_h-c_min_l)*rand(1,1);
%         initial value of agent h
        u_t=C_t-c_min_h;
        a_h=1-(1-r_h)*(1/t_max_h)^(1/beta_h);
        C_t1=c_min_h+(1-a_h)*(c_max_h-c_min_h);
        u_t1=C_t1-c_min_h;
        data(1,1)=u_t;
        data(1,2)=u_t1;
        data_value(1,1)=C_t;
        data_value(1,2)=C_t1;
        % 更新报价
        
        C_t=C_t1;
        t=2;
        
        while (t<(t_max_l+1) && t<(t_max_h+1)&& u_t<u_t1)
            if mod(t,2)==0
%                 buyer
                u_t=c_max_l-C_t;
                a_l=r_l+(1-r_l)*(t/t_max_l)^(1/beta_l);
                C_t1=c_min_l+a_l*(c_max_l-c_min_l);
                u_t1=c_max_l-C_t1;
                data(t,1)=u_t;
                data(t,2)=u_t1;
                data_value(t,1)=C_t;
                data_value(t,2)=C_t1;
                C_t=C_t1;
                t=t+1;
            else
%                 seller
                u_t=C_t-c_min_h;
                a_h=r_h+(1-r_h)*(t/t_max_h)^(1/beta_h);
                C_t1=c_min_h+(1-a_h)*(c_max_h-c_min_h);
                u_t1=C_t1-c_min_h;
                data(t,1)=u_t;
                data(t,2)=u_t1;
                data_value(t,1)=C_t;
                data_value(t,2)=C_t1;
                C_t=C_t1;
                t=t+1;
            end
        end
        negotiation_times_for_knowledge(i,1)=size(data_value,1);    
        final_negotiation_value=data_value(size(data_value,1),1);
        improve_value=(final_negotiation_value-c_max_h)/c_max_h;
        
        for j =1:size(middle,1)
            middle{j,4}=(1+improve_value)*c_result(j,1);
        end

        each_student_knowledege_matrix{i,2}=middle;
    end
    middle_each_student_knowledege_matrix=each_student_knowledege_matrix{i,2,1};
    middle_ability_factor=ability_factor(i,1);
    %计算适应学生能力的题目难度
    for j=1:size(collect_factor_difficulty_matrix,1)
        if middle_ability_factor>=collect_factor_difficulty_matrix(j,1)&&middle_ability_factor<collect_factor_difficulty_matrix(j,2)
            middle_reference_difficulty=collect_factor_difficulty_matrix(j,3);
        end
    end
    a=1.505;
    b=-0.774;
    mood_factor_now=-3+6*rand(1,1);
    mood_factor_now=a/(1+exp(-mood_factor_now))+b;
    middle_reference_difficulty=(1+mood_factor_now)*middle_reference_difficulty/compute_max_value;
    %通过negotiation来获取每一道题对应的难度
    c_max_h=1;
    c_min_l=middle_reference_difficulty;
    %|教师最大期望-教师最低期望| 大
%     c_min_h=(c_max_h+c_min_l)/2-(c_max_h-c_min_l)/4;
    %|教师最大期望-教师最低期望| 小
    c_min_h=(c_max_h+c_min_l)/2-(c_max_h-c_min_l)/4;
% 学生大
%     c_max_l=(c_max_h+c_min_l)/2+(c_max_h-c_min_l)/4;
    %学生小
    c_max_l=(c_max_h+c_min_l)/2+(c_max_h-c_min_l)/4;
    
    data=[];
    data_value=[];
    C_t=c_min_l+(c_min_h-c_min_l)*rand(1,1);
    % initial value of agent h
    u_t=C_t-c_min_h;
    a_h=1-(1-r_h)*(1/t_max_h)^(1/beta_h);
    C_t1=c_min_h+(1-a_h)*(c_max_h-c_min_h);
    u_t1=C_t1-c_min_h;
    data(1,1)=u_t;
    data(1,2)=u_t1;
    data_value(1,1)=C_t;
    data_value(1,2)=C_t1;
    % 更新报价
    C_t=C_t1;
    t=2;
    
    while (t<(t_max_l+1) && t<(t_max_h+1)&& u_t<u_t1)
        if mod(t,2)==0
            %         buyer
            u_t=c_max_l-C_t;
            a_l=r_l+(1-r_l)*(t/t_max_l)^(1/beta_l);
            C_t1=c_min_l+a_l*(c_max_l-c_min_l);
            u_t1=c_max_l-C_t1;
            data(t,1)=u_t;
            data(t,2)=u_t1;
            data_value(t,1)=C_t;
            data_value(t,2)=C_t1;
            C_t=C_t1;
            t=t+1;
        else
            %         seller
            u_t=C_t-c_min_h;
            a_h=r_h+(1-r_h)*(t/t_max_h)^(1/beta_h);
            C_t1=c_min_h+(1-a_h)*(c_max_h-c_min_h);
            u_t1=C_t1-c_min_h;
            data(t,1)=u_t;
            data(t,2)=u_t1;
            data_value(t,1)=C_t;
            data_value(t,2)=C_t1;
            C_t=C_t1;
            t=t+1;
        end
    end
    negotiation_times(i,1)=size(data_value,1);
    final_negotiation_value=data_value(size(data_value,1),1);
    middle_reference_difficulty=final_negotiation_value*compute_max_value;
    difficulty_matrix(i,1)=final_negotiation_value;
    final_collective_question={};
    %依据each_student_knowledge_matrix and difficulty_matrix
    %推荐作业通过content-based algorithm
    %相似矩阵
    similarity_matrix=[];
    %数据预处理
    %推荐的知识点
     if size(each_student_knowledege_matrix{i,2},1)==0
        final_recommandation_matrix={};
    else
        recommended_knowledge_item=each_student_knowledege_matrix{i,2};
        recommended_knowledge_item=recommended_knowledge_item(:,1);
        for j=1:size(recommended_knowledge_item,1)
            recommended_knowledge_item{j,1}=char(recommended_knowledge_item{j,1});
        end
        middle=string(zeros(size(recommended_knowledge_item,1),1));
        for j=1:size(recommended_knowledge_item,1)
            middle(j,1)=cell2mat(recommended_knowledge_item(j,1));
        end
        recommended_knowledge_item=middle;
        %推荐的难度
        recommended_difficulty=difficulty_matrix(i,1);
        
        %计算知识点的相似性
        for j=1:size(collect_knowledge_point,1)
            count=0;
            question_knowledge_item=collect_knowledge_point(j,2:4);
            for number_1=1:size(question_knowledge_item,2)
                for number_2=1:size(recommended_knowledge_item,1)
                    if strcmp (question_knowledge_item{1,number_1},recommended_knowledge_item(number_2,1))==1
                        count=count+1;
                        break
                    end
                end
            end
            similarity_matrix(j,1)=count/size(recommended_knowledge_item,1);
        end
        %计算知识点和推荐的相似度
        %欧几里得距离
        middle_collect_knowledge_point=collect_knowledge_point;
        for j=1:size(middle_collect_knowledge_point,1)
            distance=(1-similarity_matrix(j,1))^2+(recommended_difficulty-middle_collect_knowledge_point{j,6})^2;
            distance=sqrt(distance);
            middle_collect_knowledge_point{j,7}=distance;
        end
        %曼哈顿距离
        for j=1:size(middle_collect_knowledge_point,1)
            distance=abs((1-similarity_matrix(j,1)))+abs(recommended_difficulty-middle_collect_knowledge_point{j,6});
            middle_collect_knowledge_point{j,8}=distance;
        end
        %Cosine Similarity
        for j=1:size(middle_collect_knowledge_point,1)
            %X代表的是baseline，Y代表每道题目的数值,维度1代表知识点，维度2代表难度
            X=[1,recommended_difficulty];
            Y=[similarity_matrix(j,1),middle_collect_knowledge_point{j,6}];
            distance=dot(X,Y)/(norm(X)*norm(Y));
            middle_collect_knowledge_point{j,9}=1-distance;
        end
        %{
        %Pearson相关系数
        for j=1:size(middle_collect_knowledge_point,1)
            X=[1,recommended_difficulty]';
            Y=[similarity_matrix(j,1),middle_collect_knowledge_point{j,6}]';
            fenzi = sum(X .* Y) - (sum(X) * sum(Y)) / length(X);
            fenmu = sqrt((sum(X .^2) - sum(X)^2 / length(X)) * (sum(Y .^2) - sum(Y)^2 / length(X)));
            distance = fenzi / fenmu; 
            middle_collect_knowledge_point{j,10}=distance;
        end
        %}
        %按照欧几里得距离推荐作业
        a2=cell2mat(middle_collect_knowledge_point(:,7));
        [~,ind] = sort(a2);
        middle_collect_knowledge_point_Euclid = middle_collect_knowledge_point(ind,:);
        final_each_recommandation_collection_question_Euclid{i,1}=middle_collect_knowledge_point_Euclid(1:10,1:6);
        %按照曼哈顿距离推荐作业
        a2=cell2mat(middle_collect_knowledge_point(:,8));
        [~,ind] = sort(a2);
        middle_collect_knowledge_point_Manhattan = middle_collect_knowledge_point(ind,:);
        final_each_recommandation_collection_question_Manhattan{i,1}=middle_collect_knowledge_point_Manhattan(1:10,1:6);
        %按照cos距离推荐作业
        a2=cell2mat(middle_collect_knowledge_point(:,9));
        [~,ind] = sort(a2);
        middle_collect_knowledge_point_cos = middle_collect_knowledge_point(ind,:);
        final_each_recommandation_collection_question_cos{i,1}=middle_collect_knowledge_point_cos(1:10,1:6);
        %计算推荐难度和知识点难度之间的差距，1列欧几里得距离，2列曼哈顿距离，3cos距离
        middle=final_each_recommandation_collection_question_Euclid{i,1};
        mean_difficulty=mean(cell2mat(middle(:,6)));
        difficulty_error(i,1)=abs(mean_difficulty-recommended_difficulty)/recommended_difficulty;
        middle=final_each_recommandation_collection_question_Manhattan{i,1};
        mean_difficulty=mean(cell2mat(middle(:,6)));
        difficulty_error(i,2)=abs(mean_difficulty-recommended_difficulty)/recommended_difficulty;
        middle=final_each_recommandation_collection_question_cos{i,1};
        mean_difficulty=mean(cell2mat(middle(:,6)));
        difficulty_error(i,3)=abs(mean_difficulty-recommended_difficulty)/recommended_difficulty;
        %计算知识点准确度，1列欧几里得距离，2列曼哈顿距离，3cos距离
        middle=final_each_recommandation_collection_question_Euclid{i,1};
        middle=middle(:,2:4);
        middle_knowledge_item=string(zeros(10,3));
        for j=1:size(middle_knowledge_item,1)
            for k=1:size(middle_knowledge_item,2)
                middle_knowledge_item(j,k)=cell2mat(middle(j,k));
            end
        end
        middle_knowledge_item=unique(middle_knowledge_item);
        middle_knowledge_item = middle_knowledge_item(~ismember(middle_knowledge_item,'空'));
        count=0;
        for j=1:size(middle_knowledge_item,1)
            for k=1:size(recommended_knowledge_item,1)
                if strcmp (middle_knowledge_item(j,1),recommended_knowledge_item(k,1))==1
                    count=count+1;
                end
            end
        end
        knowledge_error(i,1)=count/size(recommended_knowledge_item,1);
        knowledge_coverage(i,1)=size(middle_knowledge_item,1)/size(total_knowledge_item,1);
        %曼哈顿距离
        middle=final_each_recommandation_collection_question_Manhattan{i,1};
        middle=middle(:,2:4);
        middle_knowledge_item=string(zeros(10,3));
        for j=1:size(middle_knowledge_item,1)
            for k=1:size(middle_knowledge_item,2)
                middle_knowledge_item(j,k)=cell2mat(middle(j,k));
            end
        end
        middle_knowledge_item=unique(middle_knowledge_item);
        middle_knowledge_item = middle_knowledge_item(~ismember(middle_knowledge_item,'空'));
        count=0;
        for j=1:size(middle_knowledge_item,1)
            for k=1:size(recommended_knowledge_item,1)
                if strcmp (middle_knowledge_item(j,1),recommended_knowledge_item(k,1))==1
                    count=count+1;
                end
            end
        end
        knowledge_error(i,2)=count/size(recommended_knowledge_item,1);
        knowledge_coverage(i,2)=size(middle_knowledge_item,1)/size(total_knowledge_item,1);
        %cos距离
        middle=final_each_recommandation_collection_question_cos{i,1};
        middle=middle(:,2:4);
        middle_knowledge_item=string(zeros(10,3));
        for j=1:size(middle_knowledge_item,1)
            for k=1:size(middle_knowledge_item,2)
                middle_knowledge_item(j,k)=cell2mat(middle(j,k));
            end
        end
        middle_knowledge_item=unique(middle_knowledge_item);
        middle_knowledge_item = middle_knowledge_item(~ismember(middle_knowledge_item,'空'));
        count=0;
        for j=1:size(middle_knowledge_item,1)
            for k=1:size(recommended_knowledge_item,1)
                if strcmp (middle_knowledge_item(j,1),recommended_knowledge_item(k,1))==1
                    count=count+1;
                end
            end
        end
        knowledge_error(i,3)=count/size(recommended_knowledge_item,1);
        knowledge_coverage(i,3)=size(middle_knowledge_item,1)/size(total_knowledge_item,1);
     end
end
% % 
% load number_correct_questions
% number_correct_questions=size_1;

 
%  
%  for i=1:size(number_correct_questions,1)
%      number_correct_questions(i,1)=number_correct_questions(i,1)/10;
%  end
 
%  each_student_knowledege_matrix=each_student_knowledege_matrix(:,2:11);
%  count_matrix=[];
%  for i=1:size(each_student_knowledege_matrix,1)
%      for j=1:size(each_student_knowledege_matrix,2)
%          middle_middle=each_student_knowledege_matrix{i,j,1};
%          count_matrix(i,j)=size(middle_middle,1);
%      end
%  end

clear
clc
middle={};
difficulty_matrix=[];
change_iteration_num=10;
% deadline for agent
t_max_l=1000;
% deadline for agent h
t_max_h=1000;
possiblity_threshold=0.5;
negotiation_times=[];
knowledge_error=[];
knowledge_coverage=[];

load final_each_recommandation_collection_question
load each_student_knowledege_matrix
% 每个问题包含的知识点和难度
load collect_knowledge_point
load collect_factor_difficulty_matrix
compute_max_value=collect_factor_difficulty_matrix(12,3);

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
        
        
        c_min_h=(c_max_h+c_min_l)/2-(c_max_h-c_min_l)/8;
        c_max_l=(c_max_h+c_min_l)/2+(c_max_h-c_min_l)/4;
        
        r_l=0.3+0.1*rand(1,1);%c_max_l
        r_h=0.3+0.1*rand(1,1); %c_min_h
        %beta_l=1.5+2.5*rand(1,1);
        beta_l=10+10*rand(1,1);
        beta_h=10+10*rand(1,1);
        
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
    middle_reference_difficulty=middle_reference_difficulty/compute_max_value;
    %通过negotiation来获取每一道题对应的难度
    c_max_h=1;
    c_min_l=middle_reference_difficulty;
    %|教师最大期望-教师最低期望| 大
%     c_min_h=(c_max_h+c_min_l)/2-(c_max_h-c_min_l)/4;
    %|教师最大期望-教师最低期望| 小
    c_min_h=(c_max_h+c_min_l)/2-(c_max_h-c_min_l)/8;
% 学生大
%     c_max_l=(c_max_h+c_min_l)/2+(c_max_h-c_min_l)/4;
    %学生小
    c_max_l=(c_max_h+c_min_l)/2+(c_max_h-c_min_l)/4;
%     r_l=0.3+0.1*rand(1,1);%c_max_l
%     r_h=0.3+0.1*rand(1,1); %c_min_h
%     beta_l=1.5+2.5*rand(1,1);
%     beta_h=1.5+2.5*rand(1,1);
    
    data=[];
    data_value=[];
    C_t=c_min_l;
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
    count=1;
    %按照错误率对题目进行排序
    number=[];
    for j=1:size(middle_each_student_knowledege_matrix,1)
        number(j,1)=middle_each_student_knowledege_matrix{j,4};
        middle_each_student_knowledege_matrix{j,5}='未选择';
    end
    number=sort(number,'descend');
    final_each_student_knowledge={};
    for k=1:size(number,1)
        num=number(k,1);
        for j=1:size(middle_each_student_knowledege_matrix,1)
            if middle_each_student_knowledege_matrix{j,4}==num && strcmp(middle_each_student_knowledege_matrix{j,5},'未选择')==1&& middle_each_student_knowledege_matrix{j,4}~=0
                final_each_student_knowledge(k,:)=middle_each_student_knowledege_matrix(j,:);
                middle_each_student_knowledege_matrix{j,5}='不可选择';
                break
            end
        end
    end
    if size(final_each_student_knowledge,1)==0
        each_student_knowledege_matrix{i,2}=[];
    else
        collect_knowledge_point_1=collect_knowledge_point;
        each_student_knowledege_matrix{i,2}=final_each_student_knowledge(:,1:4);
        final_collective_question={};
        count=1;
        for k=1:size(final_each_student_knowledge,1)
            each_knowledge_point=final_each_student_knowledge{k};
            for j=1:size(collect_knowledge_point_1,1)
                if (strcmp(collect_knowledge_point_1{j,2},each_knowledge_point)==1)||(strcmp(collect_knowledge_point_1{j,3},each_knowledge_point)==1)||(strcmp(collect_knowledge_point_1{j,4},each_knowledge_point)==1)
                    collect_knowledge_point_1{j,7}='可选择';
                end
            end
        end
        for j=1:size(collect_knowledge_point_1,1)
            if strcmp(collect_knowledge_point_1{j,7},'可选择')==1
                final_collective_question(count,:)=collect_knowledge_point_1(j,:);
                count=count+1;
            end
        end
        for j=1:size(final_collective_question,1)
            final_collective_question{j,8}=abs(final_collective_question{j,6}-middle_reference_difficulty);
        end
    end
    middle_each_student_knowledege_matrix=each_student_knowledege_matrix{i,2,1};
    final_recommandation_matrix={};
    each_student_knowledge=final_each_student_knowledge;
    if size(final_each_student_knowledge,1)==0
        final_recommandation_matrix={};
    else
        middle_middle=[];
        for aaa=1:size(final_collective_question,1)
            middle_middle(aaa,1)=final_collective_question{aaa,1};
            middle_middle(aaa,2)=final_collective_question{aaa,8};
            middle_matrix=final_collective_question(aaa,2:4);
            for j=1:size(each_student_knowledge,1)
                middle_value=each_student_knowledge{j,1};
                count=0;
                for k=1:size(middle_matrix,1)
                    if (strcmp(middle_matrix{k,1},middle_value)==1)||(strcmp(middle_matrix{k,2},middle_value)==1)||(strcmp(middle_matrix{k,3},middle_value)==1)
                        count=count+1;
                    end
                end
                middle_middle(aaa,j+2)=count;
            end
        end
        f=[];
        A=[];
        intcon=[];
        Aeq=[];
        for aaa=1:size(middle_middle,1)
            f(aaa)=final_collective_question{aaa,8};
            intcon(aaa)=aaa;
            Aeq(aaa)=1;
            for j=1:size(each_student_knowledge,1)
                A(j,aaa)=middle_middle(aaa,j+2);
                A(j+size(each_student_knowledge,1),aaa)=-middle_middle(aaa,j+2);
            end
        end
        beq=change_iteration_num;%5
        lb=[];
        ub=[];
        lb=zeros(size(final_collective_question,1),1);
        ub=ones(size(final_collective_question,1),1);
        
        total=0;
        for aaa=1:size(each_student_knowledge,1)
            total=total+each_student_knowledge{aaa,4};
        end
        b=[];
        for aaa=1:size(each_student_knowledge,1)
            num=floor(beq/total*each_student_knowledge{aaa,4});
            b(1,aaa)=num+1;
        end
        for aaa=size(each_student_knowledge,1)+1:2*size(each_student_knowledge,1)
            
            b(1,aaa)=-(b(1,aaa-size(each_student_knowledge,1))-1);
        end
        [x_1,fval]=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
        if isempty(fval)==1
            error(i,1)=100;
            final_recommandation_matrix=final_each_recommandation_collection_question{i-1,1};
        else
            count=1;
            error_difficulty=fval/change_iteration_num;
            error_difficulty=error_difficulty/middle_reference_difficulty;
            error(i,1)=error_difficulty;
            for aaa=1:size(x_1,1)
                if x_1(aaa,1)==1
                    final_recommandation_matrix(count,:)=final_collective_question(aaa,:);
                    count=count+1;
                end
            end
        end
    end
    final_each_recommandation_collection_question{i,2}=final_recommandation_matrix;
    %计算知识点的准确率
    %先将知识点转换为string数组
    if size(each_student_knowledege_matrix{i,2,1},1)==0
        knowledge_error(i,1)=0;
        knowledge_coverage(i,1)=0;
    else
        recommended_knowledge_item=each_student_knowledege_matrix{i,2,1};
    recommended_knowledge_item=recommended_knowledge_item(:,1);
    for j=1:size(recommended_knowledge_item,1)
        recommended_knowledge_item{j,1}=char(recommended_knowledge_item{j,1});
    end
    middle=string(zeros(size(recommended_knowledge_item,1),1));
    for j=1:size(recommended_knowledge_item,1)
        middle(j,1)=cell2mat(recommended_knowledge_item(j,1));
    end
    recommended_knowledge_item=middle;
    %计算每个难度的准确率
    middle=final_each_recommandation_collection_question{i,2};
    middle=middle(:,6);
    middle_difficulty=0;
    for j=1:size(middle,1)
        count=middle{j,1}/compute_max_value;
        middle_difficulty=middle_difficulty+count;
    end
    middle_difficulty=middle_difficulty/size(middle,1);
    difficulty_error(i,1)=abs(middle_difficulty-difficulty_matrix(i,1));
    
    %计算每个知识点的准确率
    middle=final_each_recommandation_collection_question{i,2};
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
    end
    
    %计算每个学生答对题目的数量
    possible_value=0;
    for j=1:size(final_recommandation_matrix,1)
        middle_difficulty_value=final_recommandation_matrix{j,6};
        middle_ability_value=(z+6);
        possible_value=possible_value+exp(middle_ability_value-middle_difficulty_value)/(1+exp(middle_ability_value-middle_difficulty_value));
    end
    number_correct_questions(i,2)=possible_value/10;
end

% % 
% load number_correct_questions
% number_correct_questions=size_1;
% 
%  for aaa=1:1
%      change_num=aaa+1;
%      final_correct_question=[];
%      for number_question=1:100
%         test_num=size_1(number_question,1);
%         middle_one=[1,2,3,4,5,6,7,8,9,10];
%         m=middle_one(randperm(numel(middle_one),test_num));
%         for i=1:length(m)
%             middle_one(m(i))=0;
%         end
%         for i=1:length(middle_one)
%             if middle_one(i)==0
%                 middle_one(i)=1;
%             else
%                 middle_one(i)=0;
%             end
%         end
%         final_correct_question=[final_correct_question;middle_one];
%         
%      end
%      matirx=final_correct_question;
%      
%      leaner_assignment_statement={};
%      %计算当前对于该类题目学生能力值的变化
%      for i=1:size(final_each_recommandation_collection_question,1)
%          %获取每个学生的做题情况和难度矩阵
%          learner_assignment=final_each_recommandation_collection_question{i,change_num,1};%调入生成的作业矩阵
%          difficluty_matrix=[];
%          collect_wrong_matrix=[];
%          count=1;
%          for j=1:size(learner_assignment,1)
%              if strcmp(learner_assignment(j,1),'空')==1
%                  count=count;
%              else
%                  difficluty_matrix(count,1)= learner_assignment{j,6}/6;
%                  collect_wrong_matrix(count,1)=matirx(i,count);
%                  count=count+1;
%              end
%          end
%          %计算每个学习者的学习能力
%          syms x
%          equation=1;
%          for j=1:size(collect_wrong_matrix,1)
%              u=collect_wrong_matrix(j,1);
%              difficulty=difficluty_matrix(j,1);
%              middle_equation=exp(x-difficulty)/(1+exp(x-difficulty));
%              final_equation=(middle_equation^u)*((1-middle_equation)^(1-u));
%              equation=equation*(final_equation);
%          end
%          equation=diff(equation,x);
%          fun=@fitness;
%          z=lsqnonlin(fun,0,-3,3);
%          z=round(z,2);
%          ability_factor(i,change_num)=z;
%          %     统计每个学习者的每个知识点出现的次数
%          if size(learner_assignment,1)==0
%              leaner_assignment_statement{i,1}=[];
%              compare_number_of_knowledge_point(i,1)=size(learner_assignment,1);
%          else
%              knowledge_point=learner_assignment(:,2:4);
%              knowledge_point=unique(knowledge_point);
%              count=1;
%              final_knowledge_point={};
%              for j=1:size(knowledge_point,1)
%                  if strcmp(knowledge_point(j,1),'空')==1
%                      count=count;
%                  else
%                      final_knowledge_point{count,1}=knowledge_point{j,1};
%                      count=count+1;
%                  end
%              end
%              knowledge_point=learner_assignment(:,2:4);
%              for l=1:size(final_knowledge_point,1)
%                  count_total=0;
%                  for j=1:size(knowledge_point,1)
%                      for k=1:3
%                          if strcmp(final_knowledge_point(l,1),knowledge_point(j,k))==1
%                              count_total=count_total+1;
%                          end
%                      end
%                  end
%                  final_knowledge_point{l,2}=count_total;
%              end
%              compare_number_of_knowledge_point(i,1)=size(final_knowledge_point,1);
%              %根据matrix的情况建立每个的错误的知识点矩阵
%              wrong_knowledge_point=learner_assignment(:,2:4);
%              for j=1:size(collect_wrong_matrix,1)
%                  if collect_wrong_matrix(j,1)==1
%                      for k=1:3
%                          wrong_knowledge_point{j,k}='空';
%                      end
%                  end
%              end
%              %计算每个知识点的错误率
%              for l=1:size(final_knowledge_point,1)
%                  count_total=0;
%                  for j=1:size(knowledge_point,1)
%                      for k=1:3
%                          if strcmp(final_knowledge_point(l,1),wrong_knowledge_point(j,k))==1
%                              count_total=count_total+1;
%                          end
%                      end
%                  end
%                  final_knowledge_point{l,3}=count_total;
%                  final_knowledge_point{l,4}=count_total/final_knowledge_point{l,2};
%              end
%              middle_value={};
%              middle_value=each_student_knowledege_matrix{i,change_num,1};
%              middle={};
%              count=1;
%              for j=1:size(middle_value,1)
%                  middle_middle=middle_value{j,1};
%                  for k=1:size(final_knowledge_point,1)
%                      if strcmp(middle_middle,final_knowledge_point(k,1))==1 && final_knowledge_point{k,4}~=0
%                          middle(count,:)=final_knowledge_point(k,:);
%                          count=count+1;
%                      end
%                  end
%              end
%              count=size(middle,1);
%              
%              if count~=0
%                  for j=1:size(total_knowledge_item,1)
%                      ans=total_knowledge_item(j,1);
%                      for j_1=1:size(middle,1)
%                          if strcmp(middle(j_1,1),ans)==1 && str2num(total_knowledge_item(j,2))~=0
%                              middle{j_1,5}=str2num(total_knowledge_item(j,2));
%                          elseif strcmp(middle(j_1,1),ans)==1 && str2num(total_knowledge_item(j,2))==0
%                              middle{j_1,5}=0;
%                          end
%                      end
%                      if str2num(total_knowledge_item(j,2))~=0 && strcmp(middle(j_1,1),ans)==0
%                          count=count+1;
%                          middle{count,1}=total_knowledge_item(j,1);
%                          middle{count,4}=0;
%                          middle{count,5}=str2num(total_knowledge_item(j,2));
%                      end
%                  end
%              else
%                  count=0;
%                  for j=1:size(total_knowledge_item,1)
%                     if str2num(total_knowledge_item(j,2))~=0
%                         count=count+1;
%                         middle{count,1}=total_knowledge_item(j,1);
%                         middle{count,4}=0;
%                         middle{count,5}=str2num(total_knowledge_item(j,2));
%                     end
%                  end
%              end
%              c_max_h=[];
%              c_min_l=[];
%              for j=1:size(middle,1)
%                  value1= middle{j,5};
%                  value2= middle{j,4};
%                  c_max_h(j,1)=max(value1,value2);
%                  c_min_l(j,1)=min(value1,value2);
%              end
%              
%              c_result=c_max_h;
%              
%              c_max_h=sum(c_max_h)/size(c_max_h,1);
%              c_min_l=sum(c_min_l)/size(c_min_l,1);
%              
%              
%              c_min_h=(c_max_h+c_min_l)/2-(c_max_h-c_min_l)/4;
%              c_max_l=(c_max_h+c_min_l)/2+(c_max_h-c_min_l)/4;
%              
%              r_l=c_min_l/c_max_l;
%              r_h=c_min_h/c_max_h;
%              beta_l=0.5+0.5*rand(1,1);
%              beta_h=0.5+0.5*rand(1,1);
%              
%              data=[];
%              data_value=[];
%              C_t=c_min_l+(c_min_h-c_min_l)*rand(1,1);
% %              initial value of agent h
%              u_t=C_t-c_min_h;
%              a_h=1-(1-r_h)*(1/t_max_h)^(1/beta_h);
%              C_t1=c_min_h+(1-a_h)*(c_max_h-c_min_h);
%              u_t1=C_t1-c_min_h;
%              data(1,1)=u_t;
%              data(1,2)=u_t1;
%              data_value(1,1)=C_t;
%              data_value(1,2)=C_t1;
%              % 更新报价
%              
%              C_t=C_t1;
%              t=2;
%              
%              while (t<(t_max_l+1) && t<(t_max_h+1)&& u_t<u_t1)
%                  if mod(t,2)==0
% %                      buyer
%                      u_t=c_max_l-C_t;
%                      a_l=r_l+(1-r_l)*(t/t_max_l)^(1/beta_l);
%                      C_t1=c_min_l+a_l*(c_max_l-c_min_l);
%                      u_t1=c_max_l-C_t1;
%                      data(t,1)=u_t;
%                      data(t,2)=u_t1;
%                      data_value(t,1)=C_t;
%                      data_value(t,2)=C_t1;
%                      C_t=C_t1;
%                      t=t+1;
%                  else
% %                      seller
%                      u_t=C_t-c_min_h;
%                      a_h=r_h+(1-r_h)*(t/t_max_h)^(1/beta_h);
%                      C_t1=c_min_h+(1-a_h)*(c_max_h-c_min_h);
%                      u_t1=C_t1-c_min_h;
%                      data(t,1)=u_t;
%                      data(t,2)=u_t1;
%                      data_value(t,1)=C_t;
%                      data_value(t,2)=C_t1;
%                      C_t=C_t1;
%                      t=t+1;
%                  end
%              end
%              negotiation_times_for_knowledge(i,change_num+1)=size(data_value,1);  
%              final_negotiation_value=data_value(size(data_value,1),1);
%              improve_value=(final_negotiation_value-c_max_h)/c_max_h;
%              
%              ans=0;
%              for j =1:size(middle,1)
%                  ans=ans+(1+improve_value)*c_result(j,1);
%                  middle{j,4}=(1+improve_value)*c_result(j,1);
%              end
%              
%              for j =1:size(middle,1)
%                  middle{j,4}=middle{j,4}/ans;
%              end
%              
%              each_student_knowledege_matrix{i,change_num+1}=middle;
%          end
%          middle_each_student_knowledege_matrix=each_student_knowledege_matrix{i,change_num+1,1};
%          middle_ability_factor=ability_factor(i,change_num);
%          %计算适应学生能力的题目难度
%          for j=1:size(collect_factor_difficulty_matrix,1)
%              if middle_ability_factor>=collect_factor_difficulty_matrix(j,1)&&middle_ability_factor<collect_factor_difficulty_matrix(j,2)
%                  middle_reference_difficulty=collect_factor_difficulty_matrix(j,3);
%              end
%          end
%          middle_reference_difficulty=middle_reference_difficulty/compute_max_value;
%          %通过negotiation来获取每一道题对应的难度
%          c_max_h=1;
%          c_min_l=middle_reference_difficulty;
%          
% %         c_min_h=(c_max_h+c_min_l)/2-(c_max_h-c_min_l)/4;
%         c_min_h=(c_max_h+c_min_l)/2-(c_max_h-c_min_l)/8;
%         c_max_l=(c_max_h+c_min_l)/2+(c_max_h-c_min_l)/8;
% %          c_max_l=(c_max_h+c_min_l)/2+(c_max_h-c_min_l)/4;
%          r_l=c_max_l;
%          r_h=c_min_h;
%          beta_l=0.5+0.5*rand(1,1);
%          beta_h=0.5+0.5*rand(1,1);
%          
%          data=[];
%          data_value=[];
%          C_t=c_min_l+(c_min_h-c_min_l)*rand(1,1);
%          % initial value of agent h
%          u_t=C_t-c_min_h;
%          a_h=1-(1-r_h)*(1/t_max_h)^(1/beta_h);
%          C_t1=c_min_h+(1-a_h)*(c_max_h-c_min_h);
%          u_t1=C_t1-c_min_h;
%          data(1,1)=u_t;
%          data(1,2)=u_t1;
%          data_value(1,1)=C_t;
%          data_value(1,2)=C_t1;
%          % 更新报价
%          C_t=C_t1;
%          t=2;
%          
%          while (t<(t_max_l+1) && t<(t_max_h+1)&& u_t<u_t1)
%              if mod(t,2)==0
%                  %         buyer
%                  u_t=c_max_l-C_t;
%                  a_l=r_l+(1-r_l)*(t/t_max_l)^(1/beta_l);
%                  C_t1=c_min_l+a_l*(c_max_l-c_min_l);
%                  u_t1=c_max_l-C_t1;
%                  data(t,1)=u_t;
%                  data(t,2)=u_t1;
%                  data_value(t,1)=C_t;
%                  data_value(t,2)=C_t1;
%                  C_t=C_t1;
%                  t=t+1;
%              else
%                  %         seller
%                  u_t=C_t-c_min_h;
%                  a_h=r_h+(1-r_h)*(t/t_max_h)^(1/beta_h);
%                  C_t1=c_min_h+(1-a_h)*(c_max_h-c_min_h);
%                  u_t1=C_t1-c_min_h;
%                  data(t,1)=u_t;
%                  data(t,2)=u_t1;
%                  data_value(t,1)=C_t;
%                  data_value(t,2)=C_t1;
%                  C_t=C_t1;
%                  t=t+1;
%              end
%          end
%          negotiation_times(i,change_num)=size(data_value,1);
%          final_negotiation_value=data_value(size(data_value,1),1);
% %          difficulty_matrix(i,change_num)=final_negotiation_value;
%          middle_reference_difficulty=final_negotiation_value*compute_max_value;
%          difficulty_matrix(i,change_num)=final_negotiation_value;
%          final_collective_question={};
%          count=1;
%          %按照错误率对题目进行排序
%          number=[];
%          for j=1:size(middle_each_student_knowledege_matrix,1)
%              number(j,1)=middle_each_student_knowledege_matrix{j,4};
%              middle_each_student_knowledege_matrix{j,5}='未选择';
%          end
%          number=sort(number,'descend');
%          final_each_student_knowledge={};
%          for k=1:size(number,1)
%              num=number(k,1);
%              for j=1:size(middle_each_student_knowledege_matrix,1)
%                  if middle_each_student_knowledege_matrix{j,4}==num && strcmp(middle_each_student_knowledege_matrix{j,5},'未选择')==1&& middle_each_student_knowledege_matrix{j,4}~=0
%                      final_each_student_knowledge(k,:)=middle_each_student_knowledege_matrix(j,:);
%                      middle_each_student_knowledege_matrix{j,5}='不可选择';
%                      break
%                  end
%              end
%          end
%          if size(final_each_student_knowledge,1)==0
%              each_student_knowledege_matrix{i,change_num+1}=[];
%          else
%              collect_knowledge_point_1=collect_knowledge_point;
%              each_student_knowledege_matrix{i,change_num+1}=final_each_student_knowledge(:,1:4);
%              final_collective_question={};
%              count=1;
%              for k=1:size(final_each_student_knowledge,1)
%                  each_knowledge_point=final_each_student_knowledge{k};
%                  for j=1:size(collect_knowledge_point_1,1)
%                      if (strcmp(collect_knowledge_point_1{j,2},each_knowledge_point)==1)||(strcmp(collect_knowledge_point_1{j,3},each_knowledge_point)==1)||(strcmp(collect_knowledge_point_1{j,4},each_knowledge_point)==1)
%                          collect_knowledge_point_1{j,7}='可选择';
%                      end
%                  end
%              end
%              for j=1:size(collect_knowledge_point_1,1)
%                  if strcmp(collect_knowledge_point_1{j,7},'可选择')==1
%                      final_collective_question(count,:)=collect_knowledge_point_1(j,:);
%                      count=count+1;
%                  end
%              end
%              for j=1:size(final_collective_question,1)
%                  final_collective_question{j,8}=abs(final_collective_question{j,6}-middle_reference_difficulty);
%              end
%          end
%          each_student_knowledge=final_each_student_knowledge;
%          if size(final_each_student_knowledge,1)==0
%              final_recommandation_matrix={};
%          else
%              middle_middle=[];
%              for bbb=1:size(final_collective_question,1)
%                  middle_middle(bbb,1)=final_collective_question{bbb,1};
%                  middle_middle(bbb,2)=final_collective_question{bbb,8};
%                  middle_matrix=final_collective_question(bbb,2:4);
%                  for j=1:size(each_student_knowledge,1)
%                      middle_value=each_student_knowledge{j,1};
%                      count=0;
%                      for k=1:size(middle_matrix,1)
%                          if (strcmp(middle_matrix{k,1},middle_value)==1)||(strcmp(middle_matrix{k,2},middle_value)==1)||(strcmp(middle_matrix{k,3},middle_value)==1)
%                              count=count+1;
%                          end
%                      end
%                      middle_middle(bbb,j+2)=count;
%                  end
%              end
%              f=[];
%              A=[];
%              intcon=[];
%              Aeq=[];
%              for bbb=1:size(middle_middle,1)
%                  f(bbb)=final_collective_question{bbb,8};
%                  intcon(bbb)=bbb;
%                  Aeq(bbb)=1;
%                  for j=1:size(each_student_knowledge,1)
%                      A(j,bbb)=middle_middle(bbb,j+2);
%                      A(j+size(each_student_knowledge,1),bbb)=-middle_middle(bbb,j+2);
%                  end
%              end
%              beq=change_iteration_num;
%              lb=[];
%              ub=[];
%              lb=zeros(size(final_collective_question,1),1);
%              ub=ones(size(final_collective_question,1),1);
%              
%              total=0;
%              for bbb=1:size(each_student_knowledge,1)
%                  total=total+each_student_knowledge{bbb,4};
%              end
%              b=[];
%              for bbb=1:size(each_student_knowledge,1)
%                  num=floor(beq/total*each_student_knowledge{bbb,4});
%                  b(1,bbb)=num+1;
%              end
%              for bbb=size(each_student_knowledge,1)+1:2*size(each_student_knowledge,1)
%                  
%                  b(1,bbb)=-(b(1,bbb-size(each_student_knowledge,1))-1);
%              end
%              [x_1,fval]=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
%              if isempty(fval)==1
%                  error(i,change_num)=100;
%                  final_recommandation_matrix=final_each_recommandation_collection_question{i,1};
%              else
%                  count=1;
%                  error_difficulty=fval/change_iteration_num;
%                  error_difficulty_num=error_difficulty/middle_reference_difficulty;
%                  error(i,change_num)=error_difficulty_num;
%                  for bbb=1:size(x_1,1)
%                      if x_1(bbb,1)==1
%                          final_recommandation_matrix(count,:)=final_collective_question(bbb,:);
%                          count=count+1;
%                      end
%                  end
%              end
%          end
%          final_each_recommandation_collection_question{i,change_num+1}=final_recommandation_matrix;
%          possible_value=0;
%          for j=1:size(final_recommandation_matrix,1)
%              middle_difficulty_value=final_recommandation_matrix{j,6};
%              middle_ability_value=(z+6);
%              possible_value=possible_value+exp(middle_ability_value-middle_difficulty_value)/(1+exp(middle_ability_value-middle_difficulty_value));
%          end
%          number_correct_questions(i,change_num+2)=possible_value/10;
%      end
%  end
%  
%  for i=1:size(number_correct_questions,1)
%      number_correct_questions(i,1)=number_correct_questions(i,1)/10;
%  end
%  
%  each_student_knowledege_matrix=each_student_knowledege_matrix(:,2:11);
%  count_matrix=[];
%  for i=1:size(each_student_knowledege_matrix,1)
%      for j=1:size(each_student_knowledege_matrix,2)
%          middle_middle=each_student_knowledege_matrix{i,j,1};
%          count_matrix(i,j)=size(middle_middle,1);
%      end
%  end

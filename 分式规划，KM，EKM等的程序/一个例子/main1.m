
example
test_KM
test_PSA
test_approximate_method
test_KM_my
test_iterative_method
test_change_variable_method
test_EKM
test_DA
fprintf('利用课本自带的定理9.1的程序算出的最小值为：%f \n',value_min_inter)
fprintf('利用自己写的定理9.1的程序求得的最小值为：%f \n',value_min)
fprintf('利用Primal simplex algorithm(主单纯型算法)算出的最小值为：%f \n',value_min_primal)
fprintf('利用approximate_method(邻近算法)算出的最小值为：%f \n',value_min_approximate)
fprintf('利用iterative_method(迭代算法)算出的最小值为：%f \n',value_min_iterative)
fprintf('利用change_variable_method(变量变换法)算出的最小值为：%f \n',value_min_change)
fprintf('利用EKM算出的最小值为：%f \n',EKM_min)
fprintf('利用DA算出的最小值为：%f \n',DA_min)
fprintf('\n')
fprintf('利用课本自带的定理9.1的程序算出的最大值为：%f \n',value_max_inter)
fprintf('利用自己写的定理9.1的程序求得的最大值为：%f \n',value_max)
fprintf('利用Primal simplex algorithm(主单纯型算法)算出的最大值为：%f \n',value_max_primal)
fprintf('利用approximate_method(邻近算法)算出的最大值为：%f \n',value_max_approximate)
fprintf('利用iterative_method(迭代算法)算出的最大值为：%f \n',value_max_iterative)
fprintf('利用change_variable_method(变量变换法)算出的最大值为：%f \n',value_max_change)
fprintf('利用EKM算出的最大值为：%f \n',EKM_max)
fprintf('利用DA算出的最大值为：%f \n',DA_max)
